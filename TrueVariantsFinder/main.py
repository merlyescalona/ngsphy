import argparse,csv,datetime,logging,os, sys
import numpy as np
from collections  import Counter
from MELoggingFormatter import MELoggingFormatter as mlf
################################################################################
# CONSTANTS
VERSION=0
MIN_VERSION=0
FIX_VERSION=1
PROGRAM_NAME="true.variants.finder.py"
PROGRAM_NAME_PRETTY="TrueVariantsFider"
AUTHOR="Merly Escalona <merlyescalona@uvigo.es>"
LOG_LEVEL_CHOICES=["DEBUG","INFO","WARNING","ERROR"]
################################################################################
logger=logging.getLogger('tvf')

def main(args):
    path=os.getcwd()
    startTime=datetime.datetime.now()
    endTime=None
    logging.basicConfig(format="%(asctime)s - %(levelname)s:\t%(message)s",\
    datefmt="%d/%m/%Y %I:%M:%S %p",\
    filename="{0}/{2}.{1:%Y}{1:%m}{1:%d}-{1:%H}:{1:%M}:{1:%S}.log".format(\
    path,startTime,\
    PROGRAM_NAME_PRETTY),\
    filemode='a',\
    level=logging.DEBUG)
    ch = logging.StreamHandler()
    ch.setFormatter(logging.Formatter(fmt="%(asctime)s - %(levelname)s:\t%(message)s",\
    datefmt="%d/%m/%Y %I:%M:%S %p"))
    loggerFormatter=mlf(fmt="%(asctime)s - %(levelname)s:\t%(message)s",datefmt="%d/%m/%Y %I:%M:%S %p")
    ch.setFormatter(loggerFormatter)
    ch.setLevel(args.log.upper())
    logger.addHandler(ch)
    filelistPath=os.path.abspath(args.file)
    referenceFilePath=os.path.abspath(args.reference)
    # check file with list of files exists
    if not os.path.exists(filelistPath):
        logger.error("File does not exist. Please check. Exiting.")
        sys.exit()

    # getting all the files
    fl=open(filelistPath)
    filelist=fl.readlines()
    fl.close()
    del fl

    numFiles=len(filelist)
    for index in range(0,numFiles):
        filelist[index]=filelist[index][0:-1]

    # Checking sequence length
    filepath=os.path.abspath(filelist[np.random.randint(0,numFiles,1)])
    f=open(filepath)
    seq=f.readline();seq=f.readline()
    f.close()

    if isFasta(filepath):
        lenSeqs=len(seq.strip())
        numTotalSeqs=numFiles*2
        matrix=np.chararray((numTotalSeqs,lenSeqs), itemsize=1)
        variants=dict()
    else:
        logger.error("One of the sequence files has a wrong file format. Please check. Exiting.")
        sys.exit()


    matrixIndex=0
    logger.info("Number individual files: {}".format(numFiles))
    logger.info("Sequence length: {}".format(lenSeqs))

    individualsRelation=dict()
    # Iterating over filelist
    for indexFile in range(0,len(filelist)):
        filepath=os.path.abspath(filelist[indexFile])
        f=open(filepath)
        lines=f.readlines()
        f.close()
        print(filelist[indexFile])
        for index in range(0,len(lines)):
          lines[index]=lines[index].strip()
        try:
          lines.remove("")
        except:
          pass
        numLinesFile=len(lines)
        numSeqs=len(lines)/2
        index=0
        # If I got here, i can do what i want
        indexSeqs=range(1,numLinesFile,2)
        for index in range(0,numSeqs):
            try:
                val=individualsRelation[os.path.basename(filepath)]
                individualsRelation[os.path.basename(filepath)]+=[matrixIndex]
            except:
                individualsRelation[os.path.basename(filepath)]=[matrixIndex]
            matrix[matrixIndex,:]=list(lines[indexSeqs[index]])
            matrixIndex+=1

    logger.info("Getting variant sites...")
    # Getting position with variants
    for indexCol in range(0,matrix.shape[1]):
        c=Counter(matrix[:,indexCol])
        l=np.unique(matrix[:,indexCol])
        if (len(l)>1):
            variants[str(indexCol)]=dict(c)

    logger.info("Checking reference sequence...")
    refFile=os.path.abspath(referenceFilePath)
    refF=open(refFile)
    refseq=refF.readlines()
    refF.close()
    if isFasta(refFile):
        logger.info("Getting genotypes...")
        refseq=refseq[1][0:-1]
        data=np.chararray(shape=(len(variants),9+numFiles), itemsize=50)
        pos=np.sort([int(newpos) for newpos in variants.keys()])
        dashes=["-"]*len(pos)
        data[:,0]=[os.path.splitext(os.path.basename(refFile))[0]]*len(pos)

        data[:,2]=dashes;data[:,5]=dashes;data[:,6]=dashes;data[:,7]=dashes;
        data[:,8]=["GT"]*len(pos)
        data[:,1]=pos
        ref=[];alt=[]
        for var in range(0,len(pos)):
            varPos=pos[var]
            ref+=[refseq[varPos]]
            alt+=[",".join(variants[str(varPos)].keys())]

        data[:,3]=ref;data[:,4]=alt
        for var in range(0,len(pos)):
          # var is pos
          for ind in range(0,len(individualsRelation.keys())):
            indSeqs=matrix[individualsRelation[individualsRelation.keys()[ind]]]
            varInds=indSeqs[:,str(pos[var])]
            data[var,ind+9]="/".join(varInds)
        print ref,alt
    else:
        logger.error("Reference file has wrong file format. Please check. Exiting.")
        sys.exit()

    outfile=os.path.splitext(os.path.basename(referenceFilePath))[0]
    logger.info("Writing VCF file")
    header="""##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=Integer,Description="Genotype">
    """
    indnames=[os.path.splitext(os.path.basename(ff))[0] for ff in filelist]
    # filename, file_extension = os.path.splitext('/path/to/somefile.ext')
    headerCols=["#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT"]+indnames
    # write it
    with open("{0}.vcf".format(outfile), 'w') as filevcf:
        filevcf.write(header)

    with open("{0}.vcf".format(outfile), 'a') as csvfile:
        writer = csv.writer(csvfile, delimiter="\t")
        writer.writerow(headerCols)
        [writer.writerow(r) for r in data]


#check if fasta
def isFasta(filepath):
    fastaOk=True
    f=open(filepath)
    lines=f.readlines()
    f.close()
    for index in range(0,len(lines)):
        lines[index]=lines[index].strip()
    try:
        lines.remove("")
    except:
        pass
    numLinesFile=len(lines)
    if (numLinesFile % 2 == 0): # It looks like a fasta, has even num of lines
        numSeqs=len(lines)/2
        index=0
        while (index < numSeqs and lines[index*2][0]==">"):
            index=index+1
        if not index >= numSeqs: # there must have been an error
            return False
        # If I got here, i can do what i want
        indexSeqs=range(1,numLinesFile,2)
        lenSeqs=len(lines[1])
        matrix=np.chararray((numSeqs,lenSeqs), itemsize=1)
        for index in range(0,numSeqs):
            matrix[index,:]=list(lines[indexSeqs[index]])
    else:
        return False

    return True




################################################################################
def handlingCmdArguments():
    parser = argparse.ArgumentParser(\
        prog="{0}".format(PROGRAM_NAME),\
        formatter_class=argparse.RawDescriptionHelpFormatter,\
        description='''
TrueVariantsFinder
==================

Assumptions:

- Input file is a multiple sequence alignment file.
- Sequences have the same length.
        '''
        ,\
        epilog="Version {0}.{1}.{2} (Still under development)".format(VERSION,MIN_VERSION,FIX_VERSION),\
        add_help=False
        )

    requiredGroup= parser.add_argument_group('Required arguments')
    requiredGroup.add_argument('-f','--file',metavar='<file_path>',\
        type=str,required=True,\
        help='Path to the multiple sequence alignment (MSA) file.')
    requiredGroup.add_argument('-ref','--reference',metavar='<file_path>',
        type=str,required=True,\
        help='Path to the reference sequence file path.')
    optionalGroup= parser.add_argument_group('Optional arguments')
    optionalGroup.add_argument('-l','--log',metavar='<log_level>', type=str,\
        choices=LOG_LEVEL_CHOICES, default="INFO",\
        help='Specified level of log that will be shown through the standard output. Entire log will be stored in a separate file. Values:{0}. Default: {1}. '.format(LOG_LEVEL_CHOICES,LOG_LEVEL_CHOICES[1]))
    informationGroup= parser.add_argument_group('Information arguments')
    informationGroup.add_argument('-v', '--version',\
        action='version',\
        version='{0}: Version {1}.{2}.{3}'.format(PROGRAM_NAME_PRETTY,VERSION,MIN_VERSION,FIX_VERSION),\
        help="Show program's version number and exit")
    informationGroup.add_argument('-h', '--help',\
        action='store_true',\
        help="Show this help message and exit")
    try:
        tmpArgs = parser.parse_args()
        if (tmpArgs.help): parser.print_help()
    except:
        sys.stdout.write("\n----------------------------------------------------------------------\n")
        parser.print_help()
        sys.exit()
    return tmpArgs



if __name__ == '__main__':
    cmdArgs=handlingCmdArguments()
    main(cmdArgs)
