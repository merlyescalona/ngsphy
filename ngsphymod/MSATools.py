import numpy as np
def parseMSAFile(fastapath):
    fastafile=open(fastapath, 'r')
    lines=fastafile.readlines()
    fastafile.close()
    seqdata=[]
    seqDict=dict()
    description=""; seq=""; tmp="";count=1
    for line in lines:
        if line.strip().startswith(">"):
            seqdata+=[line.strip()]
        else:
            if(seqdata[-1].startswith(">")):
                seqdata+=[line.strip()]
            else:
                seqdata[-1]+=line.strip()

    for line in seqdata:
        if not (line.strip()==''):
            if (count%2==0):
                seq=line
                try:
                    test=seqDict[tag]
                except:
                    seqDict[tag]={}

                try:
                    seqDict[tag].update({tmp[2]:{\
                        'description':description,\
                        'sequence':seq\
                    }})
                except:
                    seqDict[tag][tmp[2]]={}
                    seqDict[tag].update({tmp[2]:{\
                        'description':description,\
                        'sequence':seq\
                    }})
                seq=None
                description=None
                tmp=None
            else:
                description=line[1:len(line)]
                tmp=description.split("_")
                tag="{0}_{1}".format(tmp[0], tmp[1])
        count+=1
    return seqDict

def parseMSAFileWithDescriptions(fastapath):
    fastafile=open(fastapath, 'r')
    lines=fastafile.readlines()
    fastafile.close()
    seqDict=dict()
    description=""; seq=""; tmp="";count=1
    seqdata=[]
    for line in lines:
        if line.strip().startswith(">"):
            seqdata+=[line.strip()]
        else:
            if(seqdata[-1].startswith(">")):
                seqdata+=[line.strip()]
            else:
                seqdata[-1]+=line.strip()
    for line in seqdata:
        if not (line.strip()==''):
            if (count%2==0):
                seq=line[0:-1].strip()
                seqDict[tmp]=seq
                seq=None
                description=None
                tmp=None
            else:
                description=line.strip()
                tmp=description[1:len(description)]
        count+=1
    return seqDict

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
