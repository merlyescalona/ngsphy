#!/usr/bin/home/python
import argparse,datetime,logging,os,sys
import numpy as np
import random as rnd
import Settings as sp
from select import select

################################################################################
class Mating:
    def __init__(self, settings):
        self.appLogger=logging.getLogger('sngsw')
        self.appLogger.debug("(class Mating | __init())")

        # Number of species trees replicates/folder to work with
        self.numSpeciesTrees=0
        self.numSpeciesTreesDigits=0
        # Number of gene trees per replicate
        self.numGeneTrees=0
        # Number of fasta per replicate
        self.numFastaFiles=0
        self.numFastaFilesDigits=0
        # Checking the format of the inputted simphy directory
        if (settings.parser.get("general", "simphy_folder")[-1]=="/"):
            self.projectName=os.path.basename(settings.parser.get("general", "simphy_folder")[0:-1])
        else:
            self.projectName=os.path.basename(settings.parser.get("general", "simphy_folder"))

        self.path=os.path.abspath(settings.parser.get("general", "simphy_folder"))
        # Prefix of the datafiles that contain the FASTA sequences for
        # its corresponding gene tree.
        self.dataprefix=settings.parser.get("general","data_prefix")
        # Number of files to output per genetree
        self.numOutput=settings.parser.get("mating","output-per-individual")
        self.output=""

        # Output folder
        self.output="{0}/individuals".format(self.path)
        # Checking output path
        try:
           os.mkdir(self.output)
           self.output=os.path.abspath(self.output)
        except:
           self.appLogger.info("Output folder ({0}) exists. ".format(self.output))
        print(self.output)


    def checkArgs(self):
        self.appLogger.info("Checking SimPhy folder...")
        matingArgsOk=True
        matingArgsMessageCorrect="Settings are correct, Mating process can be run."
        matingArgsMessageWrong="Something went wrong.\n"
        # Dir exists
        simphydir=os.path.exists(self.path)
        if simphydir:
            self.appLogger.debug("SimPhy folder exists:\t{0}".format(simphydir))
        else:
            matingArgsOk=False
            matingArgsMessageWrong+="\n\tSimPhy folder does not exist."

        # List all the things in the project directory
        fileList=os.listdir(os.path.abspath(self.path))
        for index in range(0,len(fileList)):
            fileList[index]=os.path.abspath(os.path.join(self.path,fileList[index]))

        command = os.path.join(self.path,"{0}.command".format(self.projectName))
        params = os.path.join(self.path,"{0}.params".format(self.projectName))
        db = os.path.join(self.path,"{0}.db".format(self.projectName))

        self.appLogger.debug("SimPhy files (command, params, db)")
        self.appLogger.debug("{0}:\t{1}".format(os.path.basename(db),db in fileList))
        self.appLogger.debug("{0}:\t{1}".format(os.path.basename(command),command in fileList))
        self.appLogger.debug("{0}:\t{1}".format(os.path.basename(params),params in fileList))

        simphyfiles=((command in fileList) and (params in fileList) and(db in fileList))

        # check if  command, db, params files
        if not simphyfiles:
            matingArgsOk=False
            matingArgsMessageWrong+="\n\tSimPhy files do not exist."

        # check how many of them are dirs
        for item in fileList:
            baseitem=os.path.basename(item)
            if (os.path.isdir(os.path.abspath(item)) and  baseitem.isdigit()):
                self.numSpeciesTrees=self.numSpeciesTrees+1
        self.numSpeciesTreesDigits=len(str(self.numSpeciesTrees))
        # check if at least one
        self.appLogger.debug("Num species trees:\t{0}".format(self.numFastaFiles))

        if not (self.numSpeciesTrees>0):
            matingArgsOk=False
            matingArgsMessageWrong+="\n\tNo enough number (at least 1) of species trees replicate/folders:\t{0}".format(self.numSpeciesTrees>0)

        # Checking output path
        if not os.path.exists(self.output):
            self.appLogger.info("Output folder does not exist. Creating: {0} ".format(self.output))
        else:
            self.appLogger.info("Output folder ({0}) exists. ".format(self.output))

        if matingArgsOk:
            return True, matingArgsMessageCorrect
        else:
            return False, matingArgsMessageWrong
    """
    ############################################################################
    #                          Printing Configuration                          #
    ############################################################################
    """
    def print_configuration(self):
        self.appLogger.debug(\
            "\n\t{0}Configuration...{1}\n\tSimPhy project name:\t{2}\n\tSimPhy path:\t{3}\n\tOutput folder:\t{4}\n\tDataset prefix(es) (INDELible):\t{5}\n\tNumber of species trees replicates/folders:\t{6}".format(\
                mof.BOLD,mof.END,\
                self.projectName,\
                self.path,\
                self.output,\
                self.dataprefix,\
                self.numSpeciesTrees\
            )
        )

    """
    ############################################################################
    #                       Iterationg over STs                                #
    ############################################################################
    """
    def iteratingOverST(self):
        self.appLogger.debug("(class Mating | iteratingOverST() )")
        sys.exit()
        # iterate over STs
        for indexST in range(1, self.numSpeciesTrees+1):
            curReplicatePath="{0}/{1:0{2}d}/".format(self.path,indexST, self.numSpeciesTreesDigits)
            self.appLogger.info("Replicate {0}/{2} ({1})".format(indexST, curReplicatePath,self.numSpeciesTrees))
            fasta=0;geneTrees=0
            fileList=os.listdir(curReplicatePath)
            for item in fileList:
                if ("{0}_".format(self.dataprefix) in item) and (".fasta" in item):
                    fasta+=1
                if  ("g_trees" in item) and (".trees" in item):
                    geneTrees+=1

            self.appLogger.warning("Number of fasta files:\t{0}".format(fasta))
            self.numFastaFiles=fasta
            self.numFastaFilesDigits=len(str(self.numFastaFiles))
            if (fasta<1): # Do not have fasta files from the given replicate to work, I'll skip it.
                self.appLogger.warning("Replicate {0}({1}): It is not possible to do the mating for this replicate".format(indexST, curReplicatePath))
                self.appLogger.warning("There are no sequences o there is a missmatch between the prefixes and the number of sequences in the folder.")
            if (geneTrees<1): # I have at least 1 fasta file to work
                    self.ending(False,"This run may not have sense, you are trying to mate sequences, but there are no gene tree files to back that up.")
            else: # I have at least same number of files to work as gene trees (it seems that all makes sense)
                self.appLogger.info("Working directly with FASTA files")
                self.appLogger.info("Mating...")
                for indexLOC in range(1,self.numFastaFiles+1):
                    self.appLogger.debug("Number of FASTA file: {0}".format(indexLOC))
                    newFolder="{0}/{1:0{2}d}/{3:0{4}d}".format(self.output,\
                        indexST, self.numSpeciesTreesDigits,\
                        indexLOC,self.numFastaFilesDigits\
                    )
                    try:
                        os.makedirs(newFolder)
                    except:
                        self.appLogger.info("Folder {0} exists.".format(newFolder))
                    self.mate(indexST,indexLOC)


    """
    ############################################################################
    #                                  Mating                                  #
    ############################################################################
    """
    def mate(self, indexST,indexLOC):
        self.appLogger.debug("( class Mating | mate(self, indexST,indexLOC))")
        indexesDict=dict()

        self.appLogger.debug("Prefix:\t{0}".format(self.dataprefix))
        # 1. Reading the sequences of the file into a dictionary
        fastapath="{0}/{1:0{2}d}/{3}_{4:0{5}d}.fasta".format(\
            self.path,\
            indexST,\
            self.numSpeciesTreesDigits,\
            self.dataprefix,\
            indexLOC,\
            self.numFastaFilesDigits)
        fastafile=open(fastapath, 'r')
        lines=fastafile.readlines()
        fastafile.close()
        seqDict=dict()
        description=""; seq=""; tmp="";count=1
        for line in lines:
            if not (line.strip()==''):
                if (count%2==0):
                    seq=line[0:-1].strip()
                    try:
                        test=seqDict[tmp[0]]
                    except:
                        seqDict[tmp[0]]={}

                    try:
                        seqDict[tmp[0]].update({tmp[2]:{\
                            'description':description,\
                            'sequence':seq\
                        }})
                    except:
                        seqDict[tmp[0]][tmp[2]]={}
                        seqDict[tmp[0]].update({tmp[2]:{\
                            'description':description,\
                            'sequence':seq\
                        }})

                    seq=None
                    description=None
                    tmp=None
                else:
                    description=line[0:-1].strip()
                    tmp=description[1:len(description)].split("_")
                count=count+1
        """
        ####################################################################
        # Mating
        ####################################################################
        """
        species=seqDict.keys()
        pair=();subkeys=()
        # 2. Generating the indexes for the current gene tree if not existent
        #    For all the prefixes inputted.
        if indexesDict=={}:
            # If here, I'm in the first prefix, so I need to generate the
            # indexes
            test=[]; subset=set()
            self.appLogger.info("Generating indexes...")
            indId=1
            for key in species:
                self.appLogger.debug("Species {0}".format(key))
                subkeys=seqDict[key].keys()
                if not key=="0":
                    while not subkeys==[]:
                        try:
                            test=indexesDict[key]
                            self.appLogger.debug("Indexes dictionary is NOT empty")
                            for item in test:
                                subset.add(item[0])
                                subset.add(item[1])
                        except KeyError:
                            self.appLogger.debug("Indexes dictionary IS empty")
                            indexesDict[key]=[]

                        newSubkeys=set(subkeys)-set(indexesDict[key])
                        self.appLogger.debug("{0} | {1}".format(subkeys,test))
                        try:
                            pairs=rnd.sample(newSubkeys,2)
                            logging.debug("Pairs 1: {0}\t2:{1}".format(pairs[0],pairs[0]))
                            indexesDict[key]+=[(indId, pairs[0], pairs[1])]
                            subkeys.pop(subkeys.index(pairs[0]))
                            subkeys.pop(subkeys.index(pairs[1]))
                            indId+=1
                        except:
                            subkeys.pop()
                else:
                    subkeys=[]

            self.writeIndexesDictionary(indexST,indexLOC, indexesDict)
        else:
            self.appLogger.debug("Indexes dictionary is generated, just writing files for the other prefixes...")
        """
        ####################################################################
        # Getting number of individuals that will be generated
        ####################################################################
        """
        numSeqs=0
        for key in species:
            subspecies=seqDict[key].keys()
            for subkey in subspecies:
                numSeqs+=len(seqDict[key][subkey])

        numInds=np.trunc(numSeqs/2)+1
        currentInd=1
        self.appLogger.info("Writing {1} individuals from {0} number of sequences.".format(numSeqs,numInds))
        outputFolder="{0}/{1:0{2}d}/{3:0{4}d}".format(self.output,\
            indexST,self.numSpeciesTreesDigits,\
            indexLOC,self.numFastaFilesDigits\
        )

        """
        ####################################################################
        # Writing individuals into separate files
        ####################################################################
        """

        for key in species:
            seq1="";des1="";seq2="";des2="";
            if key == "0":
                # For the outgroup
                shortDesc=seqDict[key]["0"]["description"][1:len(seqDict[key]["0"]["description"])]
                des1=">{0}:{1:0{6}d}:{2:0{7}d}:{3}:{4}:{5}".format(self.projectName,\
                    indexST,indexLOC,self.dataprefix,0,shortDesc,\
                    self.numSpeciesTreesDigits,self.numFastaFilesDigits\
                    )
                seq1=seqDict[key]["0"]["sequence"]

                # Writing the outgroup -----------------------------------------
                if (self.numOutput == 1):
                    indFilename="{0}/{1}_{2:0{6}d}_{3:0{7}d}_{4}_{5}.fasta".format(\
                        outputFolder,self.projectName,indexST,indexLOC,\
                        self.dataprefix,0,\
                        self.numSpeciesTreesDigits,self.numFastaFilesDigits\
                    )
                    indFile=open(indFilename, "w")
                    indFile.write("{0}_S1\n{1}\n{2}_S2\n{3}\n".format(des1,seq1,des1,seq1))
                    indFile.close()
                if (self.numOutput == 2):
                    # Strand 1 -----------------------------------------
                    indFilename="{0}/{1}_{2:0{7}d}_{3:0{8}d}_{4}_{5}_{6}_S1.fasta".format(\
                        outputFolder,self.projectName,indexST,indexLOC,\
                        self.dataprefix,0,shortDesc,\
                        self.numSpeciesTreesDigits,self.numFastaFilesDigits\
                    )
                    indFile=open(indFilename, "w")
                    indFile.write("{0}_S1\n{1}\n".format(des1,seq1))
                    indFile.close()
                    # Strand 2 -----------------------------------------
                    indFilename="{0}/{1}_{2:0{7}d}_{3:0{8}d}_{4}_{5}_{6}_S2.fasta".format(\
                        outputFolder,self.projectName,indexST,indexLOC,\
                        self.dataprefix,0,shortDesc,\
                        self.numSpeciesTreesDigits,self.numFastaFilesDigits\
                    )
                    indFile=open(indFilename, "w")
                    indFile.write("{0}_S2\n{1}\n".format(des1,seq1))
                    indFile.close()
                if (self.numOutput == 3):
                    # General file -------------------------------------
                    indFilename="{0}/{1}_{2:0{6}d}_{3:0{7}d}_{4}_{5}.fasta".format(\
                        outputFolder,self.projectName,indexST,indexLOC,\
                        self.dataprefix,0,\
                        self.numSpeciesTreesDigits,self.numFastaFilesDigits\
                    )
                    indFile=open(indFilename, "w")
                    indFile.write("{0}_S1\n{1}\n{2}_S2\n{3}\n".format(des1,seq1,des1,seq1))
                    indFile.close()
                    # Strand 1 -----------------------------------------
                    indFilename="{0}/{1}_{2:0{7}d}_{3:0{8}d}_{4}_{5}_{6}_S1.fasta".format(\
                        outputFolder,self.projectName,indexST,indexLOC,\
                        self.dataprefix,0, shortDesc,\
                        self.numSpeciesTreesDigits,self.numFastaFilesDigits\
                    )
                    indFile=open(indFilename, "w")
                    indFile.write("{0}_S1\n{1}\n".format(des1,seq1))
                    indFile.close()
                    # Strand 2 -----------------------------------------
                    indFilename="{0}/{1}_{2:0{7}d}_{3:0{8}d}_{4}_{5}_{6}_S2.fasta".format(\
                        outputFolder,self.projectName,indexST,indexLOC,\
                        self.dataprefix,0, shortDesc,\
                        self.numSpeciesTreesDigits,self.numFastaFilesDigits\
                    )
                    indFile=open(indFilename, "w")
                    indFile.write("{0}_S2\n{1}\n".format(des1,seq1))
                    indFile.close()

            else:
                # for all species except the outgroup
                pairs=indexesDict[key]
                for pair in pairs:
                    # Extracting info from the dictionary
                    currentInd=pair[0];pair1=pair[1]; pair2=pair[2];
                    # Organizing strings
                    logging.debug("{0}|{1}".format(pair, key))
                    seq1=seqDict[key][pair1]["sequence"]
                    seq2=seqDict[key][pair2]["sequence"]
                    shortDesc1=seqDict[key][pair1]["description"][1:len(seqDict[key][pair1]["description"])]
                    des1=">{0}:{1:0{6}d}:{2:0{7}d}:{3}:{4}:{5}_S1".format(self.projectName,\
                        indexST,indexLOC,self.dataprefix,currentInd,shortDesc1,\
                        self.numSpeciesTreesDigits,self.numFastaFilesDigits\
                    )
                    shortDesc2=seqDict[key][pair2]["description"][1:len(seqDict[key][pair2]["description"])]
                    des2=">{0}:{1:0{6}d}:{2:0{7}d}:{3}:{4}:{5}_S2".format(\
                        self.projectName,indexST,indexLOC,self.dataprefix,\
                        currentInd,shortDesc2,
                        self.numSpeciesTreesDigits,self.numFastaFilesDigits\
                    )

                    if (self.numOutput == 1):
                        indFilename="{0}/{1}_{2:0{6}d}_{3:0{7}d}_{4}_{5}.fasta".format(\
                            outputFolder,self.projectName,indexST,indexLOC,\
                            self.dataprefix,currentInd,\
                            self.numSpeciesTreesDigits,self.numFastaFilesDigits\
                        )
                        indFile=open(indFilename, "w")
                        indFile.write("{0}\n{1}\n{2}\n{3}\n".format(des1,seq1,des2,seq2))
                        indFile.close()
                    if (self.numOutput == 2):
                        # Strand 1 -----------------------------------------
                        indFilename="{0}/{1}_{2:0{7}d}_{3:0{8}d}_{4}_{5}_{6}_S1.fasta".format(\
                            outputFolder,self.projectName,indexST,indexLOC,\
                            self.dataprefix,currentInd, shortDesc1,\
                            self.numSpeciesTreesDigits,self.numFastaFilesDigits\
                        )
                        indFile=open(indFilename, "w")
                        indFile.write("{0}\n{1}\n".format(des1,seq1))
                        indFile.close()
                        # Strand 2 -----------------------------------------
                        indFilename="{0}/{1}_{2:0{7}d}_{3:0{7}d}_{4}_{5}_{6}_S2.fasta".format(\
                            outputFolder,self.projectName,indexST,indexLOC,\
                            self.dataprefix,currentInd, shortDesc2,\
                            self.numSpeciesTreesDigits,self.numFastaFilesDigits\
                        )
                        indFile=open(indFilename, "w")
                        indFile.write("{0}\n{1}\n".format(des2,seq2))
                        indFile.close()
                    if (self.numOutput == 3):
                        # General file -------------------------------------
                        indFilename="{0}/{1}_{2:0{6}d}_{3:0{7}d}_{4}_{5}.fasta".format(\
                            outputFolder,self.projectName,indexST,indexLOC,\
                            self.dataprefix,currentInd,\
                            self.numSpeciesTreesDigits,self.numFastaFilesDigits\
                        )
                        indFile=open(indFilename, "w")
                        indFile.write("{0}\n{1}\n{2}\n{3}\n".format(des1,seq1,des2,seq2))
                        indFile.close()
                        # Strand 1 -----------------------------------------
                        indFilename="{0}/{1}_{2:0{7}d}_{3:0{8}d}_{4}_{5}_{6}_S1.fasta".format(\
                            outputFolder,self.projectName,indexST,indexLOC,\
                            self.dataprefix,currentInd, shortDesc1,\
                            self.numSpeciesTreesDigits,self.numFastaFilesDigits\
                        )
                        indFile=open(indFilename, "w")
                        indFile.write("{0}\n{1}\n".format(des1,seq1))
                        indFile.close()
                        # Strand 2 -----------------------------------------
                        indFilename="{0}/{1}_{2:0{7}d}_{3:0{8}d}_{4}_{5}_{6}_S2.fasta".format(\
                            outputFolder,self.projectName,indexST,indexLOC,\
                            self.dataprefix,currentInd, shortDesc2,\
                            self.numSpeciesTreesDigits,self.numFastaFilesDigits\
                        )
                        indFile=open(indFilename, "w")
                        indFile.write("{0}\n{1}\n".format(des2,seq2))
                        indFile.close()

                    del seqDict[key][pair1]
                    del seqDict[key][pair2]
                if not seqDict[key]=={}:
                    self.appLogger.warning("Number of individuals (sequences) generated per species is odd.")
                    self.appLogger.warning("Sequence {0} from species {1} will not be paired.".format( seqDict[key].keys(), key))
                    for item in seqDict[key].keys():
                        del seqDict[key][item]

        return None


    """
    #########################################################################
    # writeIndexesDictionary
    #########################################################################
    """
    def writeIndexesDictionary(self, indexST,indexLOC, indexDict):
        self.appLogger.info("Writing indexes into file...")
        indexFilename="{0}/{1}.mating".format(self.path,self.projectName)
        if not os.path.isfile(indexFilename):
            indexFile=open(indexFilename,"w")
            indexFile.write("indexST,indexLOC,indID,speciesID,mateID1,mateID2\n")
            indexFile.close()

        indexFile=open(indexFilename,"a")
        for speciesID in indexDict.keys():
            for item in indexDict[speciesID]:
                indexFile.write("{0:0{1}d},{2:0{3}d},{4},{5},{6},{7}\n".format(\
                    indexST,self.numSpeciesTreesDigits,
                    indexLOC,self.numFastaFilesDigits,\
                    item[0],speciesID,item[1],item[2]))

        indexFile.write("{0:0{1}d},{2:0{3}d},{4},{5},{6},{7}\n".format(\
            indexST,self.numSpeciesTreesDigits,
            indexLOC,self.numFastaFilesDigits,\
            0,0,0,0))
        indexFile.close()

    # def run(self):
    #     self.checkArgs()
    #     self.print_configuration()
    #     self.iteratingOverST()
