#!/usr/bin/home/python
import argparse,copy,datetime,logging,os,sys, sqlite3
import numpy as np
import random as rnd
import Settings as sp
from select import select

################################################################################
class Mating:
    def __init__(self, settings):
        self.appLogger=logging.getLogger('sngsw')
        self.appLogger.info("Mating: Run started")
        self.settings=settings
        # Number of species trees replicates/folder to work with
        self.numSpeciesTrees=0
        self.numSpeciesTreesDigits=0
        # Number of fasta per replicate
        self.numFastaFiles=0
        self.numFastaFilesDigits=0
        # Checking the format of the inputted simphy directory
        self.path=os.path.abspath(self.settings.parser.get("general", "simphy_folder"))
        if (self.settings.parser.get("general", "simphy_folder")[-1]=="/"):
            self.projectName=os.path.basename(self.settings.parser.get("general", "simphy_folder")[0:-1])
        else:
            self.projectName=os.path.basename(self.settings.parser.get("general", "simphy_folder"))

        # outgroup
        self.outgroup=False
        if (self.settings.parser.has_option("general","og")):
            self.outgroup=self.settings.parser.getboolean("general","og")
        if (self.settings.parser.has_option("general","outgroup")):
            self.outgroup=self.settings.parser.getboolean("general","outgroup")


        # Prefix of the datafiles that contain the FASTA sequences for
        # its corresponding gene tree.
        self.dataprefix=self.settings.parser.get("general","data_prefix")
        self.output=""

        self.output=os.path.abspath(self.settings.parser.get("general","output_folder_name"))
        self.appLogger.info("Generating output folder ({0}) and subfolders (individuals,mating,reads)".format(self.output))
        self.outputmating=os.path.abspath("{0}/mating".format(self.settings.parser.get("general","output_folder_name")))


    def checkArgs(self):
        self.appLogger.info("Checking SimPhy folder...")
        matingArgsMessageCorrect="Settings are correct, Mating process can be run."
        matingArgsMessageWrong="Something went wrong.\n"
        # Dir exists
        simphydir=os.path.exists(self.path)
        if simphydir:
            self.appLogger.debug("SimPhy folder exists:\t{0}".format(simphydir))
        else:
            matingArgsMessageWrong+="\n\tSimPhy folder does not exist."
            return False, matingArgsMessageWrong

        # List all the things in the project directory
        fileList=os.listdir(os.path.abspath(self.path))
        for index in range(0,len(fileList)):
            fileList[index]=os.path.abspath(os.path.join(self.path,fileList[index]))

        self.command = os.path.join(self.path,"{0}.command".format(self.projectName))
        self.params = os.path.join(self.path,"{0}.params".format(self.projectName))
        self.db = os.path.join(self.path,"{0}.db".format(self.projectName))

        self.appLogger.debug("SimPhy files (command, params, db)")
        self.appLogger.debug("{0}:\t{1}".format(os.path.basename(self.db),self.db in fileList))
        self.appLogger.debug("{0}:\t{1}".format(os.path.basename(self.command),self.command in fileList))
        self.appLogger.debug("{0}:\t{1}".format(os.path.basename(self.params),self.params in fileList))

        simphyfiles=((self.command in fileList) and (self.params in fileList) and(self.db in fileList))

        # check if  command, db, params files
        if not simphyfiles:
            matingArgsMessageWrong+="\n\tSimPhy files do not exist."
            return False, matingArgsMessageWrong

        # check how many of them are dirs
        for item in fileList:
            baseitem=os.path.basename(item)
            if (os.path.isdir(os.path.abspath(item)) and  baseitem.isdigit()):
                self.numSpeciesTrees=self.numSpeciesTrees+1
        self.numSpeciesTreesDigits=len(str(self.numSpeciesTrees))
        # check if at least one

        self.appLogger.debug("Num species trees:\t{0}".format(self.numSpeciesTrees))
        if not (self.numSpeciesTrees>0):
            matingArgsMessageWrong+="\n\tNot enough number of species tree replicates (at least 1):\t{0}".format(self.numSpeciesTrees>0)
            return False, matingArgsMessageWrong

        # Checking output path
        self.appLogger.info("Checking output folder...")
        # Checking output path
        self.outputinds="{0}/individuals".format(self.output)
        try:
            self.appLogger.info("Creating output folder: {0} ".format(self.output))
            self.output=os.path.abspath(self.output)
            os.mkdir(self.output)
        except:
           self.appLogger.debug("Output folder ({0}) exists. ".format(self.output))

        try:
            self.appLogger.info("Generated individuals/")
            os.mkdir("{0}/individuals".format(self.output))
        except:
            self.appLogger.debug("Output folder exists ({0})".format(self.output))

        try:
            self.appLogger.info("Generated reads/")
            os.mkdir("{0}/reads".format(self.output))
        except:
            self.appLogger.debug("Output folder exists ({0})".format(self.output))

        try:
            self.appLogger.info("Generated mating/")
            os.mkdir("{0}/mating".format(self.output))
        except:
            self.appLogger.debug("Output folder exists ({0})".format(self.output))


        self.filteredSts=self.stWithEvenNumberIndsPerSpecies()
        gtperstOK,message=self.getNumGTST()
        if (not gtperstOK):
            return gtperstOK,message
        self.settings.parser.set("general","filtered_ST",self.filteredSts)
        matingArgsMessageCorrect+="\n{0}".format(message)
        return True, matingArgsMessageCorrect



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

    def stWithEvenNumberIndsPerSpecies(self):
        query="select SID from Species_Trees WHERE Ind_per_sp % 2 = 0"
        con = sqlite3.connect(self.db)
        res=con.execute(query).fetchall()
        con.close()
        res=[item for sublist in res for item in sublist]
        return res


    def getNumGTST(self):
        self.numFASTAperST=np.repeat(0,self.numSpeciesTrees)
        for indexST in self.filteredSts:
            curReplicatePath="{0}/{1:0{2}d}/".format(self.path,indexST, self.numSpeciesTreesDigits)
            self.numFastaFiles=0;numGeneTrees=0
            fileList=os.listdir(curReplicatePath)
            # check composition of the current indexST folder
            for item in fileList:
                if ("{0}_".format(self.dataprefix) in item) and (".fasta" in item):
                    self.numFastaFiles+=1
                if  ("g_trees" in item) and (".trees" in item):
                    numGeneTrees+=1
            self.numFASTAperST[indexST-1]=self.numFastaFiles
            self.appLogger.warning("Number of fasta files:\t{0}".format(self.numFastaFiles))
            self.numFastaFilesDigits=len(str(self.numFastaFiles))
            if (self.numFastaFiles<1): # Do not have fasta files from the given replicate to work, I'll skip it.
                self.appLogger.warning("Replicate {0}({1}): It is not possible to do the mating for this replicate".format(indexST, curReplicatePath))
                self.appLogger.warning("There are no sequences o there is a missmatch between the prefixes and the number of sequences in the folder.")
            if (numGeneTrees<1):
                return False,"Trying to mate sequences, but there are no gene tree files to back that up. Please, finish the SimPhy run and try again afterwards."
        return True,"Got number of gene trees per species trees"

    """
    ############################################################################
    #                       Iterationg over STs                                #
    ############################################################################
    """
    def iteratingOverST(self):
        self.appLogger.debug("(class Mating | iteratingOverST() )")
        # iterate over STs
        for indexST in self.filteredSts:
            curReplicatePath="{0}/{1:0{2}d}/".format(self.outputinds,indexST, self.numSpeciesTreesDigits)
            self.appLogger.info("Replicate {0}/{2} ({1})".format(indexST, curReplicatePath,self.numSpeciesTrees))
            self.appLogger.info("Mating...")
            # generating and writing mating table
            matingTable=self.generateMatingTable(indexST)
            self.writeMatingTable(indexST,matingTable)
            # iterating over the number of gts per st
            for indexLOC in range(0,self.numFASTAperST[indexST-1]):
                self.appLogger.debug("Number of FASTA file: {0}".format(indexLOC+1))
                newFolder="{0}/{1:0{2}d}/{3:0{4}d}".format(self.outputinds,\
                    indexST, self.numSpeciesTreesDigits,\
                    indexLOC+1,self.numFastaFilesDigits\
                )
                # creating folder belonging to specific GT
                try:
                    os.makedirs(newFolder)
                except:
                    self.appLogger.debug("Folder {0} exists.".format(newFolder))
                # parsingMSA file
                seqDict=self.parseMSAFile(indexST,indexLOC+1)
                self.mateAll(indexST,matingTable,seqDict)

    def generateMatingTable(self,indexST):
        # missing outgroup
        self.appLogger.info("Connecting to the db...")
        con = sqlite3.connect(self.db)
        query="select SID, Leaves, Ind_per_sp from Species_Trees WHERE SID={0}".format(indexST)
        res=con.execute(query).fetchone()
        con.close()
        indexST=res[0];leaves=res[1];nIndsPerSp=res[2]
        # by default there are no outgroups, if there are, this
        # value will change
        nInds=leaves/2
        mates=[]
        if self.outgroup:
            mates+=[(indexST,0,0,0)]
            nInds=(leaves-1)/2
        inds=range(0,nIndsPerSp)
        species=range(1,leaves)
        self.appLogger.debug("indexST: {0} / inds:{1} ".format(indexST,inds))
        # I'm always assuming there's an outgroup
        for sp in species:
            t=copy.deepcopy(inds)
            while not t==[]:
                p1=0;p2=0
                try:
                    p1=t.pop(rnd.sample(range(0,len(t)),1)[0])
                    p2=t.pop(rnd.sample(range(0,len(t)),1)[0])
                except Exception as e:
                    break
                pair=(indexST,sp,p1,p2)
                mates+=[pair]
                self.appLogger.debug("Pair generated: {0}".format(pair))
        return mates

    def writeMatingTable(self,indexST,matingTable):
        # mating table
        # indexST,SP,ind-tip1,ind-tip2
        # If i have outgroups mate info should be in the table already
        self.appLogger.debug("Writing indexes into file...")
        indexFilename="{0}/{1}.{2:0{3}d}.mating.csv".format(self.outputmating,self.projectName, indexST, self.numSpeciesTreesDigits)
        print(indexFilename)
        if not os.path.isfile(indexFilename):
            indexFile=open(indexFilename,"w")
            indexFile.write("indexST,indID,speciesID,mateID1,mateID2\n")
            indexFile.close()
        indexFile=open(indexFilename,"a")
        for indexRow in range(0,len(matingTable)):
            indexST=indexST
            indID=indexRow
            speciesID=matingTable[indexRow][1]
            mateID1=matingTable[indexRow][2]
            mateID2=matingTable[indexRow][3]
            indexFile.write("{0:0{1}d},{2},{3},{4},{5}\n".format(\
            indexST,self.numSpeciesTreesDigits,
                indID,speciesID,mateID1,mateID2))
        indexFile.close()

    def parseMSAFile(self, indexST, indexLOC):
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
            count+=1
        return seqDict

    def mateAll(self,indexST,matingTable,sequenceDictionary):
        # Proper mating
        originalDict=copy.deepcopy(sequenceDictionary)
        species=originalDict.keys()
        numSeqs=0
        for key in species:
            subspecies=originalDict[key].keys()
            for subkey in subspecies:
                numSeqs+=len(originalDict[key][subkey])
        numInds=np.trunc(numSeqs/2)+1

        self.appLogger.debug("Writing {1} individuals from {0} number of sequences.".format(numSeqs,numInds))

        for indexLOC in range(1,self.numFastaFiles):
            outputFolder="{0}/{1:0{2}d}/{3:0{4}d}".format(self.outputinds,\
                indexST,self.numSpeciesTreesDigits,\
                indexLOC,self.numFastaFilesDigits\
            )
            self.appLogger.debug("Output folder:".format(outputFolder))
            """
            ####################################################################
            # Writing individuals into separate files
            ####################################################################
            """
            self.appLogger.error(species)
            seqDict=copy.deepcopy(originalDict)
            seq1="";des1="";seq2="";des2="";
            # for all species except the outgroup
            self.appLogger.error(outputFolder)
            for currentInd in range(0,len(matingTable)):
                # Extracting info from the dictionary
                st=str(matingTable[currentInd][0])
                sp=str(matingTable[currentInd][1])
                pair1=str(matingTable[currentInd][2])
                pair2=str(matingTable[currentInd][3])
                # Organizing strings
                self.appLogger.debug("{0}|{1}-{2}".format(sp, pair1,pair2))
                seq1=seqDict[sp][pair1]["sequence"]
                seq2=seqDict[sp][pair2]["sequence"]
                shortDesc1=seqDict[sp][pair1]["description"][1:len(seqDict[sp][pair1]["description"])]
                des1=">{0}:{1:0{6}d}:{2:0{7}d}:{3}:{4}:{5}_S1".format(self.projectName,\
                    indexST,indexLOC,self.dataprefix,currentInd,shortDesc1,\
                    self.numSpeciesTreesDigits,self.numFastaFilesDigits\
                )
                shortDesc2=seqDict[sp][pair2]["description"][1:len(seqDict[sp][pair2]["description"])]
                des2=">{0}:{1:0{6}d}:{2:0{7}d}:{3}:{4}:{5}_S2".format(\
                    self.projectName,indexST,indexLOC,self.dataprefix,\
                    currentInd,shortDesc2,
                    self.numSpeciesTreesDigits,self.numFastaFilesDigits\
                )

                indFilename="{0}/{1}_{2:0{6}d}_{3:0{7}d}_{4}_{5}.fasta".format(\
                    outputFolder,self.projectName,indexST,indexLOC,\
                    self.dataprefix,currentInd,\
                    self.numSpeciesTreesDigits,self.numFastaFilesDigits\
                )
                indFile=open(indFilename, "w")
                indFile.write("{0}\n{1}\n{2}\n{3}\n".format(des1,seq1,des2,seq2))
                indFile.close()
                if not sp== "0":
                    del seqDict[sp][pair1]
                    del seqDict[sp][pair2]
                else:
                    del seqDict["0"]
            if not seqDict[sp]=={}:
                self.appLogger.warning("Number of individuals (sequences) generated per species is odd.")
                self.appLogger.warning("Sequence {0} from species {1} will not be paired.".format( seqDict[sp].keys(), sp))
                for item in seqDict[sp].keys():
                    del seqDict[sp][item]

            return None
