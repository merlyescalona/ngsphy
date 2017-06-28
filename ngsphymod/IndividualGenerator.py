 #!/usr/bin/home/python
import argparse,copy,datetime,logging,os,sys, sqlite3
import numpy as np
import random as rnd
import ngsphymod.Settings as sp
from ngsphymod.MSATools import *
from select import select

################################################################################
class IndividualGenerator:
    appLogger=None
    settings=None
    filteredSts=None
    # Output path
    output=""
    outputinds=""
    outputtables=""
    # Number of species trees replicates/folder to work with
    numSpeciesTrees=0
    numSpeciesTreesDigits=0
    # Number of fasta per replicate
    numLociPerSpeciesTree=[]
    numLociPerSpeciesTreeDigits=[]
    # Number of individuals per replicaate (ST)
    numIndividualsPerSpeciesTree=[]
    # Prefix of the datafiles that contain the FASTA sequences for
    # its corresponding gene tree.
    filteredSts=[]
    outgrou=False

    def __init__(self, settings):
        self.appLogger=logging.getLogger('ngsphy')
        self.appLogger.info("IndividualGenerator: Run started")
        self.settings=settings
        self.output=self.settings.outputFolderPath

    def checkArgs(self):
        self.appLogger.info("Checking arguments folder...")
        matingArgsMessageCorrect="Settings are correct, Mating process can be run."
        matingArgsMessageWrong="Something went wrong.\n"
        self.generateFolderStructure()
        # need to know whether i'm working with simphy or indelible
        if (self.settings.parser.has_option("general","numspeciestrees")):
            self.numSpeciesTrees=self.settings.parser.getint("general","numspeciestrees")
        self.filteredSts=range(1,self.numSpeciesTrees+1)
        self.numSpeciesTreesDigits=len(str(self.numSpeciesTrees))
        self.numIndividualsPerSpeciesTree=[0]*self.numSpeciesTrees
        self.numLociPerSpeciesTree=[0]*self.numSpeciesTrees
        self.numLociPerSpeciesTreeDigits=[0]*self.numSpeciesTrees
        if self.settings.originMode==1:
            # checking the replicate that are going to be used
            # checking if I'll used the filtered in case there's a possibility
            # that one or many sts do not match the ploidy and number of gene copies
            if self.settings.filterSimphy:
                self.filteredSts=self.filterSTMatchingIndPerSpeciesAndPloidy(self.settings.ploidy)
            self.command = os.path.join(\
                self.settings.path,\
                self.settings.projectName,\
                "{0}.command".format(self.settings.projectName))
            self.params = os.path.join(\
                self.settings.path,\
                self.settings.projectName,\
                "{0}.params".format(self.settings.projectName))
            self.db = os.path.join(\
                self.settings.path,\
                self.settings.projectName,\
                "{0}.db".format(self.settings.projectName))
            # check that the species tree replicate folder have the correct data
            gtperstOK,message=self.checkDataWithinReplicates()
            if (not gtperstOK):
                return gtperstOK,message
            self.numLociPerSpeciesTree=self.getSimPhyNumLociPerSpeciesTree()
            self.numLociPerSpeciesTreeDigits=[len(str(a))for a in self.numLociPerSpeciesTree]

            self.settings.parser.set(\
                "general",\
                "numLociPerSpeciesTree",\
                ",".join([str(a) for a in self.numLociPerSpeciesTree]))
        elif self.settings.originMode in [2,3] :
            if not self.settings.parser.has_option("general","numLociPerSpeciesTree"):
                return False,"Information about number of loci per folder is missing. Please verify. Exiting"
            self.numLociPerSpeciesTree=[self.settings.parser.getint("general","numLociPerSpeciesTree")]
            self.numLociPerSpeciesTreeDigits=[len(str(self.numLociPerSpeciesTree[0]))]
        # elif self.settings.originMode==3:
        else:
            return False, "{0}\n\t{1}\n\t{2}".format(\
                "Individual Generator process."
                "Something is wrong with the origin of the data.",\
                "Please verify. Exiting."\
            )

        self.settings.parser.set("general","filtered_ST",",".join([str(a) for a in self.filteredSts]))
        self.outgroup=self.outgroupExists()
        return True, matingArgsMessageCorrect

    def generateFolderStructure(self):
            # Checking output path
        self.appLogger.info("Creating folder structure for individual generation...")
        self.outputinds="{0}/individuals".format(self.output)
        self.outputtables="{0}/tables".format(self.output)
        try:
            self.appLogger.info("Generated individuals/")
            os.makedirs(self.outputinds)
        except:
            self.appLogger.debug("Output folder exists ({0})".format(self.output))
        try:
            self.appLogger.info("Generated tables/")
            os.makedirs(self.outputtables)
        except:
            self.appLogger.debug("Output folder exists ({0})".format(self.outputtables))


    def printConfiguration(self):
        self.appLogger.debug(\
            "\n\t{0}Configuration...{1}\n\tSimPhy project name:\t{2}\n\tSimPhy path:\t{3}\n\tOutput folder:\t{4}\n\tDataset prefix(es) (INDELible):\t{5}\n\tNumber of species trees replicates/folders:\t{6}".format(\
                mof.BOLD,mof.END,\
                self.settings.projectName,\
                self.settings.path,\
                self.output,\
                self.settings.dataPrefix,\
                self.numSpeciesTrees\
            )
        )

    def getSimPhyNumLociPerSpeciesTree(self):
        query="select N_Loci from Species_Trees"
        con = sqlite3.connect(self.db)
        res=con.execute(query).fetchall()
        con.close()
        res=[item for sublist in res for item in sublist]
        return res

    def filterSTMatchingIndPerSpeciesAndPloidy(self, ploidy):
        query="select SID from Species_Trees"
        if not (ploidy == 1):
            query="select SID from Species_Trees WHERE Ind_per_sp % {0} = 0".format(ploidy)
        con = sqlite3.connect(self.db)
        res=con.execute(query).fetchall()
        con.close()
        res=[item for sublist in res for item in sublist]
        return res

    def outgroupExists(self):
        if self.settings.originMode==1:
            for line in open(self.command,"r"):
                if "-so" in line or "-SO" in line:
                    self.settings.parser.set("general","outgroup","on")
                    return True
        if self.settings.originMode in [2,3]:
            value=self.settings.parser.getboolean("general","outgroup")
            return value
        self.settings.parser.set("general","outgroup","off")
        return False


    def checkDataWithinReplicates(self):
        for indexST in self.filteredSts:
            curReplicatePath="{0}/{1}/{2:0{3}d}/".format(\
                self.settings.path,\
                self.settings.projectName,\
                indexST,\
                self.numSpeciesTreesDigits)
            numFastaFiles=0;numGeneTrees=0
            fileList=os.listdir(curReplicatePath)
            # check composition of the current indexST folder
            for item in fileList:
                if ("{0}_".format(self.settings.dataPrefix) in item) and (".fasta" in item):
                    numFastaFiles+=1
                if  ("g_trees" in item) and (".trees" in item):
                    numGeneTrees+=1

            self.numLociPerSpeciesTree[indexST-1]=numFastaFiles
            self.appLogger.warning(\
                "Number of fasta files:\t{0}".format(numFastaFiles))
            self.numLociPerSpeciesTreeDigits[indexST-1]=len(str(numFastaFiles))
            if (numFastaFiles<1):
                # Do not have fasta files from the given replicate to work, I'll skip it.
                self.appLogger.warning("Replicate {0}({1}): It is not possible to do the mating for this replicate".format(indexST, curReplicatePath))
                self.appLogger.warning("There are no sequences o there is a missmatch between the prefixes and the number of sequences in the folder.")
                return False, "Please verify. Exiting."
            if (numGeneTrees<1):
                return False,"Trying to mate sequences, but there are no gene tree files to back that up. Please, finish the SimPhy run and try again afterwards."
        return True,"Got number of gene trees per species trees"

    """
    ############################################################################
    #                       Iterationg over STs                                #
    ############################################################################
    """
    def iteratingOverST(self):
        self.appLogger.debug("Ploidy: {0}".format(self.settings.ploidy))
        if (self.settings.ploidy==1):
            self.iterationHaploid()
        else:
            self.iterationPolyploid()
        self.settings.parser.set(\
            "general",\
            "numIndividualsPerSpeciesTree",\
            ",".join([str(a) for a in self.numIndividualsPerSpeciesTree]))

    def iterationPolyploid(self):
        for indexST in self.filteredSts:
            curReplicatePath="{0}/{1:0{2}d}/".format(self.outputinds,indexST, self.numSpeciesTreesDigits)
            self.appLogger.info("Replicate {0}/{2} ({1})".format(indexST, curReplicatePath,self.numSpeciesTrees))
            # iterating over the number of gts per st
            self.appLogger.info("Generating individuals...")
            matingTable=self.generateMatingTable(indexST)
            self.writeMatingTable(indexST,matingTable)

            for indexLOC in range(1,self.numLociPerSpeciesTree[indexST-1]+1):
                self.appLogger.debug("Number of FASTA file: {0}".format(indexLOC))
                # parsingMSA file
                fastapath="{0}/{1}/{2:0{3}d}/{4}_{5:0{6}d}.fasta".format(\
                    self.settings.path,\
                    self.settings.projectName,\
                    indexST,\
                    self.numSpeciesTreesDigits,\
                    self.settings.dataPrefix,\
                    indexLOC,\
                    self.numLociPerSpeciesTreeDigits[indexST-1])
                seqDict=parseMSAFile(fastapath)
                self.checkFileForIndels(seqDict)
                outputFolder="{0}/{1:0{2}d}/{3:0{4}d}".format(self.outputinds,\
                    indexST,self.numSpeciesTreesDigits,\
                    indexLOC,self.numLociPerSpeciesTreeDigits[indexST-1]\
                )
                try:
                    self.appLogger.debug("Output folder: {0}".format(outputFolder))
                    os.makedirs(outputFolder)
                except OSError as err:
                    print("OS error: {0}".format(err))
                    self.appLogger.debug("Folder {0} exists.".format(outputFolder))
                # generating and writing mating table
                self.mate(indexST,indexLOC,matingTable,seqDict)



    def iterationHaploid(self):
        for indexST in self.filteredSts:
            curReplicatePath="{0}/{1:0{2}d}/".format(self.outputinds,indexST, self.numSpeciesTreesDigits)
            self.appLogger.info("Replicate {0}/{2} ({1})".format(indexST, curReplicatePath,self.numSpeciesTrees))
            # iterating over the number of gts per st
            self.appLogger.info("Generating individuals...")
            individualTable=self.generateindividualTable(indexST)
            self.writeIndividualTable(indexST,individualTable)
            for indexLOC in range(1,self.numLociPerSpeciesTree[indexST-1]+1):
                self.appLogger.debug("Number of FASTA file: {0}".format(indexLOC))
                # parsingMSA file
                self.appLogger.debug("Using ST={0}, LOC={1}".format(indexST,indexLOC))
                fastapath="{0}/{1}/{2:0{3}d}/{4}_{5:0{6}d}.fasta".format(\
                    self.settings.path,\
                    self.settings.projectName,\
                    indexST,\
                    self.numSpeciesTreesDigits,\
                    self.settings.dataPrefix,\
                    indexLOC,\
                    self.numLociPerSpeciesTreeDigits[indexST-1])
                seqDict=parseMSAFileWithDescriptions(fastapath)
                self.checkFileForIndels(seqDict)
                outputFolder="{0}/{1:0{2}d}/{3:0{4}d}".format(self.outputinds,\
                    indexST,self.numSpeciesTreesDigits,\
                    indexLOC,self.numLociPerSpeciesTreeDigits[indexST-1]\
                )
                try:
                    self.appLogger.debug("Output folder: {0}".format(outputFolder))
                    os.makedirs(outputFolder)
                except OSError as err:
                    print("OS error: {0}".format(err))
                    self.appLogger.debug("Folder {0} exists.".format(outputFolder))
                # generating and writing mating table
                self.generateIndividuals(indexST,indexLOC,individualTable,seqDict)

    def generateindividualTable(self,indexST):
        # get first gt of the st to get the descriptions
        indexLOC=1
        self.appLogger.debug("Using ST={0}, LOC={1}".format(indexST,indexLOC))
        fastapath="{0}/{1}/{2:0{3}d}/{4}_{5:0{6}d}_TRUE.fasta".format(\
            self.settings.path,\
            self.settings.projectName,\
            indexST,\
            self.numSpeciesTreesDigits,\
            self.settings.dataPrefix,\
            indexLOC,\
            self.numLociPerSpeciesTreeDigits[indexST-1])

        descriptions=parseMSAFileWithDescriptions(fastapath).keys()
        descriptions.sort()
        table=[(item, descriptions[item]) for item in range(0,len(descriptions)) ]
        self.numIndividualsPerSpeciesTree[indexST-1]=len(table)
        return table

    def generateIndividuals(self,indexST,indexLOC,individualTable,seqDict):
        self.appLogger.debug("{0}/{1} - {2}".format(indexLOC,self.numLociPerSpeciesTree[indexST-1],self.numLociPerSpeciesTreeDigits[indexST-1]))
        outputFolder="{0}/{1:0{2}d}/{3:0{4}d}".format(self.outputinds,\
            indexST,self.numSpeciesTreesDigits,
            indexLOC,self.numLociPerSpeciesTreeDigits[indexST-1]\
        )
        for currentInd in range(0,len(individualTable)):
            # Extracting info from the dictionary
            indID=str(individualTable[currentInd][0])
            description=str(individualTable[currentInd][1])
            # Organizing strings
            # if description=="6_0_0":
                # print indID,outputFolder
            seq=seqDict[description]
            des=">{0}:{1:0{6}d}:{2:0{7}d}:{3}:{4}:{5}".format(self.settings.projectName,\
                indexST,indexLOC,self.settings.dataPrefix,indID,description[1:len(description)],\
                self.numSpeciesTreesDigits,self.numLociPerSpeciesTreeDigits[indexST-1]\
            )

            indFilename="{0}/{1}_{2:0{3}d}_{4:0{5}d}_{6}_{7}.fasta".format(\
                outputFolder,\
                self.settings.projectName,\
                indexST,\
                self.numSpeciesTreesDigits,\
                indexLOC,\
                self.numLociPerSpeciesTreeDigits[indexST-1],\
                self.settings.dataPrefix,\
                indID\
            )
            indFile=open(indFilename, "w")
            indFile.write("{0}\n{1}\n".format(des,seq))
            indFile.close()

    def writeIndividualTable(self,indexST,individualTable):
        self.appLogger.debug("Writing table")
        # mating table
        # indexST,SP,ind-tip1,ind-tip2
        self.appLogger.debug("Writing indexes into file...")
        indexFilename="{0}/{1}.{2:0{3}d}.individuals.csv".format(\
            self.outputtables,\
            self.settings.projectName,\
            indexST,\
            self.numSpeciesTreesDigits)
        if not os.path.isfile(indexFilename):
            indexFile=open(indexFilename,"w")
            indexFile.write("indexST,indID,speciesID,locusID,geneID\n")
            indexFile.close()
        indexFile=open(indexFilename,"a")
        for indexRow in range(0,len(individualTable)):
            indID=individualTable[indexRow][0]
            seqDescription=individualTable[indexRow][1]
            speciesID=seqDescription.strip().split("_")[0]
            locusID=seqDescription.strip().split("_")[1]
            geneID=seqDescription.strip().split("_")[2]
            indexFile.write("{0:0{1}d},{2},{3},{4},{5}\n".format(\
                indexST,\
                self.numSpeciesTreesDigits,\
                indID,\
                speciesID,\
                locusID,\
                geneID))
        indexFile.close()


    # def generateMatingTableFromSequenceDescription(self,indexST):
    def generateMatingTable(self,indexST):
        self.appLogger.debug("Mating table")
        filename=os.path.join(\
            self.settings.path,self.settings.projectName,\
            "{0:0{1}d}".format(\
                indexST,\
                self.numSpeciesTreesDigits),\
            "{0}_{1:0{2}d}.fasta".format(\
                self.settings.dataPrefix,\
                1,\
                self.numLociPerSpeciesTreeDigits[indexST-1]))
        print filename
        f=open(filename,"r")
        lines=f.readlines()
        f.close()
        leaves=[ item.strip()[1:] for item in lines if item.startswith(">")]
        leavesSplit=[ item.split("_") for item in leaves]
        leavesDict=dict()
        for tip in leavesSplit:
            geneFamily="_".join(tip[0:2])
            try:
                val=leavesDict[geneFamily]
                leavesDict[geneFamily]+=1
            except:
                leavesDict[geneFamily]=1
        # Till here i have information about the number of tips per gene family
        nInds=0;numLeaves=0
        for item in leavesDict:
            numLeaves+=leavesDict[item]
        if self.settings.ploidy==2:
            nInds=numLeaves/2
        mates=[]
        if self.outgroup:
            mates+=[(indexST,0,0,0,0)]
            nInds=(numLeaves-1)/2
        for geneFamily in leavesDict:
            t=range(0,leavesDict[geneFamily])
            while not t==[]:
                p1=0;p2=0
                try:
                    p1=t.pop(rnd.sample(range(0,len(t)),1)[0])
                    p2=t.pop(rnd.sample(range(0,len(t)),1)[0])
                except Exception as e:
                    break
                sp=geneFamily.split("_")[0]
                lt=geneFamily.split("_")[1]
                pair=(indexST,sp,lt,p1,p2)
                mates+=[pair]
                self.appLogger.debug("Pair generated: {0}".format(pair))
        self.numIndividualsPerSpeciesTree[indexST-1]=len(mates)
        return mates



    def generateMatingTableFromDB(self,indexST):
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
        self.numIndividualsPerSpeciesTree[indexST-1]=len(mates)
        return mates

    def writeMatingTable(self,indexST,matingTable):
        # mating table
        # indexST,SP,ind-tip1,ind-tip2
        # If i have outgroups mate info should be in the table already
        self.appLogger.debug("Writing indexes into file...")
        indexFilename="{0}/{1}.{2:0{3}d}.individuals.csv".format(self.outputtables,self.settings.projectName, indexST, self.numSpeciesTreesDigits)
        self.appLogger.debug(indexFilename)
        if not os.path.isfile(indexFilename):
            indexFile=open(indexFilename,"w")
            indexFile.write("indexST,indID,speciesID,locusID,mateID1,mateID2\n")
            indexFile.close()
        indexFile=open(indexFilename,"a")
        for indexRow in range(0,len(matingTable)):
            indID=indexRow
            speciesID=matingTable[indexRow][1]
            locusID=matingTable[indexRow][2]
            mateID1=matingTable[indexRow][3]
            mateID2=matingTable[indexRow][4]
            indexFile.write("{0:0{1}d},{2},{3},{4},{5},{6}\n".format(\
                indexST,\
                self.numSpeciesTreesDigits,\
                indID,\
                speciesID,\
                locusID,\
                mateID1,\
                mateID2))
        indexFile.close()

    def mate(self,indexST,indexLOC,matingTable,sequenceDictionary):
        # Proper mating
        seqDict=copy.deepcopy(sequenceDictionary)
        species=seqDict.keys()
        numSeqs=0
        for key in species:
            subspecies=seqDict[key].keys()
            for subkey in subspecies:
                numSeqs+=len(seqDict[key][subkey])
        numInds=np.trunc(numSeqs/2)+1
        self.appLogger.debug("Writing {1} individuals from {0} number of sequences.".format(numSeqs,numInds))
        self.appLogger.debug("{0}/{1} - {2}".format(indexLOC,self.numLociPerSpeciesTree[indexST-1],self.numLociPerSpeciesTreeDigits[indexST-1]))
        outputFolder="{0}/{1:0{2}d}/{3:0{4}d}".format(self.outputinds,\
            indexST,self.numSpeciesTreesDigits,\
            indexLOC,self.numLociPerSpeciesTreeDigits[indexST-1]\
        )

        seq1="";des1="";seq2="";des2="";
        # for all species except the outgroup
        for currentInd in range(0,len(matingTable)):
            # Extracting info from the dictionary
            st=str(matingTable[currentInd][0])
            sp=str(matingTable[currentInd][1])
            lt=str(matingTable[currentInd][2])
            pair1=str(matingTable[currentInd][3])
            pair2=str(matingTable[currentInd][4])
            tag="{0}_{1}".format(sp,lt)
            # Organizing strings
            self.appLogger.debug("{0}|{1}-{2}".format(sp, pair1,pair2))
            seq1=seqDict[tag][pair1]["sequence"]
            seq2=seqDict[tag][pair2]["sequence"]
            shortDesc1=seqDict[tag][pair1]["description"][1:len(seqDict[tag][pair1]["description"])]
            des1=">{0}:{1:0{6}d}:{2:0{7}d}:{3}:{4}:{5}_S1".format(self.settings.projectName,\
                indexST,indexLOC,self.settings.dataPrefix,currentInd,shortDesc1,\
                self.numSpeciesTreesDigits,self.numLociPerSpeciesTreeDigits[indexST-1]\
            )
            shortDesc2=seqDict[tag][pair2]["description"][1:len(seqDict[tag][pair2]["description"])]
            des2=">{0}:{1:0{6}d}:{2:0{7}d}:{3}:{4}:{5}_S2".format(\
                self.settings.projectName,indexST,indexLOC,self.settings.dataPrefix,\
                currentInd,shortDesc2,
                self.numLociPerSpeciesTree[indexST-1],self.numLociPerSpeciesTreeDigits[indexST-1]\
            )

            indFilename="{0}/{1}_{2:0{3}d}_{4:0{5}d}_{6}_{7}.fasta".format(\
                outputFolder,\
                self.settings.projectName,\
                indexST,\
                self.numSpeciesTreesDigits,\
                indexLOC,\
                self.numLociPerSpeciesTreeDigits[indexST-1],\
                self.settings.dataPrefix,\
                currentInd\
            )
            indFile=open(indFilename, "w")
            indFile.write("{0}\n{1}\n{2}\n{3}\n".format(des1,seq1,des2,seq2))
            indFile.close()

            del seqDict[tag][pair1]
            if not sp== "0": del seqDict[tag][pair2]

        if not seqDict[tag]=={}:
            self.appLogger.warning("Number of individuals (sequences) generated per species is odd.")
            self.appLogger.warning("Sequence {0} from species {1} will not be paired.".format( seqDict[tag].keys(), sp))
            for item in seqDict[tag].keys():
                del seqDict[sp][item]

        return None

    def getReferencesSequences(self, indexSP,indexTip):
        #print("getReferencesSequences: {0}.{1}".format(indexSP,indexTip))
        for indexST in self.filteredSts:
            self.appLogger.debug("indexST: {0}".format(indexST))
            for indexLOC in range(1,self.numLociPerSpeciesTree[indexST-1]+1):
                #print("indexLOC: {0}".format(indexLOC))
                fastapath="{0}/{1}/{2:0{3}d}/{4}_{5:0{6}d}.fasta".format(\
                    self.settings.path,\
                    self.settings.projectName,\
                    indexST,\
                    self.numSpeciesTreesDigits,\
                    self.settings.dataPrefix,\
                    indexLOC,\
                    self.numLociPerSpeciesTreeDigits[indexST-1])
                msafileDict=parseMSAFile(fastapath)
                self.checkFileForIndels(msafileDict)
                #print(msafileDict[str(indexSP)][str(indexTip)]['description'])
                #print(msafileDict[str(indexSP)][str(indexTip)]['sequence'])
                f=open("{0}/references/REFS_ST_{1:0{2}d}.fasta".format(self.output,indexST,self.numSpeciesTreesDigits),'a')
                f.write(">ST_{0:0{1}d}_LOC_{2:0{3}d}\n{4}\n".format(\
                    indexST,\
                    self.numSpeciesTreesDigits,\
                    indexLOC,\
                    self.numLociPerSpeciesTreeDigits[indexST-1],\
                    msafileDict[str(indexSP)][str(indexTip)]['sequence']))


    def checkFileForIndels(self,msa):
        keys=msa.keys()
        item=""; status=True;message=""
        if self.settings.ploidy==1:
            for item in keys:
                if "-" in msa[item]:
                    break
            if "-" in msa[item]:
                status=False
        if self.settings.ploidy==2:
            for item in keys:
                if "-" in msa[item]["sequence"]:
                    break
            if "-" in msa[item]["sequence"]:
                status=False
        # only changeing the value if dataset has no indels and the value I'm about to change
        # is different from the value it has
        if not self.settings.indels  and (self.settings.indels==not status):
            self.settings.indels=not status
