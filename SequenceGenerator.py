#!/usr/bin/home/python
import argparse,copy,datetime,ete2,logging,os,re,sys, threading, subprocess
import numpy as np
import random as rnd
import Settings as sp
from MSATools import *
from select import select

class SequenceGenerator:
    appLoger=None
    settings=None

    ngsphyindeliblecontrol=""
    indeliblecontrol=""
    newickFile=""
    projectName=""
    path=""
    output=""

    evolve=[]
    partition=[]

    numGTs=0
    numGTDigits=0

    def __init__(self,settings):
        self.appLogger=logging.getLogger('ngsphy')
        self.appLogger.debug('INDELible run')
        self.settings=settings
        self.ngsphyindeliblecontrol=self.settings.ngsphyindeliblecontrol
        self.newickFile=self.settings.newickFile
        self.projectName=self.settings.projectName
        self.path=self.settings.path
        self.output=os.path.join(self.path,self.projectName,"1")
        self.indeliblecontrol=os.path.join(self.path,self.projectName,"1","control.txt")

    def run(self):
        self.generateFolderStructure()
        # check naming of the leaves
        namingStatus, namingMessage=self.checkTreeLabels()
        # adding extra information on presence of outgroup
        self.addOutgroupInfoToSettings()

        if (namingStatus):
            # check ploidy and tree correspondance
            ploidyStatus,ploidyMessage=self.checkPloidyTreeRelation()
            if (ploidyStatus):
                self.appLogger.debug("Out tree-ploidy relation")
                self.writeIndelibleControlFile()
                runStatus,runMessage=self.runIndelible()
                if not runStatus:
                    return False,runMessage
            else:
                return False, ploidyMessage
        else:
            return False, namingMessage
        return True, "Run finished"

    def generateFolderStructure(self):
        self.appLogger.info("Creating folder structure for INDELible run")
        # create pporject folder
        try:
            os.makedirs(os.path.join(self.path,self.projectName))
            self.appLogger.info("Generating project folder ({0})".format(os.path.join(self.path,self.projectName)))
        except:
            self.appLogger.debug("Project folder exists ({0})".format(os.path.join(self.path,self.projectName)))
        # create data folder
        try:
            os.makedirs(os.path.join(self.path,self.projectName,"1"))
            self.appLogger.info("Generating data folder ({0})".format(os.path.join(self.path,self.projectName,"1")))
        except:
            self.appLogger.debug("Data folder exists ({0})".format(os.path.join(self.path,self.projectName,"1")))

    def checkTreeLabels(self):
        self.appLogger.debug("Checking labels")
        messageCorrect="Labels of the tree are correct"
        messageWrong="INDELible control file - Something's wrong!\n\t"
        try:
            tree=ete2.Tree(self.newickFile)
        except ete2.parser.newick.NewickError as nerror:
            return False, "There is a problem with the newick file\n{}\nPlease verify. Exiting.".format(nerror)

        leaves=[ item.get_leaf_names()[0] for item in tree.get_leaves()]
        pattern = re.compile("^[0-9]+_[0-9]+_[0-9]+$")
        for item in leaves:
            if not pattern.match(item):
                break
        if not pattern.match(item):
            messageWrong+="{0}\n\t{1}\n\t{2}\n\t{3}\n\t{4}".format(\
            "Labels chosen for the tips of the tree are not correct.",\
            "Labels should follow this pattern: SpeciesID_LocusID_IndividualID",\
            "Where SpeciesID,LocusID,IndividualID are numbers and bigger than 0.",\
            "Outgroup, if present, should be represented as 0_0_0",\
            "Please verify. Exiting."\
            )
            return False,messageWrong
        return True, messageCorrect

    def addOutgroupInfoToSettings(self):
        self.appLogger.debug("Outgroup")
        try:
            tree=ete2.Tree(self.newickFile)
        except ete2.parser.newick.NewickError as nerror:
            return False, "There is a problem with the newick file\n{}\nPlease verify. Exiting.".format(nerror)

        leaves=[ item.get_leaf_names()[0] for item in tree.get_leaves()]

        for item in leaves:
            if item == "0_0_0":
                break
        if item=="0_0_0":
            self.settings.parser.set("general","outgroup","on")
        else:
            self.settings.parser.set("general","outgroup","off")

    def checkPloidyTreeRelation(self):
        self.appLogger.debug("Checking ploidy - num tips relation")
        messageCorrect="Ploidy and number of gene copies per gene family match properly."
        messageWrong="INDELible control file - Something's wrong!\n\t"
        tree=ete2.Tree(self.newickFile)
        leaves=[ item.get_leaf_names()[0] for item in tree.get_leaves()]
        leavesSplit=[ item.split("_") for item in leaves]
        leavesDict=dict()
        for tip in leavesSplit:
            geneFamily="_".join(tip[0:2])
            try:
                val=leavesDict[geneFamily]
                leavesDict[geneFamily]+=1
            except:
                leavesDict[geneFamily]=1
        for item in leavesDict:
            if (leavesDict[item] % self.settings.ploidy) != 0:
                messageWrong+="{0}\n\t{1}".format(\
                "The number of gene copies within one of the gene families does not match the ploidy selected for this run.",\
                "Please verify. Exiting."\
                )
                return False,messageWrong
        return True, messageCorrect

    def writeIndelibleControlFile(self):
        self.appLogger.debug("Writing new control file")
        f=open(self.ngsphyindeliblecontrol,"r")
        lines=f.readlines()
        f.close()
        newlines=copy.copy(lines)
        newlines.reverse()
        controllines=[]
        modelname=""
        while len(newlines)>0:
            line=newlines.pop()
            if "[NGSPHYEVOLVE]" in line:
                self.evolve=line.split() # ill get 2 elems + label
            if "[NGSPHYPARTITION]" in line:
                self.partition=line.split() # ill get 3 elems + label
            if "[MODEL]" in line:
                modelname=line.split()[1]

        newlines=copy.copy(lines)
        # for item in newlines:
            # print item
        for item in newlines:
            if item.strip().startswith("[NGSPHY"):
                break
            controllines+=[item.strip("\n")]

        f=open(self.newickFile)
        newicklines=f.readlines()
        f.close()
        newicktree=[ item.strip() for item in newicklines if item.strip()!=""]
        newicktree="".join(newicktree)
        if newicktree[-1]!=";":
            newicktree+=";"
        controllines+=["{0} {1} {2}".format(\
            "[TREE]",\
            "ngsphytree",\
            newicktree
        )]
        controllines+=["{0} {1} [{2} {3} {4}]".format(\
            "[PARTITIONS]",\
            "ngsphypartition",\
            "ngsphytree",\
            modelname,\
            self.partition[3]
        )]

        numGTs=int(self.evolve[1])
        numGTDigits=len(str(numGTs))
        controllines+=["[EVOLVE]"]
        for indexGT in range(1, numGTs+1):
            controllines+=["  {0} 1 {1}".format(\
                "ngsphypartition",\
                "{0}_{1:0{2}d}".format(self.evolve[2],indexGT, numGTDigits)
            )]

        # full control file, missing checking settings of output and fastaextension
        fastaoutput="\t[output] FASTA"
        fastaoutputext="\t[fastaextension] fasta"
        output=[]; outputext=[]
        settings=False
        for indexControl in range(0, len(controllines)):
            data=controllines[indexControl].strip()
            if data=="[SETTINGS]":
                settings=True
            if data.startswith("[output]"):
                ss=data.split()[1]
                output+=[ss.upper()]
            if data.startswith("[fastaextension]"):
                outputext+=[indexControl]
        if not settings:
            controllines.insert(1,"[SETTINGS]")
        if (not "FASTA" in output) and (len(outputext) ==0):
            controllines.insert(2,"  [output] FASTA")
            controllines.insert(3,"  [fastaextension] fasta")
        elif (not "FASTA" in output):
            controllines.insert(2,"  [output] FASTA")
        elif (len(outputext) ==0):
            controllines.insert(2,"  [fastaextension] fasta")

        # write controllines to file
        f=open(self.indeliblecontrol,"w")
        for item in controllines:
            f.write("{}\n".format(item))
        f.close()

    def runIndelible(self):
        self.appLogger.debug("Running...")
        self.settings.parser.set("general","numLociPerST",self.evolve[1])
        self.settings.parser.set("general", "filtered_ST", "1")
        self.settings.parser.set("general", "number_ST", "1")
        self.settings.parser.set("general", "data_prefix",self.evolve[2])
        try:
            self.appLogger.info("Waiting for INDELible process to finish. This may take a while...")
            t = threading.Thread(target=self.indelibleLauncher())
            t.start()
            t.join()
        except RuntimeError as verror:
            return   False, verror
        return True, "INDELible's run has finished."

    def indelibleLauncher(self):
        indelibleMessage="INDELible run has finished";proc=""
        lines=[]
        try:
            self.appLogger.info("Moving to {0}".format(self.output))
            os.chdir(self.output)
            self.appLogger.info("Running INDELible")
            proc = subprocess.check_output("indelible",stderr=subprocess.STDOUT)
            cpuTime = [line.split(":")[1].split()[0] for line in proc.split('\n') if "* Block" in line]
            numGTDigits=len(str(len(cpuTime)))
            for item in range(1,len(cpuTime)):
                indexGT=item
                cpu=cpuTime[(item-1)]
                output="{0}_{1:0{2}d}".format(self.evolve[2],item,numGTDigits )
                lines+=[indexGT,cpu,output]
        except OSError as error:
            indelibleMessage="Output folder does not exist. Please verify. Exiting."
            self.appLogger.error(indelibleMessage)
            raise RuntimeError(indelibleMessage)
        except subprocess.CalledProcessError as error:
            indelibleMessage="\nINDELible execution error. "+\
            "\n------------------------------------------------------------------------"+\
            "\n{0}".format(error)+\
            "\n------------------------------------------------------------------------"+\
            "{0}".format(error.output)+\
            "\n------------------------------------------------------------------------"+\
            "\nFor more information about this error please run the following commands separately:\n"+\
            "\n\tcd {0}\n\tindelible\n".format(self.output)
            raise RuntimeError(indelibleMessage)
        self.printRunningInfo(line)

    def printRunningInfo(self, lines):
        outputFile="{0}/{1}.indelible.info".format(
            self.output,\
            self.projectName
        )
        f=open(outputFile,"w")
        f.write("indexGT,cpuTime,outputFilePrefix\n")
        for item in lines:
            f.write(
                str(item[0])+","+\
                str(item[1])+","+\
                item[2]+"\n"
            )
        f.close()
        self.appLogger.info("File with timings of the INDELible run can be find on: {0}".format(outputFile))

    def cleanUpControlFile(self):
        f=open(self.ngsphyindeliblecontrol,"r")
        lines=f.readLines()
        f.close()
        newlines=[ item.strip() for item in lines if item.strip()!=""]
        for item in newlines:
            if item.startswith("// ") or item.startswith("//"):
                newlines.remove(item)
        for index in range(0,len(newlines)):
            item=newlines[index]
            if item.find("//")>-1:
                item=item[0:item.find("//")]
                newlines[index]=item
        for item in newlines:
            if item=="": newlines.remove(item)
        startBracket=[];endBracket=[]
        for index in range(0,len(newlines)):
            item=newlines[index]
            if item.find("/*")>-1:
                startBracket+=[index]
        for index in range(0,len(newlines)):
            item=newlines[index]
            if item.find("*/")>-1:
                endBracket+=[(index+1)]
        linestoremove=[]
        for index in range(0,len(startBracket)):
            linestoremove+=newlines[startBracket[index]:endBracket[index]]
        for item in linestoremove:
            newlines.remove(item)
        # until here my control file is out of comments
        return newlines
