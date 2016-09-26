import argparse,datetime,logging,os,sys
import numpy as np
import random as rnd
import Settings as sp
from subprocess import call

"""
art_illumina -sam   -1 $profilePath/csNGSProfile_hiseq2500_1.txt
-2 $profilePath/csNGSProfile_hiseq2500_2.txt
-f 100 -l 150 -p  -m 250 -s 50 -rs $RANDOM -ss HS25 -i \$SEED -o \${INPUTBASE}_R

"""
class NGSReadsART:
    SHORT_NAMES=["sf" ,"dp" ,"1","2","amp","c","d","ef" ,"f","h","i",\
                "ir","ir2","dr","dr2","k","l","m","mp","M","nf","na",\
                "o","p","q","qU","qs","qL","qs2","rs","s","sam","sp","ss"]

    LONG_NAMES=["simphy_folder","data_prefix","qprof1","qprof2",\
                "amplicon","rcount ","id","errfree","fcov","help",\
                "in","insRate","insRate2","delRate","delRate2","maxIndel",\
                "len","mflen","matepair","cigarM","maskN","noALN","out",\
                "paired","quiet","maxQ","qShift","minQ","qShift2","rndSeed",\
                "sdev","samout","sepProf","seqSys"]

    def __init__(self,settings):
        self.appLogger=logging.getLogger('sngsw')
        self.appLogger.debug("(class NGSReadsART | __init__())")
        self.path=os.path.abspath(settings.parser.get("general", "simphy_folder"))
        if (settings.parser.get("general", "simphy_folder")[-1]=="/"):
            self.projectName=os.path.basename(settings.parser.get("general", "simphy_folder")[0:-1])
        else:
            self.projectName=os.path.basename(settings.parser.get("general", "simphy_folder"))

        self.prefix=settigs.parser.get("general","data_prefix")
        self.params=""
        dash=""
        for p in settings.parser.items("ngs-reads-art"):
            if (p[0] in self.SHORT_NAMES): dash="-"
            if (p[0] in self.LONG_NAMES): dash="--"
            if(p[1].lower() in ["true","false"]):
                self.params+=" {0}{1}".format(dash,p[0])
            else:
                self.params+=" {0}{1} {2}".format(dash,p[0],p[1])
        self.matingDict = csv.DictReader(\
            "{0}/{1}.mating".format(\
                self.path,\
                self.projectName\
            )\
        )

    def prepFolderStructure(self):
        for row in self.matingDict:
            # indexST,indexLOC,indID,speciesID,mateID1,mateID2
            folder="{0}/reads/{1}/{2}/".format(\
                self.path,\
                row['indexST'],\
                row['indexLOC']\
            )
            self.appLogger.info("Generating folder structure")

            try: os.makedirs(folder)
            except:  self.appLogger.debug("Folder ({0}) exists.".format(folder))

    def run(self):
        ngsMessageCorrect=""
        ngsMessageWrong=""
         # need to add the INPUT parameter
         # need to add the OUTPUT parameter
        command="art_illumina {0}".format(self.params)
        # know number of cores to send processes to
        numProcessors= open('/proc/cpuinfo').read().count('processor\t:')
        for row in self.matingDict:
            # indexST,indexLOC,indID,speciesID,mateID1,mateID2
            inputFile="{0}/individuals/{1}/{2}/{3}_{1}_{2}_{4}.fasta".format(\
                self.path,\
                row['indexST'],\
                row['indexLOC'],\
                self.projectName,\
                self.prefix,\
                row['indID']\
            )
            # By default I'm generating reads from the diploid individual.
            # This means, from a multiple (2) sequence fasta file.
            outputFile="{0}/reads/{1}/{2}/{3}_{1}_{2}_{4}_R".format(\
                self.path,\
                row['indexST'],\
                row['indexLOC'],\
                self.projectName,\
                self.prefix,\
                row['indID']\
            )

        # wait to finish call
        # subprocess.call("some-program")
        # parallel processing
        # p = subprocess.Popen("some-program")
        return True, message
        #call(self.params, stdout=open('{}/art_log'.format(self.output), 'w'))
