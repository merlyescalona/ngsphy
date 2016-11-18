import argparse,csv,datetime,logging,os,subprocess,sys
import numpy as np
import random as rnd
import Settings as sp

class NGSReadsARTIllumina:
    SHORT_NAMES=["sf" ,"dp" ,"1","2","amp","c","d","ef" ,"f","h","i",\
                "ir","ir2","dr","dr2","k","l","m","mp","M","nf","na",\
                "o","p","q","qU","qs","qL","qs2","rs","s","sam","sp","ss"]
    LONG_NAMES=["simphy_folder","data_prefix","qprof1","qprof2",\
                "amplicon","rcount","id","errfree","fcov","help",\
                "in","insRate","insRate2","delRate","delRate2","maxIndel",\
                "len","mflen","matepair","cigarM","maskN","noALN","out",\
                "paired","quiet","maxQ","qShift","minQ","qShift2","rndSeed",\
                "sdev","samout","sepProf","seqSys"]
    dLONG_NAMES={i.lower():i for i in LONG_NAMES}
    dSHORT_NAMES={i.lower():i for i in SHORT_NAMES}

    def __init__(self,settings):
        self.appLogger=logging.getLogger('sngsw')
        self.appLogger.info('NGS read simulation: ART run started.')
        simphy=os.path.abspath(settings.parser.get("general", "simphy_folder"))
        self.output=os.path.abspath(\
            os.path.join(\
                simphy,\
                settings.parser.get("general", "output_folder_name")\
            )\
        )
        self.filteredST=settings.parser.get("general", "filtered_ST")
        sys.exit()
        if (simphy[-1]=="/"):
            self.projectName=os.path.basename(simphy)[0:-1]
        else:
            self.projectName=os.path.basename(simphy)

        self.prefix=settings.parser.get("general","data_prefix")
        self.params=[];
        dash=""; par=[]
        settingsParams=settings.parser.items("ngs-reads-art")
        for p in settingsParams:
            if (p[0] in self.dSHORT_NAMES.keys()): dash="-"
            if (p[0] in self.dLONG_NAMES.keys()): dash="--"
            # to be sure that i am getting the right parameter names
            if (dash=="-"): par=[self.dSHORT_NAMES[p[0]]]
            if (dash=="--"):    par=[self.dLONG_NAMES[p[0]]];
            par+=[p[1]]
            if(par[1].lower() in ["true","false"]):
                self.params+=["{0}{1}".format(dash,par[0])]
            else:
                self.params+=["{0}{1}".format(dash,par[0]),par[1]]

        self.numFiles=0
        try:
            os.mkdirs("{0}/reads".format(self.output))
            self.appLogger.info("Generating output folder ({0}/reads)".fomat(self.output))
        except:
            self.appLogger.debug("Output folder exists ({0}/reads)".fomat(self.output))

        try:
            os.mkdirs("{0}/scripts".format(self.output))
            self.appLogger.info("Generating output folder ({0}/scripts)".fomat(self.output))
        except:
            self.appLogger.debug("Output folder exists ({0}/scripts)".fomat(self.output))

    def writeSeedFile(self):
        seedfile=open("{0}/scripts/{1}.seedfile.txt".format(\
            self.output,\
            self.projectName
        ))
        for indexST in self.filteredST:
            csvfile=open("{0}/{1}.{2:0{3}d}.mating.csv".format(\
                self.output,\
                self.projectName,\
                indexST,\
                numberSTDigits
            ))
            # Generation of folder structure
            d = csv.DictReader(csvfile)
            self.matingDict = [row for row in d]
            csvfile.close()
            for row in self.matingDict:
                # indexST,indexLOC,indID,speciesID,mateID1,mateID2
                folder="{0}/reads/{1}/{2}/".format(\
                    self.output,\
                    row['indexST'],\
                    row['indexLOC']\
                )

            for row in self.matingDict:
                # indexST,indexLOC,indID,speciesID,mateID1,mateID2
                inputFile="{0}/individuals/{1}/{2}/{3}_{1}_{2}_{4}_{5}.fasta".format(\
                    self.output,\
                    row['indexST'],\
                    row['indexLOC'],\
                    self.projectName,\
                    self.prefix,\
                    row['indID']\
                )
                # This means, from a multiple (2) sequence fasta file.
                outputFile="{0}/reads/{1}/{2}/{3}_{1}_{2}_{4}_{5}_R".format(\
                    self.output,\
                    row['indexST'],\
                    row['indexLOC'],\
                    self.projectName,\
                    self.prefix,\
                    row['indID']\
                )
                seedfile.write("{0}\t{1}\n".format(inputFile,outputFile))
                self.numFiles+=1
        seedfile.close()
        self.appLogger.info("Seed file written...")

    def writeSGEScript(self):
        jobfile=open("{0}/scripts/{1}.job.sge.sh".format(\
            self.output,\
            self.projectName
        ))
        seedfile="{0}/scripts/{1}.seedfile.txt".format(\
            self.output,\
            self.projectName
        )

        inputFile="$inputfile"
        outputFile="$outputfile"
        callParams=["art_illumina"]+self.params+["--in", inputFile,"--out",outputFile]
        header="""#!/bin/bash
# SGE submission options
#$ -l num_proc=1         # number of processors to use
#$ -l h_rt=00:10:00      # Set 10 mins  - Average amount of time for up to 1000bp
#$ -t 1-{0}              # Number of jobs/files that will be treated
#$ -N art.sims           # A name for the job

inputfile=$(awk 'NR==$SGE_TASK_ID{print $1}' {1})
outputfile=$(awk 'NR==$SGE_TASK_ID{print $2}' {1})
        """.format(self.numFiles,seedfile)
        footer="".format()
        jobfile.write(header)
        jobfile.write(" ".join(callParams))
        jobfile.write(footer)
        jobfile.close()
        self.appLogger.info("SGE Job file written ({0})...".format(jobfile))

    def writeSLURMScript(self):
        jobfile=open("{0}/scripts/{1}.job.slurm.sh".format(\
            self.output,\
            self.projectName
        ))
        seedfile="{0}/scripts/{1}.seedfile.txt".format(\
            self.output,\
            self.projectName
        )
        inputFile="$inputfile"
        outputFile="$outputfile"
        callParams=["art_illumina"]+self.params+["--in", inputFile,"--out",outputFile]

        header="""#!/bin/sh
#SBATCH -n 1
#SBATCH --cpus-per-task 1
#SBATCH -t 00:10:00
#SBATCH --mem 4G
#SBATCH --array=1-{0}

inputfile=$(awk 'NR==$SLURM_ARRAY_TASK_ID{print $1}' {1})
outputfile=$(awk 'NR==$SLURM_ARRAY_TASK_ID{print $2}' {1})

""".format(self.numFiles, seedfile)
        footer="".format()
        jobfile.write(header)
        jobfile.write(" ".join(callParams))
        jobfile.write(footer)
        jobfile.close()
        self.appLogger.info("SLURM Job file written ({0})...".format(jobfile))


    def writeBashScript(self):
        bashfile=open("{0}/scripts/{1}.sh".format(\
            self.output,\
            self.projectName
        ))
        for indexST in self.filteredST:
            csvfile=open("{0}/{1}.{2:0{3}d}.mating.csv".format(\
                self.output,\
                self.projectName,\
                indexST,\
                numberSTDigits
            ))
            # Generation of folder structure
            d = csv.DictReader(csvfile)
            self.matingDict = [row for row in d]
            csvfile.close()
            for row in self.matingDict:
                # indexST,indexLOC,indID,speciesID,mateID1,mateID2
                folder="{0}/reads/{1}/{2}/".format(\
                    self.output,\
                    row['indexST'],\
                    row['indexLOC']\
                )

            for row in self.matingDict:
                # indexST,indexLOC,indID,speciesID,mateID1,mateID2
                inputFile="{0}/individuals/{1}/{2}/{3}_{1}_{2}_{4}_{5}.fasta".format(\
                    self.output,\
                    row['indexST'],\
                    row['indexLOC'],\
                    self.projectName,\
                    self.prefix,\
                    row['indID']\
                )
                # This means, from a multiple (2) sequence fasta file.
                outputFile="{0}/reads/{1}/{2}/{3}_{1}_{2}_{4}_{5}_R".format(\
                    self.output,\
                    row['indexST'],\
                    row['indexLOC'],\
                    self.projectName,\
                    self.prefix,\
                    row['indID']\
                )
                # Call to ART
                callParams=["art_illumina"]+self.params+["--in", inputFile,"--out",outputFile]
                # self.params+=["--in ",inputFile,"--out",outputFile]
                bashfile.write(" ".join(callParams))
                bashfile.write("\n")
        bashfile.close()
        self.appLogger.info("Bash script written...")

    def run(self):
        ngsMessageCorrect="ART Finished succesfully"
        ngsMessageWrong="Ops! Something went wrong.\n\t"

        run=self.parser.getboolean("execution","run")
        if (run):
            # I have to iterate over the sts, now that i have more than on
            numberSTDigits=len(str(np.max(selffilteredST)))
            for indexST in self.filteredST:
                csvfile=open("{0}/{1}.{2:0{3}d}.mating.csv".format(\
                    self.output,\
                    self.projectName,\
                    indexST,\
                    numberSTDigits
                ))
                # Generation of folder structure
                d = csv.DictReader(csvfile)
                self.matingDict = [row for row in d]
                csvfile.close()
                self.appLogger.info("Generating folder structure")
                for row in self.matingDict:
                    # indexST,indexLOC,indID,speciesID,mateID1,mateID2
                    folder="{0}/reads/{1}/{2}/".format(\
                        self.output,\
                        row['indexST'],\
                        row['indexLOC']\
                    )
                    try:
                        os.makedirs(folder)
                    except:
                        self.appLogger.debug("Folder ({0}) exists.".format(folder))

                for row in self.matingDict:
                    # indexST,indexLOC,indID,speciesID,mateID1,mateID2
                    inputFile="{0}/individuals/{1}/{2}/{3}_{1}_{2}_{4}_{5}.fasta".format(\
                        self.output,\
                        row['indexST'],\
                        row['indexLOC'],\
                        self.projectName,\
                        self.prefix,\
                        row['indID']\
                    )
                    # This means, from a multiple (2) sequence fasta file.
                    outputFile="{0}/reads/{1}/{2}/{3}_{1}_{2}_{4}_{5}_R".format(\
                        self.output,\
                        row['indexST'],\
                        row['indexLOC'],\
                        self.projectName,\
                        self.prefix,\
                        row['indID']\
                    )
                    # Call to ART
                    callParams=["art_illumina"]+self.params+["--in", inputFile,"--out",outputFile]
                    # self.params+=["--in ",inputFile,"--out",outputFile]
                    self.appLogger.debug(" ".join(callParams))

                    proc=""
                    try:
                        proc = subprocess.check_output(callParams,stderr=subprocess.STDOUT)
                    except subprocess.CalledProcessError as error:
                        ngsMessageWrong+="\n------------------------------------------------------------------------\n\n"+\
                        "{}".format(error.output)+\
                        "\n\n------------------------------------------------------------------------"+\
                        "\n\nFor more information about this error please check the log file.\n"+\
                        "You can also run the 'art' command separately.\n\n"+\
                        "art_illumina command used:\n==========================\n"+\
                        "{}\n\n".format(" ".join(callParams))
                        return False, matingArgsMessageWrong

                    cpuTime = [line for line in proc.split('\n') if "CPU" in line][0].split(":")[1]
                    seed = [line for line in proc.split('\n') if "seed" in line][0].split(":")[1]
                    #print simType,cpuTime,seed
            else:
                environment=self.parser.get("execution","environment")
                if (environment=="sge"):
                    self.writeSGEScript()
                if (environment=="slurm"):
                    self.writeSLURMScript()
                else:
                    self.writeBashScript()
        return True,ngsMessageCorrect
