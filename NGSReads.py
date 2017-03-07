import argparse,csv,datetime,logging,os,subprocess,sys,threading, time
from NGSPhyDistribution import NGSPhyDistribution as ngsphydistro
import numpy as np
import random as rnd
import Settings as sp

class RunningInfo(object):
    def __init__(self):
        # self.appLogger=logging.getLogger('sngsw')
        self.lock = threading.Lock()
        self.value = []
    def addLine(self, line):
        # self.appLogger.debug('Waiting for lock')
        self.lock.acquire()
        try:
            # self.appLogger.debug('Acquired lock')
            self.value += [line]
        finally:
            # self.appLogger.debug('Released lock')
            self.lock.release()

class NGSReadsARTIllumina:
    __SHORT_NAMES=["sf" ,"dp","ploidy","ofn","1","2","amp","c","d","ef","f","h","i",\
                "ir","ir2","dr","dr2","k","l","m","mp","M","nf","na",\
                "o","p","q","qU","qs","qL","qs2","rs","s","sam","sp","ss"]
    __LONG_NAMES=["simphy_folder","data_prefix","output_folder_name","ploidy","qprof1","qprof2",\
                "amplicon","rcount","id","errfree","fcov","help",\
                "in","insRate","insRate2","delRate","delRate2","maxIndel",\
                "len","mflen","matepair","cigarM","maskN","noALN","out",\
                "paired","quiet","maxQ","qShift","minQ","qShift2","rndSeed",\
                "sdev","samout","sepProf","seqSys"]
    __dLONG_NAMES={i.lower():i for i in __LONG_NAMES}
    __dSHORT_NAMES={i.lower():i for i in __SHORT_NAMES}

    params=[]
    commands=[]     # Init of all the commands that will be generated
    numFiles=0

    # file path related variables
    ploidyName="individuals"
    readsFolder=""
    scriptsFolder=""
    coverageFolder=""

    # Coverage Distribution variables
    experimentCoverageDistro=None
    individualCoverageDistro=None
    locusCoverageDistro=None

    def __init__(self,settings):
        self.appLogger=logging.getLogger('ngsphy')
        self.appLogger.info('NGS read simulation: ART run started.')
        self.settings=settings
        simphy=os.path.abspath(self.settings.parser.get("general", "simphy_folder"))
        self.output=os.path.abspath(\
            os.path.join(\
                simphy,\
                self.settings.parser.get("general", "output_folder_name")))
        self.filteredST=self.settings.parser.get("general", "filtered_ST")
        self.filteredST=[ int(numST) for numST in self.filteredST.split(",")]
        numSTs=self.settings.parser.getint("general","number_ST")
        self.numberSTDigits=len(str(numSTs))
        self.numLociPerST=[int(numST) for numST in self.settings.parser.get("general","numLociPerST").split(",")]
        self.indexLOCDigits=len(str(max(self.numLociPerST)))
        self.projectName=os.path.basename(simphy)
        if (simphy[-1]=="/"):
            self.projectName=os.path.basename(simphy)[0:-1]
        self.prefix=self.settings.parser.get("general","data_prefix")

        dash=""; par=[]
        settingsParams=self.settings.parser.items("ngs-reads-art")
        for p in settingsParams:
            if (p[0] in self.__dSHORT_NAMES.keys()): dash="-"
            if (p[0] in self.__dLONG_NAMES.keys()): dash="--"
            # to be sure that i am getting the right parameter names
            if (dash=="-"):
                if (p[0]=="m" and ((par[1].lower() in ["true","false","on","off"]) or (par[1] in [0,1]))):
                    par=["M"]
                elif (p[0]=="m"):
                    par=["m"]
                else:
                    par=[self.__dSHORT_NAMES[p[0]]]
            if (dash=="--"):
                par=[self.__dLONG_NAMES[p[0]]]

            par+=[p[1]]

            if((par[1].lower() in ["true","false","on","off"]) or (par[1] in [0,1])):
                self.params+=["{0}{1}".format(dash,par[0])]
            else:
                self.params+=["{0}{1}".format(dash,par[0]),par[1]]

        # generating specific folder structure
        self.generateFolderStructure()
        # Coverage related options
        self.experimentCoverageDistro=ngsphydistro(0,\
            self.settings.parser.get("coverage","experimentCoverage"))
        if (self.settings.parser.has_option("coverage","individualCoverage")):
            self.individualCoverageDistro=ngsphydistro(1,\
                self.settings.parser.get("coverage","individualCoverage"),\
                True)
        if (self.settings.parser.has_option("coverage","locusCoverage")):
            self.locusCoverageDistro=ngsphydistro(2,\
                self.settings.parser.get("coverage","locusCoverage"),\
                True)

        self.runningInfo=RunningInfo()


    def generateFolderStructure(self):
        self.appLogger.info("Creating folder structure for [ngs-reads-art]")
        self.readsFolder="{0}/reads".format(self.output)
        self.scriptsFolder="{0}/scripts".format(self.output)
        self.coverageFolder="{0}/coverage".format(self.output)

        try:
            os.makedirs(self.readsFolder)
            self.appLogger.info("Generating output folder ({0})".format(self.readsFolder))
        except:
            self.appLogger.debug("Output folder exists ({0})".format(self.readsFolder))

        try:
            os.makedirs(self.scriptsFolder)
            self.appLogger.info("Generating output folder ({0})".format(self.scriptsFolder))
        except:
            self.appLogger.debug("Output folder exists ({0})".format(self.scriptsFolder))

        try:
            os.makedirs(self.coverageFolder)
            self.appLogger.info("Generating output folder ({0})".format(self.coverageFolder))
        except:
            self.appLogger.debug("Output folder exists ({0})".format(self.coverageFolder))

    def writeSeedFile(self):
        seedfile="{0}/scripts/{1}.seedfile.txt".format(\
            self.output,\
            self.projectName
        )
        sf=open(seedfile,"w")
        for indexST in self.filteredST:
            for indexST in self.filteredST:
                csvfile=open("{0}/tables/{1}.{2:0{3}d}.{4}.csv".format(\
                    self.output,\
                    self.projectName,\
                    indexST,\
                    self.numberSTDigits,\
                    self.ploidyName
                ))
                # Generation of folder structure
                d = csv.DictReader(csvfile)
                self.matingDict = [row for row in d]
                csvfile.close()
                for indexLOC in range(1,self.numLociPerST[indexST-1]+1):
                    for row in self.matingDict:
                        # indexST,indexLOC,indID,speciesID,mateID1,mateID2
                        inputFile="{0}/individuals/{1}/{2:0{3}d}/{4}_{1}_{2:0{3}d}_{5}_{6}.fasta".format(\
                            self.output,\
                            row['indexST'],\
                            indexLOC,\
                            self.indexLOCDigits,
                            self.projectName,\
                            self.prefix,\
                            row['indID']\
                        )
                        # This means, from a multiple (2) sequence fasta file.
                        outputFile="{0}/reads/{1}/{2:0{3}d}/{4}_{1}_{2:0{3}d}_{5}_{6}_R".format(\
                            self.output,\
                            row['indexST'],\
                            indexLOC,\
                            self.indexLOCDigits,
                            self.projectName,\
                            self.prefix,\
                            row['indID']\
                        )
                        sf.write("{0}\t{1}\n".format(inputFile,outputFile))
                self.numFiles+=1
        sf.close()
        self.appLogger.info("Seed file written...")

    def writeSGEScript(self):
        jobfile="{0}/scripts/{1}.job.sge.sh".format(\
            self.output,\
            self.projectName
        )
        j=open(jobfile,"w")
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

inputfile=$(awk 'NR==$SGE_TASK_ID{{print $1}}' {1})
outputfile=$(awk 'NR==$SGE_TASK_ID{{print $2}}' {1})\n""".format(self.numFiles,seedfile)
        j.write(header)
        j.write(" ".join(callParams))
        footer="".format()
        j.write(footer)
        j.close()
        self.appLogger.info("SGE Job file written ({0})...".format(jobfile))

    def writeSLURMScript(self):
        jobfile="{0}/scripts/{1}.job.slurm.sh".format(\
            self.output,\
            self.projectName
        )
        j=open(jobfile,"w")
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

inputfile=$(awk 'NR==$SLURM_ARRAY_TASK_ID{{print $1}}' {1})
outputfile=$(awk 'NR==$SLURM_ARRAY_TASK_ID{{print $2}}' {1})

""".format(self.numFiles, seedfile)
        footer="".format()
        jobfile.write(header)
        jobfile.write(" ".join(callParams))
        jobfile.write(footer)
        jobfile.close()
        self.appLogger.info("SLURM Job file written ({0})...".format(jobfile))

    def computeCoverageMatrix(self, nInds, nLoci,expCov, indCov, locCov):
    	# coverage matrix per ST - row:indv - col:loci
    	# each cov introduced as parameter is a NGSPhyDistribution
    	self.appLogger.debug("Computing coverage matrix")
        covMatrix=np.zeros(shape=(nInds, nLoci))
        for loc in range(0,nLoci):
            for ind in range(0,nInds):
                expc=0;indc=0;locc=0
                expc=expCov.value()
                coverage=expc
                if (indCov):
                    indCov.updateValue(expc)
                    indc=indCov.value()
                    coverage=indc
                if (locCov):
                    locCov.updateValue(indc)
                    locc=locCov.value()
                    coverage=locc
                covMatrix[ind][loc]=coverage
        return covMatrix

    def writeCoverageMatrixIntoFile(self, coverageMatrix, filename):
        filepath=os.path.abspath(filename)
        with open(filepath, 'w') as csvfile:
            writer = csv.writer(csvfile)
            [writer.writerow(r) for r in coverageMatrix]

    def getCommands(self):
        for indexST in self.filteredST:
            csvfile=open("{0}/tables/{1}.{2:0{3}d}.{4}.csv".format(\
                self.output,\
                self.projectName,\
                indexST,\
                self.numberSTDigits,\
                self.ploidyName
            ))
            # Generation of folder structure
            d = csv.DictReader(csvfile)
            self.matingDict = [row for row in d]
            csvfile.close()
            nInds=len(self.matingDict)
            nLoci=self.numLociPerST[indexST-1]
            coverageMatrix=self.computeCoverageMatrix(nInds,nLoci,\
                self.experimentCoverageDistro,\
                self.individualCoverageDistro,\
                self.locusCoverageDistro)
            coverageMatrixFilename="{0}/coverage/{1}.{2:0{3}d}.csv".format(\
                self.output,\
                self.projectName,\
                indexST,\
                self.numberSTDigits\
            )
            self.writeCoverageMatrixIntoFile(coverageMatrix,coverageMatrixFilename)
            for indexLOC in range(1,self.numLociPerST[indexST-1]+1):
                for row in self.matingDict:
                    # indexST,indexLOC,indID,speciesID,mateID1,mateID2
                    inputFile="{0}/individuals/{1:0{2}d}/{3:0{4}d}/{5}_{1}_{3:0{4}d}_{6}_{7}.fasta".format(\
                        self.output,\
                        int(row['indexST']),\
                        self.numberSTDigits,\
                        indexLOC,\
                        self.indexLOCDigits,
                        self.projectName,\
                        self.prefix,\
                        int(row['indID']))
                    # This means, from a multiple (2) sequence fasta file.
                    outputFile="{0}/reads/{1:0{2}d}/{3:0{4}d}/{5}_{1}_{3:0{4}d}_{6}_{7}_R".format(\
                        self.output,\
                        int(row['indexST']),\
                        self.numberSTDigits,\
                        indexLOC,\
                        self.indexLOCDigits,
                        self.projectName,\
                        self.prefix,\
                        int(row['indID'])\
                    )
                    coverage=coverageMatrix[int(row['indID'])][indexLOC-1]
                    # Call to ART
                    callParams=["art_illumina"]+self.params+["--fcov",str(coverage),"--in", inputFile,"--out",outputFile]
                    # self.params+=["--in ",inputFile,"--out",outputFile]
                    # print(callParams)
                    self.commands+=[[row['indexST'],indexLOC,row['indID'],inputFile, outputFile]+callParams]

        self.appLogger.info("Commands have been generated...")

    def writeBashScript(self):
        bashfile="{0}/scripts/{1}.sh".format(\
            self.output,\
            self.projectName
        )
        j=open(bashfile,"w")
        for item in self.commands:
            c=item[5:len(item)]
            j.write(" ".join(c))
            j.write("\n")
        j.close()
        self.appLogger.info("Bash script written...")

    def commandLauncher(self, command):
        ngsMessage="";proc=""
        try:
            proc = subprocess.check_output(command[5:],stderr=subprocess.STDOUT)
            cpuTime = [line for line in proc.split('\n') if "CPU" in line][0].split(":")[1]
            seed = [line for line in proc.split('\n') if "seed" in line][0].split(":")[1]
            # line=[(command[0:3],cpuTime,seed, command[4])]
            line=command[0:3]+[cpuTime,seed,command[4]]
            # ngsMessage="Command: {0} - Finished succesfully.".format(" ".join(command[5:]))
        except subprocess.CalledProcessError as error:
            ngsMessage="{}".format(error.output)+\
            "\n\n------------------------------------------------------------------------"+\
            "\n\nFor more information about this error please run the 'art' command separately.\n"+\
            "art_illumina command used:\n==========================\n"+\
            "{}\n\n".format(" ".join(command[5:]))
            self.appLogger.error(ngsMessage)
            line=command[0:3]+["-","-",command[4]]
            raise RuntimeError("\nART execution error. Please verify. Exciting.")
        self.runningInfo.addLine(line)


    def run(self):
        status=True
        environment=self.settings.parser.get("execution","environment")
        # Generating commands
        try:
            self.getCommands()
            if (environment=="sge"):
                self.appLogger.info("Environment SGE. Writing scripts")
                self.writeSeedFile()
                self.writeSGEScript()
            elif (environment=="slurm"):
                self.appLogger.info("Environment SLURM. Writing scripts")
                self.writeSeedFile()
                self.writeSLURMScript()
            else:
                message="ART run has finished succesfully."
                self.appLogger.info("Environment BASH. Writing scripts")
                self.writeBashScript()
                run=self.settings.parser.getboolean("execution","run")
                if (run):
                    self.generateFolderStructureNGS()
                    curr=0
                    self.appLogger.info("Running...")
                    while (curr < len(self.commands)):
                        progress=(curr*100)/len(self.commands)
                        sys.stdout.write("Progress {0:02.1f} %\r".format(progress))
                        sys.stdout.flush()
                        # print("Progress {0:02.1f} %\r".format(prog))
                        command=self.commands[curr]
                        t = threading.Thread(target=self.commandLauncher(command))
                        t.start()
                        curr=curr+1
                        while (threading.activeCount()-1)==self.settings.numThreads:
                            time.sleep(0.1)

                    self.printRunningInfo()

        except ValueError as verror:
            status=False
            message="\nDistributon parameter error:\n\t{}\nPlease verify. Exciting.".format(verror)
        except RuntimeError as rte:
            status=False
            message=rte
        return status, message

    def generateFolderStructureNGS(self):
        # iterating over commands to create folders
        folders=set([])
        for command in self.commands:
            infile=os.path.dirname(command[4])
            folders.add(infile)
        for item in folders:
            try:
                os.makedirs(item)
            except:
                self.appLogger.debug("Output folder exists ({0})".format(item))

    def printRunningInfo(self):
        outputFile="{0}/{1}.info".format(
            self.output,\
            self.projectName
        )
        f=open(outputFile,"w")
        f.write("indexST,indexLOC,indID,inputFile,cpuTime,seed,outputFilePrefix\n")
        for item in self.runningInfo.value:
            f.write(
                str(item[0])+","+\
                str(item[1])+","+\
                str(item[2])+","+\
                str(item[3])+","+\
                str(item[4])+","+\
                item[5]+"\n"
            )
        f.close()
        self.appLogger.info("File with timings of the ART run can be find on: {0}".format(outputFile))
