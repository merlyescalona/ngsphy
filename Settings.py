import argparse,datetime,logging,os,sys
import numpy as np
from MELoggingFormatter import MELoggingFormatter as mlf
from NGSPhyDistribution import NGSPhyDistribution as ngsphydistro
if (sys.version_info[0:2]<(3,0)):
    import ConfigParser as cp
elif (sys.version_info>=(3,0)):
    import configparser as cp


INDELIBLE="indelible"
SIMPHY="simphy"

class Settings:
    # General
    ploidy=1
    projectName=""
    path=""
    simphy=False
    indelible=False
    dataprefix=""
    outputFolderName=""

    # simphy data origin
    simphyPath=""
    filter_simphy=False
    # indelible data origin
    ngsphyindeliblecontrol=""
    newickFile=""
    evolve=""
    partition=""

    #readcount
    readcount=False
    seqerror=0
    reference=""
    # ngsart
    ngsart=False

    def __init__(self,filename):
        # If I've got this far, then filename is a correct file
        self.path=os.path.abspath(filename)
        self.appLogger=logging.getLogger('ngsphy')
        self.appLogger.debug("(class Settings) __init__()")
        # default settings can be established.
        self.parser=cp.SafeConfigParser()
        self.parser.read(self.path)

    def checkArgs(self):
        allGood=True
        parserMessageCorrect="All parameters are correct."
        parserMessageWrong="Settings - Problem found! "
        statusData,messageData=self.checkSectionData(\
            parserMessageCorrect,parserMessageWrong
        )
        if not statusData:
            return statusData,messageData
        statusGeneral,messageGeneral= self.checkSectionGeneral(parserMessageCorrect,parserMessageWrong)
        if(statusGeneral):

            if self.parser.has_section("ngs-reads-art"):
                statusNGSArt,messageNGSArt=self.checkSectionNGSReadsArt(\
                    parserMessageCorrect,parserMessageWrong)
                if (statusNGSArt):
                    ## Next check
                        statusCoverage,messageCoverage=self.checkSectionCoverage(parserMessageCorrect,parserMessageWrong)
                        if not statusCoverage:
                            return statusCoverage, messageCoverage
                else:
                    return statusNGSArt,messageNGSArt
                    # Exit here
                if self.parser.has_section("read-count"):
                    self.parser.remove_section("read-count")
                    self.appLogger.warning("[read-count] section is incompatible with [ngs-reads-art]. Omiting this section.")
            else:
                self.ngsart=False
                self.appLogger.info("Settings: No NGS generation section available")
                # readcount
                if self.parser.has_section("read-count"):
                    statusRC,messageRC=self.checkSectionReadCount(parserMessageCorrect,parserMessageWrong)
                    if statusRC:
                        statusCoverage,messageCoverage=self.checkSectionCoverage(parserMessageCorrect,parserMessageWrong)
                        if not statusCoverage:
                            return statusCoverage, messageCoverage
                    else:
                        return statusRC,messageRC
                else:
                    self.readcount=False
                    self.appLogger.info("Settings: No read-count generation section available")

                if not (self.parser.has_section("read-count") or self.parser.has_section("ngs-reads-art")):
                    if (self.parser.has_section("coverage")):self.parser.remove_section("coverage")
                    self.appLogger.info("[coverage] section is not needed when [ngs-reads-art] nor [read-count] section are available. Omiting this section.")

            self.checkSectionExecution(parserMessageCorrect,parserMessageWrong)
        else:
            return statusGeneral,messageGeneral
            # Exit here
        self.appLogger.info(self.formatSettingsMessage())
        return allGood,parserMessageCorrect

    def checkSectionGeneral(self,parserMessageCorrect,parserMessageWrong):
        # Check GENERAL SECTION
        if not self.parser.has_section("general"):
            parserMessageWrong+="\n\t[general] section missing and required. Please verify. Exiting."
            return False, parserMessageWrong
        # CHECKING GENERAL PARAMETERS
        # data prefix
        if not (self.parser.has_option("general","data_prefix") or self.parser.has_option("general","dp")):
            parserMessageWrong+="\n\t<data_prefix | dp> field is missing. This prefix correponds to the name of the file sequences that are going to be processed. Exiting."
            return False, parserMessageWrong
        else:
            if (self.parser.has_option("general","dp")):
                value=self.parser.get("general","dp")
                self.parser.set("general","data_prefix",value)
                self.parser.remove_option("general","dp")
            if (self.parser.has_option("general","data_prefix")):
                self.dataprefix=self.parser.get("general","data_prefix")
        # checking ploidy for the output data
        if (not self.parser.has_option("general","ploidy")):
            self.ploidy=1
        else:
            p=self.parser.getint("general","ploidy")
            if (p>0 and p<=2):  self.ploidy=p
            elif (p<0): self.ploidy=1
            else:   self.ploidy=2

        # project information
        if (self.parser.has_option("general","project_name")):
            self.projectName=self.parser.get("general","project_name")
        else:
            parserMessageWrong+="\n\t Project Name option is needed. Please verify. Exiting."
            return False, parserMessageWrong

        if (self.parser.has_option("general","path")):
            self.path=os.path.abspath(self.parser.get("general","path"))
        else:
            parserMessageWrong+="\n\t Working path option is needed. Please verify. Exiting."
            return False, parserMessageWrong

        # Checking output folder information
        currentRun=""
        if (self.parser.has_option("general","ofn")):
            currentRun=self.parser.get("general","ofn")
            self.parser.set("general","output_folder_name",currentRun)
            self.parser.remove_option("general","ofn")
        if(self.parser.has_option("general","output_folder_name")):
            currentRun=self.parser.get("general","output_folder_name")
        else:
            currentRun="output"


        dataOrigin=self.parser.get("data","origin").lower()
        if (dataOrigin==INDELIBLE):
            try:
                os.makedirs("{0}/{1}".format(self.path,self.projectName))
            except:
                self.appLogger.warning("This folder exists.")

        if not os.path.exists("{0}/{1}/".format(self.path,self.projectName)):
            parserMessageWrong+="\n\tProject path: {0} does not exist.\n\tPlease verify. Exiting.".format(\
                "{0}/{1}/".format(self.path,self.projectName))
            return False, parserMessageWrong
        if os.path.exists("{0}/{1}/{2}/".format(self.path,self.projectName,currentRun)):
            listdir=os.listdir("{0}/{1}".format(self.path,self.projectName))
            counter=0
            for item in listdir:
                if currentRun in item:
                    counter+=1
            if not counter == 0:
                currentRun="output_{0}".format(counter+1)

        self.outputFolderName="{0}/{1}/{2}/".format(self.path,self.projectName,currentRun)
        print self.outputFolderName
        self.parser.set("general","output_folder_name","{0}/{1}/{2}".format(self.path,self.projectName,currentRun))
        return True,parserMessageCorrect

    def checkSectionData(self,parserMessageCorrect,parserMessageWrong):
        # DATA
        if (self.parser.has_section("data")):
            # simphy
            if (self.parser.has_option("data","origin") and self.parser.get("data","origin").lower()==SIMPHY):
                self.simphy=True
                self.indelible=False
                if (os.path.exists(os.path.join(self.path,self.projectName)) and os.path.isdir(os.path.join(self.path,self.projectName))):
                    self.appLogger.debug("SimPhy project folder exists")
                else:
                    parserMessageWrong+="\n\tSimPhy project folder does not exist, or the given path does not belong to a directory. Exiting."
                    return False, parserMessageWrong
                # checking flag for filtering species tree replicates
                if self.parser.has_option("data","filter_simphy"):
                    self.filter_simphy=True
                else:
                    self.filter_simphy=False
                # Remove options that do not belong to this execution
                if (self.parser.has_option("data","indelible_control")):
                    self.parser.remove_option("data","indelible_control")
                if (self.parser.has_option("data","newick_file")):
                    self.parser.remove_option("data","newick_file")

                # check if simphy is valid projectName
                status,message=self.checkSimPhyProjectValid()
                if not status:
                    return status, message
            # INDELible FILE
            elif (self.parser.has_option("data","origin") and self.parser.get("data","origin").lower()==INDELIBLE):
                self.indelible=True
                self.simphy=False
                # parameter is set up, now check if folder exist
                if (self.parser.has_option("data","indelible_control")):
                    self.ngsphyindeliblecontrol=os.path.abspath(\
                        self.parser.get("data","indelible_control"))
                if (os.path.exists(self.ngsphyindeliblecontrol) and os.path.isfile(self.ngsphyindeliblecontrol)):
                    self.appLogger.debug("INDELible control file exists")
                else:
                    parserMessageWrong+="\n\tINDELible control file does not exist. Exiting."
                    return False, parserMessageWrong
                if (self.parser.has_option("data","newick_file")):
                    self.newickFile=os.path.abspath(self.parser.get("data","newick_file").strip())
                    filesOk=(os.path.exists(self.newickFile) and os.path.isfile(self.newickFile))
                    if not filesOk:
                        parserMessageWrong+="\n\t The Newick file does not exist ({0}). Please verify. Exiting.".format(\
                            fileToCheck
                        )
                        return False, parserMessageWrong
                if self.parser.has_option("data","filter_simphy"):
                    self.parser.remove("data","filter_simphy")
                # checking indelible program
                stream = os.popen('which indelible').read()[0:-1]
                self.appLogger.info("Checking dependencies...")
                if not stream:
                    parserMessageWrong+="\n\t indelible not found. Program either not installed or not in your current path. Please verify the installation. Exiting."
                    return False, parserMessageWrong
            else:
                parserMessageWrong+="\n\t There is no data origin. Please verify. Exiting."
                return False, parserMessageWrong
        else:
            parserMessageWrong+="\n\t There is no [data] section. Please verify. Exiting."
            return False, parserMessageWrong
        return True, parserMessageCorrect



    def checkSectionNGSReadsArt(self,parserMessageCorrect,parserMessageWrong):
        ########################################################################
        # BLOCK: NGS-READS-ART
        ########################################################################
        # Checking art parameters.
        self.ngsart=True
        self.appLogger.info("NGS-reads-ART option selected.")
        # checking program dependencies
        stream = os.popen('which art_illumina').read()[0:-1]
        self.appLogger.info("Checking dependencies...")
        if stream:
            self.appLogger.info("art_illumina - Found running in: {}".format(stream))
            # Coverage parameters
            if self.parser.has_option("ngs-reads-art","o"):self.parser.remove_option("ngs-reads-art","o")
            if self.parser.has_option("ngs-reads-art","out"):self.parser.remove_option("ngs-reads-art","out")
            if self.parser.has_option("ngs-reads-art","i"):self.parser.remove_option("ngs-reads-art","i")
            if self.parser.has_option("ngs-reads-art","in"):self.parser.remove_option("ngs-reads-art","in")
            self.appLogger.warning("Removing I/O options. Be aware: I/O naming is auto-generated from SimPhy and Mating parameters.")
            # Coverage parameters
            if (self.parser.has_option("ngs-reads-art","fcov")): self.parser.remove_option("ngs-reads-art","fcov")
            if (self.parser.has_option("ngs-reads-art","f")): self.parser.remove_option("ngs-reads-art","f")
            if (self.parser.has_option("ngs-reads-art","rcount")): self.parser.remove_option("ngs-reads-art","rcount")
            if (self.parser.has_option("ngs-reads-art","c")): self.parser.remove_option("ngs-reads-art","c")

            self.appLogger.warning("Removing ART coverage options. Coverage is calculated with the [coverage] section (experimentCoverage and individualCoverage options).")
        else:
            parserMessageWrong+="art_illumina not found. Program either not installed or not in your current path. Please verify the installation. Exiting."
            return False, parserMessageWrong
        return True, parserMessageCorrect

    def checkSectionReadCount(self,parserMessageCorrect,parserMessageWrong):
        ########################################################################
        # BLOCK: READ COUNT
        ########################################################################
        # experimentCoverage /expCov
        # individualCoverage /indCov
        message=parserMessageCorrect
        if (self.parser.has_section("read-count")):
            self.readcount=True
            if not self.parser.has_option("read-count", "error"):
                self.appLogger.warning("[read-count] section. Sequencing error rate for this run is being considered as 0.")
                self.parser.set("read-count", "error","0")
            if not self.parser.has_option("read-count","reference"):
                self.appLogger.warning("[read-count] section. Using default references.")
                self.parser.set("read-count", "reference","None")
        else:
            # No read-count section
            self.readcount=False
            message="[read-count] section. Not available."

        return self.readcount,message

    ########################################################################
    # BLOCK: Coverage
    ########################################################################
    def checkSectionCoverage(self,parserMessageCorrect,parserMessageWrong):
        message=parserMessageCorrect
        expCov=None;indCov=None;locCov=None;
        if(self.parser.has_section("coverage")):
            if (self.parser.has_option("coverage","expCov")):
                value=self.parser.get("coverage","expCov")
                self.parser.set("coverage","experimentCoverage",value.lower())
                self.parser.remove_option("coverage","expCov")
            elif (self.parser.has_option("coverage","experimentCoverage")):
                value=self.parser.get("coverage","experimentCoverage")
                self.parser.set("coverage","experimentCoverage",value.lower())
            else:
                # parsear distribution
                parserMessageWrong+="Coverage section | Experiment Coverage distribution variable is required. Please verify. Exiting."
                return False,parserMessageWrong
            expCov=ngsphydistro(0,self.parser.get("coverage","experimentCoverage"), None,False)
            check,mess=expCov.coverageDistroCheck()
            if not (check):
                parserMessageWrong+=mess
                return check,parserMessageWrong
            # If i got here I have EXPERIMET COVERAGE DISTRIBUTION
            if (self.parser.has_option("coverage","individualCoverage") or (self.parser.has_option("coverage","indCov"))):
                if (self.parser.has_option("coverage","indCov")):
                    value=self.parser.get("coverage","indCov")
                    self.parser.set("coverage","individualCoverage",value.lower())
                    self.parser.remove_option("coverage","indCov")
                if (self.parser.has_option("coverage","individualCoverage")):
                    value=self.parser.get("coverage","individualCoverage")
                    self.parser.set("coverage","individualCoverage",value.lower())
                    print value
                    print expCov.value()
                    indCov=ngsphydistro(1,self.parser.get("coverage","individualCoverage"), expCov,False)
                check,mess=indCov.coverageDistroCheck()
                if not (check):
                    parserMessageWrong+=mess
                    return check,parserMessageWrong
                distro=None
            else:
                if (self.parser.has_option("coverage","locusCoverage") or (self.parser.has_option("coverage","locCov"))):
                    if (self.parser.has_option("coverage","locCov")):
                        self.parser.remove_option("coverage","locCov")
                    if (self.parser.has_option("coverage","locusCoverage")):
                        self.parser.remove_option("coverage","locusCoverage")
                    self.appLogger.warning("Locus-wise coverage option is being removed because it depends on Individual-wise coverage option, and is missing.")

            # If i got here I have EXPERIMENT COVERAGE DISTRIBUTION,
            # if i hae EXPERIMENT COVERAGE but NO INDIVIDUAL COVERAGE, then' there's no problem, i wont have locCoverage either
            # it will have been removed
            # I can keep going with the coverage loc validation
            if (self.parser.has_option("coverage","locusCoverage") or (self.parser.has_option("coverage","locCov"))):
                if (self.parser.has_option("coverage","locCov")):
                    value=self.parser.get("coverage","locCov")
                    self.parser.set("coverage","locusCoverage",value.lower())
                    self.parser.remove_option("coverage","locCov")
                if (self.parser.has_option("coverage","locusCoverage")):
                    value=self.parser.get("coverage","locusCoverage")
                    self.parser.set("coverage","locusCoverage",value.lower())
                locCov=ngsphydistro(2,self.parser.get("coverage","locusCoverage"), indCov,False)
                check,mess=locCov.coverageDistroCheck()
                if not (check):
                    parserMessageWrong+=mess
                    return check,parserMessageWrong
            else:
                if (self.parser.has_option("coverage","individualCoverage")):
                    self.appLogger.info("Using Experiment- and Individual-wise coverage.")
                else:
                    self.appLogger.info("Using only Experiment-wise coverage.")
        else:
            # No coverage section
            message="Settings: Coverage section | When using [ngs-reads-art] or [read-count] section. Coverage is required. Please verify. Exiting."
            return False,message

        return True, message

    def checkSectionExecution(self,parserMessageCorrect,parserMessageWrong):
        ########################################################################
        # BLOCK: Execution
        ########################################################################
        if not self.parser.has_section("execution"):
            self.appLogger.warning("Settings - Execution block: This block has been automatically generated.")
            self.parser.add_section("execution")
            self.parser.set("execution", "environment","bash")
            self.parser.set("execution", "run","off")
            self.parser.set("execution", "threads","1")
        else:
            ####################################################################
            # OPTION: Environment
            if (self.parser.has_option("execution","env")):
                # got the short name
                value=self.parser.get("execution","env")
                self.parser.set("execution","environment",value.lower())
                self.parser.remove_option("execution","environment")
            elif (self.parser.has_option("execution","environment")):
                # got the long name, make sure it is lowercase and within the options
                value=self.parser.get("execution","environment")
                if (value in ["sge","slurm","bash"]):
                    self.parser.set("execution","environment",value.lower())
                    if (value in ["sge","slurm"]):
                        self.parser.set("execution", "run","off")
                else:
                    message="Settings: Execution block | Evironment variable is incorrect or unavailable. Please check the settings file and rerun. Exiting."
                    return False,message
            else:
                # got no environment
                self.parser.set("execution", "environment","bash")
            ####################################################################
            # OPTION: RUN
            if (self.parser.has_option("execution","run")):
                try:
                    value=self.parser.getboolean("execution","run")
                except Exception as e:
                    self.appLogger.warning("Settings - Execution block: Run automatically set up to OFF.")
                    self.parser.set("execution","run","off")
            else:
                self.appLogger.warning("Settings - Execution block: Run automatically set up to OFF.")
                self.parser.set("execution","run","off")
            ####################################################################
            # OPTION: threads
            if (self.parser.has_option("execution","threads")):
                try:
                    self.numThreads=self.parser.getint("execution","threads")
                except Exception as e:
                    self.appLogger.warning("Settings - Execution block: Threads automatically set up to 1.")
                    self.parser.set("execution","threads","1")
                    self.numThreads=1
            else:
                self.numThreads=1
                self.appLogger.warning("Settings - Execution block: Threads automatically set up to 1.")
                self.parser.set("execution","threads","1")


    def checkIndelibleControlFile(self,parserMessageCorrect,parserMessageWrong):
        f=open(self.ngsphyindeliblecontrol,"r")
        lines=f.readLines()
        f.close()
        newlines=[ item.strip() for item in lines if item.strip()!=""]
        # check for NGSPHY blocks
        model=0; partition=0; evolve=0
        for index in range(0,len(newlines)):
            item=newlines[index]
            if item =="[MODEL]" and model != 0:
                if model!=0:
                    parserMessageWrong+="Validating INDELible Control file: {0}\n\t{1}".format(\
                        "There is more than one [MODEL] block.",\
                        "Please verify. Exiting"\
                        )
                    return False, parserMessageWrong
                model=index
            if item =="[NGSPHYPARTITION]" and partition != 0:
                if partition!=0:
                    parserMessageWrong+="Validating INDELible Control file: {0}\n\t{1}".format(\
                        "There is more than one [NGSPHYPARTITION] block.",\
                        "Please verify. Exiting"\
                        )
                    return False, parserMessageWrong
                partition=index
            if item =="[NGSPHYEVOLVE]" and evolve != 0:
                if evolve!=0:
                    parserMessageWrong+="Validating INDELible Control file: {0}\n\t{1}".format(\
                        "There is more than one [NGSPHYEVOLVE] block.",\
                        "Please verify. Exiting"\
                    )
                    return False, parserMessageWrong
                evolve=index

        # if i have gotten here and did not return before
        # i can check if the blocks match :D
        # 1) tree filename matches partition parameter and partition has the right
        # number of parameters
        self.partition=newlines[partition].split()
        if len(self.partition)!=4:
            parserMessageWrong+="Validating INDELible Control file: {0}\n\t{1}\n\t{2}".format(\
            "[NGSPHYPARTITION] block has the wrong number of parameters.",\
            "[NGSPHYPARTITION] <tree_filename_basename> <model_name> <sequence_length>"
            "Please verify. Exiting"
            )
            return False, parserMessageWrong
        # check tree corresponds to newick inputbasename
        newickBasename,_=os.path.splitext(os.path.basename(self.newickFiles))
        if not self.partition[1]==newickBasename:
            parserMessageWrong+="Validating INDELible Control file: {0}\n\t{1}\n\t{2}\n\t{3}".format(\
            "[NGSPHYPARTITION] block, tree name does not correspond with the Newick File introduced.",\
            "Remember! Newick filename: newick.tree.",\
            "[NGSPHYPARTITION] newick model1 200"
            "Please verify. Exiting"
            )
            return False, parserMessageWrong

        modelname=newlines[model].strip().split()[1]
        if not self.partition[2]==modelname:
            parserMessageWrong+="Validating INDELible Control file: {0}\n\t{1}\n\t{2}\n\t{3}".format(\
            "[NGSPHYPARTITION] block, model name does not correspond to the model defined.",\
            "[MODEL] modelname\n...",\
            "[NGSPHYPARTITION] tree modelname 200"
            "Please verify. Exiting"
            )
            return False, parserMessageWrong

        if not self.partition[3].isdigit():
            parserMessageWrong+="Validating INDELible Control file: {0}\n\t{1}".format(\
            "[NGSPHYPARTITION] block, sequence length is not valid.",\
            "Please verify. Exiting"
            )
            return False, parserMessageWrong

        # 2) evolve has the correct number of parameter
        self.evolve=newlines[evolve].split()
        if len(self.evolve)!=3:
            parserMessageWrong+="Validating INDELible Control file: {0}\n\t{1}\n\t{2}".format(\
            "[NGSPHYEVOLVE] block has the wrong number of parameters.",\
            "[NGSPHYEVOLVE] <number_sequence_per_gene_tree> <prefix_output_filename>"
            "Please verify. Exiting"
            )
            return False, parserMessageWrong

        if not self.evolve[1].isdigit():
            parserMessageWrong+="Validating INDELible Control file: {0}\n\t{1}".format(\
            "[NGSPHYEVOLVE] block, number of replicates is not valid.",\
            "Please verify. Exiting"
            )
            return False, parserMessageWrong

        if int(self.evolve[1])== 0:
            parserMessageWrong+="Validating INDELible Control file: {0}\n\t{1}\n\t{2}".format(\
            "[NGSPHYEVOLVE] block, number of replicates is not valid.",\
            "Number of replicates must be greater than 0 (nReplicates > 0).",\
            "Please verify. Exiting"
            )
            return False, parserMessageWrong
        # need newick files and control file partitions
        return True, parserMessageCorrect

    def checkSimPhyProjectValid(self):
        matingArgsMessageWrong="SimPhy project is not valid."
        matingArgsMessageCorrect="SimPhy project is valid."
        # List all the things in the project directory
        fileList=os.listdir(os.path.join(self.path,self.projectName))
        for index in range(0,len(fileList)):
            fileList[index]=os.path.abspath(os.path.join(self.path,self.projectName,fileList[index]))

        command = os.path.join(\
            self.path,\
            self.projectName,\
            "{0}.command".format(self.projectName))
        params = os.path.join(\
            self.path,\
            self.projectName,\
            "{0}.params".format(self.projectName))
        db = os.path.join(\
            self.path,\
            self.projectName,\
            "{0}.db".format(self.projectName))

        self.appLogger.debug("SimPhy files (command, params, db)")
        self.appLogger.debug("{0}:\t{1}".format(\
            os.path.basename(db),db in fileList))
        self.appLogger.debug("{0}:\t{1}".format(\
            os.path.basename(command),command in fileList))
        self.appLogger.debug("{0}:\t{1}".format(\
            os.path.basename(params),params in fileList))

        simphyfiles=((command in fileList) and (params in fileList) and(db in fileList))
        # check if  command, db, params files
        if not simphyfiles:
            matingArgsMessageWrong+="\n\tSimPhy files do not exist. Please verify. Exiting."
            return False, matingArgsMessageWrong
        # check how many of them are dirs
        for item in fileList:
            baseitem=os.path.basename(item)
            if (os.path.isdir(os.path.abspath(item)) and  baseitem.isdigit()):
                numSpeciesTrees=numSpeciesTrees+1

        numSpeciesTreesDigits=len(str(numSpeciesTrees))
        self.settings.parser.set("general","numSpeciesTrees",str(numSpeciesTrees))
        # check if at least one
        self.appLogger.debug("Num species trees:\t{0}".format(numSpeciesTrees))
        if not (numSpeciesTrees>0):
            matingArgsMessageWrong+="\n\tNot enough number of species tree replicates (at least 1 required):\t{0}\n\t Please verify. Exiting.".format(self.numSpeciesTrees>0)
            return False, matingArgsMessageWrong

        matingArgsMessageCorrect+="\n{0}".format(message)
        return True, matingArgsMessageCorrect

    def formatSettingsMessage(self):
        message="Settings:\n"
        sections=self.parser.sections()
        for sec in sections:
            message+="\t{0}\n".format(sec)
            items=self.parser.items(sec)
            for param in items:
                message+="\t\t{0}\t:\t{1}\n".format(param[0],param[1])
        return message
