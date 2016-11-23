import argparse,datetime,logging,os,path,sys
import numpy as np
from MELoggingFormatter import MELoggingFormatter as mlf
if (sys.version_info[0:2]<(3,0)):
    import ConfigParser as cp
elif (sys.version_info>=(3,0)):
    import configparser as cp

class Settings:
  def __init__(self,filename):
    # If I've got this far, then filename is a correct file
    self.path=os.path.abspath(filename)
    self.appLogger=logging.getLogger('sngsw')
    self.appLogger.debug("(class Settings) __init__()")
    # default settings can be established.
    self.parser=cp.SafeConfigParser()
    self.parser.read(self.path)

  def checkArgs(self):
    allGood=True
    parserMessageCorrect="All parameters are correct."
    parserMessageWrong="Please verify the parameters in the settings file."
    # checking general parameters
    if not (self.parser.has_option("general","data_prefix") or self.parser.has_option("general","dp")):
        parserMessageWrong+="\n\t<data_prefix|dp> field is missing. This prefix correponds to the name of the sequences that are going to be manipulated. Exiting."
        return False, parserMessageWrong

    if not (self.parser.has_option("general","simphy_folder") or
        self.parser.has_option("general","sf")):
        # check if option simphy_folder exist in sections
        parserMessageWrong+="\n\t<simphy_folder> field is missing."
        return False, parserMessageWrong
    else:
        # parameter is set up, now check if folder exist
        path=""
        if (self.parser.has_option("general","simphy_folder")):
            path=os.path.abspath(self.parser.get("general","simphy_folder"))
        if (self.parser.has_option("general","sf")):
            path=os.path.abspath(self.parser.get("general","sf"))
            self.parser.set("general","simphy_folder",path)
            self.parser.remove_option("general","sf")

        if (os.path.exists(path) and os.path.isdir(path)):
            self.appLogger.debug("SimPhy project folder exists")
        else:
            parserMessageWrong+="\n\tSimPhy project folder does not exist, or the given path does not belong to a directory. Exiting."
            return False, parserMessageWrong
        # checking ploidy for the output data
        if (not self.parser.has_option("general","ploidy")):
            self.ploidy=1
        else:
            p=self.parser.getint("general","ploidy")
            if (p>0 and p<=2):
                self.ploidy=p
            elif (p<0):
                self.ploidy=1
            else:
                self.ploidy=2


        # Checking output folder information
        currentRun=""
        if(self.parser.has_option("general","output_folder_name")):
            currentRun=self.parser.get("general","output_folder_name")
        elif (self.parser.has_option("general","ofn")):
            currentRun=self.parser.get("general","output_folder_name")
            self.parser.set("general","output_folder_",currentRun)
            self.parser.remove_option("general","ofn")
        else:
            currentRun="output"

        if os.path.exists("{0}/{1}".format(path,currentRun)):
            listdir=os.listdir(path)
            counter=0
            for item in listdir:
                if currentRun in item:
                    counter+=1
            if not counter == 0:
                currentRun="output_{0}".format(counter+1)
        self.parser.set("general","output_folder_name","{0}/{1}".format(path,currentRun))

        # Checking art parameters.
        if not self.parser.has_section("ngs-reads-art"):
            self.ngsart=False
            self.appLogger.info("No NGS generation section available")
        else:
            self.ngsart=True
            # checking program dependencies
            stream = os.popen('which art_illumina').read()[0:-1]
            self.appLogger.info("Checking dependencies...")
            if stream:
                self.appLogger.info("art_illumina - Found running in: {}".format(stream))
                if self.parser.has_option("ngs-reads-art","o"):self.parser.remove_option("ngs-reads-art","o")
                if self.parser.has_option("ngs-reads-art","out"):self.parser.remove_option("ngs-reads-art","out")
                if self.parser.has_option("ngs-reads-art","i"):self.parser.remove_option("ngs-reads-art","i")
                if self.parser.has_option("ngs-reads-art","in"):self.parser.remove_option("ngs-reads-art","in")
                self.appLogger.warning("Removing I/O options. Be aware: I/O naming is auto-generated from SimPhy and Mating parameters.")
            else:
                parserMessageWrong+="Exiting. art_illumina not found. Program either not installed or not in your current path. Please verify the installation. Exiting. "
                return False, parserMessageWrong

        # Execution parameters
        if not self.parser.has_section("execution"):
            self.parser.add_section("execution")
            self.parser.set("execution", "environment","bash")
            self.parser.set("execution", "run","off")
        else:
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
                    self.parser.set("execution","environment","bash")
                    self.parser.set("execution", "run","off")
            else:
                # got no environment
                self.parser.set("execution", "environment","bash")

            if (self.parser.has_option("execution","run")):
                try:
                    value=self.parser.getboolean("execution","run")
                except Exception as e:
                    self.appLogger.debug("Exception: - {0}".format(e))
                    self.parser.set("execution","run","off")


    self.appLogger.info(self.formatSettingsMessage())
    return True, parserMessageCorrect

  def formatSettingsMessage(self):
    message="Settings:\n"
    sections=self.parser.sections()
    for sec in sections:
        message+="\t{0}\n".format(sec)
        items=self.parser.items(sec)
        for param in items:
            message+="\t\t{0}\t:\t{1}\n".format(param[0],param[1])
    return message

  #
  # def writeSettingsExample1LongNames(self):
  #   self.appLogger.info("Writing settings into file (Example1 - Long names): {0}".format(\
  #       self.path))
  #   parser=cp.RawConfigParser()
  #   parser.add_section("general")
  #   parser.add_section("ngs-reads-art")
  #   parser.set("general","simphy_folder","test")
  #   parser.set("general","data_prefix","data")
  #   parser.set("general","outgroup","true")
  #   parser.set("ngs-reads-art","amplicon","true")
  #   parser.set("ngs-reads-art","rcount ",100)
  #   parser.set("ngs-reads-art","id","iddefault")
  #   parser.set("ngs-reads-art","errfree","false")
  #   parser.set("ngs-reads-art","len",150)
  #   parser.set("ngs-reads-art","mflen",250)
  #   parser.set("ngs-reads-art","paired","true")
  #   parser.set("ngs-reads-art","quiet","true")
  #   parser.set("ngs-reads-art","sdev",50)
  #   parser.set("ngs-reads-art","samout","true")
  #   parser.set("ngs-reads-art","seqSys","HS25")
  #   with open(self.path, 'wb') as configfile:
  #       parser.write(configfile)
