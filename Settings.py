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
    self.parser=cp.ConfigParser()
    self.parser.read(self.path)

  def checkArgs(self):
    allGood=True
    parserMessageCorrect="All parameters are correct."
    parserMessageWrong="Please verify the parameters in the settings file."
    # checking general parameters
    if not (self.parser.has_option("general","data_prefix") or self.parser.has_option("general","dp")):
        self.parser.set("general","data_prefix","data")
        self.appLogger.info("Setting default data prefix due to lacking of entry in the settings file.")

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
            self.appLogger.debug("Correct SimPhy project folder")
        else:
            parserMessageWrong+="\n\tSimPhy project folder does not exist, or the given path does not belong to a directory."
            return False, parserMessageWrong

        # Checking output folder information
        currentRun=""
        if(self.parser.has_option("general","output_folder")):
            currentRun=self.parser.get("general","output_folder")
        elif (self.parser.has_option("general","of")):
            currentRun=self.parser.get("general","output_folder")
            self.parser.set("general","output_folder",currentRun)
            self.parser.remove_option("general","of")
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

        self.parser.set("general","output_folder","{0}/{1}".format(path,currentRun))

        # Checking mating parameters.
        if not self.parser.has_section("mating"):
            # can run
            parserMessageWrong+="\n\tNo MATING section. Stopping run."
            return False, parserMessageWrong
        elif not (self.parser.has_option("mating","opi") or self.parser.has_option("mating","output-per-individual")):
            # I have mating section but no option
            self.parser.set("mating", "output-per-individual",1)
        else:
            # I have mating section and option
            # check which option I have, make sure it is stored in the long_name option
            # if short exist, remove, also, checking num_output_per_individual ranges
            num_output_per_individual=1
            if self.parser.has_option("mating","output-per-individual"):
                num_output_per_individual=self.parser.getint("mating","output-per-individual")
            if self.parser.has_option("mating","opi"):
                num_output_per_individual=self.parser.getint("mating","opi")
                if (num_output_per_individual > 3): num_output_per_individual=3
                if (num_output_per_individual < 1): num_output_per_individual=1
                self.parser.set("mating","output-per-individual",num_output_per_individual)
                self.parser.remove_option("mating","opi")

        # Checking art parameters.
        if not self.parser.has_section("ngs-reads-art"):
            # if no section, no run, if the section exist, I will leave it
            # to art to complain if parameters are wrong
            parserMessageWrong+="\n\tNo ART section. Stopping run."
            return False, parserMessageWrong
        else:
            if self.parser.has_option("ngs-reads-art","o"):self.parser.remove_option("ngs-reads-art","o")
            if self.parser.has_option("ngs-reads-art","out"):self.parser.remove_option("ngs-reads-art","out")
            if self.parser.has_option("ngs-reads-art","i"):self.parser.remove_option("ngs-reads-art","i")
            if self.parser.has_option("ngs-reads-art","in"):self.parser.remove_option("ngs-reads-art","in")
            self.appLogger.warning("Removing I/O options. Be aware: I/O naming is auto-generated from SimPhy and Mating parameters.")

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


  def writeSettingsExample1LongNames(self):
    self.appLogger.info("Writing settings into file (Example1 - Long names): {0}".format(\
        self.path))
    parser=cp.RawConfigParser()
    parser.add_section("general")
    parser.add_section("mating")
    parser.add_section("ngs-reads-art")
    parser.set("general","simphy_folder","test")
    parser.set("general","data_prefix","data")
    parser.set("mating","output-per-individual",1)
    parser.set("ngs-reads-art","amplicon","true")
    parser.set("ngs-reads-art","rcount ",100)
    parser.set("ngs-reads-art","id","iddefault")
    parser.set("ngs-reads-art","errfree","false")
    parser.set("ngs-reads-art","len",150)
    parser.set("ngs-reads-art","mflen",250)
    parser.set("ngs-reads-art","paired","true")
    parser.set("ngs-reads-art","quiet","true")
    parser.set("ngs-reads-art","sdev",50)
    parser.set("ngs-reads-art","samout","true")
    parser.set("ngs-reads-art","seqSys","HS25")
    with open(self.path, 'wb') as configfile:
        parser.write(configfile)

  def writeSettingsExample1ShortNames(self):
    self.appLogger.info("Writing settings into file (Example1 - Short names): {0}".format(\
      self.path))
    parser=cp.RawConfigParser()
    parser.add_section("general")
    parser.add_section("mating")
    parser.add_section("ngs-reads-art")
    parser.set("general","sf" ,"test")
    parser.set("general","dp" ,"data")
    parser.set("mating","opi",1)
    parser.set("ngs-reads-art","amp","true")
    parser.set("ngs-reads-art","c",100)
    parser.set("ngs-reads-art","d","iddefault")
    parser.set("ngs-reads-art","ef" ,"false")
    parser.set("ngs-reads-art","l",150)
    parser.set("ngs-reads-art","m",250)
    parser.set("ngs-reads-art","p","true")
    parser.set("ngs-reads-art","q","true")
    parser.set("ngs-reads-art","s",50)
    parser.set("ngs-reads-art","sam","true")
    parser.set("ngs-reads-art","ss","HS25")
    with open(self.path, 'wb') as configfile:
        parser.write(configfile)

  def writeSettingsExample2LongNames(self):
    self.appLogger.info("Writing settings into file (Example2 - Long names): {0}".format(\
        self.path))
    parser=cp.RawConfigParser()
    parser.add_section("general")
    parser.add_section("mating")
    parser.add_section("ngs-reads-art")
    parser.set("general","simphy_folder","test")
    parser.set("general","data_prefix","data")
    parser.set("mating","output-per-individual",1)
    parser.set("ngs-reads-art","qprof1","profileR1.txt")
    parser.set("ngs-reads-art","qprof2","profileR2.txt")
    parser.set("ngs-reads-art","amplicon","true")
    parser.set("ngs-reads-art","rcount ",100)
    parser.set("ngs-reads-art","id","iddefault")
    parser.set("ngs-reads-art","errfree","false")
    parser.set("ngs-reads-art","len",150)
    parser.set("ngs-reads-art","mflen",250)
    parser.set("ngs-reads-art","paired","true")
    parser.set("ngs-reads-art","quiet","true")
    parser.set("ngs-reads-art","sdev",50)
    parser.set("ngs-reads-art","samout","true")
    with open(self.path, 'wb') as configfile:
        parser.write(configfile)

  def writeSettingsExample2ShortNames(self):
    self.appLogger.info("Writing settings into file (Example2 - Short names): {0}".format(\
        self.path))
    parser=cp.RawConfigParser()
    parser.add_section("general")
    parser.add_section("mating")
    parser.add_section("ngs-reads-art")
    parser.set("general","sf","test")
    parser.set("general","dp","data")
    parser.set("mating","opi",1)
    parser.set("ngs-reads-art","1","profileR1.txt")
    parser.set("ngs-reads-art","2","profileR2.txt")
    parser.set("ngs-reads-art","amp","true")
    parser.set("ngs-reads-art","c",100)
    parser.set("ngs-reads-art","d","iddefault")
    parser.set("ngs-reads-art","ef","false")
    parser.set("ngs-reads-art","l",150)
    parser.set("ngs-reads-art","m",250)
    parser.set("ngs-reads-art","p","true")
    parser.set("ngs-reads-art","q","true")
    parser.set("ngs-reads-art","s",50)
    parser.set("ngs-reads-art","sam","true")
    with open(self.path, 'wb') as configfile:
        parser.write(configfile)
