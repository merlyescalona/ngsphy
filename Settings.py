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
            self.appLogger.debug("SimPhy project folder exists")
        else:
            parserMessageWrong+="\n\tSimPhy project folder does not exist, or the given path does not belong to a directory."
            return False, parserMessageWrong

        # Checking output folder information
        currentRun=""
        if(self.parser.has_option("general","output_folder_name")):
            currentRun=self.parser.get("general","output_folder_name")
        elif (self.parser.has_option("general","ofn")):
            currentRun=self.parser.get("general","output_folder_name")
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

        self.parser.set("general","output_folder_name","{0}/{1}".format(path,currentRun))
        # if there is no outgroup info no problem, mating will be done accordingly
        if self.parser.has_option("general", "og"):
            value=self.parser.getboolean("general","og")
            self.parser.set("general","outgroup")
            self.parser.remove_option("general","og",value)


        # Checking art parameters.
        if not self.parser.has_section("ngs-reads-art"):
            self.ngsart=False
            # I can have or not the section. If exist check, "else"
            # otherwise, i have to check if I have something to do, that being,
            # having the mating section, if not I'm just exiting the program
            self.appLogger.info("No NGS generation section available")
            if not self.mating:
                parserMessageWrong+="\n\tNo mating, nor NGS generation section. Exiting."
                return False, parserMessageWrong
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
    parser.add_section("ngs-reads-art")
    parser.set("general","simphy_folder","test")
    parser.set("general","data_prefix","data")
    parser.set("general","outgroup","true")
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
    parser.add_section("ngs-reads-art")
    parser.set("general","sf" ,"test")
    parser.set("general","og","true")
    parser.set("general","dp" ,"data")
    parser.set("ngs-reads-art","amp","true")
    parser.set("ngs-reads-art","c",100)
    parser.set("ngs-reads-art","d","iddefault")
    parser.set("ngs-reads-art","ef")
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
    parser.add_section("ngs-reads-art")
    parser.set("general","simphy_folder","test")
    parser.set("general","data_prefix","data")
    parser.set("general","outgroup","true")
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
    parser.add_section("ngs-reads-art")
    parser.set("general","sf","test")
    parser.set("general","dp","data")
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
