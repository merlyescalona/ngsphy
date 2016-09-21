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
    self.settings=cp.ConfigParser()
    self.settings.read(self.path)


  def checkArgs(self):
    allGood=True
    parserMessageCorrect="All parameters are correct."
    parserMessageWrong="Please verify the parameters in the settings file."
    # checking general parameters
    # check if option simphy_folder exist in sections
    # check if folder exist

    # if (!os.path.exist(self.settings.get("secti")))

    if (allGood):
        self.appLogger.info(self.formatSettingsMessage())
        return True, parserMessageCorrect
    else: return False, parserMessageCorrect

  def formatSettingsMessage(self):
    message="Settings:\n"
    sections=self.settings.sections()
    for sec in sections:
        message+="\t{0}\n".format(sec)
        items=self.settings.items(sec)
        for param in items:
            message+="\t\t{0}\t:\t{1}\n".format(param[0],param[1])
    return message


  def writeSettingsExample1LongNames(self):
    self.appLogger.info("Writing settings into file (Example1 - Long names): {0}".format(\
        self.path))
    parser=cp.RawConfigParser()
    parser.add_section("general")
    parser.add_section("ngs-reads-art")
    parser.set("general","simphy_folder","simphyfolder")
    parser.set("general","data_prefix","data")
    parser.set("ngs-reads-art","amplicon","true")
    parser.set("ngs-reads-art","rcount ",100)
    parser.set("ngs-reads-art","id","iddefault")
    parser.set("ngs-reads-art","errfree","false")
    parser.set("ngs-reads-art","in","input.fq")
    parser.set("ngs-reads-art","len",150)
    parser.set("ngs-reads-art","mflen",250)
    parser.set("ngs-reads-art","out","output")
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
    parser.set("general","sf" ,"simphyfolder")
    parser.set("general","dp" ,"data")
    parser.set("ngs-reads-art","amp","true")
    parser.set("ngs-reads-art","c",100)
    parser.set("ngs-reads-art","d","iddefault")
    parser.set("ngs-reads-art","ef" ,"false")
    parser.set("ngs-reads-art","i","input.fq")
    parser.set("ngs-reads-art","l",150)
    parser.set("ngs-reads-art","m",250)
    parser.set("ngs-reads-art","o","output")
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
    parser.set("general","simphy_folder","simphyfolder")
    parser.set("general","data_prefix","data")
    parser.set("ngs-reads-art","qprof1","profileR1.txt")
    parser.set("ngs-reads-art","qprof2","profileR2.txt")
    parser.set("ngs-reads-art","amplicon","true")
    parser.set("ngs-reads-art","rcount ",100)
    parser.set("ngs-reads-art","id","iddefault")
    parser.set("ngs-reads-art","errfree","false")
    parser.set("ngs-reads-art","in","input.fq")
    parser.set("ngs-reads-art","len",150)
    parser.set("ngs-reads-art","mflen",250)
    parser.set("ngs-reads-art","out","output")
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

    parser.set("general","sf","simphyfolder")
    parser.set("general","dp","data")
    parser.set("ngs-reads-art","1","profileR1.txt")
    parser.set("ngs-reads-art","2","profileR2.txt")
    parser.set("ngs-reads-art","amp","true")
    parser.set("ngs-reads-art","c",100)
    parser.set("ngs-reads-art","d","iddefault")
    parser.set("ngs-reads-art","ef","false")
    parser.set("ngs-reads-art","i","input.fq")
    parser.set("ngs-reads-art","l",150)
    parser.set("ngs-reads-art","m",250)
    parser.set("ngs-reads-art","o","output")
    parser.set("ngs-reads-art","p","true")
    parser.set("ngs-reads-art","q","true")
    parser.set("ngs-reads-art","s",50)
    parser.set("ngs-reads-art","sam","true")
    with open(self.path, 'wb') as configfile:
        parser.write(configfile)
