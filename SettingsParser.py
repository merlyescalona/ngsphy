import argparse,datetime,logging,os,path,sys
import numpy as np
from MELoggingFormatter import MELoggingFormatter as mlf

class SettingsParser:
  def __init__(self,filename):
      # If I've got this far, then filename is a correct file
      self.path=os.path.abspath(filename)
      self.appLogger=logging.getLogger('sngsw')
      self.appLogger.debug("(SettingsParser) init ...")

  def checkArgs(self):
    parserMessageCorrect="All parameters are correct."
    parserMessageWrong="Please verify the parameters in the settings file."
    if (allGood): return True, parserMessageCorrect
    else: return False, parserMessageCorrect
