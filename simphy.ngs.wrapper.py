import argparse,datetime,logging,os,sys
import numpy as np
import random as rnd
import MELoggingFormatter as mlf
import SettingsParser as sp
from select import select

################################################################################
# CONSTANTS
VERSION=1
MIN_VERSION=1
FIX_VERSION=0
PROGRAM_NAME="simphy.ngs.wrapper.py"
AUTHOR="Merly Escalona <merlyescalona@uvigo.es>"
################################################################################
# python mating.py -p <prefixES_sequence_filename> -SF <simphy_path>
# python mating.py -p data -sf /media/merly/ME-Conus/conusSim/csTest1/
################################################################################
class SimPhyNGSWrapper:
    def __init__(self,args):
        # I'll check the arguments, and only if, they are correct I'll call
        # any of the correpondent process.
        self.path=os.getcwd()
        self.startTime=datetime.datetime.now()
        self.endTime=None

        self.appLogger=logging.getLogger('sngsw')
        logging.basicConfig(format="%(asctime)s - %(levelname)s:\t%(message)s",\
            datefmt="%d/%m/%Y %I:%M:%S %p",\
            filename="{0}/SimPhyNGSWrapper.{1:%Y}{1:%m}{1:%d}-{1:%H}:{1:%M}:{1:%S}.log".format(\
                self.path,self.startTime),\
            filemode='a',\
            level=logging.DEBUG)
        ch = logging.StreamHandler()
        ch.setFormatter(logging.Formatter(fmt="%(asctime)s - %(levelname)s:\t%(message)s",\
            datefmt="%d/%m/%Y %I:%M:%S %p"))
        ch.setLevel(args.log)
        self.appLogger.addHandler(ch)
        self.appLogger.log(logging.INFO,"Starting")
        self.settingsFile="./settings.txt"
        if (args.settings):
            self.settingsFile=os.path.basename(args.settings)
        if (not os.path.exists(self.settingsFile)):
            self.ending(False,"Settings folder ({0}) does not exist. ".format(self.settingsFile))
        else:
            self.appLogger.debug("Starting process")
            self.run()

    def run(self):
        self.settings=sp.SettingsParser(self.settingsFile)
        good,message=self.settings.checkArgs()
        if (good):
            self.appLogger(message)
            # TODO: here I do something
            self.appLogger.debug("Parse")
        else:   self.ending(good,message) # did not pass the parser reqs.

    def log(self, level, message):
        if level==logging.DEBUG:    self.appLogger.debug(message)
        if level==logging.INFO: self.appLogger.info(message)
        if level==logging.WARNING:  self.appLogger.warning(message)
        if level==logging.ERROR:    self.appLogger.error(message)
        return None

    def getLogLevel(self,level):
        if level==mlf.LOG_LEVEL_CHOICES[0]: loggingLevel=logging.DEBUG
        elif level==mlf.LOG_LEVEL_CHOICES[1]:   loggingLevel=logging.INFO
        elif level==mlf.LOG_LEVEL_CHOICES[2]:   loggingLevel=logging.WARNING
        else:   loggingLevel=logging.ERROR
        return loggingLevel

    def ending(self, good, message):
        # good: Whether there's a good ending or not (error)
        if not good:    self.appLogger.error(message)
        else:   self.appLogger.info(message)
        self.endTime=datetime.datetime.now()
        self.appLogger.info("Elapsed time:\t{0}".format(self.endTime-self.startTime))
        self.appLogger.info("Ending at:\t{0}".format(self.endTime.strftime("%a, %b %d %Y. %I:%M:%S %p")))
        sys.exit()



def handlingParameters():
    parser = argparse.ArgumentParser(\
        prog="{0} (v.{1}.{2}.{3})".format(PROGRAM_NAME,VERSION,MIN_VERSION,FIX_VERSION),\
        formatter_class=argparse.RawDescriptionHelpFormatter,\
        description=\
            '''
SimPhy NGS wrapper
==================
This is a plugin for  SimPhy (https://github.com/adamallo/SimPhy) A comprehensive simulator of gene family
evolution. SimPhy NGS Wrapper generates diploid individuals from the sequences
generated of a SimPhy project, and afterwards generates reads from such
individuals with a next-generation sequencing simulator, ART.

For more information about usage and installation please go to the README.md file or
to the wiki page https://gitlab.com/merlyescalona/simphy-ngs-wrapper/wikis/home

            ''',\
        epilog="Version {0}.{1}.{2} (Still under development)".format(VERSION,MIN_VERSION,FIX_VERSION),\
        add_help=False
        )

    optionalGroup= parser.add_argument_group('Optional arguments')
    optionalGroup.add_argument('-s','--settings',metavar='<settings_file_path>', type=str,\
    help='Path to the settings file.')
    optionalGroup.add_argument('-l','--log',metavar='<log_level>', type=str,\
        choices=mlf.LOG_LEVEL_CHOICES, default="INFO",\
        help='Specified level of log that will be shown through the standard output. Entire log will be stored in a separate file. Values:{0}. Default: {1}. '.format(mlf.LOG_LEVEL_CHOICES,mlf.LOG_LEVEL_CHOICES[1]))
    informationGroup= parser.add_argument_group('Information arguments')
    informationGroup.add_argument('-v', '--version',\
        action='version',\
        version='{0}: Version {1}.{2}.{3}'.format(PROGRAM_NAME,VERSION,MIN_VERSION,FIX_VERSION),\
        help="Show program's version number and exit")
    informationGroup.add_argument('-h', '--help',\
        action='store_true',\
        help="Show this help message and exit")
    try:
        tmpArgs = parser.parse_args()
    except:
        sys.stdout.write("\n----------------------------------------------------------------------\n")
        parser.print_help()
        sys.exit()
    return tmpArgs

"""
###############################   MAIN   #######################################
"""
if __name__=="__main__":
    try:
        cmdArgs = handlingParameters()
        print(cmdArgs)
        prog = SimPhyNGSWrapper(cmdArgs)
        prog.run()
        prog.ending(True,"Run has finished correctly.")
    except KeyboardInterrupt:
        sys.stdout.write("{0}{1}\nInterrupted!{2}\nPlease run again for the expected outcome.\n".format(mof.BOLD,mof.DARKCYAN, mof.END))
        sys.exit()
