import argparse,datetime,logging,os,subprocess,sys
import numpy as np
import random as rnd
import Mating as mat
import Settings as sp
import NGSReads as ngs
from MELoggingFormatter import MELoggingFormatter as mlf
from select import select

################################################################################
# CONSTANTS
VERSION=1
MIN_VERSION=1
FIX_VERSION=0
PROGRAM_NAME="ngsphy.py"
AUTHOR="Merly Escalona <merlyescalona@uvigo.es>"
LOG_LEVEL_CHOICES=["DEBUG","INFO","WARNING","ERROR"]
################################################################################

class NGSphy:
    def __init__(self,args):
        self.path=os.getcwd()
        self.startTime=datetime.datetime.now()
        self.endTime=None

        self.appLogger=logging.getLogger('sngsw')
        logging.basicConfig(format="%(asctime)s - %(levelname)s (%(module)s|%(funcName)s:%(lineno)d):\t%(message)s",\
            datefmt="%d/%m/%Y %I:%M:%S %p",\
            filename="{0}/{2}.{1:%Y}{1:%m}{1:%d}-{1:%H}:{1:%M}:{1:%S}.log".format(\
                self.path,self.startTime,PROGRAM_NAME[0:-3]),\
            filemode='a',\
            level=logging.DEBUG)
        ch = logging.StreamHandler()
        loggerFormatter=mlf(fmt="%(asctime)s - %(levelname)s:\t%(message)s",datefmt="%d/%m/%Y %I:%M:%S %p")
        ch.setFormatter(loggerFormatter)
        ch.setLevel(args.log.upper())
        self.appLogger.addHandler(ch)
        self.appLogger.info("Starting")
        self.settingsFile=""
        if (args.settings):
            self.settingsFile=os.path.abspath(args.settings)
        else:
            self.settingsFile=os.path.abspath("./settings.txt")

    def run(self):
        # checking existence of settings file
        self.appLogger.info("Checking settings...")
        if (not os.path.exists(self.settingsFile)):
            self.ending(False,"Settings file ({0}) does not exist. Exiting. ".format(self.settingsFile))
        else:
            # settings file exist, go ahead and run
            self.appLogger.debug("Starting process")
            self.settings=sp.Settings(self.settingsFile)
            settingsOk,settingsMessage=self.settings.checkArgs()
            if (settingsOk):
                # Settings exist and are ok.
                # Generate Individuals (plody independency)
                self.mating=mat.Mating(self.settings)
                matingOk,matingMessage=self.mating.checkArgs()
                if (matingOk):
                    self.mating.iteratingOverST()
                else:
                    self.ending(matingOk,matingMessage) # did not pass the parser reqs.

                # After this I'll have generated the individuals and the folder structure
                if self.settings.ngsart:
                    # Doing NGS
                    self.ngs=ngs.NGSReadsARTIllumina(self.settings)
                    self.ngs.run()
                    self.appLogger.info("NGS read simulation process finished. Check log fo status.")

                else:
                    self.appLogger.info("NGS read simulation is not being done.")
            else:
                self.ending(settingsOk,settingsMessage) # did not pass the parser reqs.

    def log(self, level, message):
        if level==logging.DEBUG:    self.appLogger.debug(message)
        if level==logging.INFO: self.appLogger.info(message)
        if level==logging.WARNING:  self.appLogger.warning(message)
        if level==logging.ERROR:    self.appLogger.error(message)
        return None

    def getLogLevel(self,level):
        if level.upper()==LOG_LEVEL_CHOICES[0]: loggingLevel=logging.DEBUG
        elif level.upper()==LOG_LEVEL_CHOICES[1]:   loggingLevel=logging.INFO
        elif level.upper()==LOG_LEVEL_CHOICES[2]:   loggingLevel=logging.WARNING
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



def handlingCmdArguments():
    parser = argparse.ArgumentParser(\
        prog="{0} (v.{1}.{2}.{3})".format(PROGRAM_NAME,VERSION,MIN_VERSION,FIX_VERSION),\
        formatter_class=argparse.RawDescriptionHelpFormatter,\
        description=\
            '''
NGSphy
==================
NGSphy is a tool designed specifically as an addendum to SimPhy (https://github.com/adamallo/SimPhy)
- a phylogenomic simulator of gene, locus and species trees that considers incomplete lineage sorting,
gene duplication and loss and horizontal gene transfer - which is able to use SimPhy's output in order
to produce Illumina NGS data from haploid/diploid individuals.

For more information about usage and installation please go to the README.md file or
to the wiki page https://gitlab.com/merlyescalona/ngsphy/wikis/home

            ''',\
        epilog="Version {0}.{1}.{2} (Still under development)".format(VERSION,MIN_VERSION,FIX_VERSION),\
        add_help=False
        )

    optionalGroup= parser.add_argument_group('Optional arguments')
    optionalGroup.add_argument('-s','--settings',metavar='<settings_file_path>', type=str,\
    help='Path to the settings file.')
    optionalGroup.add_argument('-l','--log',metavar='<log_level>', type=str,\
        choices=LOG_LEVEL_CHOICES, default="INFO",\
        help='Specified level of log that will be shown through the standard output. Entire log will be stored in a separate file. Values:{0}. Default: {1}. '.format(LOG_LEVEL_CHOICES,LOG_LEVEL_CHOICES[1]))
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
        if (tmpArgs.help): parser.print_help()
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
        cmdArgs = handlingCmdArguments()
        print(cmdArgs)
        prog = NGSphy(cmdArgs)
        prog.run()
        prog.ending(True,"Run has finished correctly.")
    except KeyboardInterrupt:
        sys.stdout.write("{0}{1}\nInterrupted!{2}\nPlease run again for the expected outcome.\n".format("\033[91m","\033[1m","\033[0m"))
        sys.exit()
