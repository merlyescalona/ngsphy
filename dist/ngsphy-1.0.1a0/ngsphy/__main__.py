#!/usr/bin/env python
import argparse,datetime,logging, ngsphy, os,sys
import loggingformatter as lf
################################################################################
# CONSTANTS
VERSION=1
MIN_VERSION=0
FIX_VERSION=0
PROGRAM_NAME="ngsphy.py"
AUTHOR="Merly Escalona <merlyescalona@uvigo.es>"
INSTITUTION="University of Vigo, Spain."
LOG_LEVEL_CHOICES=["DEBUG","INFO","WARNING","ERROR"]
LINE="--------------------------------------------------------------------------------"
################################################################################
# Logger initialization
APPLOGGER=logging.getLogger('ngsphy')
ch = logging.StreamHandler()
loggerFormatter=lf.MELoggingFormatter(fmt="%(asctime)s - %(levelname)s:\t%(message)s",datefmt="%d/%m/%Y %I:%M:%S %p")
ch.setFormatter(loggerFormatter)
APPLOGGER.addHandler(ch)
################################################################################
def createLogFile():
	logging.basicConfig(format="%(asctime)s - %(levelname)s (%(module)s|%(funcName)s:%(lineno)d):\t%(message)s",\
		datefmt="%d/%m/%Y %I:%M:%S %p",\
		filename="{0}/{2}.{1:%Y}{1:%m}{1:%d}-{1:%H}:{1:%M}:{1:%S}.log".format(\
			os.getcwd(),datetime.datetime.now(),PROGRAM_NAME[0:-3].upper()),\
		filemode='a',\
		level=logging.DEBUG)

def handlingCmdArguments():
	"""
	handlingCmdArguments
	--------------------

	This function configurates the ArgumentParser with the specific details of
	this programs.

	Does not take parameters.
	"""
	parser = argparse.ArgumentParser(\
		prog="{0} (v.{1}.{2}.{3})".format(PROGRAM_NAME,VERSION,MIN_VERSION,FIX_VERSION),\
		formatter_class=argparse.RawDescriptionHelpFormatter,\
		description=\
			'''\033[1m
================================================================================
  NGSphy
================================================================================
\033[0m
NGSphy is a Python open-source tool for the genome-wide simulation of NGS data
(read counts or Illumina reads) obtained from thousands of gene families evolving
under a common species tree, with multiple haploid and/or diploid individuals per
species, where sequencing coverage (depth) heterogeneity can vary among
individuals and loci, including off-target loci and phylogenetic decay effects.

For more information about usage and installation please go to the README file
or to the wiki page https://gihub.com/merlyescalona/ngsphy/wiki/
			''',\
		epilog="Developed by:\n{0}\n{1}\n\nVersion:\t{2}.{3}.{4} (Under development)\n{5}\n".format(AUTHOR,INSTITUTION, VERSION,MIN_VERSION,FIX_VERSION, LINE),\
		add_help=False
		)
	optionalGroup= parser.add_argument_group("{0}Optional arguments{1}".format("\033[1m","\033[0m"))
	optionalGroup.add_argument('-s','--settings',metavar='<settings_file_path>', type=str,\
	help='Path to the settings file.')
	optionalGroup.add_argument('-l','--log',metavar='<log_level>', type=str,\
		choices=LOG_LEVEL_CHOICES, default="INFO",\
		help='Specified level of log that will be shown through the standard output. Entire log will be stored in a separate file. Values:{0}. Default: {1}. '.format(LOG_LEVEL_CHOICES,LOG_LEVEL_CHOICES[1]))
	informationGroup= parser.add_argument_group("{0}Information arguments{1}".format("\033[1m","\033[0m"))
	informationGroup.add_argument('-v', '--version',\
		action='version',\
		version='{0}: Version {1}.{2}.{3}'.format(PROGRAM_NAME,VERSION,MIN_VERSION,FIX_VERSION),\
		help="Show program's version number and exit")
	informationGroup.add_argument('-h', '--help',\
		action='store_true',\
		help="Show this help message and exit")

	status=True; message=""
	try:
		tmpArgs = parser.parse_args()
		if (tmpArgs.help): parser.print_help()
		sys.stdout("\n{}\n".format(LINE))
		# Checks if the log level debug is chosen it will also print the debug file
		if tmpArgs.log==LOG_LEVEL_CHOICES[0]:
			createLogFile()
	except:
		message="{0}\n{1}\n{2}\n".format(\
			"Something happened while parsing the arguments.",\
			 "Please verify. Exiting.", LINE)
		APPLOGGER.error(message)
		parser.print_help()
		sys.exit(-1)
	if not tmpArgs.settings and not os.path.exists(os.path.abspath(os.path.join(os.getcwd(),"settings.txt"))):
		parser.print_help()
	return tmpArgs

###############################   MAIN   #######################################
def main():
	"""
	main()
	--------------------

	Entry point of the NGSphy program.
	"""
	try:
		cmdArgs = handlingCmdArguments()
		prog = ngsphy.NGSphy(cmdArgs)
		ch.setLevel(cmdArgs.log.upper())
		prog.run()
	except ngsphy.NGSphyException as ex:
		if ex.expression: sys.exit()
		else: sys.exit(-1)
	except KeyboardInterrupt:
		APPLOGGER.error("{0}{1}\nProgram has been interrupted.{2}\nPlease run again for the expected outcome.\n{3}\n".format("\033[91m","\033[1m","\033[0m",LINE))
		sys.exit(-1)

if __name__=="__main__":
	main()
