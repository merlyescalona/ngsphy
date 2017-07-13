#!/usr/bin/env python
import argparse, sys
from ngsphy import ngsphy
################################################################################
# CONSTANTS
VERSION=1
MIN_VERSION=0
FIX_VERSION=0
PROGRAM_NAME="ngsphy.py"
AUTHOR="Merly Escalona <merlyescalona@uvigo.es>"
LOG_LEVEL_CHOICES=["DEBUG","INFO","WARNING","ERROR"]
################################################################################
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

###############################   MAIN   #######################################
if __name__=="__main__":
# def main():
	try:
		cmdArgs = handlingCmdArguments()
		# print(cmdArgs)
		prog = ngsphy.NGSphy(cmdArgs)
		prog.run()
		prog.ending(True,"Run has finished correctly.")

	except KeyboardInterrupt:
		sys.stdout.write("{0}{1}\nProgram has been interrupted!{2}\nPlease run again for the expected outcome.\n".format("\033[91m","\033[1m","\033[0m"))
		sys.exit()
