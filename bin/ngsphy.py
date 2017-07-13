#!/usr/bin/env python
# -*- coding: utf-8 -*-
import argparse,datetime,logging,os,threading,sys
import numpy as np
import random as rnd
import ngsphy.individual as ig
import ngsphy.coverage
import ngsphy.sequence as sg
import ngsphy.settings as sp
import ngsphy.reads as ngs
import ngsphy.readcounts as rc
import ngsphy.rerooter as rr
import ngsphy.loggingformatter as lf
from select import select

################################################################################
# CONSTANTS
VERSION=1
MIN_VERSION=0
FIX_VERSION=0
PROGRAM_NAME="ngsphy.py"
AUTHOR="Merly Escalona <merlyescalona@uvigo.es>"
LOG_LEVEL_CHOICES=["DEBUG","INFO","WARNING","ERROR"]
################################################################################

class NGSphy:
	appLogger=None
	path=os.getcwd()
	startTime=None
	endTime=None
	settingsFile=""

	def __init__(self,args):
		self.startTime=datetime.datetime.now()
		self.appLogger=logging.getLogger('ngsphy')
		logging.basicConfig(format="%(asctime)s - %(levelname)s (%(module)s|%(funcName)s:%(lineno)d):\t%(message)s",\
			datefmt="%d/%m/%Y %I:%M:%S %p",\
			filename="{0}/{2}.{1:%Y}{1:%m}{1:%d}-{1:%H}:{1:%M}:{1:%S}.log".format(\
				self.path,self.startTime,PROGRAM_NAME[0:-3].upper()),\
			filemode='a',\
			level=logging.DEBUG)
		ch = logging.StreamHandler()
		loggerFormatter=lf.MELoggingFormatter(fmt="%(asctime)s - %(levelname)s:\t%(message)s",datefmt="%d/%m/%Y %I:%M:%S %p")
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
		try:
			self.appLogger.info("Checking settings...")
			if (not os.path.exists(self.settingsFile)):
				self.ending(False,"Settings file ({0}) does not exist. Exiting. ".format(self.settingsFile))
			else:
				# settings file exist, go ahead and run
				self.appLogger.debug("Starting process")
				self.settings=sp.Settings(self.settingsFile)
				settingsOk,settingsMessage=self.settings.checkArgs()
				# Exit heere if not correnct parser of the settings
				if not settingsOk: self.ending(settingsOk,settingsMessage)
				# Generate BASIC folder structure - output folder
				self.generateFolderStructure()
				# Settings exist and are ok.
				# reroot-tree
				if self.settings.originMode>2:
					# reroot tree
					rerooter=rr.Rerooter(self.settings)
					status, message=rerooter.run()
					if not status: self.ending(status, message)
					rerooter.writeTree()
					self.settings.newickFilePath=rerooter.outputFilePath

				# Generate sequences
				if self.settings.originMode>1:
					self.appLogger.info("Running sequence generator")
					self.seqGenerator=sg.SequenceGenerator(self.settings)
					indelibleStatus,indelibleMessage=self.seqGenerator.run()
					if not (indelibleStatus): self.ending(indelibleStatus, indelibleMessage)

				if self.settings.originMode>0:
					# Generate Individuals (plody independency)
					self.indGenerator=ig.IndividualGenerator(self.settings)
					matingOk,matingMessage=self.indGenerator.checkArgs()
					self.settings.indels=not self.indGenerator.checkFilesForIndels()
					if not matingOk: self.ending(matingOk,matingMessage)
					self.indGenerator.iteratingOverST()

				if self.settings.ngsmode>0:
					self.appLogger.debug("Checking for Coverage. ")
					covGenerator=coverage.CoverageMatrixGenerator(self.settings)
					status,message=covGenerator.calculate()
					if not status: self.ending(status,message)
					# Have to calculate coverage

				# After this I'll have generated the individuals and folder structure
				if self.settings.ngsmode==1:
					self.appLogger.info("NGS Illumina reads - ART mode")
					self.ngs=ngs.ARTIllumina(self.settings)
					status, message=self.ngs.run()
					if not status: self.ending(status,message)
					self.appLogger.info("NGS read simulation process finished. Check log fo status.")
				elif self.settings.ngsmode==2:
					self.appLogger.info("Read counts mode")
					# self.appLogger.info("NGS read simulation is not being made.")
					# If i have read count folder structure must change - i need reference folder
					# print("ngsphy.py - Read count")
					if self.settings.indels:
						self.ending(False,"{0}\n\t{1}\n\t{2}".format(\
							"Read Counts does not support INDELs (for now)",\
							"Check the output folder. Data has been generated.",\
							"Exiting."\
						) )
					print "HA"
					self.readcount=rc.ReadCounts(self.settings)
					status, message=self.readcount.run()
					if not status:
						self.ending(status,message)
				else:
					self.appLogger.info("NGS simulation is not being made.")
		except Exception as ex:
			exc_type, exc_obj, exc_tb = sys.exc_info()
			fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
			message="ngsphy (main): Something is wrong.\n\t{0}\n\t{1}\t{2}\t{3}\n\tPlease verify. Exiting.".format(\
				ex,exc_type,\
				fname, exc_tb.tb_lineno)
			return False, message

		return True,""

	def generateFolderStructure(self):
		self.appLogger.info("Creating basic folder structure.")
		try:
			self.appLogger.info("Creating output folder: {0} ".format(self.settings.outputFolderPath))
			os.makedirs(self.settings.outputFolderPath)
		except:
			self.appLogger.debug("Output folder ({0}) exists. ".format(self.settings.outputFolderPath))

	def ending(self, good, message):
		self.endTime=datetime.datetime.now()
		# good: Whether there's a good ending or not (error)
		if not good:
			self.appLogger.error(message)
			self.appLogger.error("Elapsed time (ETA):\t{0}".format(self.endTime-self.startTime))
			self.appLogger.error("Ending at:\t{0}".format(self.endTime.strftime("%a, %b %d %Y. %I:%M:%S %p")))
		else:
			self.appLogger.info(message)
			self.appLogger.info("Elapsed time (ETA):\t{0}".format(self.endTime-self.startTime))
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

###############################   MAIN   #######################################
if __name__=="__main__":
# def main():
	try:
		cmdArgs = handlingCmdArguments()
		# print(cmdArgs)
		prog = NGSphy(cmdArgs)
		prog.run()
		prog.ending(True,"Run has finished correctly.")

	except KeyboardInterrupt:
		sys.stdout.write("{0}{1}\nProgram has been interrupted!{2}\nPlease run again for the expected outcome.\n".format("\033[91m","\033[1m","\033[0m"))
		sys.exit()
