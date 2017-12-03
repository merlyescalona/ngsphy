#!/usr/bin/env python
# -*- coding: utf-8 -*-
import argparse,copy,datetime,dendropy,logging,os,re,sys, multiprocessing,msatools,subprocess,platform,stat
import numpy as np
import random as rnd
import settings as sp
from msatools import *
from select import select

try:
	from collections import Counter
except ImportError:
	from counter import Counter

class SequenceGenerator:
	"""
	Class for the generation of genome sequences from gene trees
	----------------------------------------------------------------------------
	Attributes:
	- appLogger: logger to store status of the process flow
	- settings: Settings object withh all the program parameters
	- newIndelibleAncestralSequence: path of the reference.txt file that will be
	 used to call the sequence simulator INDELible.
	- newIndelibleControlFilePath: path of the control file that will be used to call
	the sequence simulator INDELible
	- output: path where the output of the INDELible execution will be stored
	- evolve: default parameters of the [EVOLVE] section of the INDELible control
	file
	- partition: parameters of the [PARTITION] section of the INDELible control
	file
	"""
	appLoger=None
	settings=None

	ancestralFreq=dict({"A":0.25, "C":0.25, "G":0.25,"T":0.25})
	newIndelibleAncestralSequence=""
	newIndelibleControlFilePath=""
	geneTreeFile=""

	evolve=[1,"ngsphydata_1"]
	partition=[]

	def __init__(self,settings):
		self.appLogger=logging.getLogger("ngsphy")
		self.appLogger.debug('INDELible run')
		self.settings=settings
		self.settings.alignmentsFolderPath=os.path.join(self.settings.alignmentsFolderPath,"1")
		self.newIndelibleControlFilePath=os.path.join(\
			self.settings.alignmentsFolderPath,"control.txt")
		self.newIndelibleAncestralSequence=os.path.join(\
			self.settings.alignmentsFolderPath,"ancestral.fasta")
		if self.settings.inputmode==3: # INDELible + reference
			self.geneTreeFile=os.path.join(\
				self.settings.alignmentsFolderPath,\
				"ngsphy.tree"\
			)
		self.generateFolderStructure()

	def run(self):
		"""
		Process flow for the generation of genome sequences from a gene tree
		and an evolutionary model
		"""
		self.writeIndelibleControlFile()

		runStatus,runMessage=self.runIndelible()
		if not runStatus: return False,runMessage
		return True, "Run finished"

	def generateFolderStructure(self):
		"""
		Generation of a folder structure for this process.
		"""
		self.appLogger.info("Creating folder structure for INDELible run")
		try:
			os.makedirs(os.path.join(self.settings.alignmentsFolderPath))
			self.appLogger.info("Generating data folder ({0})".format(\
				os.path.join(self.settings.alignmentsFolderPath)))
		except:
			self.appLogger.debug("Data folder exists ({0})".format(\
				os.path.join(self.settings.alignmentsFolderPath)))


	def copyAncestralSequenceToOutputFolder(self):
		"""
		In order to generate genome sequences all the required files must be
		in the same folder where INDELible is going to be ran. Hence, the need
		of copying the given reference file to the directory where data will
		be stored.
		-----------------------------------------------------------------------
		Returns:
		- boolean. Indicates the status of the process.
		"""
		# making sure there's only one sequence, and only one sequence will be written to the
		# reference.fasta file
		# that sequence will be the first from the file if there are more than 1 sequence
		status=True; message=""
		self.appLogger.debug("Copying reference sequence file ")
		self.appLogger.info("Copying reference sequence file to: {}".format(\
			self.newIndelibleAncestralSequence))
		description=""
		try:
			with open(self.settings.ancestralSequenceFilePath, "r") as f:
				description=f.readline().strip()
		except Exception as ex:
			message="\n\t{0}\n\t{1}\n\t{2}\n\t{3}\n".format(\
				"I/O problem.",\
				ex,
				"Stopped while reading the ancestral sequence file.",\
				"Please verify and rerun. Exiting."
			)
			status=False
			return status, message
		description=description[1:len(description)]
		referenceDict=msatools.parseMSAFileWithDescriptions(self.settings.ancestralSequenceFilePath)
		reference=referenceDict[description]
		self.getAncestralSequenceBaseFrequencies(reference)
		try:
			fout=open(self.newIndelibleAncestralSequence,"w")
			fout.write(">ngsphypartition\n{}\n".format(reference))
			fout.close()
		except Exception as ex:
			message="\n\t{0}\n\t{1}\n\t{2}\n\t{3}\n".format(\
				"I/O problem.",\
				ex,\
				"Stopped while copying the ancestral sequence file.",\
				"Please verify and rerun. Exiting."
			)
			status=False
			return status, message
		return status, message

	def getAncestralSequenceBaseFrequencies(self,reference):
		counter=Counter(reference)
		totalBases=len(reference)
		self.ancestralFreq["A"]=float(int(counter["A"])/float(totalBases))
		self.ancestralFreq["C"]=float(int(counter["C"])/float(totalBases))
		self.ancestralFreq["G"]=float(int(counter["G"])/float(totalBases))
		self.ancestralFreq["T"]=float(int(counter["T"])/float(totalBases))
		self.appLogger.info("Extracting base frequencies: {}".format(self.ancestralFreq))

	def writeIndelibleControlFile(self):
		"""
		Writes the modified INDELible control file into the appropriate
		directory to be able, afterwards, to run INDELible properly.
		"""
		self.appLogger.debug("Writing new control file")
		self.appLogger.debug("Given INDELible control file: ".format(\
			self.settings.indelibleControlFile))
		lines=None
		with open(self.settings.indelibleControlFile,"r") as f:
			lines=[line.strip() for line in f ]
		newlines=copy.copy(lines)
		newlines.reverse()
		controllines=[]
		modelname=""
		while len(newlines)>0:
			line=newlines.pop()
			if "[NGSPHYPARTITION]" in line:
				self.partition=line.split() # ill get 3 elems + label
			if "[MODEL]" in line:
				modelname=line.split()[1]
		del newlines
		newlines=copy.copy(lines)
		for item in newlines:
			if item.strip().startswith("[NGSPHY"):
				break
			controllines+=[item.strip("\n")]

		statesFreqPresent=-1
		# BASE FREQUENCIES - [statefreq] T C A G
		stateFreqEntry=["[statefreq]\t{0} {1} {2} {3}".format(\
			self.ancestralFreq["T"],\
			self.ancestralFreq["C"],\
			self.ancestralFreq["A"],\
			self.ancestralFreq["G"]\
			)]
		for index in range(0,len(controllines)):
			if controllines[index].strip().startswith("[statefreq]"):
				statesFreqPresent=index
				break
		if (statesFreqPresent < 0):
			controllines+=stateFreqEntry

		with open(self.settings.geneTreeFile) as f:
			geneTreeFile=[line.strip() for line in f if not line.strip() == ""]
		geneTreeFile="".join(geneTreeFile)
		geneTreeFile=geneTreeFile.replace("'","")
		# print(geneTreeFile)
		if geneTreeFile[-1]!=";":
			geneTreeFile+=";"
		controllines+=["{0} {1} {2}".format(\
			"[TREE]",\
			"ngsphytree",\
			geneTreeFile
		)]
		controllines+=["{0} {1} [{2} {3} {4}]".format(\
			"[PARTITIONS]",\
			"ngsphypartition",\
			"ngsphytree",\
			modelname,\
			self.partition[3]
		)]

		controllines+=["[EVOLVE]"]
		controllines+=[" {0} {1} {2}".format(\
			"ngsphypartition",\
			self.evolve[0],\
			self.evolve[1]\
		)]

		# full control file, missing checking settings of output and fastaextension
		fastaoutput="[output] FASTA"
		fastaoutputext="[fastaextension] fasta"
		output=[]; outputext=[]
		settings=False
		for indexControl in range(0, len(controllines)):
			data=controllines[indexControl].strip()
			if data=="[SETTINGS]":
				settings=True
			if data.startswith("[output]"):
				ss=data.split()[1]
				output+=[ss.upper()]
			if data.startswith("[fastaextension]"):
				outputext+=[indexControl]
		if not settings:
			controllines.insert(1,"[SETTINGS]")
		if (not "FASTA" in output) and (len(outputext) ==0):
			controllines.insert(2,"[output] FASTA")
			controllines.insert(3,"[fastaextension] fasta")
		elif (not "FASTA" in output):
			controllines.insert(2,"[output] FASTA")
		elif (len(outputext) ==0):
			controllines.insert(2,"[fastaextension] fasta")
		# check whether there is thr [randomseed] options
		# if exists keep, else use the program one
		randomseed=False
		for indexControl in range(0,len(controllines)):
			data=controllines[indexControl].strip()
			if data.startswith("[randomseed]"): randomseed=True
		if not randomseed:
			controllines.insert(2,"[randomseed] {0}".format(self.settings.seed))
		# write controllines to file
		f=open(self.newIndelibleControlFilePath,"w")
		for item in controllines:
			f.write("{}\n".format(item))
		f.close
		st = os.stat(self.newIndelibleControlFilePath)
		os.chmod(self.newIndelibleControlFilePath, st.st_mode | stat.S_IEXEC | stat.S_IXGRP | stat.S_IXOTH)

	def runIndelible(self):
		"""
		Launches a thread with the INDELible command.
		------------------------------------------------------------------------
		Returns:
		- boolean. Indicates status of the process
		"""
		self.appLogger.debug("Running...")
		self.settings.parser.set("general","numLociPerReplicate",str(1))
		self.settings.parser.set("general", "filtered_replicates", str(1))
		self.settings.parser.set("general", "simphy_data_prefix",self.evolve[1])
		try:
			self.appLogger.info("Waiting for INDELible process to finish. This may take a while...")
			t = multiprocessing.Process(target=self.indelibleLauncher())
			t.start()
			t.join()
			self.appLogger.info("INDELible's run has finished.")
		except RuntimeError as verror:
			return   False, verror
		except Exception as ex:
			return   False, ex
		return True, "INDELible's run has finished."

	def indelibleLauncher(self):
		"""
		Generates a subprocess that handles the calling to INDELible
		"""
		indelibleMessage="INDELible run has finished";proc=""
		lines=[]
		try:
			cwd=os.getcwd()
			self.appLogger.info("Running INDELible")
			proc=""
			self.appLogger.info("Moving to {0}".format(self.settings.alignmentsFolderPath))
			os.chdir(self.settings.alignmentsFolderPath)
			self.appLogger.info("Running INDELible")
			proc = subprocess.check_output([self.settings.programCommand,'control.txt'],stderr=subprocess.STDOUT)
			if (sys.version_info>=(3,0)): proc=str(proc, "utf-8")
			self.appLogger.info("Moving back to working directory: {0}".format(cwd))
			os.chdir(cwd)
			cpuTime = [line.split(":")[1].split()[0] for line in proc.split('\n') if "* Block" in line]
			for item in range(1,len(cpuTime)):
				cpu=cpuTime[(item-1)]
				output="{0}_{1}".format(self.evolve[1],item )
				lines+=[item,cpu,output]
		except subprocess.CalledProcessError as error:
			indelibleMessage="\nINDELible execution error. "+\
			"\n------------------------------------------------------------------------"+\
			"\n{0}".format(error)+\
			"\n------------------------------------------------------------------------"+\
			"{0}".format(error.output)+\
			"\n------------------------------------------------------------------------"+\
			"\nFor more information about this error please run the following commands separately:\n"+\
			"\n\tcd {0}\n\tindelible\n".format(self.settings.alignmentsFolderPath)
			raise RuntimeError(indelibleMessage)
		except Exception as ex:
			exc_type, exc_obj, exc_tb = sys.exc_info()
			fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
			message="{1} | File: {2} - Line:{3}\n\t{0}".format(\
			    ex,exc_type,\
			    fname, exc_tb.tb_lineno\
			)
			raise RuntimeError(message)
		if (self.settings.runningTimes): self.writeRunningInfoIntoFile(lines)

	def writeRunningInfoIntoFile(self, lines):
		"""
		Writes the information about timing into a file.
		"""
		outputFile=os.path.join(
			self.settings.alignmentsFolderPath,\
			"{}.indelible.time".format(self.settings.projectName)
		)
		f=open(outputFile,"w")
		f.write("indexGT,cpuTime,outputFilePrefix\n")
		for item in lines:
			f.write(
				str(item[0])+","+\
				str(item[1])+","+\
				item[2]+"\n"
			)
		f.close()
		self.appLogger.info("File with timings of the INDELible run can be find on: {0}".format(outputFile))
