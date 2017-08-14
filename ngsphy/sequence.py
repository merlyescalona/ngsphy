#!/usr/bin/home/python
import argparse,copy,datetime,dendropy,logging,os,re,sys, multiprocessing, subprocess
import numpy as np
import random as rnd
import settings as sp
from msatools import *
from select import select

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

	newIndelibleAncestralSequence=""
	newIndelibleControlFilePath=""
	geneTreeFile=""

	evolve=[1,"ngsphydata_1"]
	partition=[]

	programCommand=""


	def __init__(self,settings):
		self.appLogger=logging.getLogger('ngsphy')
		self.appLogger.debug('INDELible run')
		self.settings=settings
		self.settings.alignmentsFolderPath=os.path.join(self.settings.basepath,"1")
		self.newIndelibleControlFilePath=os.path.join(\
			self.settings.alignmentsFolderPath,"1","control.txt")
		self.newIndelibleAncestralSequence=os.path.join(\
			self.settings.alignmentsFolderPath,"1","ancestral.fasta")

		if self.settings.inputmode ==1:
			self.programCommand="indelible"
		if self.settings.inputmode in [2,3]:
			self.programCommand="indelible-ngsphy"
		if self.settings.inputmode==3: # INDELible + reference
			self.geneTreeFile=os.path.join(\
				self.settings.alignmentsFolderPath,\
				"ngsphy.tree"\
			)


	def run(self):
		"""
		Process flow for the generation of genome sequences from a gene tree
		and an evolutionary model
		"""
		self.generateFolderStructure()
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
			os.makedirs(os.path.join(self.settings.alignmentsFolderPath,"1"))
			self.appLogger.info("Generating data folder ({0})".format(\
				os.path.join(self.settings.alignmentsFolderPath,"1")))
		except:
			self.appLogger.debug("Data folder exists ({0})".format(\
				os.path.join(self.settings.alignmentsFolderPath,"1")))


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
		self.appLogger.info("Copying reference sequence file to: {}".format(\
			self.newIndelibleAncestralSequence))
		fin=open(self.settings.ancestralSequenceFilePath,"r")
		referenceLines=fin.readlines()
		fin.close()
		if len(referenceLines) > 1:
			fout=open(self.newIndelibleAncestralSequence,"w")
			fout.write(">ngsphypartition\n")
			for index in range(1,len(referenceLines)):
				item=referenceLines[index].strip()
				if not item.startswith(">"):
					fout.write("{}".format(item))
				else:
					break
			if item.startswith(">"):
				self.appLogger.warning("Ancestral sequence file has more than one sequence. Only one is needed. Selecting the first one.")
			fout.write("\n")
			fout.close()
		else:
			status=False
			message="\n\t{0}\n\t{1}\n\t{2}\n\t{3}".format(\
				"Ancestral sequence file given has a problem.",
				"Please verify",\
				"Exiting"
			)
		return status, message


	def writeIndelibleControlFile(self):
		"""
		Writes the modified INDELible control file into the appropriate
		directory to be able, afterwards, to run INDELible properly.
		"""
		self.appLogger.debug("Writing new control file")
		self.appLogger.debug("Given INDELible control file: ".format(\
			self.settings.ngsphyIndelibleControlFilePath))
		f=open(self.settings.ngsphyIndelibleControlFilePath,"r")
		lines=f.readlines()
		f.close()
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

		newlines=copy.copy(lines)
		for item in newlines:
			if item.strip().startswith("[NGSPHY"):
				break
			controllines+=[item.strip("\n")]

		f=open(self.settings.geneTreeFile)
		newicklines=f.readlines()
		f.close()
		geneTreeFile=[ item.strip() for item in newicklines if item.strip()!=""]
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
		fastaoutput="\t[output] FASTA"
		fastaoutputext="\t[fastaextension] fasta"
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
			controllines.insert(2,"  [output] FASTA")
			controllines.insert(3,"  [fastaextension] fasta")
		elif (not "FASTA" in output):
			controllines.insert(2,"  [output] FASTA")
		elif (len(outputext) ==0):
			controllines.insert(2,"  [fastaextension] fasta")

		# write controllines to file
		f=open(self.newIndelibleControlFilePath,"w")
		for item in controllines:
			f.write("{}\n".format(item))
		f.close()

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
			self.appLogger.info("Moving to {0}".format(self.settings.alignmentsFolderPath))
			os.chdir(self.settings.alignmentsFolderPath)
			self.appLogger.info("Running INDELible")
			proc=""
			if (sys.version_info[0:2]<(3,0)):
				proc =subprocess.check_output(self.programCommand,stderr=subprocess.STDOUT)
			elif (sys.version_info>=(3,0)):
				 proc = str(subprocess.check_output(self.programCommand,stderr=subprocess.STDOUT),'utf-8')
			cpuTime = [line.split(":")[1].split()[0] for line in proc.split('\n') if "* Block" in line]
			for item in range(1,len(cpuTime)):
				cpu=cpuTime[(item-1)]
				output="{0}_{1}".format(self.evolve[1],item )
				lines+=[item,cpu,output]
		except OSError as error:
			indelibleMessage="Output folder does not exist. Please verify. Exiting."
			self.appLogger.error(indelibleMessage)
			raise RuntimeError(indelibleMessage)
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
		if (self.runnning_times): self.writeRunningInfoIntoFile(lines)

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
