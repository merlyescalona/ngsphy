#!/usr/bin/home/python
import argparse,copy,datetime,dendropy,logging,os,re,sys, multiprocessing, subprocess
import numpy as np
import random as rnd
import settings as sp
from msatools import *
from select import select

class SequenceGenerator:
	appLoger=None
	settings=None

	newIndelibleReferenceSequence=""
	newIndelibleFilePath=""
	output=""

	evolve=[1,"ngsphydata_1"]
	partition=[]

	numGTs=0
	numGTDigits=0

	def __init__(self,settings):
		self.appLogger=logging.getLogger('ngsphy')
		self.appLogger.debug('INDELible run')
		self.settings=settings
		self.output=os.path.join(self.settings.path,self.settings.projectName,"1")
		self.newIndelibleFilePath=os.path.join(self.settings.path,self.settings.projectName,"1","control.txt")
		self.newIndelibleReferenceSequence=os.path.join(self.settings.path,self.settings.projectName,"1","reference.fasta")


	def run(self):
		self.generateFolderStructure()
		# check naming of the leaves
		# adding extra information on presence of outgroup
		self.addOutgroupInfoToSettings()
		# check ploidy and tree correspondance
		ploidyStatus,ploidyMessage=self.checkPloidyTreeRelation()
		if not ploidyStatus: return False, ploidyMessage
		self.appLogger.debug("Out tree-ploidy relation")
		self.writeIndelibleControlFile()
		statusRefSeq,messageRefSeq=self.copyReferenceSequenceToOutputFolder()

		if not statusRefSeq: return False,messageRefSeq
		runStatus,runMessage=self.runIndelible()
		if not runStatus: return False,runMessage
		return True, "Run finished"

	def generateFolderStructure(self):
		self.appLogger.info("Creating folder structure for INDELible run")
		try:
			os.makedirs(os.path.join(self.settings.path,self.settings.projectName))
			self.appLogger.info("Generating project folder ({0})".format(os.path.join(self.settings.path,self.settings.projectName)))
		except:
			self.appLogger.debug("Project folder exists ({0})".format(os.path.join(self.settings.path,self.settings.projectName)))
		# create data folder
		try:
			os.makedirs(os.path.join(self.settings.path,self.settings.projectName,"1"))
			self.appLogger.info("Generating data folder ({0})".format(os.path.join(self.settings.path,self.settings.projectName,"1")))
		except:
			self.appLogger.debug("Data folder exists ({0})".format(os.path.join(self.settings.path,self.settings.projectName,"1")))


	def addOutgroupInfoToSettings(self):
		self.appLogger.debug("Outgroup")
		tree=dendropy.Tree.get(path=self.settings.newickFilePath, schema="newick",preserve_underscores=True)
		leaves=[ node.taxon.label for node in tree.leaf_node_iter()]
		item=""
		for item in leaves:
			if item == "0_0_0": break
		if item=="0_0_0":
			self.settings.outgroup=True
			self.settings.parser.set("general","outgroup","on")
		else:
			self.settings.parser.set("general","outgroup","off")

	def checkPloidyTreeRelation(self):
		self.appLogger.debug("Checking ploidy - num tips relation")
		messageCorrect="Ploidy and number of gene copies per gene family match properly."
		messageWrong="INDELible control file - Something's wrong!\n\t"
		tree=dendropy.Tree.get(path=self.settings.newickFilePath, schema="newick",preserve_underscores=True)
		leaves=[]
		if self.settings.outgroup:
			leaves=[ node.taxon.label for node in tree.leaf_node_iter() if not node.taxon.label =="0_0_0"]
		else:
			leaves=[ node.taxon.label for node in tree.leaf_node_iter()]
		leavesSplit=[ item.split("_") for item in leaves]
		leavesDict=dict()
		for tip in leavesSplit:
			geneFamily="_".join(tip[0:2])
			try:
				val=leavesDict[geneFamily]
				leavesDict[geneFamily]+=1
			except:
				leavesDict[geneFamily]=1
		for item in leavesDict:
			self.appLogger.debug("{0} - {1}".format(item,leavesDict[item]))
			if not leavesDict[item] % self.settings.ploidy == 0:
				messageWrong+="{0}\n\t{1}".format(\
				"The number of gene copies within one of the gene families does not match the ploidy selected for this run.",\
				"Please verify. Exiting."\
				)
				return False,messageWrong
		return True, messageCorrect

	def copyReferenceSequenceToOutputFolder(self):
		# making sure there's only one sequence, and only one sequence will be written to the
		# reference.fasta file
		# that sequence will be the first from the file if there are more than 1 sequence
		status=True; message=""
		self.appLogger.info("Copying reference sequence file to: {}".format(self.newIndelibleReferenceSequence))
		fin=open(self.settings.referenceSequenceFilePath,"r")
		referenceLines=fin.readlines()
		fin.close()
		if len(referenceLines) > 1:
			fout=open(self.newIndelibleReferenceSequence,"w")
			fout.write(">ngsphypartition\n")
			for index in range(1,len(referenceLines)):
				item=referenceLines[index].strip()
				if not item.startswith(">"):
					fout.write("{}".format(item))
				else:
					break
			if item.startswith(">"):
				self.appLogger.warning("Reference sequence file has more than one sequence. Only one is needed. Selecting the first one.")
			fout.write("\n")
			fout.close()
		else:
			status=False
			message="\n\t{0}\n\t{1}\n\t{2}\n\t{3}".format(\
				"Reference sequence file given has a problem.",
				"Please verify",\
				"Exiting"
			)
		return status, message


	def writeIndelibleControlFile(self):
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

		f=open(self.settings.newickFilePath)
		newicklines=f.readlines()
		f.close()
		newicktree=[ item.strip() for item in newicklines if item.strip()!=""]
		newicktree="".join(newicktree)
		newicktree=newicktree.replace("'","")
		# print(newicktree)
		if newicktree[-1]!=";":
			newicktree+=";"
		controllines+=["{0} {1} {2}".format(\
			"[TREE]",\
			"ngsphytree",\
			newicktree
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
		f=open(self.newIndelibleFilePath,"w")
		for item in controllines:
			f.write("{}\n".format(item))
		f.close()

	def runIndelible(self):
		self.appLogger.debug("Running...")
		self.settings.parser.set("general","numLociPerSpeciesTree","1")
		self.settings.parser.set("general", "filtered_ST", "1")
		self.settings.parser.set("general", "data_prefix",self.evolve[1])
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
		indelibleMessage="INDELible run has finished";proc=""
		lines=[]
		try:
			self.appLogger.info("Moving to {0}".format(self.output))
			os.chdir(self.output)
			self.appLogger.info("Running INDELible")
			proc=""
			if (sys.version_info[0:2]<(3,0)):
				proc =subprocess.check_output("indelible",stderr=subprocess.STDOUT)
			elif (sys.version_info>=(3,0)):
				 proc = str(subprocess.check_output("indelible",stderr=subprocess.STDOUT),'utf-8')
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
			"\n\tcd {0}\n\tindelible\n".format(self.output)
			raise RuntimeError(indelibleMessage)
		self.printRunningInfo(lines)

	def printRunningInfo(self, lines):
		outputFile="{0}/{1}.indelible.info".format(
			self.settings.outputFolderPath,\
			self.settings.projectName
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
