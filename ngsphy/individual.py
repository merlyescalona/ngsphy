#!/usr/bin/env python
# -*- coding: utf-8 -*-
import copy,dendropy,glob,logging,os,sys, sqlite3
import numpy as np
import random as rnd
import settings as sp
from msatools import *
from select import select
from ngsphyexceptions import *
################################################################################
class IndividualAssignment:
	"""
	Class for the generation of individuals
	----------------------------------------------------------------------------
	Attributes:
	- appLogger: logger to store status of the process flow
	- settings: Settings object withh all the program parameters
	relation to the generated individuals
	- numReplicates: number of species trees.
	- numReplicateDigits: number of digits needed to represent numReplicates.
	- numLociPerReplicate: number of loci per species tree.
	- numLociPeReplicateDigits: number of digits needed to represent.
	- numIndividualsPerReplicate: number of individuals per species tree.
	- filteredST: identifier of the species trees that will be used.
	"""
	appLogger=None
	settings=None
	output=""
	numReplicates=0
	numReplicateDigits=0
	numLociPerReplicate=[]
	numLociPerReplicateDigits=[]
	numIndividualsPerReplicate=[]
	filteredReplicates=[]
	indelsPresence=False

	def __init__(self, settings):
		self.appLogger=logging.getLogger("ngsphy")
		self.appLogger.info("IndividualGenerator: Run started")
		self.settings=settings

	def checkArgs(self):
		"""
		Validation of the values needed for this process.
		- Input according origin mode
		"""
		self.appLogger.info("Checking arguments folder...")
		self.generateFolderStructure()
		message=""
		# need to know whether i'm working with simphy or indelible
		self.numReplicates=self.settings.numReplicates
		self.filteredReplicates=range(1,self.numReplicates+1)
		self.numReplicateDigits=len(str(self.numReplicates))
		self.numIndividualsPerReplicate=[0]*len(self.filteredReplicates)
		self.numLociPerReplicate=[0]*self.numReplicates
		self.numLociPerReplicateDigits=[0]*self.numReplicates
		########################################################################
		if self.settings.inputmode < 4 :
			# this is like this, because for this to work it is necessary that
			# sequences had been generated, and it is in the SequenceGenerator
			# that the numLociPerReplicate option is set up.
			status,message=self.checkPloidyTreeRelation()
			if not status: return status,message
			if not self.settings.parser.has_option("general","numlociperreplicate"):
				self=numLociPerReplicate=[1]
			else:
				self.numLociPerReplicate=[\
					self.settings.parser.getint("general","numlociperreplicate")]
			self.numLociPerReplicateDigits=[len(str(self.numLociPerReplicate[0]))]
		elif self.settings.inputmode==4:
			####################################################################
			# checking the replicate that are going to be used
			# checking if I'll used the filtered in case there's a possibility
			# that one or many sts do not match the ploidy and number of gene copies
			self.command = os.path.join(\
				self.settings.basepath,\
				"{0}.command".format(self.settings.projectName))
			self.params = os.path.join(\
				self.settings.basepath,\
				"{0}.params".format(self.settings.projectName))
			self.db = os.path.join(\
				self.settings.basepath,\
				"{0}.db".format(self.settings.projectName))

			self.numLociPerReplicate=self.getSimPhyNumLociPerSpeciesTree()
			self.numLociPerReplicateDigits=[len(str(a))for a in self.numLociPerReplicate]
			# check ploidy matches given data
			if self.settings.simphyFilter:
				self.appLogger.info("Filtering replicates...")
				self.filteredReplicates=self.filterReplicatesMatchingIndPerSpeciesAndPloidy(self.settings.ploidy)
				self.numLociPerReplicate=self.getSimPhyNumLociPerSpeciesTreeFiltered(self.filteredReplicates)
				self.numLociPerReplicateDigits=[len(str(a))for a in self.numLociPerReplicate]
			self.appLogger.info("Checking filtered replicates...")
			status,message=self.checkPloidySimPhyData()
			if not status: return status,message
 			# check that the species tree replicate folder have the correct data
			self.appLogger.info("Checking data within replicates...")
			gtperstOK,message=self.checkDataWithinReplicates()
			if (not gtperstOK): return gtperstOK,message
			self.settings.parser.set(\
				"general",\
				"numlociperreplicate",\
				",".join([str(a) for a in self.numLociPerReplicate]))
			self.printSimPhyConfiguration()
		else:
			return False, "{0}\n\t{1}\n\t{2}".format(\
				"Individual Assignment process."
				"Something is wrong with the input.",\
				"Please verify. Exiting."\
			)
		self.settings.parser.set("general","filtered_replicates",",".join([str(a) for a in self.filteredReplicates]))
		return True, message

	def generateFolderStructure(self):
		"""
		Generation of general folder structure needed for the individual generation
		"""
		self.appLogger.info("Creating folder structure for individual generation...")
		try:
			self.appLogger.info("Generating output folder {}".format(self.settings.individualsFolderPath))
			os.makedirs(self.settings.individualsFolderPath)
		except:
			self.appLogger.debug("Output folder exists ({0})".format(self.settings.outputFolderPath))
		try:
			self.appLogger.info("Generating output folder {}".format(self.settings.tablesFolderPath))
			os.makedirs(self.settings.tablesFolderPath)
		except:
			self.appLogger.debug("Output folder exists ({0})".format(self.settings.tablesFolderPath))


	def printSimPhyConfiguration(self):
		"""
		Print the configuration of the simphy project
		"""
		self.appLogger.debug(\
			"\n\tConfiguration...\n\t{0}\t{1}\n\t{2}\t{3}\n\t{4}\t{5}\n\t{6}\t{7}\n\t{8}\t{9}".format(\
				"SimPhy path:",\
				self.settings.basepath,\
				"SimPhy project name:",\
				self.settings.projectName,\
				"Output folder:",\
				self.settings.outputFolderPath,\
				"Dataset prefix(es) (INDELible):",\
				self.settings.simphyDataPrefix,\
				"Number of species trees replicates/folders:",\
				self.numReplicates\
			)
		)

	def getSimPhyNumLociPerSpeciesTree(self):
		"""
		Retrieves information of number of loci per species from the SimPhy
		database.
		------------------------------------------------------------------------
		Returns: a list with the number of loci per species tree replicate
		"""
		query="select N_Loci from Species_Trees"
		con = sqlite3.connect(self.db)
		res=con.execute(query).fetchall()
		con.close()
		res=[item for sublist in res for item in sublist]
		return res

	def getSimPhyNumLociPerSpeciesTreeFiltered(self, indices):
		"""
		Retrieves information of number of loci per species from the SimPhy
		database.
		------------------------------------------------------------------------
		Returns: a list with the number of loci per species tree replicate
		"""
		query="select N_Loci from Species_Trees where SID in ({0})".format(",".join([str(i) for i in indices]))
		con = sqlite3.connect(self.db)
		res=con.execute(query).fetchall()
		con.close()
		res=[item for sublist in res for item in sublist]
		return res

	def filterReplicatesMatchingIndPerSpeciesAndPloidy(self, ploidy):
		"""
		Identifies and filters the species tree replicates that  SimPhy
		database.
		------------------------------------------------------------------------
		Returns: a list with the number of loci per species tree replicate
		"""
		query="select SID from Species_Trees"
		if not (ploidy == 1):
			query="select SID from Species_Trees WHERE Ind_per_sp % {0} = 0".format(ploidy)
		con = sqlite3.connect(self.db)
		res=con.execute(query).fetchall()
		con.close()
		res=[item for sublist in res for item in sublist]
		return res

	def checkPloidySimPhyData(self):
		"""
		Checks whether the ploidy given matches the information in the
		SimPhy database
		------------------------------------------------------------------------
		Returns: a list with the number of loci per species tree replicate
		"""
		if not self.filteredReplicates == []:
			query="select Ind_per_sp from Species_Trees WHERE SID in ({})".format(",".join([str(i) for i in self.filteredReplicates]))
			con = sqlite3.connect(self.db)
			res=con.execute(query).fetchall()
			con.close()
			res=[item for sublist in res for item in sublist]
			status=True; message=""
			for item in res:
				if not item % self.settings.ploidy == 0:
					status=False
					message="{0}\n\t{1}\n\t{2}\n\t{3}".format(\
						"There has been a problem with the ploidy.",\
						"There is at least one species tree replicate, which number ",\
						"of individuals does not match the ploidy specified",\
						"Please verify. Exiting."
					)
		else:
			status=False
			message="{0}\n\t{1}\n\t{2}".format(\
				"There has been a problem with the ploidy.",\
				"No replicate satisfies the conditions.",\
				"Please verify. Exiting."
			)
		return status, message

	def checkPloidyTreeRelation(self):
		"""
		Checks whether the ploidy given matches the information given within
		the given tree.
		"""
		self.appLogger.debug("Checking ploidy - num tips relation")
		status=True
		message="Ploidy and number of gene copies per gene family match properly."
		tree=dendropy.Tree.get(\
			path=self.settings.geneTreeFile,\
		 	schema="newick",\
			preserve_underscores=True)
		leaves=[node.taxon.label for node in tree.leaf_node_iter()]
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
				status=False
				message+="\n\t{0}\n\t{1}\n\t{2}".format(\
				"INDELible control file - Something's wrong!",\
				"The number of gene copies within one of the gene families does not match the ploidy selected for this run.",\
				"Please verify. Exiting."\
				)

		return status, message

	def checkDataWithinReplicates(self):
		""""
		Checks the data files within the  replicates.
		Existence of fasta files.
		"""
		self.appLogger.info("Checking replicate data: ReplicateID - currentReplicate/numberOfWorkingReplicates [numberOfFiles]")
		# for indexREP in self.filteredReplicates:
		for index in range(0,len(self.filteredReplicates)):
			curReplicatePath=os.path.join(\
				self.settings.basepath,\
				"{0:0{1}d}".format(\
					self.filteredReplicates[index],\
					self.numReplicateDigits\
				)\
			)
			numFastaFiles=0;numGeneTrees=0
			# check composition of the current indexREP folder
			numFastaFiles=len(glob.glob("{0}/{1}_*_TRUE*".format(curReplicatePath, self.settings.simphyDataPrefix)))
			numGeneTrees=len(glob.glob("{0}/g_trees*.trees".format(curReplicatePath, self.settings.simphyDataPrefix)))
			self.numLociPerReplicate[index]=numFastaFiles
			self.appLogger.info("ReplicateID {0} - {1}/{2} [{3}]".format(\
				self.filteredReplicates[index],\
				index+1,\
			 	len(self.filteredReplicates),\
				numFastaFiles))
			self.numLociPerReplicateDigits[index]=len(str(numFastaFiles))
			if (numFastaFiles<1):
				# Do not have fasta files from the given replicate to work, I'll skip it.
				self.appLogger.warning(\
					"Replicate {0}({1}): It is not possible to do the mating for this replicate".format(self.filteredReplicates[index], curReplicatePath))
				return False, "Please verify. Exiting."
			if (numGeneTrees<1):
				return False,"Trying to mate sequences, but there are no gene tree files to back that up. Please, finish the SimPhy run and try again afterwards."
			self.appLogger.debug("numGeneTrees ({}) == numFastaFiles ({})? ".format(numGeneTrees,numFastaFiles))
			if not numGeneTrees == numFastaFiles:
				return False, "There are no sequences o there is a missmatch between the number of trees and the number of sequences in the folder."
		return True,"Got number of gene trees per species trees"

	def iteratingOverReplicates(self):
		"""
		Iterates over the species tree replicates and generates individuals
		according to the selected ploidy
		"""
		self.appLogger.debug("Ploidy: {0}".format(self.settings.ploidy))
		if (self.settings.ploidy==1):
			self.iterationHaploid()
		else:
			self.iterationPolyploid()
		self.settings.parser.set("general","numindividualsperreplicate",\
			",".join([str(a) for a in self.numIndividualsPerReplicate]))

	def iterationPolyploid(self):
		"""
		Iterates over the species tree replicates.
		Within each species tree, iterates over the gene trees, generates
		the "mating table" as well as the file with the individuals's sequences.
		"""
		self.appLogger.info("Generating individuals: replicateID - currentReplicate/numberOfWorkingReplicates [numLoci] (path)...")
		for index in range(0,len(self.filteredReplicates)):
			curReplicatePath=os.path.join(\
				self.settings.individualsFolderPath,\
				"REPLICATE_{0:0{1}d}".format(\
					self.filteredReplicates[index],\
					self.numReplicateDigits\
				)\
			)
			self.appLogger.info("ReplicateID {0} - {1}/{2} [{3}] ({4}) ".format(\
				self.filteredReplicates[index],\
				index+1,\
				len(self.filteredReplicates),\
				self.numLociPerReplicate[index],\
				curReplicatePath\
			))
			# iterating over the number of gts per st
			self.appLogger.debug("Generating individual table...")
			matingTable=self.generateMatingTable(self.filteredReplicates[index])
			self.numIndividualsPerReplicate[index]=len(matingTable)
			self.writeMatingTable(self.filteredReplicates[index],matingTable)
			for indexLOC in range(1,self.numLociPerReplicate[index]+1):
				# parsingMSA file
				if(self.settings.inputmode<4):
					fastapath=os.path.join(
						self.settings.basepath,\
						"REPLICATE_{0:0{1}d}".format(self.filteredReplicates[index],self.numReplicateDigits),\
						"{0}_{1:0{2}d}_TRUE.fasta".format(\
							self.settings.simphyDataPrefix,\
							indexLOC,\
							self.numLociPerReplicateDigits[index]\
						)\
					)
				else:
					fastapath=os.path.join(
						self.settings.basepath,\
						"{0:0{1}d}".format(self.filteredReplicates[index],self.numReplicateDigits),\
						"{0}_{1:0{2}d}_TRUE.fasta".format(\
							self.settings.simphyDataPrefix,\
							indexLOC,\
							self.numLociPerReplicateDigits[index]\
						)\
					)
				seqDict=parseMSAFile(fastapath)
				outputFolder=os.path.join(
					self.settings.individualsFolderPath,\
					"REPLICATE_{0:0{1}d}".format(\
						self.filteredReplicates[index],\
						self.numReplicateDigits),\
					"LOCUS_{0:0{1}d}".format(\
						indexLOC,\
						self.numLociPerReplicateDigits[index])\
				)
				try:
					self.appLogger.debug("Output folder: {0}".format(outputFolder))
					os.makedirs(outputFolder)
				except OSError as err:
					self.appLogger.warning("OS error: {0}".format(err))
					self.appLogger.debug("Folder {0} exists.".format(outputFolder))
				# generating and writing mating table
				self.mate(self.filteredReplicates[index],indexLOC,matingTable,seqDict)

	def iterationHaploid(self):
		"""
		Iterates over the species tree replicates.
		Within each species tree, iterates over the gene trees, generates
		the "relation table" as well as the file with the individuals's sequences.
		"""
		self.appLogger.info("Generating individuals: replicateID - currentReplicate/numberOfWorkingReplicates [numLoci] (path)...")
		# for indexREP in self.filteredReplicates:
		for index in range(0, len(self.filteredReplicates)):
			curReplicatePath=os.path.join(\
				self.settings.individualsFolderPath,
				"REPLICATE_{0:0{1}d}".format(\
					self.filteredReplicates[index],\
				 	self.numReplicateDigits)
			)
			self.appLogger.info("ReplicateID {0} - \t{1}/{2} [{3}] ({4}) ".format(\
				self.filteredReplicates[index],\
				index+1,\
				len(self.filteredReplicates),\
				self.numLociPerReplicate[index],\
				curReplicatePath\
			))
			# iterating over the number of gts per st
			individualTable=None
			self.appLogger.info("Generating individual table")
			individualTable=self.generateIndividualTable(self.filteredReplicates[index])
			self.numIndividualsPerReplicate[index]=len(individualTable)
			self.writeIndividualTable(self.filteredReplicates[index],individualTable)
			for indexLOC in range(1,self.numLociPerReplicate[index]+1):
				# parsingMSA file
				self.appLogger.debug("Using REP={0}, LOC={1}".format(self.filteredReplicates[index],indexLOC))
				if(self.settings.inputmode<4):
					fastapath=os.path.join(\
						self.settings.basepath,\
						"REPLICATE_{0:0{1}d}".format(self.filteredReplicates[index],self.numReplicateDigits),\
						"{0}_{1:0{2}d}.fasta".format(\
							self.settings.simphyDataPrefix,\
							indexLOC,\
							self.numLociPerReplicateDigits[index]
						)\
					)
				else:
					fastapath=os.path.join(\
						self.settings.basepath,\
						"{0:0{1}d}".format(self.filteredReplicates[index],self.numReplicateDigits),\
						"{0}_{1:0{2}d}.fasta".format(\
							self.settings.simphyDataPrefix,\
							indexLOC,\
							self.numLociPerReplicateDigits[index]
						)\
					)
				seqDict=parseMSAFileWithDescriptions(fastapath)
				################################################################
				# if self.settings.inputmode==3:
				# 	# need to modify the sequence of the anchor-root
				# 	seqs=None
				# 	with open(self.settings.ancestralSequenceFilePath, 'r') as f:
				# 		seqs=[line.strip() for line in f if not line.startswith("#")]
				# 	seqDict[self.settings.anchorTipLabel]=seqs[0]
				################################################################
				outputFolder=os.path.join(\
					self.settings.individualsFolderPath,\
					"REPLICATE_{0:0{1}d}".format(\
						self.filteredReplicates[index], self.numReplicateDigits),\
					"LOCUS_{0:0{1}d}".format(\
						indexLOC,\
						self.numLociPerReplicateDigits[index])\
				)
				try:
					self.appLogger.debug("Output folder: {0}".format(outputFolder))
					os.makedirs(outputFolder)
				except OSError as err:
					self.appLogger.warning("OS error: {0}".format(err))
					self.appLogger.debug("Folder {0} exists.".format(outputFolder))
				# generating and writing mating table
				self.generateIndividuals(self.filteredReplicates[index],indexLOC,individualTable,seqDict)

	def generateIndividualTable(self,indexREP):
		"""
		Generates the table that stores the relationship between the original
		sequnce description and the final individual identifier.
		------------------------------------------------------------------------
		Parameters:
		- indexREP: identifier of the Species tree replicate
		Returns:
		- A table/matrix/list of lists with the relation between the orignal
		sequence and the corresponding individual identifier.
		"""
		# get first gt of the st to get the descriptions
		indexLOC=1;
		index=self.filteredReplicates.index(indexREP)
		self.appLogger.info("ReplicateID {0} - {1}/{2} [Locus {3}]".format(\
			indexREP,
			index,\
			len(self.filteredReplicates),\
			indexLOC\
			))

		if(self.settings.inputmode<4):
			fastapath=os.path.join(\
				self.settings.basepath,\
				"REPLICATE_{0:0{1}d}".format(\
					self.filteredReplicates[index],\
					self.numReplicateDigits),\
				"{0}_{1:0{2}d}_TRUE.fasta".format(\
					self.settings.simphyDataPrefix,\
					indexLOC,\
					self.numLociPerReplicateDigits[index]\
				)\
			)
		else:
			fastapath=os.path.join(\
				self.settings.basepath,\
				"{0:0{1}d}".format(\
					self.filteredReplicates[index],\
					self.numReplicateDigits),\
				"{0}_{1:0{2}d}_TRUE.fasta".format(\
					self.settings.simphyDataPrefix,\
					indexLOC,\
					self.numLociPerReplicateDigits[index]\
				)\
			)

		descriptions=parseMSAFileWithDescriptions(fastapath).keys()
		descriptions.sort()
		table=[(item,descriptions[item]) for item in range(0,len(descriptions))]
		return table

	def generateIndividuals(self,indexREP,indexLOC,individualTable,seqDict):
		"""
		Once the table is generated and the corresponding file with the sequences
		of the specific gene tree is parsed, it is possible to generate the file
		with the sequences that correspond to an individual.
		------------------------------------------------------------------------
		Parameters:
		- indexREP: identifier of the Species tree replicate
		- indexLOC: identifier of the gene tree
		- individualTable: table with the relation between sequences and individuals.
		- seqDict: dictionary with the sequences parsed from the gene tree fasta file.

		Generates:
		- A set of files, as many as individuals were described in the table.
		"""
		index=self.filteredReplicates.index(indexREP)
		self.appLogger.debug(\
			"{0}/{1} - {2}".format(\
				indexLOC,\
				self.numLociPerReplicate[index],\
				self.numLociPerReplicateDigits[index]\
			)\
		)
		outputFolder=os.path.join(\
			self.settings.individualsFolderPath,\
			"REPLICATE_{0:0{1}d}".format(\
				self.filteredReplicates[index],\
				self.numReplicateDigits),\
			"LOCUS_{0:0{1}d}".format(\
				indexLOC,\
				self.numLociPerReplicateDigits[index])
		)
		for currentInd in range(0,len(individualTable)):
			# Extracting info from the dictionary
			indID=str(individualTable[currentInd][0])
			description=str(individualTable[currentInd][1])
			# Organizing strings
			# if description=="6_0_0":
				# print indID,outputFolder
			seq=seqDict[description]
			des=">{0}:{1:0{2}d}:{3:0{4}d}:{5}:{6}:{7}".format(
				self.settings.projectName,\
				self.filteredReplicates[index],\
				self.numReplicateDigits,\
				indexLOC,\
				self.numLociPerReplicateDigits[index],\
				self.settings.simphyDataPrefix,\
				indID,\
				description[1:len(description)]\
			)
			indFilename=os.path.join(\
				outputFolder,\
				"{0}_{1:0{2}d}_{3:0{4}d}_{5}_{6}.fasta".format(\
					self.settings.projectName,\
					self.filteredReplicates[index],\
					self.numReplicateDigits,\
					indexLOC,\
					self.numLociPerReplicateDigits[index],\
					self.settings.simphyDataPrefix,\
					indID\
				)\
			)
			indFile=open(indFilename, "w")
			indFile.write("{0}\n{1}\n".format(des,seq))
			indFile.close()

	def writeIndividualTable(self,indexREP,individualTable):
		"""
		Writes into a file the table (individualTable) with the relation between
		sequences and individuals identifiers for the specific indexREP species
		tree replicate
		------------------------------------------------------------------------
		Parameters:
		- indexREP: identifier of the Species tree replicate
		- individualTable: table with the relation between sequences and individuals.
		Generates:
		- File with the table for the indexREP species tree replicate
		"""
		# indexREP,SP,ind-tip1,ind-tip2
		index=self.filteredReplicates.index(indexREP)
		self.appLogger.debug("Writing indexes into file...")
		indexFilename=os.path.join(\
			self.settings.tablesFolderPath,\
			"{0}.{1:0{2}d}.individuals.csv".format(\
				self.settings.projectName,\
				self.filteredReplicates[index],\
				self.numReplicateDigits\
			)\
		)
		if not os.path.isfile(indexFilename):
			indexFile=open(indexFilename,"w")
			indexFile.write("repID,indID,spID,locID,geneID\n")
			indexFile.close()
		indexFile=open(indexFilename,"a")
		for indexRow in range(0,len(individualTable)):
			indID=individualTable[indexRow][0]
			seqDescription=individualTable[indexRow][1]
			speciesID=seqDescription.strip().split("_")[0]
			locusID=seqDescription.strip().split("_")[1]
			geneID=seqDescription.strip().split("_")[2]
			indexFile.write("{0:0{1}d},{2},{3},{4},{5}\n".format(\
				self.filteredReplicates[index],\
				self.numReplicateDigits,\
				indID,\
				speciesID,\
				locusID,\
				geneID))
		indexFile.close()


	# def generateMatingTableFromSequenceDescription(self,indexREP):
	def generateMatingTable(self,indexREP):
		"""
		Generates the "mating" table which stores the relation between the
		sequences that generates an individual.
		This table is generated from the sequence file of the first gene tree
		of the specific indexREP
		------------------------------------------------------------------------
		Parameters:
		- indexREP: identifier of the Species tree replicate
		Returns:
		- the mating table
		"""
		index=self.filteredReplicates.index(indexREP)
		self.appLogger.debug("Mating table")
		if(self.settings.inputmode<4):
			filename=os.path.join(\
				self.settings.basepath,\
				"REPLICATE_{0:0{1}d}".format(\
					self.filteredReplicates[index],\
					self.numReplicateDigits),\
				"{0}_{1:0{2}d}.fasta".format(\
					self.settings.simphyDataPrefix,\
					1,\
					self.numLociPerReplicateDigits[index]))
		else:
			filename=os.path.join(\
				self.settings.basepath,\
				"{0:0{1}d}".format(\
					self.filteredReplicates[index],\
					self.numReplicateDigits),\
				"{0}_{1:0{2}d}.fasta".format(\
					self.settings.simphyDataPrefix,\
					1,\
					self.numLociPerReplicateDigits[index]))
		self.appLogger.debug("Reading file: {0}".format(filename))
		leaves=None
		with open(filename,"r") as f:
			leaves=[ line.strip()[1:] for line in f if line.startswith(">")]
		leavesSplit=[ item.split("_") for item in leaves if not item == self.settings.anchorTipLabel]
		leavesDict=dict()
		for tip in leavesSplit:
			geneFamily="_".join(tip[0:2])
			try:
				val=leavesDict[geneFamily]
				leavesDict[geneFamily]+=1
			except:
				leavesDict[geneFamily]=1
		# Till here i have information about the number of tips per gene family
		numLeaves=0
		for item in leavesDict:
			numLeaves+=leavesDict[item]
		mates=[]
		for geneFamily in leavesDict:
			if leavesDict[geneFamily]==1:
				sp=geneFamily.split("_")[0]
				lt=geneFamily.split("_")[1]
				pair=(self.filteredReplicates[index],sp,lt,0,0)
				mates+=[pair]
				self.appLogger.debug("Pair generated: {0}".format(pair))
			else:
				t=range(0,leavesDict[geneFamily])
				while not t==[]:
					p1=0;p2=0
					try:
						p1=t.pop(rnd.sample(range(0,len(t)),1)[0])
						p2=t.pop(rnd.sample(range(0,len(t)),1)[0])
					except Exception as e:
						break
					sp=geneFamily.split("_")[0]
					lt=geneFamily.split("_")[1]
					pair=(self.filteredReplicates[index],sp,lt,p1,p2)
					mates+=[pair]
					self.appLogger.debug("Pair generated: {0}".format(pair))
		return mates

	def writeMatingTable(self,indexREP,matingTable):
		"""
		Writes into a file the table (individualTable) with the relation between
		sequences and individuals identifiers for the specific indexREP species
		tree replicate
		------------------------------------------------------------------------
		Parameters:
		- indexREP: identifier of the Species tree replicate
		- matingTable: table with the relation between sequences and the generated
		individuals
		Generates:
		- File with the table for the indexREP species tree replicate
		"""
		index=self.filteredReplicates.index(indexREP)
		self.appLogger.debug("Writing indexes into file...")
		indexFilename=os.path.join(\
			self.settings.tablesFolderPath,\
			"{0}.{1:0{2}d}.individuals.csv".format(\
				self.settings.projectName,\
				self.filteredReplicates[index],\
				self.numReplicateDigits\
			)\
		)
		self.appLogger.debug(indexFilename)
		if not os.path.isfile(indexFilename):
			indexFile=open(indexFilename,"w")
			indexFile.write("repID,indID,spID,locID,mateID1,mateID2\n")
			indexFile.close()
		indexFile=open(indexFilename,"a")
		for indexRow in range(0,len(matingTable)):
			indID=indexRow
			speciesID=matingTable[indexRow][1]
			locusID=matingTable[indexRow][2]
			mateID1=matingTable[indexRow][3]
			mateID2=matingTable[indexRow][4]
			indexFile.write("{0:0{1}d},{2},{3},{4},{5},{6}\n".format(\
				self.filteredReplicates[index],\
				self.numReplicateDigits,\
				indID,\
				speciesID,\
				locusID,\
				mateID1,\
				mateID2))
		indexFile.close()

	def mate(self,indexREP,indexLOC,matingTable,sequenceDictionary):
		"""
		Once the table is generated and the corresponding file with the sequences
		of the specific gene tree is parsed, it is possible to generate the file
		with the sequences that correspond to an individual.
		------------------------------------------------------------------------
		Parameters:
		- indexREP: identifier of the Species tree replicate
		- indexLOC: identifier of the gene tree
		- matingTable: table with the relation between sequences and individuals.
		- seqDict: dictionary with the sequences parsed from the gene tree fasta file.
		Generates:
		- A set of files, as many as individuals were described in the table.
		"""
		index=self.filteredReplicates.index(indexREP)
		# Proper mating
		seqDict=copy.deepcopy(sequenceDictionary)
		species=seqDict.keys()
		numSeqs=0
		for key in species:
			subspecies=seqDict[key].keys()
			for subkey in subspecies:
				numSeqs+=len(seqDict[key][subkey])
		numInds=np.trunc(numSeqs/2)+1
		self.appLogger.debug("ReplicateID {0} -  {1}/{2} |  Locus {3:0{5}d}/{4:0{5}d}. Writing {6} individuals from {7} sequecens.".format(\
			self.filteredReplicates[index],\
			index+1,\
			len(self.filteredReplicates),\
			indexLOC,\
			self.numLociPerReplicate[index],\
			self.numLociPerReplicateDigits[index],\
			numSeqs,\
			numInds\
			)\
		)
		outputFolder=os.path.join(\
			self.settings.individualsFolderPath,\
			"REPLICATE_{0:0{1}d}".format(\
				self.filteredReplicates[index],\
				self.numReplicateDigits),\
			"LOCUS_{0:0{1}d}".format(\
				indexLOC,\
				self.numLociPerReplicateDigits[index])\
		)
		seq1="";des1="";seq2="";des2="";
		for currentInd in range(0,len(matingTable)):
			# Extracting info from the dictionary
			st=str(matingTable[currentInd][0])
			sp=str(matingTable[currentInd][1])
			lt=str(matingTable[currentInd][2])
			pair1=str(matingTable[currentInd][3])
			pair2=str(matingTable[currentInd][4])
			tag="{0}_{1}".format(sp,lt)
			# Organizing strings
			self.appLogger.debug("{0}|{1}-{2}".format(sp, pair1,pair2))
			seq1=seqDict[tag][pair1]["sequence"]
			seq2=seqDict[tag][pair2]["sequence"]
			shortDesc1=seqDict[tag][pair1]["description"]
			des1=">{0}:{1:0{2}d}:{3:0{4}d}:{5}:{6}:{7}_S1".format(\
				self.settings.projectName,\
				self.filteredReplicates[index],\
				self.numReplicateDigits,\
				indexLOC,\
				self.numLociPerReplicateDigits[index],\
				self.settings.simphyDataPrefix,\
				currentInd,shortDesc1\
			)
			shortDesc2=seqDict[tag][pair2]["description"]
			des2=">{0}:{1:0{2}d}:{3:0{4}d}:{5}:{6}:{7}_S2".format(\
				self.settings.projectName,\
				self.filteredReplicates[index],\
				self.numReplicateDigits,\
				indexLOC,\
				self.numLociPerReplicateDigits[index],\
				self.settings.simphyDataPrefix,\
				currentInd,\
				shortDesc2\
			)
			indFilename=os.path.join(\
				outputFolder,\
				"{0}_{1:0{2}d}_{3:0{4}d}_{5}_{6}.fasta".format(\
					self.settings.projectName,\
					self.filteredReplicates[index],\
					self.numReplicateDigits,\
					indexLOC,\
					self.numLociPerReplicateDigits[index],\
					self.settings.simphyDataPrefix,\
					currentInd\
				)\
			)
			indFile=open(indFilename, "w")
			indFile.write("{0}\n{1}\n{2}\n{3}\n".format(des1,seq1,des2,seq2))
			indFile.close()
			del seqDict[tag][pair1]
			if not sp== "0": del seqDict[tag][pair2]

		for tag in seqDict.keys():
			if not seqDict[tag]=={}:
				self.appLogger.warning("Number of individuals (sequences) generated per species is odd.")
				self.appLogger.warning("Sequence {0} from species {1} will not be paired.".format(\
					seqDict[tag].keys(), sp))
				for item in seqDict[tag].keys():
					del seqDict[sp][item]

		return None

	def checkFilesForIndels(self):
		"""
		Checks the FASTA files within the species tree replicates
		for indels. Since the RC process won't be ran if they appear in the
		dataset.
		------------------------------------------------------------------------
		Returns:
		- a boolean:
			TRUE if it does not have indel, meaning dataset is OK for processing.
			FALSE if it has indel, meaning dataset is INCORRECT for processing.

		"""
		self.appLogger.info("Checking for indels...")
		indelsList=[];self.indelsPresence=False;message=""
		endStatus=True
		try:
			for index in range(0,len(self.filteredReplicates)):
				self.appLogger.info("ReplicateID {0} - {1}/{2} [{3}] ".format(\
					self.filteredReplicates[index],\
					index+1,\
					len(self.filteredReplicates),\
					self.numLociPerReplicate[index]
				))
				genomicdata=""
				for indexLOC in range(1,self.numLociPerReplicate[index]+1):
					self.appLogger.debug("indexREP: {0}\tlocID: {1}".format(\
						self.filteredReplicates[index], indexLOC))
					if(self.settings.inputmode<4):
						fastapath=os.path.join(
							self.settings.basepath,\
							"REPLICATE_{0:0{1}d}".format(self.filteredReplicates[index],self.numReplicateDigits),\
							"{0}_{1:0{2}d}_TRUE.fasta".format(\
								self.settings.simphyDataPrefix,\
								indexLOC,\
								self.numLociPerReplicateDigits[index]\
							)\
						)
					else:
						fastapath=os.path.join(
							self.settings.basepath,\
							"{0:0{1}d}".format(self.filteredReplicates[index],self.numReplicateDigits),\
							"{0}_{1:0{2}d}_TRUE.fasta".format(\
								self.settings.simphyDataPrefix,\
								indexLOC,\
								self.numLociPerReplicateDigits[index]\
							)\
						)
					genomicdata=None
					with open(fastapath, "r") as fasta:
						genomicdata="".join([line.strip() for line in fasta if not line.startswith(">")])
					dashpos=int(genomicdata.find("-"))
					if dashpos >= 0: indelsList+=[True] #indel
					else: indelsList+=[False] # ok
			if sum(indelsList) > 0: self.indelsPresence=True
		except IOError as e:
			endStatus=False
			message="{}\t{} {} {}\n\t{}\n\t{}\n\t{}\n\t{}\n\n\t{}\n\t{}\n\t{}\n\t{}".format(\
				"[IOERROR]",\
				"File",\
				fastapath,\
				"does not exist",\
				"Input data should be in FASTA format. Please verify your dataset.",\
				"--------------------------------------------------------------------------------",\
				"If using input modes 1,2 or 3, make sure the control file is properly structured.",\
				"If using input mode 4, make sure the INDELible control file for your SimPhy simulations contains the settings lines:",\
				"[SETTINGS]",\
				"  [fastaextension] fasta",\
				"  [output] FASTA",\
				"--------------------------------------------------------------------------------"\
				)
		return endStatus, message
