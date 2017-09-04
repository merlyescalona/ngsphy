#!/usr/bin/env python
# -*- coding: utf-8 -*-
import argparse,copy,csv,datetime,logging,os,subprocess,sys,multiprocessing,threading,time,warnings
import numpy as np
import random as rnd
import settings as sp
from msatools import *
from coverage import NGSPhyDistribution as ngsphydistro

try:
	from collections import Counter
except ImportError:
	from counter import Counter

def getScoreSingle(data):
	"""
	Calulates a single Phred score.
	----------------------------------------------------------------------------
	Parameters:
	- data: the data that will be processed
	Returns:
	- the Phred score of the given value
	"""
	self.appLogger.debug("getScoreSingle(data)")
	if data!=0:
		return -10*np.log10(data)
	else:
		return 0
def getScoreMatrix(data):
	"""
	Calulates the Phred score of a matrix.
	----------------------------------------------------------------------------
	Parameters:
	- data: the data that will be processed
	Returns:
	- the Phred score of the given values
	"""
	self.appLogger.debug("getScoreMatrix(data)")
	value=np.copy(data)
	with warnings.catch_warnings():
		warnings.filterwarnings('error')
		try:
			value=-10*np.log10(data)
		except:
			pass
		infs=np.inf==value
		value[infs]=0
	return value

class RunningInfo:
	"""
	Class for the storage of running time information Reads Counts Class
	Separated to be able to be used by several threads.
	"""
	def __init__(self, filename):
		self.filename=filename
		# self.appLogger=logging.getLogger('sngsw')
		self.lock = threading.Lock()
		self.lock.acquire()
		f=open(self.filename,"w")
		f.write("indexREP,indexLOC,indID,inputFile,cpuTime,seed,outputFilePrefix\n")
		f.close()
		self.lock.release()


	def addLine(self, line):
		"""
		Adds a line to the running information estructure
		"""
		# self.appLogger.debug('Waiting for lock')
		self.lock.acquire()
		try:
			# self.appLogger.debug('Acquired lock')
			f=open(self.filename,"a")
			f.write(
				str(line[0])+","+\
				str(line[1])+","+\
				str(line[2])+","+\
				str(line[3])+","+\
				str(line[4])+","+\
				str(line[5])+","+\
				line[6]+"\n"
			)
			f.close()
		finally:
			# self.appLogger.debug('Released lock')
			self.lock.release()

class ReadCounts:
	"""
	Class to generate READ COUNTS of the given dataset
	----------------------------------------------------------------------------
	Atributes:
	- __NUCLEOTIDES: constant array with possible nucleotides characters
	- appLogger: logger to store status of the process flow
	- settings: Settings object withh all the program parameters
	- runningInfo: varible to store timing information of the processes.
	Structure to be shared across threads.
	- coverageFolderPath: path of the coverage matrix distribution.
	- settings.readsFolderPath: path where read count data will be stored.
	- readcountNoErrorFolderPath: path where true read count data will be stored.
	- readcountErrorFolderPath: path where sampled read count data will be
	stored.
	- numReplicates: number of species tree.
	- numLociPerReplicates: list with the number of loci per species tree
	- numIndividualsPerReplicate: list with the number of individuals per
	species tree.
	- numReplicateDigits: number of digits that represent numReplicates
	- numLociPerReplicateDigits: list of the number of digits that represent
	numLociPerReplicates.
	- numIndividualsPerReplicateDigits: list of the number of digits that represent
	numIndividualsPerReplicate.
	- filteredReplicates: identifiers of the species tree replicates that will be used.
	"""

	__NUCLEOTIDES=["A","C","G","T"]
	appLogger=None
	settings=None
	runningInfo=None
	# path related variables
	readcountNoErrorFolderPath=""
	readcountErrorFolderPath=""

	numReplicates=None
	numLociPerReplicates=[]
	numIndividualsPerReplicate=[]
	numReplicateDigits=[]
	numLociPerReplicateDigits=[]
	numIndividualsPerReplicateDigits=[]
	filteredReplicates=[]

	def __init__(self,settings):
		self.appLogger=logging.getLogger('ngsphy')
		self.appLogger.info('Read counts.')
		self.settings=settings

		self.numReplicates=self.settings.parser.getint("general", "numreplicates")
		cc=self.settings.parser.get("general", "numLociPerReplicates").strip().split(",")
		self.numLociPerReplicates=[ int(item) for item in cc if not item ==""]
		cc=self.settings.parser.get("general", "numIndividualsPerReplicate").strip().split(",")
		self.numIndividualsPerReplicate=[ int(item) for item in cc if not item == ""]
		self.numReplicateDigits=len(str(self.numReplicates))
		self.numLociPerReplicateDigits=[len(str(item)) for item in self.numLociPerReplicates]
		self.numIndividualsPerReplicateDigits=[len(str(item )) for item in self.numIndividualsPerReplicate]
		self.filteredReplicates=[int(item) for item in self.settings.parser.get("general", "filtered_replicates").strip().split(",")]

		if(self.settings.runningTimes):
			# Generating folder structure
			infoFile=os.path.join(\
				self.settings.outputFolderPath,\
				"{0}.info".format(self.settings.projectName)\
			)
			self.runningInfo=RunningInfo(infoFile)
			self.appLogger.info("File with timings of the [ngs-read-counts] per loci can be find on:\t{0}".format(infoFile))

		# naming variables for shorter codeline
		self.readcountNoErrorFolderPath=os.path.join(\
			self.settings.readsFolderPath,"no_error"\
		)
		self.readcountErrorFolderPath=os.path.join(\
			self.settings.readsFolderPath,"with_error"\
		)
		self.generateFolderStructureGeneral()

	def generateFolderStructureGeneral(self):
		"""
		Generates basic folder structure needed for Read Count option.
		This includes the folders:
		- reads
		- reads/true
		- reads/sampled
		- refereces
		"""
		self.appLogger.debug("generateFolderStructureGeneral(self)")
		self.appLogger.info("Creating folder structure for [ngs-read-counts]")
		try:
			os.makedirs(self.settings.readsFolderPath)
			self.appLogger.info("Generating output folder ({0})".format(\
				self.settings.readsFolderPath))
		except:
			self.appLogger.debug("Output folder exists ({0})".format(\
				self.settings.readsFolderPath))
		try:
			os.makedirs(self.readcountNoErrorFolderPath)
			self.appLogger.info("Generating output folder ({0})".format(\
				self.readcountNoErrorFolderPath))
		except:
			self.appLogger.debug("Output folder exists ({0})".format(\
				self.readcountNoErrorFolderPath))
		try:
			os.makedirs(self.readcountErrorFolderPath)
			self.appLogger.info("Generating output folder ({0})".format(\
				self.readcountErrorFolderPath))
		except:
			self.appLogger.debug("Output folder exists ({0})".format(\
				self.readcountErrorFolderPath))
		try:
			os.makedirs(self.refAllelesFolderPath)
			self.appLogger.info("Generating output folder ({0})".format(\
				self.refAllelesFolderPath))
		except:
			self.appLogger.debug("Output folder exists ({0})".format(\
				self.refAllelesFolderPath))

	def generateFolderStructureDetail(self):
		"""
		Generates detailed folder structure needed for Read Count option.
		This includes:
		- Replicate folders  (as many as species tree replicates)
		- true/sampled folder per replicate.
		- These folders will contain the VCF files (true and sampled)
		"""
		self.appLogger.debug("generateFolderStructureDetail(self)")
		for i in range(0,len(self.filteredReplicates)):
			indexREP=self.filteredReplicates[i]
			self.appLogger.debug("{0}|{1}".format(i, self.filteredReplicates[i]))
			noErrorFolder=os.path.join(\
				self.readcountNoErrorFolderPath,\
				"{0:0{1}d}".format(indexREP,self.numReplicateDigits)\
			)
			errorFolder=os.path.join(\
				self.readcountErrorFolderPath,\
				"{0:0{1}d}".format(indexREP,self.numReplicateDigits)\
			)
			ref_alleles=os.path.join(\
				self.refAllelesFolderPath,\
				"{0:0{1}d}".format(indexREP,self.numReplicateDigits)\
			)

			try:
				os.makedirs(noErrorFolder)
				self.appLogger.info("Generating output folder ({0})".format(noErrorFolder))
			except:
				self.appLogger.debug("Output folder exists ({0})".format(noErrorFolder))
			try:
				os.makedirs(errorFolder)
				self.appLogger.info("Generating output folder ({0})".format(errorFolder))
			except:
				self.appLogger.debug("Output folder exists ({0})".format(errorFolder))
			try:
				os.makedirs(ref_alleles)
				self.appLogger.info("Generating output folder ({0})".format(ref_alleles))
			except:
				self.appLogger.debug("Output folder exists ({0})".format(ref_alleles))


	def retrieveCoverageMatrix(self, indexREP):
		"""
		Reads coverage matrix for a specific species tree identifier
		------------------------------------------------------------------------
		Parameters:
		- indexREP:specific species tree identifier
		Returns:
		- coverage Matrix
		"""
		self.appLogger.debug("Retrieving coverage matrix")
		# coverage matrix per ST - row:indv - col:loci
		# each cov introduced as parameter is a NGSPhyDistribution
		coverageMatrixFilename=os.path.join(\
			self.settings.coverageFolderPath,\
			"{0}.{1:0{2}d}.coverage.csv".format(\
				self.settings.projectName,\
				indexREP,\
				self.numReplicateDigits\
			)
		)
		coverageMatrix=np.zeros(\
			shape=(\
				self.numIndividualsPerReplicate[indexREP-1],\
				self.numLociPerReplicates[indexREP-1]\
			)\
		)
		firstLine=NOError; counter=0
		with open(coverageMatrixFilename, 'rb') as csvfile:
			coveragereader = csv.reader(csvfile, delimiter=',')
			for row in coveragereader:
				if firstLine:
					firstLine=False
				else:
					coverageMatrix[counter,]=[float(row[index]) for index in range(1,len(row))]
					counter+=1
				if not counter < self.numIndividualsPerReplicate[indexREP-1]: break
		return coverageMatrix


	def parseReferenceAllelesList(self, filename):
		"""
		Used to parse referenceList, file with format: STID,SPID,TIPID
		-----------------------------------------------------------------------
		Parameters:
		- filename: path of the referenceList file.
		There's only ONE file with the relation of the ref_alleles
		If "None" inputted (file is missing) then reference by default is 1_0_0
		for all species tree replicates.
		Returns:
		- output: list. each element of the list is a triplet
				(REPLICATEID,SPID,LOCID,GENEID)
		"""
		self.appLogger.debug("parseReferenceAllelesList(self, filename)")
		referenceList=[]
		if filename:
			# There's a file
			filepath=os.path.abspath(filename)
			f=open(filepath,"r")
			lines=f.readlines()
			f.close()
			lines=[line.strip().split(",") for line in lines]
			lines=sorted(lines, key=lambda x:x[1])
			skipped=False
			message="Parsing reference list. "+\
				"A default reference has been introduced.\n"+\
				"Replicate index:"

			index=0; itemIndex=0
			while (index < len(lines)) and (itemIndex < self.numReplicates):
				d=lines[index]
				if  (int(d[0]) == (itemIndex+1)) and (not "" in d):
					try:
						referenceList+=[(d[0],d[1],d[2],d[3])]
					except IndexError as ie:
						skipped=NOError
					index+=1
				else: skipped=NOError

				if skipped:
					self.appLogger.warning("{0}{1}".format(message,(index+1)))
					referenceList+=[(index+1,1,0,0)]
					skipped=False
				itemIndex+=1

			if (itemIndex < self.numReplicates):
			    for i in range(self.numReplicates):
			        print([i+1,1,0,0])
		else:
			# iterate get a list with same decription for all species trees
			for item in range(0,self.numReplicates):
				referenceList+=[(item+1,1,0,0)]
		return referenceList

	def extractNOErrorVariantsPositions(self, filename):
		"""
		Extract true variant positions from the MSA file used as input
		------------------------------------------------------------------------
		Parameters:
		- filename: MSA fasta file from where to extract the variable positions
		Returns:
		- dictionary. indices=variable sites, content=variable nucleotide set
		"""
		self.appLogger.debug("extractNOErrorVariantsPositions(self, filename)")
		filepath=os.path.abspath(filename)
		lines=[];variants=dict();seqDescriptions=[]
		numTotalSeqs=0;lenSeq=0; matrix=None
		self.appLogger.debug(\
			"Extracting variable positions from: {0}".format(\
			filepath
		))
		# Checking sequence length
		if isFasta(filepath):
			f=open(filepath)
			lines=f.readlines()
			f.close()
			seq=lines[1] # lines[0] will be a description
			lenSeqs=len(seq.strip())
		else:
			message="{0}\n{1}\n{2}\n".format(\
				"File has a wrong file format.",\
				filepath,\
				"Please check. Exiting."
			)
			raise TypeError(message)
			# If I get here, function run is over - goes back to  run()

		for item in lines:
			if item.startswith(">"): numTotalSeqs+=1
		matrix=np.chararray((numTotalSeqs,lenSeqs), itemsize=1)
		# Cleaning strings - removing empty lines and removing "\n"
		for index in range(0,len(lines)):
			lines[index]=lines[index].strip()

		try:
		  lines.remove("")
		except: # may raise an exception if list empty (x not in list)
			pass

		numLinesFile=len(lines);index=0
		indexSeqs=range(1,numLinesFile,2)

		for index in range(0,numTotalSeqs):
			seqDescriptions+=[lines[indexSeqs[index]-1]]
			matrix[index,:]=list(lines[indexSeqs[index]])

		for indexCol in range(0,matrix.shape[1]):
			c=Counter(matrix[:,indexCol])
			l=np.unique(matrix[:,indexCol])
			if (len(l)>1):
				variants[str(indexCol)]=c.keys()
		return variants

	def parseIndividualRelationFile(self,filename):
		"""
		Parses "Individual-description relation" file per species tree
		in order to obtain information on how sequences from INDELible
		are related. Generated within NGSPhy.
		------------------------------------------------------------------------
		Parameters:
		- filename: path for the Individual-description relation file
		Returns:
		- if ploidy=1:	returns a dict. key: indID,content: dict(indexREP,speciesID,locusID, geneID)
		- if ploidy=2:	returns a dict. key: indID,content: dict(indexREP,speciesID,locusID, mateID1,mateID2)
		"""
		self.appLogger.debug("parseIndividualRelationFile(self,filename)")
		individuals=dict()
		if (self.settings.ploidy>0 and self.settings.ploidy<=2):
			csvfile=open(os.path.abspath(filename))
			d = csv.DictReader(csvfile)
			if (self.settings.ploidy==1):
				for row in d:
					individuals[str(row["indID"])] = {\
						"indexREP":row["indexREP"],\
						"speciesID":row["speciesID"],\
						"locusID":row["locusID"],\
						"geneID":row["geneID"]}
			if (self.settings.ploidy==2):
				# indexREP,indID,speciesID,mateID1,mateID2
				for row in d:
					individuals[str(row["indID"])] = {\
						"indexREP":row["indexREP"],\
						"speciesID":row["speciesID"],\
						"locusID":row["locusID"],\
						"mateID1":row["mateID1"],\
						"mateID2":row["mateID2"]}
			csvfile.close()
		else:
			# There has been a verification in Settings class, but just in case.
			message="\n\t{0}\n\t{1}".format(\
				"Ploidy assigned is incorrect.",\
				"Please verify. Exciting."
			)
			raise ValueError(message)
		return individuals

	def getDepthCoveragePerIndividual(self, numVarSites,startingCoverage):
		"""
		Compute coverage per individuals
		since the value for SNPS is fixed to a Negative Binomial distribution
		I iterate over the number of variants for this specific locus and
		since everythime i hit play! (called value()) I'll sample a new value
		from the previously set distribution, i'll then get a random value from a
		negative bionmial distribution with mean max-coverage
		------------------------------------------------------------------------
		Parameters:
		- numVariableSites: number of variable sites
		- startingCoverage: starting value to calculate coverage
		Returns:
		- list of values that correspond to the coverage per each snp of the MSA
		"""
		self.appLogger.debug(\
			"{} - ({},{})".format(\
				"getDepthCoveragePerIndividual(self, numVarSites,startingCoverage)",\
				numVarSites,\
				startingCoverage\
			)\
		)
		distro=ngsphydistro("nb",[startingCoverage,startingCoverage])
		DP=distro.value(numVarSites)
		return DP

	def getHaploidIndividualSequence(self,msa,ind):
		"""
		Extract individuals "ind" sequence(s) from MSA dictionary
		------------------------------------------------------------------------
		Parameters:
		- msa: dictionary
		- ind: dict(indexREP,seqDEscription)
		Returns:
		- sequence of the individual
		"""
		# ind: [indID,indexREP,seqDescription]
		self.appLogger.debug("getHaploidIndividualSequence(self,msa,ind)")
		seqSize=len(msa["{0}_{1}".format(str(1),str(0))][str(0)]['sequence'])
		fullInd=None; speciesID=None; tipID=None; tmp=None
		fullInd=np.chararray(shape=(1,seqSize), itemsize=1)
		speciesID=ind["speciesID"].strip()
		locusID=ind["locusID"].strip()
		tipID=ind["geneID"].strip()
		tmp=list(msa["{0}_{1}".format(str(speciesID), str(locusID))][str(tipID)]['sequence'])
		fullInd=[item for item in tmp]
		return fullInd

	def getDiploidIndividualSequence(self,msa,ind):
		"""
		Extract individuals "ind" sequence(s) from MSA dictionary
		Parameters:
		- msa: dictionary
		- ind: dict(indexREP,speciesID, mateID1,mateID2)
		Returns:
		- matrix: 2 x seqLength - representing the sequence of the individual
		"""
		# ind:[indID, indexREP,speciesID, mateID1,mateID2]
		self.appLogger.debug("getDiploidIndividualSequence(self,msa,ind)")
		seqSize=len(msa["1_0"][str(0)]['sequence'])
		fullInd=None;  tmp=None
		speciesID=ind["speciesID"];
		locusID=ind["locusID"];
		tipID1=ind["mateID1"];
		tipID2=ind["mateID2"];
		geneFamily="{0}_{1}".format(speciesID,locusID)
		fullInd=np.chararray(shape=(2,seqSize), itemsize=1)
		tmp=list(msa[geneFamily][str(tipID1)]['sequence'])
		fullInd[0,:]=[item for item in tmp]
		tmp=list(msa[geneFamily][str(tipID2)]['sequence'])
		fullInd[1,:]=[item for item in tmp]
		return fullInd


	def computeHaploid(self,indexREP,indexGT,msa,individuals,refAllelesFilePath,\
		referenceSeqFull,variableSites,DP):
		"""
		compute the READ COUNT for the specic ploidy
		------------------------------------------------------------------------
		Parameters:
		- indexREP: species tree id
		- indexGT: locus/gene tree id
		- msa: parsed multiple sequence alignent file (dictionary)
		- individuals: parsed individual-description relation file (dictionary)
		- refAllelesFilePath: where the reference is written.
		- referenceSeqFull: reference sequence as is
		- variableSites: dictionary. keys: variable positions. content: variable nucleotide set
		- DP: variable sites coverage
		Returns:
		- get a file written in the STreplicate folder (true/sampled)
		"""
		self.appLogger.debug("computeHaploid(...)")
		nInds=len(individuals)
		nVarSites=len(variableSites.keys())
		# Getting the proper indices of the variable sites
		variableSitesPositionIndices=np.sort([int(pos) for pos in variableSites.keys()])
		alt=dict();  altInds=dict()
		rcsNOError=dict();rcsWErrors=dict()
		HTGeneral=dict();HLGeneral=dict();ADGeneral=dict()
		HTGeneralWErrors=dict();HLGeneralWErrors=dict();ADGeneralWErrors=dict()
		for pos in variableSites.keys():
			alt[pos]=[]; altInds[pos]=[]
		for pos in variableSitesPositionIndices:
			alt[str(pos)]=list(set(variableSites[str(pos)])-set([referenceSeqFull[pos]]))
		for index in range(0,nInds):
			indexIND=str(index)
			HTGeneral[indexIND]=[];HLGeneral[indexIND]=[];ADGeneral[indexIND]=[]
			HTGeneralWErrors[indexIND]=[];
			HLGeneralWErrors[indexIND]=[];ADGeneralWErrors[indexIND]=[]
			rcsNOError[indexIND]=[];rcsWErrors[indexIND]=[]
		########################################################
		# TRUE
		########################################################
		# iterate over individuals

		self.appLogger.debug("NOError")
		altInd=copy.copy(alt)
		for indIndex in range(0,nInds):
			indexIND=str(indIndex)
			self.appLogger.debug(\
				"\tReplicate tree: {0} - Locus {1} |  Individual {2} ({3}/{4})".format(\
					indexREP,\
					indexGT,\
					indIndex,\
					(indIndex+1),\
					nInds)\
				)
			ind=individuals[indexIND]
			individualSeq=self.getHaploidIndividualSequence(msa,ind)
			# individualSeq[202]+="**"
			# AD per individual different if true/sampledInd
			ADNOErrors,ADWErrors=self.getAllelicDepthPerHaploidIndividual(\
				individualSeq,variableSitesPositionIndices,DP[indIndex])
			rcNOError,rcWErrors=self.getReadCountPerIndividual(\
				ADNOErrors,ADWErrors,variableSitesPositionIndices)
			altInd=self.getAltUpdatedPerIndividual(\
				referenceSeqFull,altInd,ADWErrors)
			HTNOError=self.gettingHaplotype(\
				referenceSeqFull,individualSeq,alt,variableSitesPositionIndices)
			HLNOError=self.haplotypeLikehood(\
				rcNOError,variableSitesPositionIndices,0)
			rcsNOError[indexIND]=rcNOError
			rcsWErrors[indexIND]=rcWErrors
			ADGeneral[indexIND]=ADNOErrors
			ADGeneralWErrors[indexIND]=ADWErrors
			HTGeneral[indexIND]=HTNOError
			HLGeneral[indexIND]=HLNOError

		# here corresponds to the INDEXST that is going to be written
		self.writeVCFFile(\
			indexREP,indexGT,\
			refAllelesFilePath, referenceSeqFull,\
			alt,variableSitesPositionIndices,\
			HTGeneral,HLGeneral,\
			ADGeneral,DP,NOError)

		########################################################
		# With errors
		########################################################
		self.appLogger.debug("With errors")
		for indIndex in range(0,nInds):
			indexIND=str(indIndex)
			self.appLogger.info("Individual {0} ({1}/{2})".format(indIndex,(indIndex+1),nInds))
			ind=individuals[indexIND]
			individualSeq=self.getHaploidIndividualSequence(msa,ind)
			HTGeneralWErrors[indexIND]=self.gettingHaplotype(\
				referenceSeqFull,individualSeq,\
				altInd, variableSitesPositionIndices)
			HLGeneralWErrors[indexIND]=self.haplotypeLikehood(\
				rcsWErrors[indexIND],\
				variableSitesPositionIndices,\
				self.settings.seqerror)

		self.writeVCFFile(
			indexREP,indexGT,\
			refAllelesFilePath, referenceSeqFull,\
			altInd,variableSitesPositionIndices,\
			HTGeneralWErrors,HLGeneralWErrors,\
			ADGeneralWErrors,DP,False)

	def haplotypeLikehood(self,variantsRC,variableSitesPositionIndices,error):
		"""
		computes log10-scaled haplotype likelihood per variable site
		------------------------------------------------------------------------
		Parameters:
		- variantsRC: read count pairs of the variable positions
		- variableSitesPositionIndices: indices of the variable sites position
		- error: sequence error
		Returns:
		- GL: dictionary.
		- keys: positons.
		- content: log10-scaled likelihood for the specific variable site
		"""
		self.appLogger.debug("haplotypeLikehood(self,variantsRC,observed,variableSites,error)")
		nVariants=len(variableSitesPositionIndices)
		HL=np.ones(shape=(4,nVariants), dtype=np.float)
		error=float(error)
		for indexVar in range(0,nVariants):
			indexReads=str(variableSitesPositionIndices[indexVar])
			reads=variantsRC[indexReads]
			# possible haplotypes
			for indexNuc in range(0,4):
				nuc=self.__NUCLEOTIDES[indexNuc]
				# read to be analyzed
				for b in reads:
					if (b.upper()==nuc.upper()):
						HL[indexNuc,indexVar]=HL[indexNuc,indexVar]*(1-error)
					else:
						HL[indexNuc,indexVar]=HL[indexNuc,indexVar]*(error/3)
				# print indexVar, "b: ",b, " A:",nuc,",".join(HL[:,indexVar].astype(str))
		value=np.copy(HL)
		with warnings.catch_warnings():
			warnings.filterwarnings('error')
			try:
				value=np.log10(HL)
			except:
				pass
			infs=np.inf==value
			value[infs]=-300
		return value

	def gettingHaplotype(self,ref,seq,alt, variableSitesPositionIndices):
		"""
		get haplotype for the given ref-seq-alt triplet
		------------------------------------------------------------------------
		Parameters:
		- ref: reference sequence
		- seq: individual sequnece
		- alt: dictionary with alternative alleles per variable position
		- variableSitesPositionIndices: variableSites
		Returns:
		- GT: dictionary. keys: positons. content: haplotype for the specific variable site
		"""
		self.appLogger.debug("gettingHaplotype(self,ref,seq,alt, variableSitesPositionIndices)")
		GT=dict()
		# init variantsRC
		for indexVar in variableSitesPositionIndices:
			GT[str(indexVar)]=0

		for indexVar in variableSitesPositionIndices:
			refPos=ref[indexVar]
			seqPos=seq[indexVar]
			altValues=alt[str(indexVar)]
			iv=str(indexVar)
			for index in range(0,len(altValues)):
				altToCompare=altValues[index]
				if (index==0) and (altToCompare==seqPos):
					GT[iv]=1
				if (index==1) and (altToCompare==seqPos):
					GT[iv]=2
				if (index==2) and (altToCompare==seqPos):
					GT[iv]=3
		return copy.copy(GT)

	def getAltUpdatedPerIndividual(self,ref,alt,AD):
		"""
		update general alternative alleles list
		------------------------------------------------------------------------
		Parameters:
		- referenceSeq: reference sequence
		- alt: dictionary. current alternative allele list per variant
		- AD: allelic depth matrix
		Return:
		- newAlt: dictionary. keys: positions. content: alt alleles corresponding to that position
		"""
		self.appLogger.debug("getAltUpdatedPerIndividual(self,ref,alt,AD)")
		# update alt values
		# With errors - ALT is a DICT
		# altUpdate is going to be a dict too
		altUpdated=dict()
		for item in alt.keys(): altUpdated[item]=[]
		sortedAltKeys=np.sort([int(item) for item in alt.keys()])
		for index in range(0,len(sortedAltKeys)):
			possibleNucsAlt=[]
			if (AD[0,index]>0): possibleNucsAlt+=["A"]
			if (AD[1,index]>0): possibleNucsAlt+=["C"]
			if (AD[2,index]>0): possibleNucsAlt+=["G"]
			if (AD[3,index]>0): possibleNucsAlt+=["T"]
			newAlt=set(\
					alt[str(sortedAltKeys[index])] + \
					possibleNucsAlt\
				) - \
				set(ref[sortedAltKeys[index]])
			altUpdated[str(sortedAltKeys[index])]+=list(np.sort(list(newAlt)))
		return altUpdated

	def getReadCountPerIndividual(self,ADNOErrors,ADWErrors, variableSitesPositionIndices):
		"""
		generates list of read counts per variant per individual
		------------------------------------------------------------------------
		Parameters:
		- ADNOErrors: NOError allelic depth
		- ADWErrors: With errors allelic depth
		- variableSitesPositionIndices: position of the variable sites
		Return:
		- RCNOError, RCWErrors.
		"""
		self.appLogger.debug("getReadCountPerIndividual(self,ADNOErrors,ADWErrors, variableSitesPositionIndices)")
		nVariants=len(variableSitesPositionIndices)
		rcNOError=dict();rcWErrors=dict()
		# init variantsRC
		for indexVar in range(0,nVariants):
			indexConverted=str(variableSitesPositionIndices[indexVar])
			rcNOError[indexConverted]=[]
			rcWErrors[indexConverted]=[]
		# getting read count
		for indexVar in range(0,nVariants):
			for indexNuc in range(0,4):
				with warnings.catch_warnings():
					warnings.filterwarnings("ignore", category=np.VisibleDeprecationWarning)
					nuc=[self.__NUCLEOTIDES[indexNuc]]
					indexConverted=str(variableSitesPositionIndices[indexVar])
					val=ADNOErrors[indexNuc,indexVar]
					rcNOError[indexConverted]+=nuc*val
					val=ADWErrors[indexNuc,indexVar]
					rcWErrors[indexConverted]+=nuc*val
		return rcNOError, rcWErrors

	def getAllelicDepthPerHaploidIndividual(self,individualSeq,variableSitesPositionIndices,DP):
		"""
		Generate allelic depth per individual per site
		------------------------------------------------------------------------
		Parameters:
		- individualSequence: sequence of the individual
		- variableSitesPositionIndices: position of the variable sites
		- DP: coverage of the variable sites
		Returns:
		- ADNOErrors, ADWErrors.
		"""
		self.appLogger.debug("getAllelicDepthPerHaploidIndividual(self,fullInd,variableSites,DP)")
		nVariants=len(variableSitesPositionIndices)
		ADNOErrors=np.zeros(shape=(4,nVariants), dtype=np.int)
		ADWErrors=np.zeros(shape=(4,nVariants), dtype=np.int)
		self.appLogger.debug("Starting to iterate over the variants")
		# There was a problem in this block related to low coverage, error
		# and possible substitutions.
		for indexVar in range(0,nVariants):
			finalRC=[];indNucs=[]
			# getting general coverage per position
			posCoverage=DP[indexVar]
			indexSeq=variableSitesPositionIndices[indexVar]
			indNucs=individualSeq[indexSeq]
			finalRC=[indNucs]*(posCoverage)
			counter=Counter(finalRC)
			# TRUE READ count
			ADNOErrors[0,indexVar]=counter["A"];ADNOErrors[1,indexVar]=counter["C"]
			ADNOErrors[2,indexVar]=counter["G"];ADNOErrors[3,indexVar]=counter["T"]
			# SAMPLE READ COUNT - need to know error distribution
			errorDistro=ngsphydistro("b",[posCoverage,self.settings.seqerror])
			errorD=errorDistro.value(1)[0]
			errorPositions=[]
			# need to know possible nucleotides to substitute my position with error
			possibleNucs=list(set(self.__NUCLEOTIDES)-set([indNucs]))
			# I have some positions (at least 1) that is an error
			# errorD= array with coded error nucleotides that will be modified
			if (errorD>0):
				errorChoices=np.random.choice(possibleNucs, int(errorD), replace=NOError)
				maxAvailablePositions=posCoverage
				if not ((posCoverage % 2) == 0): maxAvailablePositions=posCoverage-1
				if maxAvailablePositions == 0:  maxAvailablePositions=posCoverage

				errorPositions=np.random.randint(maxAvailablePositions,size=int(errorD))
				for item in range(0,len(errorPositions)):
					posToChange=errorPositions[item]
					finalRC[posToChange]=errorChoices[item]
			counter=Counter(finalRC)
			# if any of the nucleotides does not have a counter retrieving it will be 0
			ADWErrors[0,indexVar]=counter["A"];ADWErrors[1,indexVar]=counter["C"]
			ADWErrors[2,indexVar]=counter["G"];ADWErrors[3,indexVar]=counter["T"]
		return ADNOErrors,ADWErrors

	def computeDiploid(self,indexREP,indexGT,msa,individuals,refAllelesFilePath,referenceSeqFull,variableSites,DP):
		"""
		compute the READ COUNT for the specic ploidy
		------------------------------------------------------------------------
		Parameters:
		- indexREP: species tree id
		- indexGT: locus/gene tree id
		- msa: parsed multiple sequence alignent file (dictionary)
		- individuals: parsed individual-description relation file (dictionary)
		- referenceSeqFull: reference sequence as is
		- variableSites: dictionary. keys: variable positions. content: variable nucleotide set
		- DP: variable sites coverage
		Returns:
		- get a file written in the STreplicate folder (true/sampled)
		"""
		self.appLogger.debug("computeDiploid(...)")
		nInds=len(individuals)
		nVarSites=len(variableSites.keys())
		# Getting the proper indices of the variable sites
		variableSitesPositionIndices=np.sort([int(pos) for pos in variableSites.keys()])
		alt=dict();  altInds=dict()
		rcsNOError=dict();rcsWErrors=dict()
		GTGeneral=dict();GLGeneral=dict();ADGeneral=dict()
		GTGeneralWErrors=dict();GLGeneralWErrors=dict();ADGeneralWErrors=dict()
		for pos in variableSites.keys():
			alt[pos]=[]; altInds[pos]=[]
		for pos in variableSitesPositionIndices:
			alt[str(pos)]=list(set(variableSites[str(pos)])-set([referenceSeqFull[pos]]))
		for index in range(0,nInds):
			indexIND=str(index)
			GTGeneral[indexIND]=[];GLGeneral[indexIND]=[];ADGeneral[indexIND]=[]
			GTGeneralWErrors[indexIND]=[];
			GLGeneralWErrors[indexIND]=[];ADGeneralWErrors[indexIND]=[]
			rcsNOError[indexIND]=[];rcsWErrors[indexIND]=[]
		########################################################
		# DIPLOID TRUE
		########################################################
		# iterate over individuals
		self.appLogger.debug("NOError")
		altInd=copy.copy(alt)
		for indIndex in range(0,nInds):
			GTNOError=dict()
			indexIND=str(indIndex)
			self.appLogger.debug("\tSpecies tree: {0} - Locus {1} |  Individual {2} ({3}/{4})".format(\
				indexREP, indexGT,
				indIndex,(indIndex+1),nInds))
			ind=individuals[indexIND]
			individualSeq=self.getDiploidIndividualSequence(msa,ind)
			ADNOErrors,ADWErrors=self.getAllelicDepthPerDiploidIndividual(\
				individualSeq,variableSitesPositionIndices,DP[indIndex])
			rcNOError,rcWErrors=self.getReadCountPerIndividual(\
				ADNOErrors,ADWErrors,variableSitesPositionIndices)
			altInd=self.getAltUpdatedPerIndividual(\
				referenceSeqFull,altInd,ADWErrors)
			GTNOErrorS0=self.gettingHaplotype(
				referenceSeqFull,individualSeq[0,:],\
				alt, variableSitesPositionIndices
			)
			GTNOErrorS1=self.gettingHaplotype(
				referenceSeqFull,individualSeq[1,:],\
				alt, variableSitesPositionIndices
			)
			for item in variableSitesPositionIndices:
				GTNOError[str(item)]=[GTNOErrorS0[str(item)],GTNOErrorS1[str(item)]]
			GLNOError=self.genotypeLikehood(\
				rcNOError,variableSitesPositionIndices,0)
			rcsNOError[indexIND]=rcNOError
			rcsWErrors[indexIND]=rcWErrors
			ADGeneral[indexIND]=ADNOErrors
			ADGeneralWErrors[indexIND]=ADWErrors
			GTGeneral[indexIND]=GTNOError
			GLGeneral[indexIND]=GLNOError

		# here corresponds to the INDEXST that is going to be written
		self.writeVCFFile(\
			indexREP,indexGT,\
			refAllelesFilePath, referenceSeqFull,\
			alt,variableSitesPositionIndices,\
			GTGeneral,GLGeneral,\
			ADGeneral,DP,NOError)

		########################################################
		# DIPLOID SAMPLED
		########################################################
		self.appLogger.debug("With errors")
		for indIndex in range(0,nInds):
			GTWErrors=dict()
			indexIND=str(indIndex)
			self.appLogger.debug("\tSpecies tree: {0} - Locus {1} |  Individual {2} ({3}/{4})".format(\
				indexREP, indexGT,
				indIndex,(indIndex+1),nInds))
			ind=individuals[indexIND]
			individualSeq=self.getDiploidIndividualSequence(msa,ind)
			GTWErrorsS0=self.gettingHaplotype(
				referenceSeqFull,individualSeq[0,:],\
				altInd, variableSitesPositionIndices
			)
			GTWErrorsS1=self.gettingHaplotype(
				referenceSeqFull,individualSeq[1,:],\
				altInd, variableSitesPositionIndices
			)
			for item in variableSitesPositionIndices:
				GTWErrors[str(item)]=[GTWErrorsS0[str(item)],GTWErrorsS1[str(item)]]
			GLWErrors=self.genotypeLikehood(\
				rcsWErrors[indexIND],\
				variableSitesPositionIndices,\
				self.settings.seqerror)
			GTGeneralWErrors[indexIND]=GTWErrors
			GLGeneralWErrors[indexIND]=GLWErrors

		self.writeVCFFile(
			indexREP,indexGT,\
			refAllelesFilePath, referenceSeqFull,\
			altInd,variableSitesPositionIndices,\
			GTGeneralWErrors,GLGeneralWErrors,\
			ADGeneralWErrors,DP,False)

	def genotypeLikehood(self, variantsRC,variableSitesPositionIndices,error):
		"""
		compute the genotype likelihoods of a specific ngs-read-counts set
		------------------------------------------------------------------------
		Parameters:
		- variantsRC:
		- variableSitesPositionIndices:
		- error
		Returns:
		- Genotype likelihoods for the given variable positions and error
		"""
		self.appLogger.debug("genotypeLikehood(self, variantsRC,variableSitesPositionIndices,error)")
		nVariants=len(variableSitesPositionIndices)
		GL=dict();GLout=dict() # Initialization of genotypeLikehood variable
		possibleGenotypes=self.getAllPossibleGenotypes()
		for indexVar in variableSitesPositionIndices:
			GL[str(indexVar)]={}
			GLout[str(indexVar)]={}
			for pg in possibleGenotypes:
				gen="".join(pg)
				GL[str(indexVar)][gen]=1
				GLout[str(indexVar)][gen]=1
		# getting likelihoods

		for indexVar in variableSitesPositionIndices:
			for pg in possibleGenotypes:
				gen="".join(pg)
				reads=variantsRC[str(indexVar)]
				for b in reads:
					A1=pg[0];A2=pg[1]
					valAl1=0; valA2=0
					if (b==A1):
						valA1=0.5*(1-error)
					else:
						valA1=0.5*(error/3)
					if (b==A2):
						valA2=0.5*(1-error)
					else:
						valA2=0.5*(error/3)
					GL[str(indexVar)][gen]=\
						GL[str(indexVar)][gen]*(valA1+valA2)
		for indexVar in variableSitesPositionIndices:
			for pg in possibleGenotypes:
				gen="".join(pg)
				value=GL[str(indexVar)][gen]
				if value==0:
					GLout[str(indexVar)][gen]=-300
				else:
					GLout[str(indexVar)][gen]=np.log10(value)
		return GL

	def getAllPossibleGenotypes(self):
		"""
		Generates all the possible genotypes
		------------------------------------------------------------------------
		returns a dictionary. key: positions. content [A1,A2]
		"""
		self.appLogger.debug("getAllPossibleGenotypes(self)")
		possibleGenotypes=[]
		for pos1 in self.__NUCLEOTIDES:
			for pos2 in self.__NUCLEOTIDES:
				possibleGenotypes+=[[pos1,pos2]]
		return possibleGenotypes

	def getPossibleGenotypesPerVariableSite(self,ref,alt, variableSitesPositionIndices):
		"""
		get the genotypes for a specfic position
		------------------------------------------------------------------------
		intput:
			ref
			alt
			variableSitesPositionIndices
		output:
			returns a dictionary. key: positions. content [A1,A2]
		"""
		self.appLogger.debug("getPossibleGenotypesPerVariableSite(self,ref,alt, variableSites)")
		self.appLogger.debug("Getting genotypes")
		possibleGenotypes=dict();elems=dict()
		# Initialization
		for varSite in variableSitesPositionIndices:
			possibleGenotypes[str(varSite)]=[]
			elems[str(varSite)]=[]

		for varSite in variableSitesPositionIndices:
			elems[str(varSite)]=[ref[varSite]]+alt[str(varSite)]

		for indexAlt in elems.keys():
			altPos=elems[indexAlt]
			d=[]
			for pos1 in altPos:
				for pos2 in altPos:
					if not (set([pos1,pos2]) in d):
						d+=[set([pos1,pos2])]
						possibleGenotypes[str(indexAlt)]+=["".join([pos1,pos2])]

		return possibleGenotypes

	def genotypeOrder(self,alleles):
		ordered=[]
		for a in range(0,len(alleles)):
			for b in range(0,(a+1)):
				ordered+=[[alleles[a],alleles[b]]]
		return ordered


	def getAllelicDepthPerDiploidIndividual(self,individualSeq,variableSitesPositionIndices,DP):
		"""
		Generate allelic depth per individual per site
		------------------------------------------------------------------------
		intput:
			individualSequence: sequence of the individual
			variableSitesPositionIndices: position of the variable sites
			DP: coverage of the variable sites
		output:
			ADNOErrors, ADWErrors.
		"""
		self.appLogger.debug("getAllelicDepthPerDiploidIndividual(self,fullInd,variableSitesPositionIndices,DP)")
		nVariants=len(variableSitesPositionIndices)
		ADNOErrors=np.zeros(shape=(4,nVariants), dtype=np.int)
		ADWErrors=np.zeros(shape=(4,nVariants), dtype=np.int)
		self.appLogger.debug("Starting to iterate over the variants")
		# There was a problem in this block related to low coverage, error
		# and possible substitutions.
		for indexVar in range(0,nVariants):
			finalRC=[];indNucs=[]
			# getting general coverage per position
			posCoverage=DP[indexVar]
			indexSeq=variableSitesPositionIndices[indexVar]
			indNucStrand1=individualSeq[0,indexSeq]
			indNucStrand2=individualSeq[1,indexSeq]
			diploidDistro=ngsphydistro("b",[posCoverage,0.5])
			diploidCoverage=diploidDistro.value(1)[0]
			finalRC=[indNucStrand1]*(diploidCoverage)
			finalRC+=[indNucStrand2]*(posCoverage-diploidCoverage)
			counter=Counter(finalRC)
			# TRUE READ count
			ADNOErrors[0,indexVar]=counter["A"];ADNOErrors[1,indexVar]=counter["C"]
			ADNOErrors[2,indexVar]=counter["G"];ADNOErrors[3,indexVar]=counter["T"]
			# SAMPLE READ COUNT - need to know error distribution
			errorDistro=ngsphydistro("b",[posCoverage,self.settings.seqerror])
			errorD=errorDistro.value(1)[0]
			errorPositions=[]
			# need to know possible nucleotides to substitute my position with error
			possibleNucs=list(set(self.__NUCLEOTIDES)-set([indNucStrand1]+[indNucStrand2]))
			# I have some positions (at least 1) that is an error
			# errorD= array with coded error nucleotides that will be modified
			if (errorD>0):
				errorChoices=np.random.choice(possibleNucs, int(errorD), replace=NOError)
				maxAvailablePositions=posCoverage
				if not ((posCoverage % 2) == 0): maxAvailablePositions=posCoverage-1
				if maxAvailablePositions == 0:  maxAvailablePositions=posCoverage

				errorPositions=np.random.randint(maxAvailablePositions,size=int(errorD))
				for item in range(0,len(errorPositions)):
					posToChange=errorPositions[item]
					finalRC[posToChange]=errorChoices[item]
			counter=Counter(finalRC)
			# if any of the nucleotides does not have a counter retrieving it will be 0
			ADWErrors[0,indexVar]=counter["A"];ADWErrors[1,indexVar]=counter["C"]
			ADWErrors[2,indexVar]=counter["G"];ADWErrors[3,indexVar]=counter["T"]
		return ADNOErrors,ADWErrors


	def formatIndividualDataForVCF(self,ref,alt,variableSitesPositionIndices,HT,HL,AD,DP):
		"""
		gets single individual data formated as a GT:GL:AD:DP VCF column
		------------------------------------------------------------------------
		intput:
			reference sequence
			alternative alleles
			variableSitesPositionIndices
			haplotype
			likelihood
			allelic depth
			read depth
		output:
			string column with information in format: GT:GL:AD:DP
		"""
		self.appLogger.debug("formatIndividualDataForVCF(self,ref,alt,variableSitesPositionIndices,HT,HL,AD,DP)")
		nVariants=len(variableSitesPositionIndices)
		nInds=len(HT.keys())
		allVariants=dict()
		possibleGenotypes=None
		if self.settings.ploidy==2:
			possibleGenotypes=self.getPossibleGenotypesPerVariableSite(ref,alt,variableSitesPositionIndices)

		for indexVAR in variableSitesPositionIndices:
			allVariants[str(indexVAR)]=[]
		# Had an error here, because the alt variable was empty (passed wrong variable name to the function)
		for indexVAR in range(0,nVariants):
			tmpInd=variableSitesPositionIndices[indexVAR]
			for indexIND in range(0,nInds):
				trueRows=None; htPerInd=None; adPerInd=None;hlPerInd=None;
				if self.settings.ploidy==1:
					valuesToCodify=[ref[tmpInd]]+alt[str(tmpInd)]
					trueRows=self.codifySequences(valuesToCodify)
					htPerInd="{}".format(HT[str(indexIND)][str(tmpInd)])
					hlPerInd=",".join(HL[str(indexIND)][trueRows,indexVAR].astype(dtype=int).astype(dtype=np.str))
					adPerInd=",".join(AD[str(indexIND)][trueRows,indexVAR].astype(dtype=int).astype(dtype=np.str))
				if self.settings.ploidy==2:
					valuesToCodify=[ref[tmpInd]]+alt[str(tmpInd)]
					trueRows=["".join(pair) for pair in self.genotypeOrder(valuesToCodify)]
					htVar=HT[str(indexIND)][str(tmpInd)]
					htPerInd="/".join([str(val) for val in htVar])
					hlPerInd=",".join(\
						[str(HL[str(indexIND)][str(tmpInd)][var]) for var in trueRows])
					alleles=self.codifySequences(valuesToCodify)
					adPerInd=",".join(AD[str(indexIND)][alleles,indexVAR].astype(dtype=int).astype(dtype=np.str))

				ind="{0}:{1}:{2}:{3}".format(\
					htPerInd,\
					hlPerInd,\
					adPerInd,\
					DP[indexIND][indexVAR])
				allVariants[str(tmpInd)]+=[ind]

		return allVariants

	def codifySequences(self,seq):
		"""
		codify nucleotidic sequence in numbers
		A=0,C=1,G=2,T=3
		------------------------------------------------------------------------
		input:
			seq: nucleotidic sequence (containing only ACGT)
		output:
			codified nucleotidic sequence in numbers
		"""
		# self.appLogger.debug("codifySequences(self,ref)")
		codedRef=[]
		for item in seq:
			if "A"==item: codedRef+=[0]
			if "C"==item: codedRef+=[1]
			if "G"==item: codedRef+=[2]
			if "T"==item: codedRef+=[3]
		return np.array(codedRef)

	def writeVCFFile(self, indexREP,indexGT,refAllelesFilePath,REF,alt,variableSitesPositionIndices,HT,HL,AD,DP,flag):
		"""
		gets single individual data formated as a GT:GL:AD:DP VCF column
		------------------------------------------------------------------------
		intput:
			indexREP
			indexGT
			refAllelesFilePath
			reference sequence
			alternative alleles
			variableSitesPositionIndices
			haplotype
			likelihood
			allelic depth
			read depth
			flag: to indicate whether is true o sampled what's being written
		output:
			string column with information in format: GT:GL:AD:DP
		"""
		self.appLogger.debug("writeVCFFile(self, indexREP,indexGT,REF,alt,variableSites,HT,HL,AD,DP,flag)")
		# flag is either true or sampled
		nInds=len(HT.keys())
		if flag:
			self.appLogger.debug("Writing VCF file (true)")
		else:
			self.appLogger.debug("Writing VCF file (sampled)")
		header="{0}\n{1}={2}\n{3}\n{4}={5}".format(\
			"##fileformat=VCFv4.0",\
			"##fileDate",\
			datetime.datetime.now(),\
			"##source=ngsphy",
			"##reference",\
			refAllelesFilePath
		)
		formatLines="{0}\n{1}\n{2}\n{3}".format(
			"##FORMAT=<ID=GT,Number=1,Type=Integer,Description=\"Genotype\">",\
			"##FORMAT=<ID=GL,Number=1,Type=Integer,Description=\"Log10 scale genotype likelihood\">",\
			"##FORMAT=<ID=AD,Number=1,Type=Integer,Description=\"Allelic Depth\">",\
			"##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Depth of coverage\">"
		)

		indnames=["Ind{0}".format(i) for i in range(0,nInds)]
		# filename, file_extension = os.path.splitext('/path/to/somefile.ext')
		headerCols=["#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT"]+indnames
		# CHROM
		numGeneTreeDigits=len(str(self.numLociPerReplicates[(indexREP-1)]))
		chromName="REP.{0:0{1}d}.GT.{2:0{3}d}".format(indexREP,\
			self.numReplicateDigits,\
			indexGT,\
			numGeneTreeDigits)
		# POS
		nVariants=len(variableSitesPositionIndices)
		POS=[(item+1) for item in variableSitesPositionIndices]
		# ID
		ID=["REP.{0:0{1}d}.GT.{2:0{3}d}.ID.{4}".format(\
			indexREP,\
			self.numReplicateDigits,\
			indexGT,
			numGeneTreeDigits,\
			ID) for ID in range(1, (nVariants+1))]
		# ALT
		ALT=[ ",".join(alt[str(pos)]) for pos in variableSitesPositionIndices ]
		# qual
		QUAL=[u'\u00B7']*nVariants
		# FILTER
		FILTER=[u'\u00B7']*nVariants
		# INFO
		INFO=[u'\u00B7']*nVariants
		# format
		FORMAT=["GT:GL:AD:DP"]*nVariants
		# extra 9 columns: #CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO,FORMAT = 9
		nLoci=self.numLociPerReplicates[indexREP-1]
		numGeneTreeDigits=len(str(nLoci))
		allVariants=self.formatIndividualDataForVCF(\
			REF,alt,variableSitesPositionIndices,HT,HL,AD,DP)
		outfile=""
		# flag true=true  - Flag false= sampled
		if flag:
			outfile=os.path.join(
				self.readcountNoErrorFolderPath,\
				"{0:0{1}d}".format(indexREP,self.numReplicateDigits),\
				"{0}_{1:0{2}d}_{3:0{4}d}_NOERROR.vcf".format(\
					self.settings.simphyDataPrefix,\
					indexREP,self.numReplicateDigits,\
					indexGT,numGeneTreeDigits\
				)\
			)
		else:
			outfile=os.path.join(
				self.readcountNoErrorFolderPath,\
				"{0:0{1}d}".format(indexREP,self.numReplicateDigits),\
				"{0}_{1:0{2}d}_{3:0{4}d}.vcf".format(\
					self.settings.simphyDataPrefix,\
					indexREP,self.numReplicateDigits,\
					indexGT,numGeneTreeDigits\
				)\
			)


		# before writing i'm getting max width of the lines written per column
		colWidths=self.getColWidhts(\
			chromName,POS,ID,\
			REF,alt,QUAL,FILTER,\
			INFO,FORMAT,\
			allVariants,variableSitesPositionIndices)
		# (sizeChrom,sizePOS,sizeID,sizeREF,sizeALT,sizeQUAL,sizeFILTER,sizeINFO,sizeFORMAT,sizeInds)
		maxLenIndName=max([len(elem) for elem in indnames])
		maxLenIndName=max(colWidths[9], maxLenIndName)
		headerWidths=[maxLenIndName]*len(indnames)
		headerWidths=colWidths+headerWidths
		headerFields=["{0:{1}s}".format(\
			headerCols[indexField],headerWidths[indexField])\
			for indexField in range(0,len(headerCols))]
		filevcf=open(outfile, 'a')
		filevcf.write("{0}\n{1}\n{2}\n".format(\
			header,\
			formatLines,\
			"\t".join(headerFields)\
			))

		for index in range(0, nVariants):
			# extra 9 columns: #CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO,FORMAT = 9
			line="{0:{1}s}\t{2:{3}s}\t{4:{5}s}\t{6:{7}s}\t{8:{9}s}\t{10:{11}s}\t{12:{13}s}\t{14:{15}s}\t{16:{17}s}\t{18}\n".format(\
				chromName,colWidths[0],\
				str(POS[index]),colWidths[1],\
				ID[index],colWidths[2],\
				REF[variableSitesPositionIndices[index]],colWidths[3],\
				",".join(alt[str(variableSitesPositionIndices[index])]),colWidths[4],\
				QUAL[index].encode("UTF-8"),colWidths[5],\
				FILTER[index].encode("UTF-8"),colWidths[6],\
				INFO[index].encode("UTF-8"),colWidths[7],\
				FORMAT[index].encode("UTF-8"),colWidths[8],\
				"\t".join(\
					["{0:{1}s}".format(indVar,maxLenIndName) for indVar in allVariants[str(variableSitesPositionIndices[index])]]\
				)
			)
			filevcf.write(line)
		filevcf.close()

	def getColWidhts(self,chromName,POS,ID,REF,ALT,QUAL,FILTER,INFO,FORMAT,allVariants, variableSitesPositionIndices):
		"""
		get string widht of the columns to format the VCF output
		------------------------------------------------------------------------
		Parameters:
		- chromName: data column
		- POS: data column
		- ID: data column
		- REF: data column
		- ALT: data column
		- QUAL: data column
		- FILTER: data column
		- INFO: data column
		- FORMAT: data column
		- allVariants: data column
		- variableSitesPositionIndices: data column
		Returns:
		- list with the leghts of the columns that are going to be written in the VCF file
		"""
		self.appLogger.debug("getColWidhts(self,chromName,POS,ID,REF,ALT,QUAL,FILTER,INFO,FORMAT,allVariants, variableSites)")
		#CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO,FORMAT
		sizeChrom=len(chromName)
		sizePOS=max([len(str(item)) for item in POS])
		sizeID=max([len(item) for item in ID])
		sizeREF=max([len(item) for item in REF])
		tmpALT=[",".join(ALT[str(var)]) for var in variableSitesPositionIndices]
		sizeALT=max([len(item) for item in tmpALT])
		sizeQUAL=max([len(item) for item in QUAL])
		sizeFILTER=max([len(item) for item in FILTER])
		sizeINFO=max([len(item) for item in INFO])
		sizeFORMAT=max([len(item) for item in FORMAT])
		tmpInds=[]
		for item in variableSitesPositionIndices:
			tmpInds+=[max([ len(elem) for elem in allVariants[str(item)]])]
		sizeInds=max(tmpInds)
		return [sizeChrom,sizePOS,sizeID,sizeREF,sizeALT,sizeQUAL,sizeFILTER,sizeINFO,sizeFORMAT,sizeInds]

	def writeReference(self,indexREP,indexGT,referenceSpeciesID,referenceLocusID,referenceTipID,referenceSeqFull):
		"""
		write reference sequence file separately
		------------------------------------------------------------------------
		Parameters:
		- indexREP: index of the species tree to which it belongs
		- indexGT: index of the gene tree to which it belongs
		- referenceSpeciesID: species ID to which the reference belongs within a MSA file
		- referenceTipID:  tip ID to which the reference belongs within a species in a MSA file
		- referenceSeqFull: nucleotidic sequence
		Returns:
		- filepath where the reference will be written.
		"""
		self.appLogger.debug(" writeReference(self,indexREP,indexGT,referenceSpeciesID,referenceTipID,referenceSeqFull):")
		refAllelesFilePath=os.path.join(\
			self.refAllelesFolderPath,\
			"{0:0{1}d}".format(indexREP,self.numReplicateDigits),\
			"{0}_REF_{1}_{2}.fasta".format(\
				self.settings.projectName,\
				indexREP,
				indexGT\
			)\
		)
		description=">REFERENCE:{0}:ST.{1:{2}}:GT.{3:{4}}:{5}_{6}_{7}".format(\
			self.settings.simphyDataPrefix,\
			indexREP,\
			self.numReplicateDigits,\
			self.numLociPerReplicates[(indexREP-1)],\
			len(str(self.numLociPerReplicates[(indexREP-1)])),\
			referenceSpeciesID,\
			referenceLocusID,\
			referenceTipID\
		)
		f=open(refAllelesFilePath,"w")
		f.write("{0}\n{1}\n".format(\
			description,\
			"".join(referenceSeqFull)
		))
		return refAllelesFilePath



	def launchCommand(self, referenceForCurrST, indexREP,indexGT, individuals,coverageMatrix):
		"""
		Launches the execution ogf the RC computation.
		------------------------------------------------------------------------
		Parameters:
		- referenceForCurrST: reference for current species tree replicate, the
		 one that  will be processed
		- indexREP: proper identifier of the species tree replicete
		- indexGT: proper identifier of the gene tree
		- individuals: number of individuals
		- coverageMatrix: coverage matrix
		"""
		try:
			self.appLogger.debug("launchCommand(...)")
			tStartTime=datetime.datetime.now()
			numIndividuals=len(individuals.keys())
			nLoci=self.numLociPerReplicates[indexREP-1]
			self.appLogger.debug("Iterating over nLoci: {0}/{1}".format(indexGT,nLoci))
			numGeneTreeDigits=len(str(nLoci))
			filepathLoc=os.path.join(\
				self.settings.path,\
				self.settings.projectName,\
				"{0:0{1}d}".format(indexREP,self.numReplicateDigits),\
				"{0}_{1:0{2}d}_TRUE.fasta".format(\
					self.settings.simphyDataPrefix,\
					indexGT,\
					numGeneTreeDigits\
				)\
			)

			# dictionary. indices=variable sites, content=variable nucleotide set
			variableSites=self.extractNOErrorVariantsPositions(filepathLoc)
			nVarSites=len(variableSites.keys())
			self.appLogger.debug(\
				"Species tree: {0}\t Locus: {1}\t - Found [{2}] variable sites.".format(\
					indexREP,\
					indexGT,\
					nVarSites))
			# Parse MSA files - dictionary[speciesID][tipID]
			#				   content={description,sequence}
			msa=parseMSAFile(filepathLoc)
			# indices to get the sequence of the REFERENCE!
			tag="{0}_{1}".format(referenceForCurrST[1],referenceForCurrST[2])
			if tag in msa.keys():
				referenceSpeciesID=str(referenceForCurrST[1])
				referenceLocusID=str(referenceForCurrST[2])
				referenceTipID=str(referenceForCurrST[3])
				referenceSeqFull=msa[tag][referenceTipID]['sequence']
				refAllelesFilePath=self.writeReference(\
					indexREP,indexGT,\
					referenceSpeciesID,referenceLocusID,referenceTipID,\
					referenceSeqFull)
				# Coverage per SNP is the same for both true/sampled dataset
				self.appLogger.info("{0} {1}\t Locus: {2}\t - {3}".format(\
					"Replicate:",\
					indexREP,\
					indexGT,\
					"Getting coverage per individual per variant"\
				))
				DP=[self.getDepthCoveragePerIndividual(\
						nVarSites,\
						coverageMatrix[index][(indexGT-1)]) \
						for index in range(0,numIndividuals)\
						]
				if self.settings.ploidy==1:
					self.computeHaploid(\
						indexREP,indexGT,\
						msa,individuals,\
						refAllelesFilePath,referenceSeqFull,\
						variableSites,DP)
				if self.settings.ploidy==2:
					self.computeDiploid(\
						indexREP,indexGT,\
						msa,individuals,\
						refAllelesFilePath,referenceSeqFull,\
						variableSites,DP)
				tEndTime=datetime.datetime.now()
				ETA=(tEndTime-tStartTime).total_seconds()
				# row['indexREP'],indexLOC,row['indID'],inputFile, outputFile]+callParams]
				if(self.settings.runningTimes):
					outfile=os.path.join(
						self.readcountNoErrorFolderPath,\
						"{0:0{1}d}".format(indexREP,self.numReplicateDigits),\
						"{0}_{1:0{2}d}_{3:0{4}d}*".format(
							self.settings.simphyDataPrefix,\
							indexREP,\
							self.numReplicateDigits,\
							indexGT,\
							numGeneTreeDigits\
						)\
					)
					line=[indexREP,indexGT,"-",filepathLoc,ETA,"-",outfile]
	 				self.runningInfo.addLine(line)
		except KeyboardInterrupt:
			message="{0}{1}Thread has been interrupted!{2}\n".format("\033[91m","\033[1m","\033[0m")
			raise KeyboardInterrupt(message)

	def run(self):
		"""
		Process flow for the generation of read counts
		------------------------------------------------------------------------
		Returns:
		- Boolean indating whether process has finished correctly or not.
		"""
		self.appLogger.debug("run(self)")
		status=NOError;	message="Read counts finished ok."
		self.appLogger.debug( "Run - read count")
		try:
			# generating folder structure for this part
			self.generateFolderStructureDetail()
			# Get list of reference sequences
			referenceList=self.parseReferenceAllelesList(self.settings.readCountsReferenceAllelesFile)
			nSTS=len(self.filteredReplicates)
			# iterate over the "iterable" species trees / filtered STs
			self.appLogger.info("Running...")
			nProcesses=sum([self.numLociPerReplicates[item-1] for item in self.filteredReplicates])
			currProcessesRunning=1
			pool=multiprocessing.Pool(self.settings.numThreads)

			filepathIndividualsRelations=[\
				"{0}/tables/{1}.{2:0{3}d}.individuals.csv".format(\
				self.settings.outputFolderPath,\
				self.settings.projectName,\
				indexREP,\
				self.numReplicateDigits)\
			for indexREP in self.filteredReplicates]

			individualsList=[self.parseIndividualRelationFile(singlefile) for singlefile in filepathIndividualsRelations]
			numIndividualsList=[len(item.keys()) for item in individualsList]
			coverageMatrices=[self.retrieveCoverageMatrix(indexREP) for indexREP in self.filteredReplicates]
			nLociList=[self.numLociPerReplicates[item-1] for item in self.filteredReplicates]
			ARGS=[(referenceList[item],\
				self.filteredReplicates[item],\
				nLociList[item],\
				individualsList[item],\
				coverageMatrices[item])\
			for item in range(0,nSTS)]

			for indexFilterST in range(0,nSTS):
				indexREP=self.filteredReplicates[indexFilterST] # Get Proper ID for ST
				nLoci=nLociList[indexFilterST]
				self.appLogger.debug("Number of individuals: {0}".format(numIndividualsList[indexFilterST]))
				self.appLogger.debug("Number of loci: {0}".format(nLoci))
				indexGT=1
				threadsToBeRan=nLoci
				while indexGT < (nLoci+1):
					# for indexGT in range(1,(nLoci+1)):
					jobs=[]
					for item in range(0,min(self.settings.numThreads, threadsToBeRan)):
						progress=((currProcessesRunning*100)*1.0)/(nProcesses+1)
						currProcessesRunning+=1
						sys.stdout.write("Progress {0:02.1f} %\r".format(progress))
						sys.stdout.flush()
						time.sleep(0.2)
						t = multiprocessing.Process(\
							target=self.launchCommand,\
							args=(\
								referenceList[indexREP-1],\
								indexREP,\
								indexGT,\
								individualsList[indexFilterST],\
								coverageMatrices[indexFilterST])\
						)
						t.daemon=NOError
						jobs.append(t)
						t.start()
						indexGT+=1
						threadsToBeRan=(nLoci+1)-indexGT

					for j in jobs:
						j.join()

		except ValueError as ve:
			# If there's a wrong ploidy inserted.
			status=False
			message="\n\t{0}\n\t{1}\n\t{2}".format(\
				"Ploidy value invalid",\
				ve,\
				"Please verify. Exiting.")
		except TypeError as te:
			# One of the files is not fasta.
			status=False
			message=message="\n\t{0}\n\t{1}\n\t{2}".format(\
				"Problem with the file type.",\
				te,\
				"Please verify. Exiting.")
		except Exception as ex:
			exc_type, exc_obj, exc_tb = sys.exc_info()
			fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
			message="\n\tUnexpected: {0} | {1} - File: {2} - Line:{3}".format(\
				ex,exc_type, fname, exc_tb.tb_lineno)
			status=False
		return status,message
