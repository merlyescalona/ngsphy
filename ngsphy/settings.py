import argparse,datetime,dendropy, logging,os,re,sys
import numpy as np
from coverage import NGSPhyDistributionParser as ngsphydistro
if (sys.version_info[0:2]<(3,0)):
	import ConfigParser as cp
elif (sys.version_info>=(3,0)):
	import configparser as cp

class Settings:
	"""This class parses and verifies that the settings file is correct."""
	# General
	originMode=1
	ploidy=1
	projectName=""
	path=""
	outputFolderName=""
	outputFolderPath=""

	dataPrefix="ngsphydata"
	# simphy data origin | origin 1
	simphyProjectPath=""
	filterSimphy=False
	# indelible data origin | origin 2
	ngsphyIndelibleControlFilePath=""
	newickFilePath=""
	evolve=""
	partition=""
	# indelible data origin + referebce | origin 3
	referenceSequenceFilePath=""
	referenceTipLabel=""

	#readcount
	seqerror=0
	readCountsReferenceFilePath=""

	# ngsart
	ngsmode=0

	# coverage
	ontarget=1
	offtarget=0
	notcaptured=0
	experiment=None
	individual=None
	locus=None
	genomicNoise=None
	phylogeneticDecay=dict()

	numSpeciesTrees=0
	numLociPerSpeciesTree=[]
	outgroup=False
	indels=False

	def __init__(self,filename):
		# If I've got this far, then filename is a correct file
		self.settingsFile=os.path.abspath(filename)
		self.appLogger=logging.getLogger('ngsphy')
		self.appLogger.debug("(class Settings) __init__()")
		# default settings can be established.
		self.parser=cp.SafeConfigParser()
		self.parser.read(self.settingsFile)

	def checkArgs(self):
		self.appLogger.debug("Check argumentss")
		allGood=True
		parserMessageCorrect="All parameters are correct."
		parserMessageWrong="Settings - Problem found! "
		# Checking general section
		statusGeneral,messageGeneral= self.checkSectionGeneral(\
			parserMessageCorrect,parserMessageWrong)

		# Exit here
		if not statusGeneral: return statusGeneral, messageGeneral
		# Checking data section
		statusData,messageData=self.checkSectionData(\
			parserMessageCorrect,parserMessageWrong)
		# Exit here
		if not statusData: return statusData, messageData
		# ----------------------------------------------------------------------
		# NGS sections
		if self.parser.has_section("ngs-reads-art"):
			self.ngsmode=1
			statusNGSArt,messageNGSArt=self.checkSectionNGSReadsArt(\
				parserMessageCorrect,parserMessageWrong)
			if (statusNGSArt):
				## Next check
					statusCoverage,messageCoverage=self.checkSectionCoverage(parserMessageCorrect,parserMessageWrong)
					# Exit here
					if not statusCoverage: return statusCoverage, messageCoverage
			else:
				# Exit here
				return statusNGSArt,messageNGSArt
			# if all parameters for NGS ART are correct I don't need the other NGS section
			# removing NGS read count
			if self.parser.has_section("ngs-read-count"):
				self.parser.remove_section("ngs-read-count")
				self.appLogger.warning("[ngs-read-count] section is incompatible with [ngs-reads-art]. Omiting this section.")
		elif self.parser.has_section("ngs-read-count"):
			self.ngsmode=2
			# readcount
			statusRC,messageRC=self.checkSectionReadCount(parserMessageCorrect,parserMessageWrong)
			if statusRC:
				statusCoverage,messageCoverage=self.checkSectionCoverage(parserMessageCorrect,parserMessageWrong)
				# Exit here
				if not statusCoverage: return statusCoverage, messageCoverage
			else:
				# Exit here
				return statusRC,messageRC
		else: # Meaning have no ngs-reads-art nor ngs-read-counts
			if (self.parser.has_section("coverage")):self.parser.remove_section("coverage")
			self.appLogger.info("[coverage] section is not needed if [ngs-reads-art] or [ngs-read-count] are not available. Omiting this section.")

		# Checking execution section and options
		self.checkSectionExecution(parserMessageCorrect,parserMessageWrong)
		# Exit here
		self.appLogger.info(self.formatSettingsMessage())
		return allGood,parserMessageCorrect

	def checkSectionGeneral(self,parserMessageCorrect,parserMessageWrong):
		self.appLogger.debug("Section General")
		# Check GENERAL SECTION
		if not self.parser.has_section("general"):
			parserMessageWrong+="\n\t{0}\n\t{1}".format(\
				"[general] section missing and required.",\
				"Please verify. Exiting."\
			)
			return False, parserMessageWrong
		# CHECKING GENERAL PARAMETERS

		# checking ploidy for the output data
		if (not self.parser.has_option("general","ploidy")):
			self.ploidy=1
		else:
			p=self.parser.getint("general","ploidy")
			if (p>0 and p<=2):  self.ploidy=p
			elif (p<0): self.ploidy=1
			else:   self.ploidy=2

		# project information
		if (self.parser.has_option("general","project_name")):
			self.projectName=self.parser.get("general","project_name")
		else:
			parserMessageWrong+="\n\t Project Name option is needed. Please verify. Exiting."
			return False, parserMessageWrong

		if (self.parser.has_option("general","path")):
			self.path=os.path.abspath(self.parser.get("general","path"))
		else:
			parserMessageWrong+="\n\t Working path option is needed. Please verify. Exiting."
			return False, parserMessageWrong

		# Checking output folder information
		if (self.parser.has_option("general","ofn")):
			self.outputFolderName=self.parser.get("general","ofn")
			self.parser.set("general","output_folder_name",self.outputFolderName)
			self.parser.remove_option("general","ofn")
		if(self.parser.has_option("general","output_folder_name")):
			self.outputFolderName=self.parser.get("general","output_folder_name")
		else:
			self.outputFolderName="output"

		try:
			os.makedirs("{0}/{1}".format(self.path,self.projectName))
		except:
			origin=0
			if self.parser.has_option("data","origin"):
				origin=self.parser.get("data","origin")
				if (origin >1): return False, "\n\t{0}\n\t{1}\n\t{2}\n\t{3}".format(\
					"Project folder exists and origin mode is in [2,3].",\
					"Features are incompatible.",\
					"A new run would overwrite the sequence data.",\
					"Please verify. Exiting."\
				)
			self.appLogger.warning("Project folder exists. New output folders will be created.")

		if not os.path.exists("{0}/{1}/".format(self.path,self.projectName)):
			parserMessageWrong+="\n\tProject path: {0} does not exist.\n\tPlease verify. Exiting.".format(\
				"{0}/{1}/".format(self.path,self.projectName))
			return False, parserMessageWrong
		if os.path.exists("{0}/{1}/{2}/".format(self.path,self.projectName,self.outputFolderName)):
			listdir=os.listdir("{0}/{1}".format(self.path,self.projectName))
			counter=0
			for item in listdir:
				if self.outputFolderName in item:
					counter+=1
			if not counter == 0: self.outputFolderName="output_{0}".format(counter+1)

		self.outputFolderPath=os.path.abspath("{0}/{1}/{2}/".format(self.path,self.projectName,self.outputFolderName))
		self.parser.set("general","output_folder_name",self.outputFolderPath)
		return True,parserMessageCorrect

	def checkSectionData(self,parserMessageCorrect,parserMessageWrong):
		self.appLogger.debug("Section Data")
		if (self.parser.has_section("data")):
			if (self.parser.has_option("data","origin")):
				try:
					self.originMode=self.parser.getint("data","origin")
				except:
					parserMessageWrong+="\n\t\n\t{0}\n\t{1}\n\t{2}".format(\
						"[data] section: Value introduced in the origin option is not valid.",\
						"Value should be an integer.",\
						"Please verify. Exiting."\
					)
					return False, parserMessageWrong
			else:
				parserMessageWrong+="\n\t\n\t{0}\n\t{1}\n\t{2}".format(\
					"[data] section: origin option value is missing. This option is required.",\
					"Please verify. Exiting."\
				)
				return False, parserMessageWrong

			####################################################################
			if not self.originMode in [1,2,3]:
				parserMessageWrong+="\n\t\n\t{0}\n\t{1}".format(\
					"[data] section: origin option value is invalid",\
					"Value is out of range. Possible range [1,3]",\
					"Please verify. Exiting."\
				)
				return False, parserMessageWrong
			if self.originMode==1:
				self.simphyProjectPath=os.path.join(self.path,self.projectName)
				if (os.path.exists(self.simphyProjectPath) and os.path.isdir(self.simphyProjectPath)):
					self.appLogger.debug("SimPhy project folder exists")
				else:
					parserMessageWrong+="\n\t\n\t{0}\n\t{1}\n\t{2}".format(\
						"[data] section: Origin mode (1) selected but invalid option.",\
						"SimPhy project folder does not exist, or the given path does not belong to a directory.",\
						"Please verify. Exiting."\
					)
					return False, parserMessageWrong
				# data prefix for the simphy sequences
				if (self.parser.has_option("data","data_prefix") or self.parser.has_option("data","dp")):
					if (self.parser.has_option("data","dp")):
						value=self.parser.get("data","dp")
						self.parser.set("data","data_prefix",value)
						self.parser.remove_option("data","dp")
					if (self.parser.has_option("data","data_prefix")):
						self.dataPrefix=self.parser.get("data","data_prefix")
				else:
					parserMessageWrong+="\n\t\n\t{0}\n\t{1}\n\t{2}".format(\
						"[data] section: Origin mode (1) selected but invalid option.",\
						"<data_prefix | dp> field is missing. This prefix correponds to the name of the file sequences that are going to be processed.",\
						"Please verify. Exiting."\
					)
					return False, parserMessageWrong
				# checking flag for filtering species tree replicates
				if self.parser.has_option("data","filter_simphy"):
					self.filterSimphy=True
				else:
					self.filterSimphy=False
				# removing options that do not match with the origin mode selected
				if (self.parser.has_option("data","indelible_control")):
					self.parser.remove_option("data","indelible_control")
				if (self.parser.has_option("data","newick_file")):
					self.parser.remove_option("data","newick_file")
				# --------------------------------------------------------------
				# check if simphy is valid projectName
				status,message=self.checkSimPhyProjectValid()
				if not status:
					return status, message
			####################################################################
			if self.originMode in [2,3]:
					# parameter is set up, now check if folder exist
				if (self.parser.has_option("data","indelible_control")):
					self.ngsphyIndelibleControlFilePath=os.path.abspath(\
						self.parser.get("data","indelible_control")
					)

					# checkin indelible_control file
					if (os.path.exists(self.ngsphyIndelibleControlFilePath) and os.path.isfile(self.ngsphyIndelibleControlFilePath)):
						self.appLogger.debug("INDELible control file exists")
					else:
						parserMessageWrong+="\n\t{0}{1}{2}\n\t{3}\n\t{4}".format(\
							"[data] section: Origin mode (",self.originMode,") selected but invalid option.",\
							"INDELible control file does not exist.",\
							"Please verify. Exiting."\
						)
						return False, parserMessageWrong


					# checkin tree file
					if (self.parser.has_option("data","newick_file")):
						self.appLogger.debug("Newick file option exists")
						self.newickFilePath=os.path.abspath(self.parser.get("data","newick_file").strip())
						filesOk=(os.path.exists(self.newickFilePath) and os.path.isfile(self.newickFilePath))
						if not filesOk:
							parserMessageWrong+="\n\t{0}{1}{2}\n\t{3}\n\t{4}".format(\
								"[data] section: Origin mode (",self.originMode,") selected but invalid option.",\
								"NEWICK file does not exist.",\
								"Please verify. Exiting."\
							)
							return False, parserMessageWrong
					else:
						parserMessageWrong+="\n\t{0}{1}{2}\n\t{3}\n\t{4}".format(\
							"[data] section: Origin mode (",self.originMode,") selected but invalid option.",\
							"INDELible control file does not exist.",\
							"Please verify. Exiting."\
						)
						return False, parserMessageWrong


					self.appLogger.debug("Checking INDELible control file")
					# checking control file format
					statusOk, message=self.checkIndelibleControlFile(parserMessageCorrect,parserMessageWrong)
					if not statusOk: return statusOk,message
					# checking tip label format
					self.appLogger.debug("Checking tree label formatting")
					statusOk, message=self.checkLabelFormatInTree()
					if not statusOk: return statusOk,message
					# checking indelible program
					stream = os.popen('which indelible').read()[0:-1]
					self.appLogger.info("Checking dependencies...")
					if stream:
						self.appLogger.info("indelible - Found running in: {}".format(stream))
					else:
						parserMessageWrong+="\n\t{0}{1}{2}\n\t{3}\n\t{4}".format(\
							"[data] section: Origin mode (",self.originMode,") selected but invalid option.",\
							"indelible not found. Program either not installed or not in the current path.",\
							"Please verify the installation. Exiting."\
						)
						return False, parserMessageWrong
					self.parser.set("general","numspeciestrees",str(1))
					self.numSpeciesTrees=1
					# removing options that do not match with the origin mode selected
					self.parser.set("data","data_prefix",self.dataPrefix)
					if self.parser.has_option("data","dp"): self.parser.remove("data","dp")
					if self.parser.has_option("data","filter_simphy"):  self.parser.remove("data","filter_simphy")
					if self.originMode ==2:
						if self.parser.has_option("data","reference_sequence_file"):  self.parser.remove_option("data","reference_sequence_file")
						if self.parser.has_option("data","reference_tip_label"):  self.parser.remove_option("data","reference_tip_label")
				else:
					parserMessageWrong+="\n\t{0}{1}{2}\n\t{3}\n\t{4}".format(\
						"[data] section: Origin mode (",self.originMode,") selected but invalid option.",\
						"INDELible  control file is required",\
						"Please verify. Exiting."\
					)
					return False, parserMessageWrong
			####################################################################
			if self.originMode==3:
				# INDELIBLE + tree + reference + tiplabel
				# referenceSequenceFilePath
				# referenceTipLabel
				if (self.parser.has_option("data","reference_sequence_file")):
					self.referenceSequenceFilePath=os.path.abspath(self.parser.get("data","reference_sequence_file").strip())
					filesOk=(os.path.exists(self.referenceSequenceFilePath) and os.path.isfile(self.referenceSequenceFilePath))
					if not filesOk:
						parserMessageWrong+="\n\t{0}{1}{2}\n\t{3}\n\t{4}".format(
							"[data] section: Origin mode (",self.originMode,") selected but invalid option.",\
							"Reference sequence file does not exist or the given path does not belong to a file.",\
							"Please verify. Exiting."\
						)
						return False, parserMessageWrong
				else:
					parserMessageWrong+="\n\t{0}{1}{2}\n\t{3}\n\t{4}".format(\
						"[data] section: Origin mode (",self.originMode,") selected but invalid option.",\
						"Reference sequence file is required",\
						"Please verify. Exiting."\
					)
					return False, parserMessageWrong

				if (self.parser.has_option("data","reference_tip_label")):
					self.referenceTipLabel=self.parser.get("data","reference_tip_label")
					# existence of tip in the given tree
					statusOk, message=self.checkReferenceTipLabelInTree()
					if not statusOk: return statusOk,message
				else:
					parserMessageWrong+="\n\t{0}{1}{2}\n\t{3}\n\t{4}".format(\
						"[data] section: Origin mode (",self.originMode,") selected but invalid option.",\
						"Reference sequence tip label is required",\
						"Please verify. Exiting."\
					)
					return False, parserMessageWrong
		else:
			parserMessageWrong+="\n\t\n\t{0}\n\t{1}\n\t{2}".format(\
				"[data] section is missing. This section is required.",\
				"Please verify. Exiting."\
			)
			return False, parserMessageWrong
		return True, parserMessageCorrect



	def checkSectionNGSReadsArt(self,parserMessageCorrect,parserMessageWrong):
		########################################################################
		# BLOCK: NGS-READS-ART
		########################################################################
		# Checking art parameters.
		self.appLogger.info("NGS-reads-ART option selected.")
		# checking program dependencies
		stream = os.popen('which art_illumina').read()[0:-1]
		self.appLogger.info("Checking dependencies...")
		if stream:
			self.appLogger.info("art_illumina - Found running in: {}".format(stream))
			# IO parameters
			if self.parser.has_option("ngs-reads-art","o"):self.parser.remove_option("ngs-reads-art","o")
			if self.parser.has_option("ngs-reads-art","out"):self.parser.remove_option("ngs-reads-art","out")
			if self.parser.has_option("ngs-reads-art","i"):self.parser.remove_option("ngs-reads-art","i")
			if self.parser.has_option("ngs-reads-art","in"):self.parser.remove_option("ngs-reads-art","in")
			self.appLogger.warning("Removing I/O options. Be aware: I/O naming is auto-generated.")
			# Coverage parameters
			if (self.parser.has_option("ngs-reads-art","fcov")): self.parser.remove_option("ngs-reads-art","fcov")
			if (self.parser.has_option("ngs-reads-art","f")): self.parser.remove_option("ngs-reads-art","f")
			if (self.parser.has_option("ngs-reads-art","rcount")): self.parser.remove_option("ngs-reads-art","rcount")
			if (self.parser.has_option("ngs-reads-art","c")): self.parser.remove_option("ngs-reads-art","c")
			self.appLogger.warning("Removing ART coverage options. Coverage is calculated with the [coverage] section parameters.")
		else:
			parserMessageWrong+="\n\t{0}\n\t{1}".format(\
				"art_illumina not found. Program either not installed or not in your current path.",\
				"Please verify the installation. Exiting."
			)
			return False, parserMessageWrong
		return True, parserMessageCorrect

	def checkSectionReadCount(self,parserMessageCorrect,parserMessageWrong):
		########################################################################
		# BLOCK: READ COUNT
		########################################################################
		message=parserMessageCorrect
		if (self.parser.has_section("ngs-read-count")):
			self.readcount=True
			if not self.parser.has_option("ngs-read-count", "error"):
				self.appLogger.warning("[ngs-read-count] section. Sequencing error rate for this run is being considered as 0.")
				self.parser.set("ngs-read-count", "error","0")
				self.seqerror=0
			else:
				self.seqerror=float(self.parser.get("ngs-read-count", "error"))

			if not self.parser.has_option("ngs-read-count","reference"):
				self.appLogger.warning("[ngs-read-count] section. Using default references.")
				self.parser.set("ngs-read-count", "reference","None")
				self.readCountsReferenceFilePath=None
			else:
				self.readCountsReferenceFilePath=os.path.abspath(self.parser.get("ngs-read-count","reference"))
		else:
			# No ngs-read-count section
			self.readcount=False
			message="[ngs-read-count] section. Not available."

		return self.readcount,message

	########################################################################
	# BLOCK: Coverage
	########################################################################
	def checkSectionCoverage(self,parserMessageCorrect,parserMessageWrong):
		self.appLogger.debug("Checking coverage")
		message=parserMessageCorrect
		expCov=None;indCov=None;locCov=None;
		if(self.parser.has_section("coverage")):
			# ------------------------------------------------------------------
			#  Targeted-sequencing related
			if self.originMode==1:
				self.appLogger.debug("Origin Mode 1")
				if (self.parser.has_option("coverage","ontarget")):
					self.appLogger.debug("Ontarget")
					self.ontarget=self.parser.getfloat("coverage","ontarget")
					if not (self.ontarget <= 1):
						return False, "\n\t{0}\n\t{1}\n\t{2}\n\t{3}".format(\
							parserMessageWrong,\
							"[coverage] section: ontarget option invalid.",\
							"Value out of range, should be within the interval [0,1]",\
							"Please verify. Exiting"\
						)
					if not self.parser.has_option("coverage","offtarget"):
						self.appLogger.debug("Ontarget only. Auto-complete")
						self.offtarget=1-self.ontarget
						self.parser.set("coverage","offtarget", str(self.offtarget))
				if (self.parser.has_option("coverage","offtarget")):
					self.appLogger.debug("offtarget")
					self.offtarget=self.parser.getfloat("coverage","offtarget")
					if not (self.offtarget <= 1):
						return False, "\n\t{0}\n\t{1}\n\t{2}\n\t{3}".format(\
							parserMessageWrong,\
							"[coverage] section: offtarget option invalid.",\
							"Value out of range, should be within the interval [0,1]",\
							"Please verify. Exiting"\
							)
					if not self.ontarget:
						self.appLogger.debug("Offtarget only. Auto-complete")
						self.ontarget=1-self.offtarget
						self.parser.set("coverage","ontarget", str(self.ontarget))

				if ((self.parser.has_option("coverage","ontarget")) and (self.parser.has_option("coverage","offtarget")) and (self.ontarget+self.offtarget > 1)):
					self.appLogger.debug("on/offtarget")
					temp=self.ontarget+self.offtarget
					self.ontarget=self.ontarget/temp
					self.offtarget=self.offtarget/temp
					self.appLogger.warning("[coverage] - on/off target parameters had been normalized.")
					# If got in here, values had been already assigned

				if (self.parser.has_option("coverage","notcaptured")):
					self.appLogger.debug("not captured")
					self.notcaptured=self.parser.getfloat("coverage","notcaptured")
					if not (self.notcaptured <= 1):
						return False, "\n\t{0}\n\t{1}\n\t{2}\n\t{3}".format(\
							parserMessageWrong,\
							"[coverage] section: notcaptured option invalid.",\
							"Value out of range, should be within the half-closed interval [0,1)",\
							"Please verify. Exiting"\
						)
			else:
				if (self.parser.has_option("coverage","ontarget")): self.parser.remove("coverage","ontarget")
				if (self.parser.has_option("coverage","offtarget")): self.parser.remove("coverage","offtarget")
				if (self.parser.has_option("coverage","notcaptured")): self.parser.remove("coverage","notcaptured")
				if (self.parser.has_option("coverage","ontarget") or self.parser.has_option("coverage","offtarget") or self.parser.has_option("coverage","notcaptured")):
					self.appLogger.warning("Coverage options - ontarget, offtarget and notcaptured - only possible with [data] origin = 1 (SimPhy project)")
			# ------------------------------------------------------------------
			#  Coverage distribution related
			if (self.parser.has_option("coverage","experiment")):
				self.appLogger.debug("Coverage experiment")
				value=self.parser.get("coverage","experiment")
				self.parser.set("coverage","experiment",value.lower())
			else:
				# parsear distribution
				parserMessageWrong+="\n\t{0}\n\t{1}".format(\
					"[coverage] section:  experiment option is required.",\
					"Please verify. Exiting."\
				)
				return False,parserMessageWrong

			self.experiment=ngsphydistro(self.parser.get("coverage","experiment"), False)
			check,mess=self.experiment.validate()
			if not (check):
				parserMessageWrong+=mess
				return check,parserMessageWrong
			# GENOMIC NOISE -------------------------------------------------------
			if (self.parser.has_option("coverage","genomic_noise")):
				value=self.parser.get("coverage","genomic_noise")
				self.parser.set("coverage","genomic_noise",value.lower())
				self.genomicNoise=ngsphydistro(self.parser.get("coverage","genomic_noise"), False)
				check,mess=self.genomicNoise.validate()
				if not (check):
					parserMessageWrong+=mess
					return check,parserMessageWrong
			# INDIVIDUAL -------------------------------------------------------
			if (self.parser.has_option("coverage","individual")):
				self.appLogger.debug("Coverage individual")
				value=self.parser.get("coverage","individual")
				self.parser.set("coverage","individual",value.lower())
				self.individual=ngsphydistro(self.parser.get("coverage","individual"), False)
				check,mess=self.individual.validate()
				if not (check):
					parserMessageWrong+=mess
					return check,parserMessageWrong
			# LOCUS ------------------------------------------------------------
			if (self.parser.has_option("coverage","locus")):
				self.appLogger.debug("Coverage locus")
				value=self.parser.get("coverage","locus")
				self.parser.set("coverage","locus",value.lower())
				self.locus=ngsphydistro(self.parser.get("coverage","locus"), False)
				check,mess=self.locus.validate()
				if not (check):
					parserMessageWrong+=mess
					return check,parserMessageWrong
			# ------------------------------------------------------------------
			#  Phylogenetic Decay
			if (self.parser.has_option("coverage","phylogenetic_decay")):
				values=self.parser.get("coverage","phylogenetic_decay").strip().split(",")
				if not (len(values)%2==0):
					return False, "\n\t{0}\n\t{1}\n\t{2}\n\t{3}".format(\
						parserMessageWrong,\
						"[coverage] section: phylogenetic_decay option invalid.",\
						"Incorrect number of values.",\
						"Please verify. Exiting"\
						)
				phyloDecayOK=True
				for item in range(0,len(values),2):
					val=float(values[item+1])
					if (val >= 0 and val <=1):
						self.phylogeneticDecay[str(values[item])]=val
					else:
						phyloDecayOK=False
				if not phyloDecayOK:
					return False, "\n\t{0}\n\t{1}\n\t{2}\n\t{3}".format(\
						parserMessageWrong,\
						"[coverage] section: phylogenetic_decay option invalid.",\
						"One or many values are out of range. Value should be in the interval [0,1].",\
						"Please verify. Exiting"\
						)
				self.appLogger.info("Applying phylogenetic coverage decay")
		else:
			# No coverage section
			message="\n\t{0}\n\t{1}\n\t{2}\n\t".format(
				"Settings: [coverage] section",\
				"When using [ngs-reads-art] or [ngs-read-count] section. Coverage is required.",\
				"Please verify. Exiting."\
			)
			return False,message

		return True, message

	def checkSectionExecution(self,parserMessageCorrect,parserMessageWrong):
		########################################################################
		# BLOCK: Execution
		########################################################################
		if not self.parser.has_section("execution"):
			self.appLogger.warning("Settings - Execution block: This block has been automatically generated.")
			self.parser.add_section("execution")
			self.parser.set("execution", "environment","bash")
			self.parser.set("execution", "run","off")
			self.parser.set("execution", "threads","1")
		else:
			####################################################################
			# OPTION: Environment
			if (self.parser.has_option("execution","env")):
				# got the short name
				value=self.parser.get("execution","env")
				self.parser.set("execution","environment",value.lower())
				self.parser.remove_option("execution","environment")
			elif (self.parser.has_option("execution","environment")):
				# got the long name, make sure it is lowercase and within the options
				value=self.parser.get("execution","environment")
				if (value in ["sge","slurm","bash"]):
					self.parser.set("execution","environment",value.lower())
					if (value in ["sge","slurm"]):
						self.parser.set("execution", "run","off")
				else:
					message="Settings: Execution block | Evironment variable is incorrect or unavailable. Please check the settings file and rerun. Exiting."
					return False,message
			else:
				# got no environment
				self.parser.set("execution", "environment","bash")
			####################################################################
			# OPTION: RUN
			if (self.parser.has_option("execution","run")):
				try:
					value=self.parser.getboolean("execution","run")
				except Exception as e:
					self.appLogger.warning("Settings - Execution block: Run automatically set up to OFF.")
					self.parser.set("execution","run","off")
			else:
				self.appLogger.warning("Settings - Execution block: Run automatically set up to OFF.")
				self.parser.set("execution","run","off")
			####################################################################
			# OPTION: threads
			if (self.parser.has_option("execution","threads")):
				try:
					self.numThreads=self.parser.getint("execution","threads")
				except Exception as e:
					self.appLogger.warning("Settings - Execution block: Threads automatically set up to 1.")
					self.parser.set("execution","threads","1")
					self.numThreads=1
			else:
				self.numThreads=1
				self.appLogger.warning("Settings - Execution block: Threads automatically set up to 1.")
				self.parser.set("execution","threads","1")


	def checkIndelibleControlFile(self,parserMessageCorrect,parserMessageWrong):
		self.appLogger.debug("checkIndelibleControlFile(self,parserMessageCorrect,parserMessageWrong)")
		f=open(self.ngsphyIndelibleControlFilePath,"r")
		lines=f.readlines()
		f.close()
		# keeping only lines with content
		newlines=[ item.strip() for item in lines if not item.strip()==""]
		# check for NGSPHY blocks
		model=0; partition=0
		for index in range(0,len(newlines)):
			if newlines[index].startswith("[MODEL]") and  model == 0:
				if model== 0:
					model=index
				else:
					parserMessageWrong+="Validating INDELible Control file: \n\t{0}\n\t{1}".format(\
						"There is more than one [MODEL] block.",\
						"Please verify. Exiting"\
						)
					return False, parserMessageWrong
			if newlines[index].startswith("[NGSPHYPARTITION]") and partition == 0:
				if partition==0:
					partition=index
				else:
					parserMessageWrong+="Validating INDELible Control file: \n\t{0}\n\t{1}".format(\
						"There is more than one [NGSPHYPARTITION] block.",\
						"Please verify. Exiting"\
						)
					return False, parserMessageWrong

		# if i have gotten here and did not return before
		# i can check if the blocks match :D
		# 1) tree filename matches partition parameter and partition has the right
		# number of parameters
		self.partition=newlines[partition].split()
		if len(self.partition)!=4:
			parserMessageWrong+="Validating INDELible Control file: \n\t{0}\n\t{1}\n\t{2}".format(\
			"[NGSPHYPARTITION] block has the wrong number of parameters.",\
			"[NGSPHYPARTITION] <tree_filename_basename> <model_name> <sequence_length>"
			"Please verify. Exiting"
			)
			return False, parserMessageWrong
		# check tree corresponds to newick inputbasename
		newickBasename,_=os.path.splitext(os.path.basename(self.newickFilePath))
		if not self.partition[1]==newickBasename:
			parserMessageWrong+="Validating INDELible Control file: \n\t{0}\n\t{1}\n\t{2}\n\t{3}".format(\
			"[NGSPHYPARTITION] block, tree name does not correspond with the Newick File introduced.",\
			"Remember! Newick filename: newick.tree.",\
			"[NGSPHYPARTITION] newick model1 200"
			"Please verify. Exiting"
			)
			return False, parserMessageWrong

		modelname=newlines[model].strip().split()[1]
		if not self.partition[2]==modelname:
			parserMessageWrong+="Validating INDELible Control file: \n\t{0}\n\t{1}\n\t{2}\n\t{3}".format(\
			"[NGSPHYPARTITION] block, model name does not correspond to the model defined.",\
			"[MODEL] modelname\n...",\
			"[NGSPHYPARTITION] tree modelname 200"
			"Please verify. Exiting"
			)
			return False, parserMessageWrong

		if not self.partition[3].isdigit():
			parserMessageWrong+="Validating INDELible Control file: \n\t{0}\n\t{1}".format(\
			"[NGSPHYPARTITION] block, sequence length is not valid.",\
			"Please verify. Exiting"
			)
			return False, parserMessageWrong

		return True, parserMessageCorrect

	def checkSimPhyProjectValid(self):
		matingArgsMessageWrong="SimPhy project is not valid."
		matingArgsMessageCorrect="SimPhy project is valid."
		# List all the things in the project directory
		fileList=os.listdir(os.path.join(self.path,self.projectName))
		for index in range(0,len(fileList)):
			fileList[index]=os.path.abspath(os.path.join(self.path,self.projectName,fileList[index]))

		command = os.path.join(\
			self.path,\
			self.projectName,\
			"{0}.command".format(self.projectName))
		params = os.path.join(\
			self.path,\
			self.projectName,\
			"{0}.params".format(self.projectName))
		db = os.path.join(\
			self.path,\
			self.projectName,\
			"{0}.db".format(self.projectName))

		self.appLogger.debug("SimPhy files (command, params, db)")
		self.appLogger.debug("{0}:\t{1}".format(\
			os.path.basename(db),db in fileList))
		self.appLogger.debug("{0}:\t{1}".format(\
			os.path.basename(command),command in fileList))
		self.appLogger.debug("{0}:\t{1}".format(\
			os.path.basename(params),params in fileList))

		simphyfiles=((command in fileList) and (params in fileList) and(db in fileList))
		# check if  command, db, params files
		if not simphyfiles:
			matingArgsMessageWrong+="\n\tSimPhy files do not exist. Please verify. Exiting."
			return False, matingArgsMessageWrong
		# check how many of them are dirs

		for item in fileList:
			baseitem=os.path.basename(item)
			if (os.path.isdir(os.path.abspath(item)) and  baseitem.isdigit()):
				self.numSpeciesTrees=self.numSpeciesTrees+1

		numSpeciesTreesDigits=len(str(self.numSpeciesTrees))
		self.parser.set("general","numspeciestrees",str(self.numSpeciesTrees))
		# check if at least one
		self.appLogger.debug("Num species trees:\t{0}".format(self.numSpeciesTrees))
		if not (self.numSpeciesTrees>0):
			matingArgsMessageWrong+="\n\tNot enough number of species tree replicates (at least 1 required):\t\n\t{0}\n\t Please verify. Exiting.".format(self.numSpeciesTrees>0)
			return False, matingArgsMessageWrong

		return True, matingArgsMessageCorrect


	def checkLabelFormatInTree(self):
		self.appLogger.debug("Checking labels")
		messageCorrect="Labels of the tree are correct"
		messageWrong="INDELible control file - Something's wrong!\n\t"
		tree=dendropy.Tree.get(path=self.newickFilePath, schema="newick",preserve_underscores=True)
		leaves=[ node.taxon.label for node in tree.leaf_node_iter()]
		item=""
		for item in leaves:
			if not  bool(re.match("^([0-9]+_[0-9]+_[0-9]+){1}",item)):
				break
		if not bool(re.match("^([0-9]+_[0-9]+_[0-9]+){1}",item)):
			messageWrong+="\n\t{0}\n\t{1}\n\t{2}\n\t{3}\n\t{4}\n\t{5}".format(\
			"Labels chosen for the tips of the tree are not correct.",\
			"Labels should follow this pattern: SpeciesID_LocusID_IndividualID",\
			"Where SpeciesID,LocusID,IndividualID are numbers.",\
			"Where SpeciesID > 0, LocusID > 0, IndividualID > 0",\
			"Except for the Outgroup, if present, should be represented as 0_0_0",\
			"Please verify. Exiting."\
			)
			return False,messageWrong
		return True, messageCorrect

	def checkReferenceTipLabelInTree(self):
		# checking if label in the given tree
		self.appLogger.debug("Checking if label in the given tree")
		messageCorrect="Labels of the tree are correct"
		messageWrong="INDELible control file - Something's wrong!\n\t"
		tree=dendropy.Tree.get(path=self.newickFilePath, schema="newick",preserve_underscores=True)
		filter = lambda taxon: True if taxon.label==self.referenceTipLabel else False
		node = tree.find_node_with_taxon(filter)
		if not node:
			messageWrong+="\n\t{0}\n\t{1}".format(\
			"Tip label for reference sequence is not in the Newick file.",\
			"Please verify. Exiting.")
			return False, messageWrong
		return True, messageCorrect

	def formatSettingsMessage(self):
		message="Settings:\n"
		sections=self.parser.sections()
		for sec in sections:
			message+="\t{0}\n".format(sec)
			items=self.parser.items(sec)
			for param in items:
				extra=""
				if param[0]=="origin":
					if int(param[1])==1: extra="(SimPhy output)"
					if int(param[1])==2: extra="(Gene tree + INDELible control file)"
					if int(param[1])==3: extra="(Gene tree + INDELible control file + Reference sequence)"
				if param[0]=="ploidy":
					if int(param[1])==1: extra="(Haploid individuals)"
					if int(param[1])==2: extra="(Diploid individuals)"
				message+="\t\t{0}\t:\t{1} {2}\n".format(param[0],param[1],extra)
		return message
