import argparse,datetime,dendropy, logging,os,re,sys
import numpy as np
from coverage import NGSPhyDistributionParser as ngsphydistro
if (sys.version_info[0:2]<(3,0)):
	import ConfigParser as cp
elif (sys.version_info>=(3,0)):
	import configparser as cp

class Settings:
	"""
	This class parses and verifies that the settings file is correct.
	----------------------------------------------------------------------------
	Attributes:

	General
	-------
	- inputmode: Represents where the input data comes from. Has 3 different values.
	(1) when using the output of a SimPhy runART. (2) when using a single gene tree
	with a random root/reference sequence and (3) when using a single gene tree with a
	specific root/reference sequence.
	- ploidy: This is to identify whether the user wants to generate haploid or
	diploid individuals from the sequences given/generated. Not all ploidies are
	accepted for all oriiginModes. Please verify the external documentation for a
	detail explanation.
	- projectName: NGSPHY or the name of the SimPhy folder.
	- path: indicates where the output folder is/will be created.
	- outputFolderName: Name of the folder that will be used to store the generated
	data.
	- outputFolderPath: absolute path of the output folder.
	- simphyDataPrefix: Prefix of the sequences. If generated within NGSphy, default is
	"ngsphydata", otherwise, this informations is parsed from the settings
	file and will have to correspond, to the prefix given to SimPhy.
	- simphyFolderPath: Related to InputMode=1. Corresponds to the absolute path
	of the SimPhy project.
	- simphyFilter:  Related to InputMode=1. Using species tree distributions.
	Allows the user to filter out from the program execution, those species tree
	replicates that do not match the given ploify to the number of individuals
	per species to be able to generate individuals according to the given ploidy.


	- indelibleControlFile: Related to InputMode in [2,3]. Path
	to the modified INDELible control file.
	- geneTreeFile: Related to InputMode=2. Path to the gene tree file in
	Newick format. Must contain a single tree. Name of the file without extension
	must match the name of the tree within the control file.
	- evolve: Related to InputMode=2. To handle evolve parameters of the control
	file.
	- partition: Related to InputMode=2. To handle partition parameters of the
	control file.


	- ancestralSequenceFilePath: Related to InputMode=3. Path to the reference
	sequene file where the root/reference sequence is.
	- anchorTipLabel: Related to InputMode=3. Label of the tip that will be
	used as reference/root.

	- readCountsError: Related to NGSMode=1. sequencing read_counts_error that will be used for the read counts.
	- readCountsReferenceAllelesFile: Related to NGSMode=1. Path to the file that contains all the
	information (ReplicateIdentifier, Species, Locus, Gene) needed to
	select a sequence and use it as a reference to generate the read counts.

	- ngsmode: Variable that identifies which ngsmode will be used.

	- ontarget: Parameter related to coverage. Only possible with inputmode=1.
	When emulating Targeted-sequencing, represents the fraction of the loci that
	will be considered as on-target.
	- offtarget: Parameter related to coverage.  Only possible with inputmode=1.
	When emulating Targeted-sequencing, represents the fraction of the loci that
	will be considered as on-target.
	- notcaptured: Parameter related to coverage.  Only possible with inputmode=1.
	When emulating Targeted-sequencing, represents the fraction of the on-target
	loci that will not be sequece, coverage == 0.
	- experiment: Parameter related to coverage. Value for the expected coverage.
	Parameterization can be based on a distribution or a fixed value.
	- individual: Parameter related to coverage. Expected coverage multiplier. Hyper
	parameter. Value used here can be based on a distribution or a fixed value.
	Value sampled will be used as the alpha shape of a Gamma distribution with
	mean = 1.
	- locus: Parameter related to coverage. Expected coverage multiplier. Hyper
	parameter. Value used here can be based on a distribution or a fixed value.
	Value sampled will be used as the alpha shape of a Gamma distribution with
	mean = 1.
	- taxon: Parameter related to coverage. This parameters is to
	represent the possible coverage decay due to the distance of the species
	to the reference sequence within the given tree.
	Is a list of pairs, where the first element is the species identifier and
	the second, is the fraction of coverage that the species will keep.

	- numReplicates: number of species tree replicates.
	- numLociPerReplicate: digits used to represente the numReplicates
	- indels: Indicates whether the sequences (given or generated) contain
	indels. Important for ngsmode=1 (ReadCounts), because this mode does not
	handle sequence with indels.
	"""
	# General
	inputmode=1
	ploidy=1
	projectName="NGSphy"
	path=""
	outputFolderName="NGSphy_output"
	outputFolderPath=""
	basepath=""


	# simphy data origin | origin 1
	simphyDataPrefix="ngsphydata"
	simphyFolderPath=""
	simphyFilter=False

	# indelible data origin | origin 2
	indelibleControlFile=""
	geneTreeFile=""
	evolve=""
	partition=""
	# indelible data origin + referebce | origin 3
	ancestralSequenceFilePath=""
	anchorTipLabel=""

	#readcount
	readCountsError=0
	readCountsReferenceAllelesFile=""

	ngsmode=0

	# coverage
	ontarget=1
	offtarget={"loci":0, "coverage":1}
	notcaptured=0
	experiment=None
	individual=None
	locus=None
	taxon=dict()

	numReplicates=0
	numLociPerReplicate=[]
	indels=False

	alignmentsFolderPath=""
	coverageFolderPath=""
	individualsFolderPath=""
	readsFolderPath=""
	refAllelesFolderPath=""
	scriptsFolderPath=""
	tablesFolderPath=""

	runART=False
	threads=1
	executionMode=1


	__NUCLEOTIDES=["A","C","G","T"]
	def __init__(self,filename):
		# If I've got this far, then filename is a correct file
		self.settingsFile=os.path.abspath(filename)
		self.appLogger=logging.getLogger('ngsphy')
		self.appLogger.debug("(class Settings) __init__()")
		# default settings can be established.
		self.parser=cp.SafeConfigParser()
		self.parser.read(self.settingsFile)

	def checkArgs(self):
		"""
		Checks the existence and validity of the settings introduced in the
		settings file per block and option.s
		------------------------------------------------------------------------
		Returns:
		- boolean, message: the status of the process and the message related to
		such status
		"""
		self.appLogger.debug("Check arguments")
		allGood=True
		parserMessageCorrect="All parameters are correct."
		parserMessageWrong="Settings - Problem found! "
		# Checking general block
		statusGeneral,messageGeneral= self.checkBlockGeneral(\
			parserMessageCorrect,parserMessageWrong)

		# Exit here
		if not statusGeneral: return statusGeneral, messageGeneral
		# Checking data block
		statusData,messageData=self.checkBlockData(\
			parserMessageCorrect,parserMessageWrong)
		# Exit here
		if not statusData: return statusData, messageData
		# ----------------------------------------------------------------------
		# NGS sections
		if self.parser.has_section("ngs-reads-art"):
			self.ngsmode=1
			statusNGSArt,messageNGSArt=self.checkBlockNGSReadsArt(\
				parserMessageCorrect,parserMessageWrong)
			if (statusNGSArt):
				## Next check
					statusCoverage,messageCoverage=self.checkBlockCoverage(parserMessageCorrect,parserMessageWrong)
					# Exit here
					if not statusCoverage: return statusCoverage, messageCoverage
			else:
				# Exit here
				return statusNGSArt,messageNGSArt
			# if all parameters for NGS ART are correct I don't need the other NGS block
			# removing NGS read count
			if self.parser.has_section("ngs-read-counts"):
				self.parser.remove_section("ngs-read-counts")
				self.appLogger.warning("[ngs-read-counts] block is incompatible with [ngs-reads-art]. Omiting this block.")
		elif self.parser.has_section("ngs-read-counts"):
			self.ngsmode=2
			# readcount
			statusRC,messageRC=self.checkBlockNGSReadCounts(parserMessageCorrect,parserMessageWrong)
			if statusRC:
				statusCoverage,messageCoverage=self.checkBlockCoverage(parserMessageCorrect,parserMessageWrong)
				# Exit here
				if not statusCoverage: return statusCoverage, messageCoverage
			else:
				# Exit here
				return statusRC,messageRC
		else: # Meaning have no ngs-reads-art nor ngs-read-countss
			if (self.parser.has_section("coverage")):self.parser.remove_section("coverage")
			self.appLogger.info("[coverage] block is not needed if [ngs-reads-art] or [ngs-read-counts] are not available. Omiting this block.")

		# Checking execution block and options
		self.checkBlockExecution(parserMessageCorrect,parserMessageWrong)
		# Exit here
		self.appLogger.info(self.formatSettingsMessage())
		if self.settings.inputmode <4:
			self.basepath=self.alignmentsFolderPath
		else:
			self.basepath=self.settings.simphyFolderPath

		return allGood,parserMessageCorrect

	def checkBlockGeneral(self,parserMessageCorrect,parserMessageWrong):
		"""
		Checks the General block of the settings, regarding general paths,
		project name and ploidy.
		------------------------------------------------------------------------
		Returns:
		- boolean, message: the status of the process and the message related to
		such status
		"""
		self.appLogger.debug("Section General")
		# Check GENERAL SECTION
		if not self.parser.has_section("general"):
			parserMessageWrong="\n\t{0}\n\t{1}".format(\
				"[general] block missing and required.",\
				"Please verify. Exiting."\
			)
			return False, parserMessageWrong
		# CHECKING GENERAL PARAMETERS
		# checking ploidy for the output data
		if (not self.parser.has_option("general","ploidy")):
			self.ploidy=1
		else:
			p=self.parser.getint("general","ploidy")
			if p in [1,2]:
				self.ploidy=p
			else:
				parserMessageWrong="\n\t{0}\n\t{1}\n\t{2}".format(\
					"Ploidy value is out of range.",\
					"Value must be in [1,2].",\
					"Please verify. Exiting."\
				)
				return False, parserMessageWrong

		if (self.parser.has_option("general","path")):
			self.path=os.path.abspath(self.parser.get("general","path"))
		else:
			parserMessageWrong="\n\t{0}\n\t{1}".format(\
				"Path option is needed.",\
				"Please verify. Exiting."
			)
			return False, parserMessageWrong

		# Checking output folder information
		if(self.parser.has_option("general","output_folder_name")):
			self.outputFolderName=self.parser.get("general","output_folder_name")

		if os.path.exists("{0}/{1}/{2}/".format(self.path,self.outputFolderName)):
			listdir=os.listdir("{0}/{1}".format(self.path))
			counter=0
			for item in listdir:
				if self.outputFolderName in item:
					counter+=1
			if not counter == 0: self.outputFolderName+="_{0}".format(counter+1)

		self.outputFolderPath=os.path.join(self.path,self.outputFolderName)
		self.alignmentsFolderPath=os.path.join(\
			self.outputFolderName,"alignments")
		self.coverageFolderPath=os.path.join(\
			self.outputFolderName,"coverage")
		self.individualsFolderPath=os.path.join(\
			self.outputFolderName,"individuals")
		self.readsFolderPath=os.path.join(\
			self.outputFolderName,"reads")
		self.refAllelesFolderPath=os.path.join(\
			self.outputFolderName,"ref_alleles")
		self.scriptsFolderPath=os.path.join(\
			self.outputFolderName,"scripts")
		self.tablesFolderPath=os.path.join(\
			self.outputFolderName,"ind_labels")

		return True,parserMessageCorrect

	def checkBlockData(self,parserMessageCorrect,parserMessageWrong):
		"""
		Checks Data block of the settings file, regarding origin and parameters
		related to each origin mode.
		------------------------------------------------------------------------
		Returns:
		- boolean, message: the status of the process and the message related to
		such status
		"""
		self.appLogger.debug("Section Data")
		if (self.parser.has_section("data")):
			if (self.parser.has_option("data","inputmode")):
				try:
					self.inputmode=self.parser.getint("data","inputmode")
				except:
					parserMessageWrong="\n\t\n\t{0}\n\t{1}\n\t{2}".format(\
						"[data] block: Value introduced in the inputmode option is not valid.",\
						"Value should be an integer.",\
						"Please verify. Exiting."\
					)
					return False, parserMessageWrong
			else:
				parserMessageWrong="\n\t\n\t{0}\n\t{1}\n\t{2}".format(\
					"[data] block: inputmode option value is missing. This option is required.",\
					"Please verify. Exiting."\
				)
				return False, parserMessageWrong

			####################################################################
			if not self.inputmode in [1,2,3,4]:
				parserMessageWrong="\n\t\n\t{0}\n\t{1}".format(\
					"[data] block: inputmode option value is invalid",\
					"Value is out of range. Possible range [1,3]",\
					"Please verify. Exiting."\
				)
				return False, parserMessageWrong
			####################################################################
			# SINGLE GENE TREE MODES GENERAL
			if self.inputmode < 4: # 1,2,3
					# parameter is set up, now check if folder exist
				if (self.parser.has_option("data","indelible_control_file")):
					self.indelibleControlFile=os.path.abspath(\
						self.parser.get("data","indelible_control_file")
					)
				else:
					parserMessageWrong="\n\t{0}{1}{2}\n\t{3}\n\t{4}".format(\
						"[data] block: Input mode (",self.inputmode,") selected but invalid option.",\
						"INDELible  control file is required",\
						"Please verify. Exiting."\
					)
					return False, parserMessageWrong
				# checkin indelible_control_file file
				if (os.path.exists(self.indelibleControlFile) and os.path.isfile(self.indelibleControlFile)):
					self.appLogger.debug("INDELible control file exists")
				else:
					parserMessageWrong="\n\t{0}{1}{2}\n\t{3}\n\t{4}".format(\
						"[data] block: Input mode (",self.inputmode,") selected but invalid option.",\
						"INDELible control file does not exist.",\
						"Please verify. Exiting."\
					)
					return False, parserMessageWrong

				# checkin tree file
				if (self.parser.has_option("data","gene_tree_file")):
					self.appLogger.debug("Gene tree file exists")
					self.geneTreeFile=os.path.abspath(self.parser.get("data","gene_tree_file").strip())
					filesOk=(os.path.exists(self.geneTreeFile) and os.path.isfile(self.geneTreeFile))
					if not filesOk:
						parserMessageWrong="\n\t{0}{1}{2}\n\t{3}\n\t{4}".format(\
							"[data] block: Input mode (",self.inputmode,") selected but invalid option.",\
							"Gene tree file does not exist.",\
							"Please verify. Exiting."\
						)
						return False, parserMessageWrong
				else:
					parserMessageWrong="\n\t{0}{1}{2}\n\t{3}\n\t{4}".format(\
						"[data] block: Input mode (",self.inputmode,") selected but invalid option.",\
						"Gene tree file does not exist.",\
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
				statusOk, message=self.checkBranchLengthsInTree()
				if not statusOk: return statusOk,message

				self.parser.set("general","numspeciestrees",str(1))
				self.numReplicates=1
				# removing options that do not match with the origin mode selected
				if self.parser.has_option("data","simphy_data_prefix"): self.parser.remove("data","simphy_data_prefix")
				if self.parser.has_option("data","simphy_folder_path"): self.parser.remove("data","simphy_folder_path")
				if self.parser.has_option("data","simphy_filter"):  self.parser.remove("data","simphy_filter")

			####################################################################
			if self.inputmode ==1:
				# checking indelible program
				stream1 = os.popen('which indelible').read()[0:-1]
				stream2 = os.popen('which indelible-ngsphy').read()[0:-1]
				self.appLogger.info("Checking dependencies...")
				if stream1:
					self.appLogger.info("indelible - Found running in: {}".format(stream))
				elif stream2:
					self.appLogger.info("indelible-ngsphy - Found running in: {}".format(stream))
				else:
					parserMessageWrong="\n\t{0}{1}{2}\n\t{3}\n\t{4}".format(\
						"[data] block: Input mode (",self.inputmode,") selected but invalid option.",\
						"indelible nor indelible-ngsphy found. Programs either not installed or not in the current path.",\
						"Please verify the installation. Exiting."\
					)
					return False, parserMessageWrong
				if self.parser.has_option("data","ancestral_sequence_file"):  self.parser.remove_option("data","ancestral_sequence_file")
				if self.parser.has_option("data","anchor_tip_label"):  self.parser.remove_option("data","anchor_tip_label")
			####################################################################
			if self.inputmode in [2,3]:
				stream = os.popen('which indelible-ngsphy').read()[0:-1]
				self.appLogger.info("Checking dependencies...")
				if stream:
					self.appLogger.info("indelible-ngsphy - Found running in: {}".format(stream))
				else:
					parserMessageWrong="\n\t{0}{1}{2}\n\t{3}\n\t{4}".format(\
						"[data] block: Input mode (",self.inputmode,") selected but invalid option.",\
						"indelible-ngsphy not found. Program either not installed or not in the current path.",\
						"Please verify the installation. Exiting."\
					)
					return False, parserMessageWrong
				if (self.parser.has_option("data","ancestral_sequence_file")):
					self.ancestralSequenceFilePath=os.path.abspath(self.parser.get("data","ancestral_sequence_file").strip())
					filesOk=(os.path.exists(self.ancestralSequenceFilePath) and os.path.isfile(self.ancestralSequenceFilePath))
					if not filesOk:
						parserMessageWrong="\n\t{0}{1}{2}\n\t{3}\n\t{4}".format(
						"[data] block: Input mode (",self.inputmode,") selected but invalid option.",\
						"Ancestral sequence file does not exist or the given path does not belong to a file.",\
						"Please verify. Exiting."\
						)
						return False, parserMessageWrong
					statusSeq,messageSeq=self.correctContentReferenceSequence()
					if not statusSeq: return statusSeq, messageSeq
				else:
					parserMessageWrong="\n\t{0}{1}{2}\n\t{3}\n\t{4}".format(\
					"[data] block: Input mode (",self.inputmode,") selected but invalid option.",\
					"Ancestral sequence file is required",\
					"Please verify. Exiting."\
					)
					return False, parserMessageWrong

			####################################################################
			if self.inputmode==3:
				if (self.parser.has_option("data","anchor_tip_label")):
					self.anchorTipLabel=self.parser.get("data","anchor_tip_label")
					statusOk, message=self.checkAnchorTipLabelInTree()
					if not statusOk: return statusOk,message
				else:
					parserMessageWrong="\n\t{0}{1}{2}\n\t{3}\n\t{4}".format(\
						"[data] block: Input mode (",self.inputmode,") selected but invalid option.",\
						"Anchor tip label is required",\
						"Please verify. Exiting."\
					)
					return False, parserMessageWrong
			####################################################################
			if self.inputmode==4:
				if (self.parser.has_option("data","simphy_folder_path")):
					self.simphyFolderPath=os.path.abspath(self.parser.get("data","simphy_folder_path"))
				else:
					parserMessageWrong="\n\t{0}{1}{2}\n\t{3}\n\t{4}".format(\
						"[data] block: Input mode (",self.inputmode,") selected but invalid option.",\
						"SimPhy folder does not exist, or the given path does not belong to a directory.",\
						"Please verify. Exiting."\
					)
					return False, parserMessageWrong
				if (os.path.exists(self.simphyFolderPath) and os.path.isdir(self.simphyFolderPath)):
					self.projectName=os.path.basename(self.simphyFolderPath)
				else:
					parserMessageWrong="\n\t{0}{1}{2}\n\t{3}\n\t{4}".format(\
						"[data] block: Input mode (",self.inputmode,") selected but invalid option.",\
						"SimPhy folder does not exist, or the given path does not belong to a directory.",\
						"Please verify. Exiting."\
					)
					return False, parserMessageWrong
				# data prefix for the simphy sequences
				if (self.parser.has_option("data","simphy_data_prefix")):
					self.simphyDataPrefix=self.parser.get("data","simphy_data_prefix")
				else:
					parserMessageWrong="\n\t{0}{1}{2}\n\t{3}\n\t{4}".format(\
						"[data] block: Input mode (",self.inputmode,") selected but invalid option.",\
						"<simphy_data_prefix> field is missing. This prefix correponds to the name of the file sequences that are going to be processed.",\
						"Please verify. Exiting."\
					)
					return False, parserMessageWrong
				# checking flag for filtering species tree replicates
				if self.parser.has_option("data","simphy_filter"):
					self.simphyFilter=True
				else:
					self.simphyFilter=False
				# removing options that do not match with the origin mode selected
				if (self.parser.has_option("data","indelible_control_file")):
					self.parser.remove_option("data","indelible_control_file")
				if (self.parser.has_option("data","gene_tree_file")):
					self.parser.remove_option("data","gene_tree_file")
				if (self.parser.has_option("data","anchor_tip_label")):
					self.parser.remove_option("data","anchor_tip_label")
				if (self.parser.has_option("data","ancestral_sequence_file")):
					self.parser.remove_option("data","ancestral_sequence_file")
				# --------------------------------------------------------------
				# check if simphy is valid projectName
				status,message=self.checkSimPhyProjectValid()
				if not status:
					return status, message
				####################################################################
		else:
			parserMessageWrong="\n\t\n\t{0}\n\t{1}\n\t{2}".format(\
				"[data] block is missing. This block is required.",\
				"Please verify. Exiting."\
			)
			return False, parserMessageWrong
		return True, parserMessageCorrect


	def checkBlockNGSReadsArt(self,parserMessageCorrect,parserMessageWrong):
		"""
		Checks NGS-reads-ART block of the settings file, regarding parameters
		of the NGS simulation with ART
		------------------------------------------------------------------------
		Returns:
		- boolean, message: the status of the process and the message related to
		such status
		"""
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
			self.appLogger.warning("Removing ART coverage options. Coverage is calculated with the [coverage] block parameters.")
		else:
			parserMessageWrong="\n\t{0}\n\t{1}".format(\
				"art_illumina not found. Program either not installed or not in your current path.",\
				"Please verify the installation. Exiting."
			)
			return False, parserMessageWrong
		return True, parserMessageCorrect

	def checkBlockNGSReadCounts(self,parserMessageCorrect,parserMessageWrong):
		"""
		Checks read counts  block of the settings file, regarding parameters
		of the read counts process.
		------------------------------------------------------------------------
		Returns:
		- boolean, message: the status of the process and the message related to
		such status
		"""
		########################################################################
		# BLOCK: READ COUNT
		########################################################################
		message=parserMessageCorrect
		if (self.parser.has_section("ngs-read-counts")):
			if not self.parser.has_option("ngs-read-counts", "read_counts_error"):
				self.appLogger.warning("[ngs-read-counts] block. Read counts error rate for this runART is being considered as 0.")
				self.parser.set("ngs-read-counts", "read_counts_error","0")
				self.readCountsError=0
			else:
				try:
					self.readCountsError=float(self.parser.get("ngs-read-counts", "read_counts_error"))
				except:
					message="\n\t{0}\n\t{1}".format(\
						"[ngs-read-counts] block."
					)
					return False, message

			if not self.parser.has_option("ngs-read-counts","reference_alleles_file"):
				self.appLogger.warning("[ngs-read-counts] block. Using default references.")
				self.parser.set("ngs-read-counts", "reference_alleles_file","None")
				self.readCountsReferenceAllelesFile=None
			else:
				self.readCountsReferenceAllelesFile=os.path.abspath(self.parser.get("ngs-read-counts","reference_alleles_file"))
				fileOk=os.path.exists(self.readCountsReferenceAllelesFile) and os.path.isfile(self.readCountsReferenceAllelesFile)
				if not fileOk:
					status=False
					message="\n\t{0}\n\t{1}\n\t{2}".format(\
						"[ngs-read-counts] block.",\
						"Reference file does not exist or path is incorrect.",\
						"Please verify. Exiting.")
		else:
			# No ngs-read-counts block
			status=False
			message="\n\t{0}\n\t{1}".format(\
				"[ngs-read-counts] block is missing.",\
				"Please verify. Exiting.")

		return status,message

	########################################################################
	# BLOCK: Coverage
	########################################################################
	def checkBlockCoverage(self,parserMessageCorrect,parserMessageWrong):
		"""
		Checks coverage block of the settings file, regarding parameters that
		model the coverage variation within the execution
		------------------------------------------------------------------------
		Returns:
		- boolean, message: the status of the process and the message related to
		such status
		"""
		self.appLogger.debug("Checking coverage")
		message=parserMessageCorrect
		expCov=None;indCov=None;locCov=None;
		if(self.parser.has_section("coverage")):
			# ------------------------------------------------------------------
			#  Targeted-sequencing related
			if self.inputmode==1:
				self.appLogger.debug("Inputmode 1")
				if (self.parser.has_option("coverage","offtarget")):
					self.appLogger.debug("offtarget")
					offTargetValues=self.parser.get("coverage","offtarget").strip().split(",")
					if not (len(offTargetValues) ==2):
						return False, "\n\t{0}\n\t{1}\n\t{2}\n\t{3}".format(\
							parserMessageWrong,\
							"[coverage] block: offtarget option invalid.",\
							"Incorrect number of parameters.",\
							"Please verify. Exiting"\
							)
					try:
						self.offtarget["loci"]=float(self.offtarget[0])
						self.offtarget["coverage"]=float(self.offtarget[1])

					except:
						return False, "\n\t{0}\n\t{1}\n\t{2}\n\t{3}".format(\
							parserMessageWrong,\
							"[coverage] block: offtarget option invalid.",\
							"One or both values are incorrect.",\
							"Please verify. Exiting"\
							)

					if (self.offtarget["loci"] < 0 and self.offtarget["loci"] > 1) and \
						(self.offtarget["coverage"] < 0 and self.offtarget["coverage"] > 1):
						return False, "\n\t{0}\n\t{1}\n\t{2}\n\t{3}".format(\
							"[coverage] block: offtarget option invalid.",\
							"One or both values are out of range",\
							"Value should be in the interval [0,1]",\
							"Please verify. Exiting"\
							)
					self.ontarget=1-self.offtarget["loci"]
					self.parser.set("coverage","ontarget", str(self.ontarget))

				if (self.parser.has_option("coverage","notcaptured")):
					self.appLogger.debug("not captured")
					self.notcaptured=self.parser.getfloat("coverage","notcaptured")
					if not (self.notcaptured <= 1):
						return False, "\n\t{0}\n\t{1}\n\t{2}\n\t{3}".format(\
							parserMessageWrong,\
							"[coverage] block: notcaptured option invalid.",\
							"Value out of range, should be within the half-closed interval [0,1)",\
							"Please verify. Exiting"\
						)
			else:
				if (self.parser.has_option("coverage","ontarget")): self.parser.remove("coverage","ontarget")
				if (self.parser.has_option("coverage","offtarget")): self.parser.remove("coverage","offtarget")
				if (self.parser.has_option("coverage","notcaptured")): self.parser.remove("coverage","notcaptured")
				if (self.parser.has_option("coverage","ontarget") or self.parser.has_option("coverage","offtarget") or self.parser.has_option("coverage","notcaptured")):
					self.appLogger.warning("Coverage options - offtarget and notcaptured - only possible with [data] inputmode = 1 (gene-tree distribution)")
			# ------------------------------------------------------------------
			#  Coverage distribution related
			if (self.parser.has_option("coverage","experiment")):
				self.appLogger.debug("Coverage experiment")
				value=self.parser.get("coverage","experiment")
				self.parser.set("coverage","experiment",value.lower())
			else:
				# parsear distribution
				parserMessageWrong="\n\t{0}\n\t{1}".format(\
					"[coverage] block:  experiment option is required.",\
					"Please verify. Exiting."\
				)
				return False,parserMessageWrong

			self.experiment=ngsphydistro(self.parser.get("coverage","experiment"), False)
			check,mess=self.experiment.validate()
			if not (check):
				parserMessageWrong=mess
				return check,parserMessageWrong

			# INDIVIDUAL -------------------------------------------------------
			if (self.parser.has_option("coverage","individual")):
				self.appLogger.debug("Coverage individual")
				value=self.parser.get("coverage","individual")
				self.parser.set("coverage","individual",value.lower())
				self.individual=ngsphydistro(self.parser.get("coverage","individual"), False)
				check,mess=self.individual.validate()
				if not (check):
					parserMessageWrong=mess
					return check,parserMessageWrong
			# LOCUS ------------------------------------------------------------
			if (self.parser.has_option("coverage","locus")):
				self.appLogger.debug("Coverage locus")
				value=self.parser.get("coverage","locus")
				self.parser.set("coverage","locus",value.lower())
				self.locus=ngsphydistro(self.parser.get("coverage","locus"), False)
				check,mess=self.locus.validate()
				if not (check):
					parserMessageWrong=mess
					return check,parserMessageWrong
			# ------------------------------------------------------------------
			#  taxon
			if (self.parser.has_option("coverage","taxon")):
				values=self.parser.get("coverage","taxon").strip().split(",")
				if not (len(values)%2==0):
					return False, "\n\t{0}\n\t{1}\n\t{2}\n\t{3}".format(\
						parserMessageWrong,\
						"[coverage] block: taxon option invalid.",\
						"Incorrect number of values.",\
						"Please verify. Exiting"\
						)
				taxonOk=True
				for item in range(0,len(values),2):
					# TODO: may need a try catch block here, string to float conversion
					try:
						val=float(values[item+1])
					except:
						return False, "\n\t{0}\n\t{1}\n\t{2}\n\t{3}".format(\
							parserMessageWrong,\
							"[coverage] block: taxon option invalid.",\
							"One or many values are invalid.",\
							"Please verify. Exiting"\
							)
					if (val >= 0 and val <=1):
						self.taxon[str(values[item])]=val
					else:
						taxonOk=False
				if not taxonOk:
					return False, "\n\t{0}\n\t{1}\n\t{2}\n\t{3}".format(\
						parserMessageWrong,\
						"[coverage] block: taxon option invalid.",\
						"One or many values are out of range. Value should be in the interval [0,1].",\
						"Please verify. Exiting"\
						)
				self.appLogger.info("Applying phylogenetic coverage decay")
		else:
			# No coverage block
			message="\n\t{0}\n\t{1}\n\t{2}\n\t".format(
				"Settings: [coverage] block",\
				"When using [ngs-reads-art] or [ngs-read-counts] block. Coverage is required.",\
				"Please verify. Exiting."\
			)
			return False,message

		return True, message

	def checkBlockExecution(self,parserMessageCorrect,parserMessageWrong):
		"""
		Checks execution block of the settings file, regarding parameters
		of the general execution of the programs, type, threading...
		------------------------------------------------------------------------
		Returns:
		- boolean, message: the status of the process and the message related to
		such status
		"""
		########################################################################
		# BLOCK: Execution
		########################################################################
		if not self.parser.has_section("execution"):
			self.appLogger.warning("Settings - Execution block: This block has been automatically generated.")
			self.parser.add_section("execution")
			self.parser.set("execution", "environment","bash")
			self.parser.set("execution", "runART","off")
			self.parser.set("execution", "threads","1")
			self.parser.set("execution", "runnint_times","off")
		else:
			####################################################################
			# OPTION: Environment
			if (self.parser.has_option("execution","environment")):
				# got the long name, make sure it is lowercase and within the options
				value=self.parser.get("execution","environment")
				if (value in ["sge","slurm","bash"]):
					self.parser.set("execution","environment",value.lower())
					if value =="bash": self.executionMode=1
					if value =="sge": self.executionMode=2
					if value =="slurm": self.executionMode=3
					if self.executionMode > 1:	self.runART=False
				else:
					message="\n\t{0}\n\t{1}\n\t{2}".format(\
						"Settings: Execution block",\
						"Environment variable is incorrect or unavailable.",\
						"Please check the settings file and rerun. Exiting."
					)
					return False,message
			else:
				message="\n\t{0}\n\t{1}\n\t{2}".format(\
					"Settings: Execution block",\
					"Environment variable is incorrect or unavailable.",\
					"Please check the settings file and rerun. Exiting."
				)
				return False,message

			####################################################################
			# OPTION: RUN
			if (self.parser.has_option("execution","runART")):
				try:
					value=self.parser.getboolean("execution","runART")
				except Exception as e:
					message="\n\t{0}\n\t{1}\n\t{2}".format(\
						"Settings: Execution block",\
						"runART variable is incorrect.",\
						"Please check the settings file and rerun. Exiting."
					)
					return False,message
			else:
				self.appLogger.warning("Settings - Execution block: Run automatically set up to OFF.")
				self.parser.set("execution","runART","off")
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
			####################################################################
			 # get running times
 			if (self.parser.has_option("execution","running_times")):
 				try:
 					self.runningTimes=self.parser.getboolean("execution","running_times")
 				except Exception as e:
 					self.runningTimes=False
 			else:
 				self.runningTimes=False
 				self.parser.set("execution","running_times","0")


	def checkIndelibleControlFile(self,parserMessageCorrect,parserMessageWrong):
		"""
		Checks indelible control file according to the modifications for this program.
		------------------------------------------------------------------------
		Returns:
		- boolean, message: the status of the process and the message related to
		such status
		"""
		self.appLogger.debug("checkIndelibleControlFile(self,parserMessageCorrect,parserMessageWrong)")
		f=open(self.indelibleControlFile,"r")
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
					parserMessageWrong="Validating INDELible Control file: \n\t{0}\n\t{1}".format(\
						"There is more than one [MODEL] block.",\
						"Please verify. Exiting"\
						)
					return False, parserMessageWrong
			if newlines[index].startswith("[NGSPHYPARTITION]") and partition == 0:
				if partition==0:
					partition=index
				else:
					parserMessageWrong="Validating INDELible Control file: \n\t{0}\n\t{1}".format(\
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
			parserMessageWrong="Validating INDELible Control file: \n\t{0}\n\t{1}\n\t{2}".format(\
			"[NGSPHYPARTITION] block has the wrong number of parameters.",\
			"[NGSPHYPARTITION] <tree_filename_basename> <model_name> <sequence_length>"
			"Please verify. Exiting"
			)
			return False, parserMessageWrong
		# check tree corresponds to newick inputbasename
		newickBasename,_=os.path.splitext(os.path.basename(self.geneTreeFile))
		if not self.partition[1]==newickBasename:
			parserMessageWrong="Validating INDELible Control file: \n\t{0}\n\t{1}\n\t{2}\n\t{3}".format(\
			"[NGSPHYPARTITION] block, tree name does not correspond with the Newick File introduced.",\
			"Remember! Newick filename: newick.tree.",\
			"[NGSPHYPARTITION] newick model1 200"
			"Please verify. Exiting"
			)
			return False, parserMessageWrong

		modelname=newlines[model].strip().split()[1]
		if not self.partition[2]==modelname:
			parserMessageWrong="Validating INDELible Control file: \n\t{0}\n\t{1}\n\t{2}\n\t{3}".format(\
			"[NGSPHYPARTITION] block, model name does not correspond to the model defined.",\
			"[MODEL] modelname\n...",\
			"[NGSPHYPARTITION] tree modelname 200"
			"Please verify. Exiting"
			)
			return False, parserMessageWrong

		if not self.partition[3].isdigit():
			parserMessageWrong="Validating INDELible Control file: \n\t{0}\n\t{1}".format(\
			"[NGSPHYPARTITION] block, sequence length is not valid.",\
			"Please verify. Exiting"
			)
			return False, parserMessageWrong

		return True, parserMessageCorrect

	def checkSimPhyProjectValid(self):
		"""
		Verifies that the data given as input from a SimPhy project is complete
		and correct
		------------------------------------------------------------------------
		Returns:
		- boolean, message: the status of the process and the message related to
		such status
		"""
		message=""
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
			message="\n\t{0}\n\t{1}".format(\
				"One of the mandatory SimPhy files does not exist.",\
				"Please verify. Exiting.")
			return False, message
		# check how many of them are dirs

		for item in fileList:
			baseitem=os.path.basename(item)
			if (os.path.isdir(os.path.abspath(item)) and  baseitem.isdigit()):
				self.numReplicates=self.numReplicates+1

		numReplicatesDigits=len(str(self.numReplicates))
		self.parser.set("general","numspeciestrees",str(self.numReplicates))
		# check if at least one
		self.appLogger.debug("Num species trees:\t{0}".format(self.numReplicates))
		if not (self.numReplicates>0):
			message="\n\t{0}\n\t{1}".format(
				"Not enough number of species tree replicates (at least 1 required)",\
				"Please verify. Exiting."\
			)
			return False, message

		return True, message


	def checkBranchLengthsInTree(self):
		"""
		Verifies existence of branch lengths in the given tree
		------------------------------------------------------------------------
		Returns:
		- boolean, message: the status of the process and the message related to
		such status
		"""
		self.appLogger.debug("Checking branch lengths")
		message=""
		tree=dendropy.Tree.get(path=self.geneTreeFile, schema="newick",preserve_underscores=True)
		leaves=[ node.taxon.label for node in tree.leaf_node_iter()]
		leafedge=None
		for leaf in leaves:
			leafedge= leaf.edge_length
			if not leafedge:
				break
		if not leafedge:
			message="\n\t{0}\n\t{1}\n\t{2}\n\t{3}\n\t{4}\n\t{5}\n\t".format(\
			"INDELible control file - Something's wrong!",\
			"Branch lengths are not specified in the gene tree file.",\
			"Please verify. Exiting."\
			)
			return False,message

		pass

	def checkLabelFormatInTree(self):
		"""
		Verifies that the data given in gene tree format follows the formatting
		rules given. Specifically, that the labels of the tree are in the
		<SpeciesID_LocusID_IndividualID> format.
		------------------------------------------------------------------------
		Returns:
		- boolean, message: the status of the process and the message related to
		such status
		"""
		self.appLogger.debug("Checking labels")
		message=""
		tree=dendropy.Tree.get(path=self.geneTreeFile, schema="newick",preserve_underscores=True)
		leaves=[ node.taxon.label for node in tree.leaf_node_iter()]
		item=""
		for item in leaves:
			if not  bool(re.match("^([1-9]+_[0-9]+_[0-9]+){1}",item)):
				break
		if not bool(re.match("^([1-9]+_[0-9]+_[0-9]+){1}",item)):
			message="\n\t{0}\n\t{1}\n\t{2}\n\t{3}\n\t{4}\n\t{5}\n\t".format(\
			"INDELible control file - Something's wrong!",\
			"Labels chosen for the tips of the tree are not correct.",\
			"Labels should follow this pattern: SpeciesID_LocusID_IndividualID",\
			"Where SpeciesID,LocusID,IndividualID are numbers.",\
			"Where SpeciesID > 1, LocusID > 0, IndividualID > 0",\
			"Please verify. Exiting."\
			)
			return False,message
		return True, message

	def checkAnchorTipLabelInTree(self):
		"""
		Verifies that the parameter given to select a root sequence is within
		the gene tree given.
		------------------------------------------------------------------------
		Returns:
		- boolean, message: the status of the process and the message related to
		such status
		"""
		# checking if label in the given tree
		self.appLogger.debug("Checking if label in the given tree")
		messageCorrect="Labels of the tree are correct"
		messageWrong="INDELible control file - Something's wrong!\n\t"
		tree=dendropy.Tree.get(path=self.geneTreeFile, schema="newick",preserve_underscores=True)
		filter = lambda taxon: True if taxon.label==self.anchorTipLabel else False
		node = tree.find_node_with_taxon(filter)
		if not node:
			messageWrong+="\n\t{0}\n\t{1}".format(\
			"Tip label for reference sequence is not in the Newick file.",\
			"Please verify. Exiting.")
			return False, messageWrong
		return True, messageCorrect

	def formatSettingsMessage(self):
		"""
		Prints all the settings.
		"""
		message="Settings:\n"
		sections=self.parser.sections()
		for sec in sections:
			message+="\t{0}\n".format(sec)
			items=self.parser.items(sec)
			for param in items:
				extra=""
				if param[0]=="inputmode":
					if int(param[1])==1: extra="(Gene-tree distribution - SimPhy output)"
					if int(param[1])==2: extra="(Single gene tree)"
					if int(param[1])==3: extra="(Single gene tree with ancestral sequence)"
				if param[0]=="ploidy":
					if int(param[1])==1: extra="(Haploid individuals)"
					if int(param[1])==2: extra="(Diploid individuals)"
				message+="\t\t{0}\t:\t{1} {2}\n".format(param[0],param[1],extra)
		return message

	def correctContentReferenceSequence(self):
		"""
		Verifies content of the reference sequence file for nucleotides only.
		------------------------------------------------------------------------
		Returns:
		- boolean, message: the status of the process and the message related to
		such status
		"""
		status=True; message=""
		referenceDict=msatools.parseMSAFileWithDescriptions(self.ancestralSequenceFilePath)
		reference=referenceDict[0]
		for item in reference:
			if not item.uppercase() in self.__NUCLEOTIDES:
				status=False
				message="\n\t{0}{1}{2}\n\t{3}\n\t{4}".format(
					"[data] block: Input mode (",self.inputmode,") selected but invalid option.",\
					"Reference sequence should be a nucleotidic sequence, but some other characters exist..",\
					"Please verify. Exiting."\
					)
				break
		return status, message
