import csv,logging,os,subprocess,sys,threading, time
import numpy as np
import random as rnd
from coverage import NGSPhyDistribution as ngsphydistro
import settings as sp

class RunningInfo(object):
	"""
	Class for the storage of running time information of ARTIllumina Class
	Separated to be able to be used by several threads.
	"""
	def __init__(self):
		# self.appLogger=logging.getLogger('sngsw')
		self.lock = threading.Lock()
		self.value = []

	def addLine(self, line):
		"""
		Adds a line to the running information estructure
		"""
		# self.appLogger.debug('Waiting for lock')
		self.lock.acquire()
		try:
			# self.appLogger.debug('Acquired lock')
			self.value += [line]
		finally:
			# self.appLogger.debug('Released lock')
			self.lock.release()

class ARTIllumina:
	"""
	Class for the generation of Illumina reads with the ART simulator.
	----------------------------------------------------------------------------
	Attributes:
	- CONSTANT: __SHORT_NAMES: all available coded-short parameters of the ART simulator
	- CONSTANT: __LONG_NAMES: all available  coded parameters of the ART simulator

	Technicallity. Parsing the simulator parameters from the settings file
	returns all of them into lower cases, in order to be able to compare them
	was necessarty to create these variable. Special cases are handled separate.

	- CONSTANT: __dLONG_NAMESHORT_NAMES: coded-short parameters of the ART simulator.
	- CONSTANT: __dLONG_NAMES: coded parameters of the ART simulator

	- params: final parameters used in the read simulation
	- commands: list of all the generated commands for the call of the ART simulator
	- numFiles: number of files that will be used for the simulation
	- individualsFolderName: pat that stores the individual FASTA files
	- coverageFolderPath: path that stores coverage matrices
	- readsFolderPath: path where reads generated will be stored
	- scriptsFolderPath: path for the folder that will store the generated scripts

	"""
	__SHORT_NAMES=["sf" ,"dp","ploidy","ofn","1","2","amp","c","d","ef","f","h","i",\
				"ir","ir2","dr","dr2","k","l","m","mp","M","nf","na",\
				"o","p","q","qU","qs","qL","qs2","rs","s","sam","sp","ss"]
	__LONG_NAMES=["simphy_folder","data_prefix","output_folder_name","ploidy","qprof1","qprof2",\
				"amplicon","rcount","id","errfree","fcov","help",\
				"in","insRate","insRate2","delRate","delRate2","maxIndel",\
				"len","mflen","matepair","cigarM","maskN","noALN","out",\
				"paired","quiet","maxQ","qShift","minQ","qShift2","rndSeed",\
				"sdev","samout","sepProf","seqSys"]
	__dLONG_NAMES={i.lower():i for i in __LONG_NAMES}
	__dSHORT_NAMES={i.lower():i for i in __SHORT_NAMES}

	params=[]
	commands=[]	 # Init of all the commands that will be generated
	numFiles=0

	# file path related variables
	individualsFolderName="individuals"
	coverageFolderPath=""
	readsFolderPath=""
	scriptsFolderPath=""


	def __init__(self,settings):
		self.appLogger=logging.getLogger('ngsphy')
		self.appLogger.info('NGS read simulation: ART run started.')
		self.settings=settings
		self.numSpeciesTrees=self.settings.parser.getint("general","numspeciestrees")
		self.numSpeciesTreesDigits=len(str(self.numSpeciesTrees))
		self.filteredST=[ int(numST) for numST in self.settings.parser.get("general", "filtered_ST").strip().split(",")]
		self.numLociPerSpeciesTree=[int(numST) for numST in self.settings.parser.get("general","numLociPerSpeciesTree").split(",")]
		self.numLociPerSpeciesTreeDigits=[ len(str(item)) for item in self.numLociPerSpeciesTree]
		cc=self.settings.parser.get("general", "numIndividualsPerSpeciesTree").strip().split(",")
		self.numIndividualsPerSpeciesTree=[ int(item) for item in cc if not item == ""]
		self.numIndividualsPerSpeciesTreeDigits=[len(str(item )) for item in self.numIndividualsPerSpeciesTree]

		dash=""; par=[]
		settingsParams=self.settings.parser.items("ngs-reads-art")
		for p in settingsParams:
			if (p[0] in self.__dSHORT_NAMES.keys()): dash="-"
			if (p[0] in self.__dLONG_NAMES.keys()): dash="--"
			# to be sure that i am getting the right parameter names
			if (dash=="-"):
				if (p[0]=="m" and ((par[1].lower() in ["true","false","on","off"]) or (par[1] in [0,1]))):
					par=["M"]
				elif (p[0]=="m"):
					par=["m"]
				else:
					par=[self.__dSHORT_NAMES[p[0]]]
			if (dash=="--"):
				par=[self.__dLONG_NAMES[p[0]]]

			par+=[p[1]]

			if((par[1].lower() in ["true","false","on","off"]) or (par[1] in [0,1])):
				self.params+=["{0}{1}".format(dash,par[0])]
			else:
				self.params+=["{0}{1}".format(dash,par[0]),par[1]]

		# generating specific folder structure
		self.generateFolderStructure()
		self.runningInfo=RunningInfo()


	def generateFolderStructure(self):
		"""
		Generation of basic output folder structure for this process
		"""
		self.appLogger.info("Creating folder structure for [ngs-reads-art]")
		self.readsFolderPath=os.path.join(self.settings.outputFolderPath,"reads")
		self.scriptsFolderPath=os.path.join(self.settings.outputFolderPath,"scripts")
		self.coverageFolderPath=os.path.join(self.settings.outputFolderPath,"coverage")
		try:
			os.makedirs(self.readsFolderPath)
			self.appLogger.info("Generating output folder ({0})".format(self.readsFolderPath))
		except:
			self.appLogger.debug("Output folder exists ({0})".format(self.readsFolderPath))

		try:
			os.makedirs(self.scriptsFolderPath)
			self.appLogger.info("Generating output folder ({0})".format(self.scriptsFolderPath))
		except:
			self.appLogger.debug("Output folder exists ({0})".format(self.scriptsFolderPath))


	def writeSeedFile(self):
		"""
		For the generation of scripts for cluster environments (Execution section,
		options run: SGE/SLUM)
		This generates the seed file for the calling of jobs for cluster environments.
		"""
		self.appLogger.debug("Start")
		seedfile="{0}/scripts/{1}.seedfile.txt".format(\
			self.settings.outputFolderPath,\
			self.settings.projectName
		)
		sf=open(seedfile,"w")
		for indexST in self.filteredST:
			for indexST in self.filteredST:
				csvfile=open("{0}/tables/{1}.{2:0{3}d}.{4}.csv".format(\
					self.settings.outputFolderPath,\
					self.settings.projectName,\
					indexST,\
					self.numSpeciesTreesDigits,\
					self.individualsFolderName
				))
				# Generation of folder structure
				d = csv.DictReader(csvfile)
				self.matingDict = [row for row in d]
				csvfile.close()
				for indexLOC in range(1,self.numLociPerSpeciesTree[indexST-1]+1):
					for row in self.matingDict:
						# indexST,indexLOC,indID,speciesID,mateID1,mateID2
						inputFile="{0}/individuals/{1}/{2:0{3}d}/{4}_{1}_{2:0{3}d}_{5}_{6}.fasta".format(\
							self.settings.outputFolderPath,\
							row['indexST'],\
							indexLOC,\
							self.numLociPerSpeciesTreeDigits[indexST-1],
							self.settings.projectName,\
							self.settings.dataPrefix,\
							row['indID']\
						)
						# This means, from a multiple (2) sequence fasta file.
						outputFile="{0}/reads/{1}/{2:0{3}d}/{4}_{1}_{2:0{3}d}_{5}_{6}_R".format(\
							self.settings.outputFolderPath,\
							row['indexST'],\
							indexLOC,\
							self.numLociPerSpeciesTreeDigits[indexST-1],
							self.settings.projectName,\
							self.settings.dataPrefix,\
							row['indID']\
						)
						sf.write("{0}\t{1}\n".format(inputFile,outputFile))
				self.numFiles+=1
		sf.close()
		self.appLogger.info("Seed file written...")

	def writeSGEScript(self):
		"""
		For the generation of scripts for cluster environments (Execution section,
		options run: SGE)
		This generates job script that will be launched in a SGE cluster.
		"""
		self.appLogger.debug("Start")
		jobfile="{0}/scripts/{1}.job.sge.sh".format(\
			self.settings.outputFolderPath,\
			self.settings.projectName
		)
		j=open(jobfile,"w")
		seedfile="{0}/scripts/{1}.seedfile.txt".format(\
			self.settings.outputFolderPath,\
			self.settings.projectName
		)

		inputFile="$inputfile"
		outputFile="$outputfile"
		callParams=["art_illumina"]+self.params+["--in", inputFile,"--out",outputFile]
		header="""#!/bin/bash
# SGE submission options
#$ -l num_proc=1		 # number of processors to use
#$ -l h_rt=00:10:00	  # Set 10 mins  - Average amount of time for up to 1000bp
#$ -t 1-{0}			  # Number of jobs/files that will be treated
#$ -N art.sims		   # A name for the job

inputfile=$(awk 'NR==$SGE_TASK_ID{{print $1}}' {1})
outputfile=$(awk 'NR==$SGE_TASK_ID{{print $2}}' {1})\n""".format(self.numFiles,seedfile)
		j.write(header)
		j.write(" ".join(callParams))
		footer="".format()
		j.write(footer)
		j.close()
		self.appLogger.info("SGE Job file written ({0})...".format(jobfile))

	def writeSLURMScript(self):
		"""
		For the generation of scripts for cluster environments (Execution section,
		options run: SLURM)
		This generates job script that will be launched in a SLURM cluster.
		"""
		self.appLogger.debug("Start")
		jobfile="{0}/scripts/{1}.job.slurm.sh".format(\
			self.settings.outputFolderPath,\
			self.settings.projectName
		)
		j=open(jobfile,"w")
		seedfile="{0}/scripts/{1}.seedfile.txt".format(\
			self.settings.outputFolderPath,\
			self.settings.projectName
		)
		inputFile="$inputfile"
		outputFile="$outputfile"
		callParams=["art_illumina"]+self.params+["--in", inputFile,"--out",outputFile]

		header="""#!/bin/sh
#SBATCH -n 1
#SBATCH --cpus-per-task 1
#SBATCH -t 00:10:00
#SBATCH --mem 4G
#SBATCH --array=1-{0}

inputfile=$(awk 'NR==$SLURM_ARRAY_TASK_ID{{print $1}}' {1})
outputfile=$(awk 'NR==$SLURM_ARRAY_TASK_ID{{print $2}}' {1})

""".format(self.numFiles, seedfile)
		footer="".format()
		jobfile.write(header)
		jobfile.write(" ".join(callParams))
		jobfile.write(footer)
		jobfile.close()
		self.appLogger.info("SLURM Job file written ({0})...".format(jobfile))

	def retrieveCoverageMatrix(self, indexST):
		"""
		Reads coverage matrix for a specific species tree identifier
		------------------------------------------------------------------------
		Parameters:
		- indexST:specific species tree identifier
		Returns:
		- coverage Matrix
		"""
		self.appLogger.debug("Retrieving coverage matrix")
		# coverage matrix per ST - row:indv - col:loci
		# each cov introduced as parameter is a NGSPhyDistribution
		coverageMatrixFilename=os.path.join(\
			self.coverageFolderPath,\
			"{0}.{1:0{2}d}.coverage.csv".format(\
				self.settings.projectName,\
				indexST,\
				self.numSpeciesTreesDigits\
			)
		)

		coverageMatrix=np.zeros(shape=(self.numIndividualsPerSpeciesTree[indexST-1], self.numLociPerSpeciesTree[indexST-1]))
		firstLine=True; counter=0
		with open(coverageMatrixFilename, 'rb') as csvfile:
			coveragereader = csv.reader(csvfile, delimiter=',')
			for row in coveragereader:
				if firstLine:
					firstLine=False
				else:
					coverageMatrix[counter,]=[float(row[index]) for index in range(1,len(row))]
					counter+=1
				if not counter < self.numIndividualsPerSpeciesTree[indexST-1]: break
		return coverageMatrix

	def getCommands(self):
		"""
		Generates command lines that will be called to generate the reads
		"""
		self.appLogger.debug("Start")
		for indexST in self.filteredST:
			csvfile=open("{0}/tables/{1}.{2:0{3}d}.{4}.csv".format(\
				self.settings.outputFolderPath,\
				self.settings.projectName,\
				indexST,\
				self.numSpeciesTreesDigits,\
				self.individualsFolderName
			))
			# Generation of folder structure
			d = csv.DictReader(csvfile)
			self.matingDict = [row for row in d]
			csvfile.close()
			nInds=len(self.matingDict)
			nLoci=self.numLociPerSpeciesTree[indexST-1]
			coverageMatrix=self.retrieveCoverageMatrix(indexST)
			for indexLOC in range(1,self.numLociPerSpeciesTree[indexST-1]+1):
				for row in self.matingDict:
					# indexST,indexLOC,indID,speciesID,mateID1,mateID2
					inputFile="{0}/individuals/{1:0{2}d}/{3:0{4}d}/{5}_{1}_{3:0{4}d}_{6}_{7}.fasta".format(\
						self.settings.outputFolderPath,\
						int(row['indexST']),\
						self.numSpeciesTreesDigits,\
						indexLOC,\
						self.numLociPerSpeciesTreeDigits[indexST-1],
						self.settings.projectName,\
						self.settings.dataPrefix,\
						int(row['indID']))
					# This means, from a multiple (2) sequence fasta file.
					outputFile="{0}/reads/{1:0{2}d}/{3:0{4}d}/{5}_{1}_{3:0{4}d}_{6}_{7}_R".format(\
						self.settings.outputFolderPath,\
						int(row['indexST']),\
						self.numSpeciesTreesDigits,\
						indexLOC,\
						self.numLociPerSpeciesTreeDigits[indexST-1],
						self.settings.projectName,\
						self.settings.dataPrefix,\
						int(row['indID'])\
					)
					coverage=coverageMatrix[int(row['indID'])][indexLOC-1]
					# Call to ART
					callParams=["art_illumina"]+self.params+["--fcov",str(coverage),"--in", inputFile,"--out",outputFile]
					# self.params+=["--in ",inputFile,"--out",outputFile]
					# print(callParams)
					self.commands+=[[row['indexST'],indexLOC,row['indID'],inputFile, outputFile]+callParams]

		self.appLogger.info("Commands have been generated...")

	def writeBashScript(self):
		"""
		Generates a bash script with all the command lines used to generate
		the reads (ART commands)
		"""
		self.appLogger.debug("Start")
		bashfile="{0}/scripts/{1}.sh".format(\
			self.settings.outputFolderPath,\
			self.settings.projectName
		)
		j=open(bashfile,"w")
		for item in self.commands:
			c=item[5:len(item)]
			j.write(" ".join(c))
			j.write("\n")
		j.close()

	def commandLauncher(self, command):
		"""
		Launches the execution of the third-party app for the simulation
		of  NGS reads.
		"""
		self.appLogger.debug("Start")
		ngsMessage="";proc=""
		try:
			proc = subprocess.check_output(command[5:],stderr=subprocess.STDOUT)
			cpuTime = [line for line in proc.split('\n') if "CPU" in line][0].split(":")[1]
			seed = [line for line in proc.split('\n') if "seed" in line][0].split(":")[1]
			# line=[(command[0:3],cpuTime,seed, command[4])]
			line=command[0:3]+[cpuTime,seed,command[4]]
			# ngsMessage="Command: {0} - Finished succesfully.".format(" ".join(command[5:]))
		except subprocess.CalledProcessError as error:
			ngsMessage="{}".format(error.output)+\
			"\n\n------------------------------------------------------------------------"+\
			"\n\nFor more information about this error please run the 'art' command separately.\n"+\
			"art_illumina command used:\n==========================\n"+\
			"{}\n\n".format(" ".join(command[5:]))
			self.appLogger.error(ngsMessage)
			line=command[0:3]+["-","-",command[4]]
			raise RuntimeError("\nART execution error. Please verify. Exciting.")
		self.runningInfo.addLine(line)


	def run(self):
		"""
		Process flow for the generation of NGS reads
		"""
		status=True
		environment=self.settings.parser.get("execution","environment")
		# Generating commands
		try:
			self.getCommands()
			if (environment=="sge"):
				self.appLogger.info("Environment SGE. Writing scripts")
				self.writeSeedFile()
				self.writeSGEScript()
			elif (environment=="slurm"):
				self.appLogger.info("Environment SLURM. Writing scripts")
				self.writeSeedFile()
				self.writeSLURMScript()
			else:
				message="ART run has finished succesfully."
				self.appLogger.info("Environment BASH. Writing scripts")
				self.writeBashScript()
				run=self.settings.parser.getboolean("execution","run")
				if (run):
					self.generateFolderStructureNGS()
					curr=0
					self.appLogger.info("Running...")
					while (curr < len(self.commands)):
						progress=(curr*100)/len(self.commands)
						sys.stdout.write("Progress {0:02.1f} %\r".format(progress))
						sys.stdout.flush()
						# print("Progress {0:02.1f} %\r".format(prog))
						command=self.commands[curr]
						t = threading.Thread(target=self.commandLauncher(command))
						t.start()
						curr=curr+1
						while (threading.activeCount()-1)==self.settings.numThreads:
							time.sleep(0.1)

					self.writeRunningInfoIntoFile()

		except ValueError as verror:
			status=False
			message="\n\t{0}\n\t{1}\n\t{2}".format(\
				"Distributon parameter error:",\
				verror,\
				"Please verify. Exciting."\
			)
		except RuntimeError as rte:
			status=False
			message="\n\t{0}\n\t{1}\n\t{2}\n\t{3}".format(\
				"Distributon parameter error.",\
				rte,\
				"Please verify. Exciting.",\
			)

		return status, message

	def generateFolderStructureNGS(self):
		"""
		Generation of detail folder structure used for this process
		"""
		# iterating over commands to create folders
		folders=set([])
		for command in self.commands:
			infile=os.path.dirname(command[4])
			folders.add(infile)
		for item in folders:
			try:
				os.makedirs(item)
			except:
				self.appLogger.debug("Output folder exists ({0})".format(item))

	def writeRunningInfoIntoFile(self):
		"""
		Writes running time information into a file
		"""
		outputFile=os.path.join(self.settings.outputFolderPath,"{0}.info".format(self.settings.projectName))
		f=open(outputFile,"w")
		f.write("indexST,indexLOC,indID,inputFile,cpuTime,seed,outputFilePrefix\n")
		for item in self.runningInfo.value:
			f.write(
				str(item[0])+","+\
				str(item[1])+","+\
				str(item[2])+","+\
				str(item[3])+","+\
				str(item[4])+","+\
				item[5]+"\n"
			)
		f.close()
		self.appLogger.info("File with timings of the ART run can be find on: {0}".format(outputFile))
