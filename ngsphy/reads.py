#!/usr/bin/env python
# -*- coding: utf-8 -*-
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

	individualsFileNameSuffix="individuals"

	def __init__(self,settings):
		self.appLogger=logging.getLogger("ngsphy")
		self.appLogger.info('NGS read simulation: ART run started.')
		self.settings=settings
		self.numReplicates=self.settings.numReplicates
		self.numReplicateDigits=len(str(self.numReplicates))
		self.filteredReplicates=[ int(numST) for numST in self.settings.parser.get("general", "filtered_replicates").strip().split(",")]
		self.numLociPerReplicate=[int(numST) for numST in self.settings.parser.get("general","numlociperreplicate").split(",")]
		self.numLociPerReplicateDigits=[ len(str(item)) for item in self.numLociPerReplicate]
		cc=self.settings.parser.get("general", "numindividualsperreplicate").strip().split(",")
		self.numIndividualsPerReplicate=[ int(item) for item in cc if not item == ""]
		self.numIndividualsPerReplicateDigits=[len(str(item )) for item in self.numIndividualsPerReplicate]

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
		self.settings.coverageFolderPath=os.path.join(self.settings.outputFolderPath,"coverage")
		self.settings.readsFolderPath=os.path.join(self.settings.outputFolderPath,"reads")
		self.settings.scriptsFolderPath=os.path.join(self.settings.outputFolderPath,"scripts")
		self.generateFolderStructure()
		if (self.settings.runningTimes): self.runningInfo=RunningInfo()


	def generateFolderStructure(self):
		"""
		Generation of basic output folder structure for this process
		"""
		self.appLogger.info("Creating folder structure for [ngs-reads-art]")
		try:
			os.makedirs(self.settings.readsFolderPath)
			self.appLogger.info("Generating output folder ({0})".format(self.settings.readsFolderPath))
		except:
			self.appLogger.debug("Output folder exists ({0})".format(self.settings.readsFolderPath))

		try:
			os.makedirs(self.settings.scriptsFolderPath)
			self.appLogger.info("Generating output folder ({0})".format(self.settings.scriptsFolderPath))
		except:
			self.appLogger.debug("Output folder exists ({0})".format(self.settings.scriptsFolderPath))


	def writeSGEScript(self):
		"""
		For the generation of scripts for cluster environments (Execution section,
		options run: SGE)
		This generates job script that will be launched in a SGE cluster.
		"""
		self.appLogger.debug("Start")
		jobfile=os.path.join(\
			self.settings.scriptsFolderPath,\
			"{0}.job.sge.sh".format(self.settings.projectName)\
		)

		j=open(jobfile,"w")
		seedfile=os.path.join(\
			self.settings.scriptsFolderPath,\
			"{0}.sh".format(self.settings.projectName)\
		)
		script="""#!/bin/bash
# SGE submission options
#$ -l num_proc=1		 # number of processors to use
#$ -l h_rt=00:10:00	  # Set 10 mins  - Average amount of time for up to 1000bp
#$ -t 1-{0}			  # Number of jobs/files that will be treated
#$ -N art.sims		   # A name for the job

command=$(awk 'NR==$SGE_TASK_ID{{print $1}}' {1})
$command
""".format(self.numFiles,seedfile)
		j.write(script)
		j.close()
		self.appLogger.info("SGE Job file written ({0})...".format(jobfile))

	def writeSLURMScript(self):
		"""
		For the generation of scripts for cluster environments (Execution section,
		options run: SLURM)
		This generates job script that will be launched in a SLURM cluster.
		"""
		self.appLogger.debug("Start")
		jobfile=os.path.join(\
			self.settings.scriptsFolderPath,\
			"{0}.job.slurm.sh".format(self.settings.projectName)
		)
		j=open(jobfile,"w")
		seedfile=os.path.join(\
			self.settings.scriptsFolderPath,\
			"{0}.sh".format(self.settings.projectName)\
		)

		script="""#!/bin/sh
#SBATCH -n 1
#SBATCH --cpus-per-task 1
#SBATCH -t 00:10:00
#SBATCH --mem 4G
#SBATCH --array=1-{0}

command=$(awk 'NR==$SLURM_ARRAY_TASK_ID{{print $1}}' {1})
$command


""".format(self.numFiles, seedfile)
		footer="".format()
		jobfile.write(script)
		jobfile.close()
		self.appLogger.info("SLURM Job file written ({0})...".format(jobfile))

	def retrieveCoverageMatrix(self, indexREP):
		"""
		Reads coverage matrix for a specific species tree identifier
		------------------------------------------------------------------------
		Parameters:
		- indexREP:specific species tree identifier
		Returns:
		- coverage Matrix
		"""
		index=self.filteredReplicates.index(indexREP)
		self.appLogger.debug("Retrieving coverage matrix")
		# coverage matrix per ST - row:indv - col:loci
		# each cov introduced as parameter is a NGSPhyDistribution
		coverageMatrixFilename=os.path.join(\
			self.settings.coverageFolderPath,\
			"{0}.{1:0{2}d}.coverage.csv".format(\
				self.settings.projectName,\
				self.filteredReplicates[index],\
				self.numReplicateDigits\
			)
		)
		coverageMatrix=np.zeros(shape=(self.numIndividualsPerReplicate[index], self.numLociPerReplicate[index]))
		firstLine=True; counter=0
		with open(coverageMatrixFilename, 'rb') as csvfile:
			coveragereader = csv.reader(csvfile, delimiter=',')
			for row in coveragereader:
				if firstLine:
					firstLine=False
				else:
					coverageMatrix[counter,]=[float(row[i]) for i in range(1,len(row))]
					counter+=1
				if not counter < self.numIndividualsPerReplicate[index]: break
		return coverageMatrix

	def getCommands(self):
		"""
		Generates command lines that will be called to generate the reads
		"""
		self.appLogger.debug("Start")
		for index in range(0, len(self.filteredReplicates)):
			csvfilename=os.path.join(\
				self.settings.tablesFolderPath,\
				"{0}.{1:0{2}d}.{3}.csv".format(\
					self.settings.projectName,\
					self.filteredReplicates[index],\
					self.numReplicateDigits,\
					self.individualsFileNameSuffix\
				)
			)
			with open(csvfilename) as csvfile:
				d = csv.DictReader(csvfile)
				self.matingDict = [row for row in d]
				csvfile.close()
				nInds=len(self.matingDict)
				nLoci=self.numLociPerReplicate[index]
				coverageMatrix=self.retrieveCoverageMatrix(self.filteredReplicates[index])
				if (self.settings.ploidy == 2):
					coverageMatrix=coverageMatrix*0.5
				for indexLOC in range(1,self.numLociPerReplicate[index]+1):
					for row in self.matingDict:
						# indexREP,indexLOC,indID,speciesID,mateID1,mateID2
						inputFile=os.path.join(\
							self.settings.individualsFolderPath,\
							"REPLICATE_{0:0{1}d}".format(\
								int(row['repID']),\
								self.numReplicateDigits
							),\
							"LOCUS_{0:0{1}d}".format(\
								indexLOC,\
								self.numLociPerReplicateDigits[index]\
							),\
							"{0}_{1:0{2}d}_{3:0{4}d}_{5}_{6}.fasta".format(\
								self.settings.projectName,\
								int(row['repID']),\
								self.numReplicateDigits,\
								indexLOC,\
								self.numLociPerReplicateDigits[index],
								self.settings.simphyDataPrefix,\
								int(row['indID'])\
							)\
						)
						# This means, from a multiple (2) sequence fasta file.
						outputFile=os.path.join(\
							self.settings.readsFolderPath,\
							"REPLICATE_{0:0{1}d}".format(\
								int(row['repID']),\
								self.numReplicateDigits
							),\
							"LOCUS_{0:0{1}d}".format(\
								indexLOC,\
								self.numLociPerReplicateDigits[index]\
							),\
							"{0}_{1:0{2}d}_{3:0{4}d}_{5}_{6}_R".format(\
								self.settings.projectName,\
								int(row['repID']),\
								self.numReplicateDigits,\
								indexLOC,\
								self.numLociPerReplicateDigits[index],
								self.settings.simphyDataPrefix,\
								int(row['indID'])\
							)\
						)
						coverage=coverageMatrix[int(row['indID'])][indexLOC-1]
						coverageParam=""
						if self.settings.coveragemodeART:
							coverageParam="--rcount"
						else:
							coverageParam="--fcov"
						# Call to ART
						callParams=["art_illumina"]+self.params+[coverageParam,str(coverage),"--in", inputFile,"--out",outputFile]
						# self.params+=["--in ",inputFile,"--out",outputFile]
						# print(callParams)
						self.commands+=[[row['repID'],indexLOC,row['indID'],inputFile, outputFile]+callParams]

		self.appLogger.info("Commands have been generated...")

	def writeBashScript(self):
		"""
		Generates a bash script with all the command lines used to generate
		the reads (ART commands)
		"""
		self.appLogger.debug("Start")
		bashfile=os.path.join(\
			self.settings.scriptsFolderPath,\
			"{0}.sh".format(self.settings.projectName)\
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
		if (self.settings.runningTimes): self.runningInfo.addLine(line)


	def run(self):
		"""
		Process flow for the generation of NGS reads
		"""
		status=True
		environment=self.settings.executionMode
		# Generating commands
		try:
			self.getCommands()
			self.generateFolderStructureNGS()
			if environment==1:
				message="ART run has finished succesfully."
				if (self.settings.runART):
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
					if (self.settings.runningTimes): self.writeRunningInfoIntoFile()
				else:
					self.appLogger.info("Environment BASH. Writing scripts")
					self.writeBashScript()
			if (environment==2):
				self.writeBashScript()
				self.appLogger.info("Environment SGE. Writing scripts")
				self.writeSeedFile()
			if (environment==3):
				self.appLogger.info("Environment SLURM. Writing scripts")
				self.writeBashScript()
				self.writeSLURMScript()
		except ValueError as verror:
			status=False
			message="\t{0}\n\t{1}\n\t{2}".format(\
				"Distribution parameter error:",\
				verror,\
				"Please verify. Exciting."\
			)
		except RuntimeError as rte:
			status=False
			message="\t{0}\n\t{1}\n\t{2}".format(\
				"Execution error.",\
				rte,\
				"Please verify. Exciting."\
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
		f.write("repID,locID,indID,inputFile,cpuTime,seed,outputFilePrefix\n")
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
