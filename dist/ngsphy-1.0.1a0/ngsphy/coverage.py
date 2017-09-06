import csv,datetime,logging,os,subprocess,sys,threading,time
import numpy as np
import random as rnd
from scipy.stats import  binom,expon,gamma,lognorm,norm,nbinom,poisson,uniform

class NGSAvailableDistributions:
	"""
	Class for the definition of the available distributions for this program.
	----
	Attributes:
	- relationNumParams: Dictionary with the pairing code:number_of_parameters.
	- distributions: Codes of all the available distributions.
	- names: Full name of the available distributions.
	"""

	relationNumParams=dict({\
		"b":2,#mean,percentages\
		"e":1,#mean\
		"f":1,\
		"g":2,\
		"g1":1,\
		"ln":2,\
		"n":2,\
		"nb":2,\
		"p":1,\
		"u":1\
		})
	distributions=relationNumParams.keys()
	names=dict({\
		"b":"Binomial",\
		"e":"Exponential",\
		"f":"Fixed",\
		"g":"Gamma",\
		"g1":"Gamma (mean=1)",\
		"ln":"Log Normal",\
		"n":"Normal",\
		"nb":"Negative Binomial",\
		"p":"Poisson",\
		"u":"Uniform"\
	})

class NGSPhyDistributionParser:
	"""
	Class for parsing the distribution "words" introduced in the settings file.
	----------------------------------------------------------------------------
	Attributes:
	- __params: to store distribution parameters
	- __name: to store distribution code name
	- __dependency: to identify whether the whole word should be taking into
	account or if the mean (or main) parameter of the distribution depends on
	previously calculated values.
	"""
	__params=[]
	__name=""
	__dependency=None

	def __init__(self,distroline,dependency):
		distroline=distroline.split(":")
		self.__dependency=dependency or False
		# i depend on another level value
		if self.__dependency:
			# Value is empty
			if (len(distroline)==1):
				# only have distribution name, meaning
				# distribution has only 1 param and it comes
				# from the upperLevelDistribution
				self.__name=distroline[0].lower()
			if (len(distroline)==2):
				# i have name and a parameter,
				# meaning, first param comes from above
				# second is the second parameter of my distro
				self.__name=distroline[0].lower()
				self.__params=distroline[1].split(",")
		else:
			self.__name=distroline[0].lower()
			if len(distroline) > 1:
				self.__params= distroline[1].split(",")


	def setParams(self,params):
		"""
		Allows to modify the values of the distribution parameters
		-----------------------------------------------------------------------
		- input: params
		"""
		self.__params=params

	def getParams(self):
		"""
		Returns distribution parameters
		"""
		return self.__params

	def printParams(self):
		"""
		Prints distribution parameters
		"""
		print(self.__params)

	def setName(self,name):
		"""
		Allows to modify the code name of the distribution
		-----------------------------------------------------------------------
		- input: name
		"""
		self.__name=name

	def getName(self):
		"""
		Returns distribution parameters
		"""
		return self.__name

	def printName(self):
		"""
		Prints distribution name
		"""
		print(self.__name)

	def asNGSPhyDistribution(self):
		"""
		Returns a NGSPhyDistribution object from the parsed distribution word.
		"""
		if not self.__dependency:
			return NGSPhyDistribution(self.__name,self.__params)
		else: return None

	def validate(self):
		"""
		Validates if the input of the distribution word correspond to its definition.
		Meaning, codename is within the available distribution code names,
		number of parameters correspond to the selected distribution and the
		parameters introduced are numbers.
		"""
		message="Coverage distribution parsing OK.\n\t"
		status=True
		if not self.__name in NGSAvailableDistributions.distributions:
			message="{0}\n\t{1}\n\t{2}".format(\
				"Problem parsing coverage distribution.",\
				"Distribution selected invalid."
				"Please verify. Exiting"\
				)
			status=False
		if (NGSAvailableDistributions.relationNumParams[self.__name]==len(self.__params)) or (NGSAvailableDistributions.relationNumParams[self.__name]-1==len(self.__params) and self.__dependency):
			# print(len(self.__params))
			for index in range(0,len(self.__params)):
				try:
					self.__params[index]=float(self.__params[index])
				except ValueError:
					message="{0}\n\t{1}\n\t{2}".format(\
					"Problem parsing coverage distribution.",\
					"Not a correct parameter value.",\
					"Please verify. Exiting"\
					)
					status=False
					break
		else:
			message="{0}\n\t{1}\n\t{2}\t{3} ({4}) - Num params: {5}\n\t{6}\t{7}\n\t{8}".format(\
				"Problem parsing coverage distribution.",\
				"Not the correct number of parameters for the selected distribution",\
				"Given: ",NGSAvailableDistributions.names[self.__name],self.__name,len(self.__params),\
				"Required: ",NGSAvailableDistributions.relationNumParams[self.__name],\
				"Please verify. Exiting"\
				)
			status=False

		return status, message

class NGSPhyDistribution:
	"""
	Class that samples values from the available distributions.
	----------------------------------------------------------------------------
	Attributes:
	- __name: Code name of the distribution
	- __params: Parameters of the distribution to be sampled.
	"""
	__name=None
	__params=[]
	def __init__(self,name,params):
		self.__name=name
		self.__params=params
		# print(self.__name)
		# print(self.__params)

	def value(self,samples=1):
		"""
		Samples number of values given from the specific distribution.
		------------------------------------------------------------------------
		- samples: number of values that will be returned.
		"""
		value=0
		try:
			for item in self.__params:
				if item==0: break
			if item==0: value=[0]*samples
			else:
				if self.__name=="b": value=self.binom(samples)
				if self.__name=="e": value=self.exponential(samples)
				if self.__name=="f": value=self.fixed(samples)
				if self.__name=="g": value=self.gamma(samples)
				if self.__name=="g1": value=self.gamma1(samples)
				if self.__name=="ln": value=self.lognormal(samples)
				if self.__name=="n": value=self.normal(samples)
				if self.__name=="nb": value=self.nbinom(samples)
				if self.__name=="p": value=self.poisson(samples)
				if self.__name=="u": value=self.uniform(samples)
		except Exception as ex:
			exc_type, exc_obj, exc_tb = sys.exc_info()
			fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
			message="\n\tUnexpected: {0} | {1} - File: {2} - Line:{3}".format(\
				ex,exc_type, fname, exc_tb.tb_lineno)
			status=False
			raise Exception(message)
			# self.appLogger.error(message)
			# sys.exit()

		return value

	def binom(self,samples):
		"""
		Sampling from a binomial distribution
		------------------------------------------------------------------------
		- samples: number of values that will be returned.
		"""
		#n=number of times tested
		#p=probability
		n=self.__params[0]
		p=self.__params[1]
		distro=binom(n,p)
		f=distro.rvs(size=samples)
		return f

	def exponential(self,samples):
		"""
		Sampling from a exponential distribution
		------------------------------------------------------------------------
		- samples: number of values that will be returned.
		"""
		f=np.random.exponential(float(self.__params[0]),size=samples)
		return value

	def fixed(self,samples):
		"""
		Fixed values
		------------------------------------------------------------------------
		- samples: number of values that will be returned.
		"""
		value=[self.__params[0]]*samples
		return value

	def gamma1(self,samples):
		"""
		Sampling from a Gamma distribution with mean 1
		------------------------------------------------------------------------
		- samples: number of values that will be returned.
		- gamma parameterized with shape and scale
		"""
		# E|x| = k.theta (alpha*theta)
		# If i want  mean=1, theta=E|x|/alpha=1/alpha
		shape=float(self.__params[0]*1.0)
		theta=float(1/(shape*1.0))
		distro=gamma(shape,theta)
		f=distro.rvs(size=samples)
		return f


	def gamma(self,samples):
		"""
		Sampling from a Gamma distribution with mean 1

		The parameterization with alpha and beta is more common in Bayesian statistics,
		where the gamma distribution is used as a conjugate prior distribution
		for various types of inverse scale (aka rate) parameters, such as the
		lambda of an exponential distribution or a Poisson distribution[4]  or for
		t hat matter, the beta of the gamma distribution itself.
		(The closely related inverse gamma distribution is used as a conjugate
		prior for scale parameters, such as the variance of a normal distribution.)
		shape, scale = 2., 2. # mean=4, std=2*sqrt(2)

		(Wikipedia: https://en.wikipedia.org/wiki/Gamma_distribution)

		s = np.random.gamma(shape, scale, 1000)
		E|x| = k.theta (alpha*theta)
		If i want a specific mean, theta=E|x|/alpha
		------------------------------------------------------------------------

		- samples: number of values that will be returned.
		"""
		shape=float(self.__params[0]*1.0)
		theta=float(self.__params[1]*1.0)
		distro=gamma(shape,theta)
		f=distro.rvs(size=samples)
		return f


	def lognormal(self,samples):
		"""
		Sampling from a Log Normal distribution

		Parameters:
		m=mean
		s=sigma - standard deviation
		------------------------------------------------------------------------
		- samples: number of values that will be returned.
		"""
		mean=float(self.__params[0]*1.0)
		sigma=float(self.__params[1]*1.0)
		distro=lognorm(mean,sigma)
		f=distro.rvs(size=samples)
		return f

	def nbinom(self,samples):
		"""
		Sampling from a Negative binomial distribution

		Parameters:
		mu= poissonon mean
		r controls the deviation from the poisson

		This makes the negative binomial distribution suitable as a robust
		alternative to the Poisson, which approaches the Poisson for large r,
		but which has larger variance than the Poisson for small r.
		------------------------------------------------------------------------
		- samples: number of values that will be returned.
		"""
		mu=float(self.__params[0])
		r=float(self.__params[1])
		p=(r*1.0)/(r+mu)
		distro=nbinom(r,p)
		f=distro.rvs(size=samples)
		return f

	def normal(self,samples):
		"""
		Sampling from a Normal distribution

		Parameters:
		mean (location)
		variance (squared scale)
		------------------------------------------------------------------------
		- samples: number of values that will be returned.
		"""
		locMean=float(self.__params[0]*1.0)
		scaleVariance=float(self.__params[1]*1.0)
		distro=norm(loc=locMean,scale=scaleVariance)
		f=distro.rvs(size=samples)
		return f

	def poisson(self,samples):
		"""
		Sampling from a Poisson distribution
		Parameters:
		mean
		------------------------------------------------------------------------
		- samples: number of values that will be returned.
		"""
		l=float(self.__params[0]*1.0)
		distro=poisson(l)
		f=distro.rvs(size=samples)
		return f

	def uniform(self,samples):
		"""
		Sampling from a Poisson distribution
		Parameters:
		meean
		------------------------------------------------------------------------
		- samples: number of values that will be returned.
		"""
		mean=float(self.__params[0]*1.0)
		distro=uniform(mean)
		f=distro.rvs(size=samples)
		return f

class CoverageMatrixGenerator:
	"""
	Class for the calculation of coverage matrices.
	Able to intoduced:
	 	- "genomic noise" for entire datasets.
		- Variation across individuals and loci.
		- on/off target effect
		- coverage decay related to distance to the reference used
	----------------------------------------------------------------------------
	Attributes:
	- OFFTARGET_COVERAGE_FRACTION: Specific for the on/off target effect.
	- appLogger: Logger object
	- settings: Settings object withh all the program parameters
	- alphashapesLocus: parameters (alpha shapes) for the Gamma distribution with
	mean 1 that will be used as multiplier for across loci coverage variation.
	- alphashapesIndividuals: parameters (alpha shapes) for the Gamma distribution
	with mean 1 that will be used as multiplier for across individuals coverage
	variation.
	- locusMultiplier: parameters to be used to sample alphashapesLocus.
	- individualsMultiplier: parameters to be used to sample alphashapesIndividuals.
	- experiment: expected coverage for a specific replicate (species tree data).
	- outputFolderPath: folder where coverage matrices will be stored.
	- numReplicate: number of species trees.
	- numLociPerReplicate: number of loci per species tree.
	- numIndividualsPerReplicate: number of individuals per species tree.
	- numReplicateDigits: number of digits needed to represent numReplicate.
	- numLociPerReplicateDigits: number of digits needed to represent.
	- filteredReplicates: identifier of the species trees that will be used.

	"""
	# OFFTARGET_COVERAGE_FRACTION=0.1#
	appLogger=None
	settings=None

	alphashapesLocus=[]
	alphashapesIndividuals=[]
	locusMultiplier=[]
	individualsMultiplier=[]

	experiment=None

	outputFolderPath=""

	numReplicate=0
	numLociPerReplicate=[]
	numIndividualsPerReplicate=[]
	numReplicateDigits=0
	numLociPerReplicateDigits=[]
	filteredReplicates=[]

	def __init__(self, settings):
		self.appLogger=logging.getLogger('ngsphy')
		self.appLogger.info("Coverage calculations")
		self.settings=settings

		self.experiment=self.settings.experiment.asNGSPhyDistribution()
		self.numReplicate=self.settings.parser.getint("general", "numreplicates")
		cc=self.settings.parser.get("general", "numLociPerReplicate").strip().split(",")
		self.numLociPerReplicate=[ int(item) for item in cc if not item ==""]
		cc=self.settings.parser.get("general", "numIndividualsPerReplicate").strip().split(",")
		self.numIndividualsPerReplicate=[ int(item) for item in cc if not item == ""]
		self.numReplicateDigits=len(str(self.numReplicate))
		self.numLociPerReplicateDigits=[len(str(item)) for item in self.numLociPerReplicate]
		self.numIndividualsPerReplicateDigits=[len(str(item )) for item in self.numIndividualsPerReplicate]
		self.filteredReplicates=[int(item) for item in self.settings.parser.get("general", "filtered_replicates").strip().split(",")]
		if (self.settings.locus):
			distro=self.settings.locus.asNGSPhyDistribution()
			self.alphashapesLocus=distro.value(self.settings.numReplicate)
			self.locusMultiplier=[NGSPhyDistribution("g1",[self.alphashapesLocus[index]]).value(self.numLociPerReplicate[index]) for index in range(0,self.numReplicate)]
		if (self.settings.individual):
			distro=self.settings.individual.asNGSPhyDistribution()
			self.alphashapesIndividuals=distro.value(self.settings.numReplicate)
			self.individualsMultiplier=[NGSPhyDistribution("g1",[self.alphashapesIndividuals[index]]).value(self.numLociPerReplicate[index]) for item in range(0,self.numReplicate)]
		self.generateFolderStructure()

	def generateFolderStructure(self):
		"""
		Generation of folder structure needed for the coverage matrices.
		"""
		self.appLogger.debug("Coverage folder structure generation...")
		try:
			os.makedirs(self.settings.coverageFolderPath)
			self.appLogger.info("Generating output folder ({0})".format(self.settings.coverageFolderPath))
		except:
			self.appLogger.debug("Output folder exists ({0})".format(self.settings.coverageFolderPath))

	def calculate(self):
		"""
		Calulation of coverage matrices
		"""
		message=""
		status=True;
		self.appLogger.debug("Coverage calculations...")
		for indexRep in self.filteredReplicates:
			nInds=self.numIndividualsPerReplicate[indexRep-1]
			nLoci=self.numLociPerReplicate[indexRep-1]
			# expectedCoverage=self.experiment.value(nInds*nLoci)
			val=self.experiment.value(1)
			# self.appLogger.debug(val)
			expectedCoverage=val*(nInds*nLoci)
			coverageMatrix=np.matrix(expectedCoverage)
			coverageMatrix.shape=[nInds,nLoci]
			# individuals + loci coverage variation
			# individuals + loci multipliers
			try:
				if self.settings.locus:
					multipliers=self.locusMultiplier[indexRep-1]
					for loc in range(0,nLoci):
						coverageMatrix[:, loc]=coverageMatrix[:,loc]*multipliers[loc]
				if self.settings.individual:
					multipliers=self.individualsMultiplier[indexRep-1]
					for ind in range(0,nInds):
						coverageMatrix[ind,]=coverageMatrix[ind,]*multipliers[ind]
			except Exception as ex:
				exc_type, exc_obj, exc_tb = sys.exc_info()
				fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
				message="\n\t{0} {1} | {2} - File: {3} - Line:{4}".format(\
					"Unexpected. Coverage computation:",
					ex,exc_type, fname, exc_tb.tb_lineno)
				return False, message
			# taxonomic variation
			# I need the file with the realtion of the individuals
			# dont car if haploid or diploid because ill only need the first 3 columns
			# which match in both cases
			if self.settings.taxon:
				individualsTableFilename=os.path.join(\
					self.settings.outputFolderPath,\
					"tables",
					"{0}.{1:0{2}d}.individuals.csv".format(\
						self.settings.projectName,\
						indexRep,\
						self.numReplicateDigits)
					)
				individualsTableFile=open(individualsTableFilename)
				d = csv.DictReader(individualsTableFile)
				individualsTable = [row for row in d]
				individualsTableFile.close()
				for row in individualsTable:
					if str(row['speciesID']) in self.settings.taxon.keys():
						coverageMatrix[ind,]=np.array(coverageMatrix[ind,])*self.settings.taxon[str(row['speciesID'])]
			#on/off target
			onTargetLoci=range(0,nLoci)
			offTargetLoci=[]
			try:
				if self.settings.offtarget["loci"] > 0:
					nSamples=int(self.settings.offtarget["loci"]*nLoci)
					offTargetLoci=np.random.choice(nLoci,nSamples,replace=F)
					onTargetLoci=set(range(0,nLoci))-set(offTargetLoci)
					for loc in offTargetLoci:
						column=np.array(np.array(coverageMatrix[:,loc])*self.offtarget["coverage"])
						coverageMatrix[:,loc]=np.transpose(column)
			except Exception as ex:
				return False, "".format(\
					"Coverage Matrix Generator",\
					"An unexpected error occurr with the on/off target parameter. ",\
					"Please verify. Exiting.")
			# not captured
			try:
				if self.settings.notcaptured:
					# random number of loci  (respecting probabilities)
					nSamples=int(self.settings.notcaptured)*(nLoci-len(offTargetLoci))
					notcaptured=np.random.choice(onTargetLoci,nSamples,replace=F)
					for loc in notcaptured:
						column=np.array(coverageMatrix[:,loc])*0
						coverageMatrix[:loc]=np.transpose(column)
			except Exception as ex:
				return False, "".format(\
					"Coverage Matrix Generator",\
					"An unexpected error occurr with the not captured parameter. ",\
					"Please verify. Exiting.")
			# print(coverageMatrix)
			self.write(coverageMatrix,indexRep)
		return status,message

	def write(self,coverageMatrix,indexRep):
		"""
		Writing into file the coverage matrix
		"""
		self.appLogger.debug("Writing coverage matrix for replicate: {0:0{1}d}".format(indexRep, self.numReplicateDigits))
		filename=os.path.join(self.settings.coverageFolderPath,"{0}.{1:0{2}d}.coverage.csv".format(\
			self.settings.projectName,\
			indexRep,\
			self.numReplicateDigits\
			)
		)
		filepath=os.path.abspath(filename)
		header=["indID"]+[ "L.{0:0{1}d}".format(loc+1,self.numLociPerReplicateDigits[indexRep-1]) for loc in range(0,self.numLociPerReplicate[indexRep-1])]
		nInds=self.numIndividualsPerReplicate[indexRep-1]
		indIds=np.transpose(np.array(range(0,nInds)))
		f=open(filepath,"w")
		f.write("{}\n".format(" ,".join(header)))
		for r in range(0,nInds):
			weird=coverageMatrix[r,].tolist()[0]
			weird=["{0:.4f}".format(float(item)) for item in weird]
			line=", ".join(weird)
			line=line.replace("[","").replace("]","")
			# print(line)
			f.write("{0} ,{1}\n".format(r,line))
		f.close()
