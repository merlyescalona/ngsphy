import csv,datetime,logging,os,subprocess,sys,threading,time
import numpy as np
import random as rnd
from scipy.stats import  binom,expon,gamma,lognorm,norm,nbinom,poisson,uniform

class NGSAvailableDistributions:
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
		self.__params=params

	def getParams(self):
		return self.__params

	def printParams(self):
		print(self.__params)

	def setName(self,name):
		self.__name=name

	def getName(self):
		return self.__name

	def printName(self):
		print(self.__name)

	def asNGSPhyDistribution(self):
		if not self.__dependency:
			return NGSPhyDistribution(self.__name,self.__params)
		else: return None

	def validate(self):
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
	__name=None
	__params=[]
	def __init__(self,name,params):
		self.__name=name
		self.__params=params
		# print(self.__name)
		# print(self.__params)

	def value(self,samples=1):
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
			message="{0}\n{1}\t{2}\t{3}".format(\
				ex,exc_type,\
				fname, exc_tb.tb_lineno)
			raise Exception(message)
			# self.appLogger.error(message)
			# sys.exit()

		return value

	def binom(self,samples):
		#n=number of times tested
		#p=probability
		n=self.__params[0]
		p=self.__params[1]
		distro=binom(n,p)
		f=distro.rvs(size=samples)
		return f

	def exponential(self,samples):
		f=np.random.exponential(float(self.__params[0]),size=samples)
		return value

	def fixed(self,samples):
		value=[self.__params[0]]*samples
		return value

	def gamma1(self,samples):
		# E|x| = k.theta (alpha*theta)
		# If i want  mean=1, theta=E|x|/alpha=1/alpha
		shape=float(self.__params[0]*1.0)
		theta=float(1/(shape*1.0))
		distro=gamma(shape,theta)
		f=distro.rvs(size=samples)
		return f


	def gamma(self,samples):
		# The parameterization with alpha and beta is more common in Bayesian statistics,
		# where the gamma distribution is used as a conjugate prior distribution
		# for various types of inverse scale (aka rate) parameters, such as the
		#lambda of an exponential distribution or a Poisson distribution[4]  or for
		#t hat matter, the beta of the gamma distribution itself.
		#(The closely related inverse gamma distribution is used as a conjugate
		# prior for scale parameters, such as the variance of a normal distribution.)
		# shape, scale = 2., 2. # mean=4, std=2*sqrt(2)
		# s = np.random.gamma(shape, scale, 1000)
		# E|x| = k.theta (alpha*theta)
		# If i want a specific mean, theta=E|x|/alpha
		shape=float(self.__params[0]*1.0)
		theta=float(self.__params[1]*1.0)
		distro=gamma(shape,theta)
		f=distro.rvs(size=samples)
		return f


	def lognormal(self,samples):
		# m=mean, s=sigma - standard deviation
		mean=float(self.__params[0]*1.0)
		sigma=float(self.__params[1]*1.0)
		distro=lognorm(mean,sigma)
		f=distro.rvs(size=samples)
		return f

	def nbinom(self,samples):
		# mu= poissonon mean
		# r controls the deviation from the poisson
		# This makes the negative binomial distribution suitable as a robust alternative to the Poisson,
		# which approaches the Poisson for large r, but which has larger variance than the Poisson for small r.
		mu=float(self.__params[0])
		r=float(self.__params[1])
		p=(r*1.0)/(r+mu)
		distro=nbinom(r,p)
		f=distro.rvs(size=samples)
		return f

	def normal(self,samples):
		# mean (location) - variance (squared scale)
		locMean=float(self.__params[0]*1.0)
		scaleVariance=float(self.__params[1]*1.0)
		distro=norm(loc=locMean,scale=scaleVariance)
		f=distro.rvs(size=samples)
		return f

	def poisson(self,samples):
		l=float(self.__params[0]*1.0)
		distro=poisson(l)
		f=distro.rvs(size=samples)
		return f

	def uniform(self,samples):
		# mean
		mean=float(self.__params[0]*1.0)
		distro=uniform(mean)
		f=distro.rvs(size=samples)
		return f

class CoverageMatrixGenerator:
	OFFTARGET_COVERAGE_FRACTION=0.1
	appLogger=None
	settings=None

	alphashapesLocus=[]
	alphashapesIndividuals=[]
	locusMultiplier=[]
	individualsMultiplier=[]

	experiment=None

	outputFolderPath=""

	numSpeciesTrees=0
	numLociPerSpeciesTree=[]
	numIndividualsPerSpeciesTree=[]
	numSpeciesTreesDigits=0
	numLociPerSpeciesTreeDigits=[]
	numIndividualsPerSpeciesTree=[]
	filteredST=[]

	def __init__(self, settings):
		self.appLogger=logging.getLogger('ngsphy')
		self.appLogger.info("Coverage calculations")
		self.settings=settings
		self.outputFolderPath=os.path.join(\
			self.settings.outputFolderPath,\
			"coverage"
		)

		self.experiment=self.settings.experiment.asNGSPhyDistribution()
		self.numSpeciesTrees=self.settings.parser.getint("general", "numspeciestrees")
		cc=self.settings.parser.get("general", "numLociPerSpeciesTree").strip().split(",")
		self.numLociPerSpeciesTree=[ int(item) for item in cc if not item ==""]
		cc=self.settings.parser.get("general", "numIndividualsPerSpeciesTree").strip().split(",")
		self.numIndividualsPerSpeciesTree=[ int(item) for item in cc if not item == ""]
		self.numSpeciesTreesDigits=len(str(self.numSpeciesTrees))
		self.numLociPerSpeciesTreeDigits=[len(str(item)) for item in self.numLociPerSpeciesTree]
		self.numIndividualsPerSpeciesTreeDigits=[len(str(item )) for item in self.numIndividualsPerSpeciesTree]
		self.filteredST=[int(item) for item in self.settings.parser.get("general", "filtered_ST").strip().split(",")]
		if (self.settings.locus):
			distro=self.settings.locus.asNGSPhyDistribution()
			self.alphashapesLocus=distro.value(self.settings.numSpeciesTrees)
			self.locusMultiplier=[NGSPhyDistribution("g1",[self.alphashapesLocus[index]]).value(self.numLociPerSpeciesTree[index]) for index in range(0,self.numSpeciesTrees)]
		if (self.settings.individual):
			distro=self.settings.individual.asNGSPhyDistribution()
			self.alphashapesIndividuals=distro.value(self.settings.numSpeciesTrees)
			self.individualsMultiplier=[NGSPhyDistribution("g1",[self.alphashapesIndividuals[index]]).value(self.numLociPerSpeciesTree[index]) for item in range(0,self.numSpeciesTrees)]
		self.generateFolderStructure()

	def generateFolderStructure(self):
		self.appLogger.debug("Coverage folder structure generation...")
		try:
			os.makedirs(self.outputFolderPath)
			self.appLogger.info("Generating output folder ({0})".format(self.outputFolderPath))
		except:
			self.appLogger.debug("Output folder exists ({0})".format(self.outputFolderPath))

	def calculate(self):
		message=""
		status=True;
		self.appLogger.debug("Coverage calculations...")
		for indexST in self.filteredST:
			nInds=self.numIndividualsPerSpeciesTree[indexST-1]
			nLoci=self.numLociPerSpeciesTree[indexST-1]
			# initialCoverage=self.experiment.value(nInds*nLoci)
			val=self.experiment.value(1)
			# self.appLogger.debug(val)
			initialCoverage=val*(nInds*nLoci)
			coverageMatrix=np.matrix(initialCoverage)
			coverageMatrix.shape=[nInds,nLoci]
			# individuals + loci coverage variation
			# individuals + loci multipliers
			try:
				# --------------------------------------------------------------
				# For genomic stochasticity
				# --------------------------------------------------------------
				if self.settings.genomicNoise:
					params=[]
					for loc in range(0,nLoci):
						for ind in range(0,nInds):
							params=[coverageMatrix[ind,loc]]+self.locus.getParams()
							coverageMatrix[ind,loc]=NGSPhyDistribution(self.locus.getName(),params).value(1)
				if self.settings.locus:
					multipliers=self.locusMultiplier[indexST-1]
					for loc in range(0,nLoci):
						coverageMatrix[:, loc]=coverageMatrix[:,loc]*multipliers[loc]
				if self.settings.individual:
					multipliers=self.individualsMultiplier[indexST-1]
					for ind in range(0,nInds):
						coverageMatrix[ind,]=coverageMatrix[ind,]*multipliers[ind]


			except Exception as ex:
				exc_type, exc_obj, exc_tb = sys.exc_info()
				fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
				message="Coverage calculation error.\n\t{0}\n\t{1}\t{2}\t{3}".format(\
					ex,exc_type,\
					fname, exc_tb.tb_lineno)
				return False, message
			# phylogenetic decay
			# I need the file with the realtion of the individuals
			# dont car if haploid or diploid because ill only need the first 3 columns
			# which match in both cases
			if self.settings.phylogeneticDecay:
				individualsTableFilename=os.path.join(\
					self.settings.outputFolderPath,\
					"tables",
					"{0}.{1:0{2}d}.individuals.csv".format(\
						self.settings.projectName,\
						indexST,\
						self.numSpeciesTreesDigits)
					)
				individualsTableFile=open(individualsTableFilename)
				d = csv.DictReader(individualsTableFile)
				individualsTable = [row for row in d]
				individualsTableFile.close()
				for row in individualsTable:
					if str(row['speciesID']) in self.settings.phylogeneticDecay.keys():
						coverageMatrix[ind,]=np.array(coverageMatrix[ind,])*self.settings.phylogeneticDecay[str(row['speciesID'])]
			#on/off target
			onTargetLoci=range(0,nLoci)
			offTargetLoci=[]
			try:
				if self.settings.offtarget > 0:
					nSamples=int(self.settings.offtarget*nLoci)
					offTargetLoci=np.random.choice(nLoci,nSamples,replace=F)
					onTargetLoci=set(range(0,nLoci))-set(offTargetLoci)
					for loc in offTargetLoci:
						column=np.array(np.array(coverageMatrix[:,loc])*self.OFFTARGET_COVERAGE_FRACTION)
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
			self.write(coverageMatrix,indexST)
		return status,message

	def write(self,coverageMatrix,indexST):
		self.appLogger.debug("Writing coverage matrix for ST replicate: {0:0{1}d}".format(indexST, self.numSpeciesTreesDigits))
		filename=os.path.join(self.outputFolderPath,"{0}.{1:0{2}d}.coverage.csv".format(\
			self.settings.projectName,\
			indexST,\
			self.numSpeciesTreesDigits\
			)
		)
		filepath=os.path.abspath(filename)
		header=["indID"]+[ "L.{0:0{1}d}".format(loc+1,self.numLociPerSpeciesTreeDigits[indexST-1]) for loc in range(0,self.numLociPerSpeciesTree[indexST-1])]
		nInds=self.numIndividualsPerSpeciesTree[indexST-1]
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
