#!/usr/bin/env python
# -*- coding: utf-8 -*-
import argparse,copy,csv,datetime,logging,os,subprocess,sys,threading,time,warnings
import numpy as np
import random as rnd
import Settings as sp
from MELoggingFormatter import MELoggingFormatter as mlf
from MSATools import *
from NGSPhyDistribution import NGSPhyDistribution as ngsphydistro

try:
    from collections import Counter
except ImportError:
    from counter import Counter

def getScoreSingle(data):
    self.appLogger.debug("getScoreSingle(data)")
    if data!=0:
        return -10*np.log10(data)
    else:
        return 0
def getScoreMatrix(data):
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

class ReadCount:
    __NUCLEOTIDES=["A","C","G","T"]
    # path related variables
    readcountFolder=""
    referencesFolder=""
    coverageFolder=""

    # Coverage Distribution variables
    experimentCoverageDistro=None
    individualCoverageDistro=None
    locusCoverageDistro=None

    seqerror=0
    refereceFilepath=None
    dataprefix=""

    output=""
    numLociPerST=[]
    numSpeciesTrees=0
    numSpeciesTreesDigits=0
    filteredST=[]
    numSts=0

    def __init__(self,settings):
        self.appLogger=logging.getLogger('ngsphy')
        self.appLogger.info('Read counts.')
        self.settings=settings
        self.path=os.path.abspath(self.settings.parser.get("general", "simphy_folder"))
        if (self.settings.parser.get("general", "simphy_folder")[-1]=="/"):
            self.projectName=os.path.basename(self.settings.parser.get("general", "simphy_folder")[0:-1])
        else:
            self.projectName=os.path.basename(self.settings.parser.get("general", "simphy_folder"))

        # This has been previously checked at Settings Class!
        self.experimentCoverageDistro=ngsphydistro(0,\
            self.settings.parser.get("coverage","experimentCoverage"))
        if (self.settings.parser.has_option("coverage","individualCoverage")):
            self.individualCoverageDistro=ngsphydistro(1,\
                self.settings.parser.get("coverage","individualCoverage"),\
                True)
        if (self.settings.parser.has_option("coverage","locusCoverage")):
            self.locusCoverageDistro=ngsphydistro(2,\
                self.settings.parser.get("coverage","locusCoverage"),\
                True)

        self.output=self.settings.parser.get("general","output_folder_name")
        self.numLociPerST=[int(numST) for numST in self.settings.parser.get("general","numLociPerST").split(",")]
        self.numSpeciesTrees=int(self.settings.parser.get("general","numSpeciesTrees"))
        self.numSpeciesTreesDigits=len(str(self.numSpeciesTrees))
        self.dataprefix=self.settings.parser.get("general","data_prefix")
        self.filteredST=self.settings.parser.get("general", "filtered_ST")
        self.filteredST=[ int(numST) for numST in self.filteredST.split(",")]
        self.numSTs=self.settings.parser.getint("general","number_ST")

        self.refereceFilepath=self.settings.parser.get("read-count","reference")
        if self.refereceFilepath=="None":
            self.refereceFilepath=None
        else:
            self.refereceFilepath=os.path.abspath(self.refereceFilepath)

        try:
            self.seqerror=float(self.settings.parser.get("read-count","error"))
        except:
            self.seqerror=0
        # Generating folder structur - basic for
        self.generateFolderStructureGeneral()

    """
    @function:
        generateFolderStructureGeneral(self)
    @description:
        Generates basic folder structure needed for Read Count option.
        This includes:
            - read-count
            - read-count/true
            - read-count/sampled
            - refereces
            - coverage
    """
    def generateFolderStructureGeneral(self):
        self.appLogger.debug("generateFolderStructureGeneral(self)")
        # Checking output path
        self.appLogger.info("Creating folder structure for [read-count]")
        self.readcountFolder="{0}/read-count".format(self.output)
        self.readcountTrueFolder="{0}/read-count/true".format(self.output)
        self.readcountSampledFolder="{0}/read-count/sampled".format(self.output)
        self.referencesFolder="{0}/refereces".format(self.output)
        self.coverageFolder="{0}/coverage".format(self.output)
        try:
            os.makedirs(self.readcountFolder)
            self.appLogger.info("Generating output folder ({0})".format(self.readcountFolder))
        except:
            self.appLogger.debug("Output folder exists ({0})".format(self.readcountFolder))
        try:
            os.makedirs(self.readcountTrueFolder)
            self.appLogger.info("Generating output folder ({0})".format(self.readcountTrueFolder))
        except:
            self.appLogger.debug("Output folder exists ({0})".format(self.readcountTrueFolder))
        try:
            os.makedirs(self.readcountSampledFolder)
            self.appLogger.info("Generating output folder ({0})".format(self.readcountSampledFolder))
        except:
            self.appLogger.debug("Output folder exists ({0})".format(self.readcountSampledFolder))
        try:
            os.makedirs(self.referencesFolder)
            self.appLogger.info("Generating output folder ({0})".format(self.referencesFolder))
        except:
            self.appLogger.debug("Output folder exists ({0})".format(self.referencesFolder))

        try:
            os.makedirs(self.coverageFolder)
            self.appLogger.info("Generating output folder ({0})".format(self.coverageFolder))
        except:
            self.appLogger.debug("Output folder exists ({0})".format(self.coverageFolder))

    """
    @function:
        generateFolderStructureDetail(self)
    @description:
        Generates detailed folder structure needed for Read Count option.
        This includes:
            - Replicate folders  (as many as species tree replicates)
            - true/sampled folder per replicate.
                - These folders will contain the VCF files (true and sampled)
    """
    def generateFolderStructureDetail(self):
        self.appLogger.debug("generateFolderStructureDetail(self)")
        for i in range(0,len(self.filteredST)):
            indexST=self.filteredST[i]
            trueFolder="{0}/{1:0{2}d}".format(\
                self.readcountTrueFolder,\
                indexST,\
                self.numSpeciesTreesDigits)
            sampledFolder="{0}/{1:0{2}d}".format(\
                self.readcountSampledFolder,\
                indexST,\
                self.numSpeciesTreesDigits)
            references="{0}/{1:0{2}d}".format(\
                self.referencesFolder,\
                indexST,\
                self.numSpeciesTreesDigits
            )
            try:
                os.makedirs(trueFolder)
                self.appLogger.info("Generating output folder ({0})".format(trueFolder))
            except:
                self.appLogger.debug("Output folder exists ({0})".format(trueFolder))
            try:
                os.makedirs(sampledFolder)
                self.appLogger.info("Generating output folder ({0})".format(sampledFolder))
            except:
                self.appLogger.debug("Output folder exists ({0})".format(sampledFolder))
            try:
                os.makedirs(references)
                self.appLogger.info("Generating output folder ({0})".format(references))
            except:
                self.appLogger.debug("Output folder exists ({0})".format(references))

    """
    @function:
        computeCoverageMatrix(self, nInds, nLoci,expCov, indCov, locCov):
    @description:
        Computes coverage matrix for a whole replicate
    @intput:
        nInds: number of individuals
        nLoci: number of loci
        expCov: NGSPhyDistribution object used for experiment-wide coverage
        indCov: NGSPhyDistribution object used for individual-specific coverage
        locCov: NGSPhyDistribution object used for locus-specific coverage
    @output:
        covMatrix: nInds (rows) x nLoci (columns)
    """
    def computeCoverageMatrix(self, nInds, nLoci,expCov, indCov, locCov):
        self.appLogger.debug("computeCoverageMatrix(self, nInds, nLoci,expCov, indCov, locCov)")
        self.appLogger.debug("Computing coverage matrix")

        # coverage matrix per ST - row:indv - col:loci
        # each cov introduced as parameter is a NGSPhyDistribution
        covMatrix=np.zeros(shape=(nInds, nLoci))
        for loc in range(0,nLoci):
            for ind in range(0,nInds):
                expc=0;indc=0;locc=0
                expc=expCov.value()
                coverage=expc
                if (indCov):
                    indCov.updateValue(expc)
                    indc=indCov.value()
                    coverage=indc
                if (locCov):
                    locCov.updateValue(indc)
                    locc=locCov.value()
                    coverage=locc
                covMatrix[ind][loc]=coverage
        return covMatrix

    """
    @function:
        writeCoverageMatrixIntoFile(self, coverageMatrix, filename):
    @description:
        Write coverageMatrix into filename
    @intput:
        coverageMatrix: computed by - computeCoverageMatrix()
        filename: path of the file where the coverageMatrix will be stored
    @result:
        file will be generated in the given apth
    """
    def writeCoverageMatrixIntoFile(self, coverageMatrix, filename):
        self.appLogger.debug("writeCoverageMatrixIntoFile(self, coverageMatrix, filename)")
        filepath=os.path.abspath(filename)
        with open(filepath, 'w') as csvfile:
            writer = csv.writer(csvfile)
            [writer.writerow(r) for r in coverageMatrix]

    """
    @function:
        parseReferenceList(self, filename)
    @description:
        Used to parse referenceList, file with format: STID,SPID,TIPID
    @intput:
        filename: path of the referenceList file.
        There's only ONE file with the relation of the references
        If "None" inputted (file is missing) then reference by default is 1_0_0
        for all species tree replicates.
    @output:
        output: list. each element of the list is a triplet (STID,SPID,TIPID)
    """
    def parseReferenceList(self, filename):
        self.appLogger.debug("parseReferenceList(self, filename)")
        referenceList=[]
        if filename:
            # There's a file
            filepath=os.path.abspath(filename)
            f=open(filepath,"r")
            lines=f.readline()
            lines=f.readlines()
            f.close()
            for index in range(0, len(lines)):
                d=lines[index].strip().split(",")
                try:
                    referenceList+=[(d[0],d[1],d[2])]
                except IndexError as ie:
                    self.appLogger.warning("Parsing reference list. "+\
                        "A default reference has been introduced.\n"+\
                        "Replicate index: {0}".format(\
                        index
                    ))
                    referenceList+=[(index+1,1,0)]
        else:
            # iterate get a list with same decription for all species trees
            for item in range(0,self.numSpeciesTrees):
                referenceList+=[(item+1,1,0)]
        return referenceList

    """
    @function:
        extractTrueVariantsPositions(self,filename)
    @description:
        Extract true variant positions from the MSA file used as input
    @intput:
        filename: MSA fasta file from where to extract the variable positions
    @output:
        dictionary. indices=variable sites, content=variable nucleotide set
    """
    def extractTrueVariantsPositions(self, filename):
        self.appLogger.debug("extractTrueVariantsPositions(self, filename)")
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
            numTotalSeqs=len(lines)/2
            matrix=np.chararray((numTotalSeqs,lenSeqs), itemsize=1)
        else:
            raise TypeError("File has a wrong file format.\n{}\nPlease check. Exiting.".format(filepath))
            # If I get here, function run is over - goes back to  run()

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

    """
    @function:
        parseIndividualRelationFile(self,filename)
    @description:
        Parses "Individual-description relation" file per species tree
        in order to obtain information on how sequences from INDELible
        are related. Generated within NGSPhy.
    @intput:
        filename: path for the Individual-description relation file
    @output:
        if ploidy=1:    returns a dict. key: indID,content: dict(indexST,seqDEscription)
        if ploidy=2:    returns a dict. key: indID,content: dict(indexST,speciesID, mateID1,mateID2)
    """
    def parseIndividualRelationFile(self,filename):
        self.appLogger.debug("parseIndividualRelationFile(self,filename)")
        individuals=dict()
        if (self.settings.ploidy>0 and self.settings.ploidy<=2):
            csvfile=open(os.path.abspath(filename))
            d = csv.DictReader(csvfile)
            if (self.settings.ploidy==1):
                for row in d:
                    individuals[str(row["indID"])] = {\
                        "indexST":row["indexST"],\
                        "seqDescription":row["seqDescription"]}
            if (self.settings.ploidy==2):
                # indexST,indID,speciesID,mateID1,mateID2
                for row in d:
                    individuals[str(row["indID"])] = {\
                        "indexST":row["indexST"],\
                        "speciesID":row["speciesID"],\
                        "mateID1":row["mateID1"],\
                        "mateID2":row["mateID2"]}
            csvfile.close()
        else:
            # There has been a verification in Settings class, but just in case.
            raise ValueError("Ploidy assigned is incorrect. Please verify. Exciting.")
        return individuals

    """
    @function:
        getDepthCoveragePerIndividual(self, numVarSites,startingCoverage):
    @description:
        Compute coverage per individuals
        since the value for SNPS is fixed to a Negative Binomial distribution
        I iterate over the number of variants for this specific locus and
        since everythime i hit play! (called value()) I'll sample a new value
        from the previously set distribution, i'll then get a random value from a
        negative bionmial distribution with mean max-coverage
    @intput:
        numVariableSites: number of variable sites
        startingCoverage: starting value to calculate coverage
    @output:
        list of values that correspond to the coverage per each snp of the MSA
    """
    def getDepthCoveragePerIndividual(self, numVarSites,startingCoverage):
        self.appLogger.debug(\
            "getDepthCoveragePerIndividual(self, numVarSites,startingCoverage) - ({},{})".format(\
                numVarSites,\
                startingCoverage
            ))
        distro=ngsphydistro(0,"nb:{0},{1}".format(startingCoverage,startingCoverage))
        DP=[ distro.value()  for index in range(0,numVarSites) ]
        return DP


    """
    @function:
        getHaploidIndividualSequence(self,msa,ind)
    @description:
        Extract individuals "ind" sequence(s) from MSA dictionary
    @intput:
        msa: dictionary
        ind: dict(indexST,seqDEscription)
    @output:
        sequence of the individual
    """
    def getHaploidIndividualSequence(self,msa,ind):
        # ind: [indID,indexST,seqDescription]
        self.appLogger.debug("getHaploidIndividualSequence(self,msa,ind)")
        seqSize=len(msa[str(1)][str(0)]['sequence'])
        fullInd=None; speciesID=None; tipID=None; tmp=None
        fullInd=np.chararray(shape=(1,seqSize), itemsize=1)
        seqDescription=ind["seqDescription"].strip().split("_")
        speciesID=seqDescription[0]
        tipID=seqDescription[2]
        tmp=list(msa[str(speciesID)][str(tipID)]['sequence'])
        fullInd=[item for item in tmp]
        return fullInd

    """
    @function:
        getDiploidIndividualSequence(self,msa,ind)
    @description:
        Extract individuals "ind" sequence(s) from MSA dictionary
    @intput:
        msa: dictionary
        ind: dict(indexST,speciesID, mateID1,mateID2)
    @output:
        matrix: 2 x seqLength - representing the sequence of the individual
    """
    def getDiploidIndividualSequence(self,msa,ind):
        # ind:[indID, indexST,speciesID, mateID1,mateID2]
        self.appLogger.debug("getDiploidIndividualSequence(self,msa,ind)")
        seqSize=len(msa[str(1)][str(0)]['sequence'])
        fullInd=None;  tmp=None
        speciesID=ind["speciesID"];
        tipID1=ind["mateID1"];
        tipID2=ind["mateID2"];
        fullInd=np.chararray(shape=(2,seqSize), itemsize=1)
        tmp=list(msa[str(speciesID)][str(tipID1)]['sequence'])
        fullInd[0,:]=[item for item in tmp]
        tmp=list(msa[str(speciesID)][str(tipID2)]['sequence'])
        fullInd[1,:]=[item for item in tmp]
        return fullInd


    """
    @function:
        computeHaploid(indexST,indexGT,msa,individuals,referenceSeqFull,variableSites,DP)
    @description:
        compute the READ COUNT for the specic ploidy
    @intput:
        indexST: species tree id
        indexGT: locus/gene tree id
        msa: parsed multiple sequence alignent file (dictionary)
        individuals: parsed individual-description relation file (dictionary)
        referenceFilepath: where the reference is written.
        referenceSeqFull: reference sequence as is
        variableSites: dictionary. keys: variable positions. content: variable nucleotide set
        DP: variable sites coverage
    @result:
        get a file written in the STreplicate folder (true/sampled)
    """
    def computeHaploid(self,indexST,indexGT,msa,individuals,referenceFilepath,referenceSeqFull,variableSites,DP):
        self.appLogger.debug("computeHaploid(msa,individuals,referenceFilepath,referenceSeqFull,variableSites,DP)")
        nInds=len(individuals)
        nVarSites=len(variableSites.keys())
        # Getting the proper indices of the variable sites
        variableSitesPositionIndices=np.sort([int(pos) for pos in variableSites.keys()])
        alt=dict();  altInds=dict()
        rcsTrue=dict();rcsSampled=dict()
        HTGeneral=dict();HLGeneral=dict();ADGeneral=dict()
        HTGeneralSampled=dict();HLGeneralSampled=dict();ADGeneralSampled=dict()
        for pos in variableSites.keys():
            alt[pos]=[]; altInds[pos]=[]
        for pos in variableSitesPositionIndices:
            alt[str(pos)]=list(set(variableSites[str(pos)])-set([referenceSeqFull[pos]]))
        for index in range(0,nInds):
            indexIND=str(index)
            HTGeneral[indexIND]=[];HLGeneral[indexIND]=[];ADGeneral[indexIND]=[]
            HTGeneralSampled[indexIND]=[];
            HLGeneralSampled[indexIND]=[];ADGeneralSampled[indexIND]=[]
            rcsTrue[indexIND]=[];rcsSampled[indexIND]=[]
        ########################################################
        # TRUE
        ########################################################
        # iterate over individuals
        self.appLogger.debug("True")
        altInd=copy.copy(alt)
        for indIndex in range(0,nInds):
            indexIND=str(indIndex)
            self.appLogger.debug("Individual {0} ({1}/{2})".format(indIndex,(indIndex+1),nInds))
            ind=individuals[indexIND]
            individualSeq=self.getHaploidIndividualSequence(msa,ind)
            # individualSeq[202]+="**"
            # AD per individual different if true/sampledInd
            ADTrue,ADSampled=self.getAllelicDepthPerHaploidIndividual(\
                individualSeq,variableSitesPositionIndices,DP[indIndex])
            rcTrue,rcSampled=self.getReadCountPerIndividual(\
                ADTrue,ADSampled,variableSitesPositionIndices)
            altInd=self.getAltUpdatedPerIndividual(\
                referenceSeqFull,altInd,ADSampled)
            HTTrue=self.gettingHaplotype(\
                referenceSeqFull,individualSeq,alt,variableSitesPositionIndices)
            HLTrue=self.haplotypeLikehood(\
                rcTrue,variableSitesPositionIndices,0)
            rcsTrue[indexIND]=rcTrue
            rcsSampled[indexIND]=rcSampled
            ADGeneral[indexIND]=ADTrue
            ADGeneralSampled[indexIND]=ADSampled
            HTGeneral[indexIND]=HTTrue
            HLGeneral[indexIND]=HLTrue

        # here corresponds to the INDEXST that is going to be written
        self.writeVCFFile(\
            indexST,indexGT,\
            referenceFilepath, referenceSeqFull,\
            alt,variableSitesPositionIndices,\
            HTGeneral,HLGeneral,\
            ADGeneral,DP,True)

        ########################################################
        # Sampled
        ########################################################
        self.appLogger.debug("Sampled")
        for indIndex in range(0,nInds):
            indexIND=str(indIndex)
            self.appLogger.debug("Individual {0} ({1}/{2})".format(indIndex,(indIndex+1),nInds))
            ind=individuals[indexIND]
            individualSeq=self.getHaploidIndividualSequence(msa,ind)
            HTSampled=self.gettingHaplotype(\
                referenceSeqFull,individualSeq,\
                altInd, variableSitesPositionIndices)
            HLSampled=self.haplotypeLikehood(\
                rcsSampled[indexIND],\
                variableSitesPositionIndices,\
                self.seqerror)
            HTGeneralSampled[indexIND]=HTSampled
            HLGeneralSampled[indexIND]=HLSampled
        self.writeVCFFile(
            indexST,indexGT,\
            referenceFilepath, referenceSeqFull,\
            altInd,variableSitesPositionIndices,\
            HTGeneralSampled,HLGeneralSampled,\
            ADGeneralSampled,DP,False)
        sys.exit()
        pass

    """
    @function:
        haplotypeLikehood(self,variantsRC,observed,variableSites,error):
    @description:
        computes log10-scaled haplotype likelihood per variable site
    @intput:
        readcount:
        variableSitesPositionIndices:
        error:
    @output:
        GL: dictionary.
            keys: positons.
            content: log10-scaled likelihood for the specific variable site
    """
    def haplotypeLikehood(self,variantsRC,variableSitesPositionIndices,error):
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
            value[infs]=0
        return value

    """
    @function:
        gettingHaplotype(self,ref,seq,alt, variableSitesPositionIndices):
    @description:
        get haplotype for the given ref-seq-alt triplet
    @intput:
        ref: reference sequence
        seq: individual sequnece
        alt: dictionary with alternative alleles per variable position
        variableSitesPositionIndices: variableSites
    @output:
        GT: dictionary. keys: positons. content: haplotype for the specific variable site
    """
    def gettingHaplotype(self,ref,seq,alt, variableSitesPositionIndices):
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
                if refPos==seqPos:
                    GT[iv]=0
                else:
                    if (index==0) and (altToCompare==seqPos):
                        GT[iv]=1
                    if (index==1) and (altToCompare==seqPos):
                        GT[iv]=2
                    if (index==2) and (altToCompare==seqPos):
                        GT[iv]=3
        return GT

    """
    @function:
        getAltUpdatedPerIndividual(self,ref,alt,AD)
    @description:
        update general alternative alleles list
    @intput:
        referenceSeq: reference sequence
        alt: dictionary. current alternative allele list per variant
        AD: allelic depth matrix
    @output:
        newAlt: dictionary. keys: positions. content: alt alleles corresponding to that position
    """
    def getAltUpdatedPerIndividual(self,ref,alt,AD):
        self.appLogger.debug("getAltUpdatedPerIndividual(self,ref,alt,AD)")
        # update alt values
        # Sampled - ALT is a DICT
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

    """
    @function:
        getReadCountPerIndividual(self,ADTrue,ADSampled, variableSitesPositionIndices)
    @description:
        generates list of read counts per variant per individual
    @intput:
        ADTrue: True allelic depth
        ADSampled: Sampled allelic depth
        variableSitesPositionIndices: position of the variable sites
    @output:
        RCTrue, RCSampled.
    """
    def getReadCountPerIndividual(self,ADTrue,ADSampled, variableSitesPositionIndices):
        self.appLogger.debug("getReadCountPerIndividual(self,ADTrue,ADSampled, variableSitesPositionIndices)")
        nVariants=len(variableSitesPositionIndices)
        rcTrue=dict();rcSampled=dict()
        # init variantsRC
        for indexVar in range(0,nVariants):
            indexConverted=str(variableSitesPositionIndices[indexVar])
            rcTrue[indexConverted]=[]
            rcSampled[indexConverted]=[]
        # getting read count
        for indexVar in range(0,nVariants):
            for indexNuc in range(0,4):
                with warnings.catch_warnings():
                    warnings.filterwarnings("ignore", category=np.VisibleDeprecationWarning)
                    nuc=[self.__NUCLEOTIDES[indexNuc]]
                    indexConverted=str(variableSitesPositionIndices[indexVar])
                    val=ADTrue[indexNuc,indexVar]
                    rcTrue[indexConverted]+=nuc*val
                    val=ADSampled[indexNuc,indexVar]
                    rcSampled[indexConverted]+=nuc*val
        return rcTrue, rcSampled

    """
    @function:
        getAllelicDepthPerHaploidIndividual(self,fullInd,variableSites,DP)
    @description:
        Generate allelic depth per individual per site
    @intput:
        individualSequence: sequence of the individual
        variableSitesPositionIndices: position of the variable sites
        DP: coverage of the variable sites
    @output:
        ADTrue, ADSampled.
    """
    def getAllelicDepthPerHaploidIndividual(self,individualSeq,variableSitesPositionIndices,DP):
        self.appLogger.debug("getAllelicDepthPerHaploidIndividual(self,fullInd,variableSites,DP)")
        nVariants=len(variableSitesPositionIndices)
        ADTrue=np.zeros(shape=(4,nVariants), dtype=np.int)
        ADSampled=np.zeros(shape=(4,nVariants), dtype=np.int)
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
            ADTrue[0,indexVar]=counter["A"];ADTrue[1,indexVar]=counter["C"]
            ADTrue[2,indexVar]=counter["G"];ADTrue[3,indexVar]=counter["T"]
            # SAMPLE READ COUNT - need to know error distribution
            errorDistro=ngsphydistro(0,"b:{0},{1}".format(posCoverage,self.seqerror))
            errorD=errorDistro.value()
            errorPositions=[]
            # need to know possible nucleotides to substitute my position with error
            possibleNucs=list(set(self.__NUCLEOTIDES)-set([indNucs]))
            # I have some positions (at least 1) that is an error
            # errorD= array with coded error nucleotides that will be modified
            if (errorD>0):
                errorChoices=np.random.choice(possibleNucs, int(errorD), replace=True)
                maxAvailablePositions=posCoverage
                if not ((posCoverage % 2) == 0): maxAvailablePositions=posCoverage-1
                if maxAvailablePositions == 0:  maxAvailablePositions=posCoverage

                errorPositions=np.random.randint(maxAvailablePositions,size=int(errorD))
                for item in range(0,len(errorPositions)):
                    posToChange=errorPositions[item]
                    finalRC[posToChange]=errorChoices[item]
            counter=Counter(finalRC)
            # if any of the nucleotides does not have a counter retrieving it will be 0
            ADSampled[0,indexVar]=counter["A"];ADSampled[1,indexVar]=counter["C"]
            ADSampled[2,indexVar]=counter["G"];ADSampled[3,indexVar]=counter["T"]
        return ADTrue,ADSampled

    """
    @function:
        computeDiploid(indexST,indexGT,msa,individuals,referenceSeqFull,variableSites,DP)
    @description:
        compute the READ COUNT for the specic ploidy
    @intput:
        indexST: species tree id
        indexGT: locus/gene tree id
        msa: parsed multiple sequence alignent file (dictionary)
        individuals: parsed individual-description relation file (dictionary)
        referenceSeqFull: reference sequence as is
        variableSites: dictionary. keys: variable positions. content: variable nucleotide set
        DP: variable sites coverage
    @result:
        get a file written in the STreplicate folder (true/sampled)
    """
    def computeDiploid(self,indexST,indexGT,msa,individuals,referenceSeqFull,variableSites,DP):
        self.appLogger.debug("computeDiploid(self,indexST,indexGT,msa,individuals,referenceSeqFull,variableSites,DP)")
        pass


    """
    @function:
        formatIndividualDataForVCF(self,ref,alt,variableSites,HT,HL,AD,DP):
    @description:
        gets single individual data formated as a GT:GL:AD:DP VCF column
    @intput:
        reference sequence
        alternative alleles
        variableSitesPositionIndices
        haplotype
        likelihood
        allelic depth
        read depth
    @output:
        string column with information in format: GT:GL:AD:DP
    """
    def formatIndividualDataForVCF(self,ref,alt,variableSitesPositionIndices,HT,HL,AD,DP):
        self.appLogger.debug("formatIndividualDataForVCF(self,ref,alt,variableSitesPositionIndices,HT,HL,AD,DP)")
        nVariants=len(variableSitesPositionIndices)
        nInds=len(HT.keys())
        allVariants=dict()
        for indexVAR in variableSitesPositionIndices:
            allVariants[str(indexVAR)]=[]
        # print "Init allvariants dict"
        # print len(HT.keys()),len(indices)
        # Had an error here, because the alt variable was empty (passed wrong variable name to the function)
        for indexVAR in range(0,nVariants):
            tmpInd=variableSitesPositionIndices[indexVAR]
            for indexIND in range(0,nInds):
                # HT is a dict, HT[indexIND] -> is HT for individual indexIND
                htPerInd=HT[str(indexIND)]
                trueRows=None
                # trueRows is the ref nucleotides
                valuesToCodify=[ref[tmpInd]]+alt[str(tmpInd)]
                trueRows=self.codifySequences(valuesToCodify)
                # print "Before truerows ",DP[indexIND][indexVAR]
                ind="{0}:{1}:{2}:{3}".format(\
                    htPerInd[str(tmpInd)],\
                    ",".join(HL[str(indexIND)][trueRows,indexVAR].astype(dtype=int).astype(dtype=np.str)),\
                    ",".join(AD[str(indexIND)][trueRows,indexVAR].astype(dtype=int).astype(dtype=np.str)),\
                    DP[indexIND][indexVAR])
                allVariants[str(tmpInd)]+=[ind]

        return allVariants

    """
        @function:
            codifySequences(self,seq):
        @description:
            codify nucleotidic sequence in numbers
            A=0,C=1,G=2,T=3
        @input:
            seq: nucleotidic sequence (containing only ACGT)
        @output:
            codified nucleotidic sequence in numbers
    """
    def codifySequences(self,seq):
        # self.appLogger.debug("codifySequences(self,ref)")
        codedRef=[]
        for item in seq:
            if "A"==item: codedRef+=[0]
            if "C"==item: codedRef+=[1]
            if "G"==item: codedRef+=[2]
            if "T"==item: codedRef+=[3]
        return np.array(codedRef)

    """
    @function:
        writeVCFFile(self, indexST,indexGT,REF,alt,variableSitesPositionIndices,HT,HL,AD,DP,flag):
    @description:
        gets single individual data formated as a GT:GL:AD:DP VCF column
    @intput:
        indexST
        indexGT
        referenceFilepath
        reference sequence
        alternative alleles
        variableSitesPositionIndices
        haplotype
        likelihood
        allelic depth
        read depth
        flag: to indicate whether is true o sampled what's being written
    @output:
        string column with information in format: GT:GL:AD:DP
    """
    def writeVCFFile(self, indexST,indexGT,referenceFilepath,REF,alt,variableSitesPositionIndices,HT,HL,AD,DP,flag):
        self.appLogger.debug("writeVCFFile(self, indexST,indexGT,REF,alt,variableSites,HT,HL,AD,DP,flag)")
        # flag is either true or sampled
        nInds=len(HT.keys())
        if flag:
            self.appLogger.info("Writing VCF file (true)")
        else:
            self.appLogger.info("Writing VCF file (sampled)")
        header="{0}\n{1}={2}\n{3}\n{4}={5}".format(\
            "##fileformat=VCFv4.0",\
            "##fileDate",\
            datetime.datetime.now(),\
            "##source=ngsphy.py",
            "##reference",\
            referenceFilepath
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
        numGeneTreeDigits=len(str(self.numLociPerST[(indexST-1)]))
        chromName="ST.{0:0{1}d}.GT.{2:0{3}d}".format(indexST,\
            self.numSpeciesTreesDigits,\
            indexGT,\
            numGeneTreeDigits)
        # POS
        nVariants=len(variableSitesPositionIndices)
        POS=[(item+1) for item in variableSitesPositionIndices]
        # ID
        ID=["ST.{0:0{1}d}.GT.{2:0{3}d}.ID.{4}".format(\
            indexST,\
            self.numSpeciesTreesDigits,\
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
        nLoci=self.numLociPerST[indexST-1]
        numGeneTreeDigits=len(str(nLoci))
        allVariants=self.formatIndividualDataForVCF(\
            REF,alt,variableSitesPositionIndices,HT,HL,AD,DP)
        outfile=""
        # flag true=true  - Flag false= sampled
        if flag:
            outfile="{0}/{2:0{3}d}/{1}_{2:0{3}d}_{4:0{5}d}_TRUE.vcf".format(\
                self.readcountTrueFolder,\
                self.dataprefix,
                indexST,\
                self.numSpeciesTreesDigits,\
                indexGT,\
                numGeneTreeDigits\
            )
        else:
            outfile="{0}/{2:0{3}d}/{1}_{2:0{3}d}_{4:0{5}d}.vcf".format(\
                self.readcountSampledFolder,\
                self.dataprefix,
                indexST,\
                self.numSpeciesTreesDigits,\
                indexGT,\
                numGeneTreeDigits\
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
        headerFields=["{0:{1}s}".format(headerCols[indexField],headerWidths[indexField]) for indexField in range(0,len(headerCols))]
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

    """
    @function:
        getColWidhts(self,chromName,POS,ID,REF,ALT,QUAL,FILTER,INFO,FORMAT,allVariants, variableSitesPositionIndices):
    @description:
        get string widht of the columns to format the VCF outpu
    @intput:
        chromName: data column
        POS: data column
        ID: data column
        REF: data column
        ALT: data column
        QUAL: data column
        FILTER: data column
        INFO: data column
        FORMAT: data column
        allVariants: data column
        variableSitesPositionIndices: data column
    @output:
        list with the leghts of the columns that are going to be written in the VCF file
    """
    def getColWidhts(self,chromName,POS,ID,REF,ALT,QUAL,FILTER,INFO,FORMAT,allVariants, variableSitesPositionIndices):
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

    """
    @function:
        writeReference(self,indexST,indexGT,referenceSpeciesID,referenceTipID,referenceSeqFull):
    @description:
        write reference sequence file separately
    @intput:
        indexST: index of the species tree to which it belongs
        indexGT: index of the gene tree to which it belongs
        referenceSpeciesID: species ID to which the reference belongs within a MSA file
        referenceTipID:  tip ID to which the reference belongs within a species in a MSA file
        referenceSeqFull: nucleotidic sequence
    @output:
        filepath where the reference will be written.
    """
    def writeReference(self,indexST,indexGT,referenceSpeciesID,referenceTipID,referenceSeqFull):
        self.appLogger.debug(" writeReference(self,indexST,indexGT,referenceSpeciesID,referenceTipID,referenceSeqFull):")
        referenceFilepath="{0}/{3:{4}}/{1}_{2}_REF_{3:{4}}_{5:{6}}.fasta".format(\
            self.referencesFolder,\
            self.projectName,\
            self.dataprefix,\
            indexST,\
            self.numSpeciesTreesDigits,\
            indexGT,\
            len(str(self.numLociPerST[indexGT-1]))\
        )
        description=">REFERENCE:{0}:ST.{1:{2}}:GT.{3:{4}}:{5}_0_{6}".format(\
            self.dataprefix,\
            indexST,\
            self.numSpeciesTreesDigits,\
            self.numLociPerST[indexGT-1],\
            len(str(self.numLociPerST[indexGT-1])),\
            referenceSpeciesID,\
            referenceTipID\
        )
        f=open(referenceFilepath,"w")
        f.write("{0}\n{1}\n".format(\
            description,\
            "".join(referenceSeqFull)
        ))
        return referenceFilepath

    def run(self):
        self.appLogger.debug("run(self)")
        status=True;    message="Read counts finished ok."
        self.appLogger.debug( "Run - read count")
        try:
            # generating folder structure for this part
            self.generateFolderStructureDetail()
            # Get list of reference sequences
            referenceList=self.parseReferenceList(self.refereceFilepath)
            # iterate over the "iterable" species trees / filtered STs
            for indexFilterST in range(0,len(self.filteredST)):
                indexST=self.filteredST[indexFilterST] # Get Proper ID for ST
                referenceForCurrST=referenceList[(indexST-1)] # get proper index for current INDEXST
                self.appLogger.debug(\
                    "Iterating over filteredST: {0} ({1}/{2})".format(\
                        indexST,\
                        (indexFilterST+1),\
                        len(self.filteredST)\
                ))
                # need the "Individual-description relation" file per species tree
                # to generate individuals
                filepathIndividualsRelation=\
                    "{0}/tables/{1}.{2:0{3}d}.individuals.csv".format(\
                        self.output,\
                        self.projectName,\
                        indexST,\
                        self.numSpeciesTreesDigits
                    )
                # Parse the file to get individual relation
                # is a dictionary
                # if ploidy=1: key: indID, content (also a dict): indexST,seqDEscription
                # if ploidy=2: key: indID, content (also a dict): indexST,speciesID, mateID1,mateID2
                individuals=self.parseIndividualRelationFile(filepathIndividualsRelation)
                # I need information for the generation of the coverage MATRIX
                # for this ST
                numIndividuals=len(individuals.keys())
                self.appLogger.debug("Number of individuals: {0}".format(numIndividuals))
                nLoci=self.numLociPerST[indexST-1]
                self.appLogger.debug("Number of loci: {0}".format(nLoci))
                coverageMatrix=self.computeCoverageMatrix(\
                    numIndividuals,nLoci,\
                    self.experimentCoverageDistro,\
                    self.individualCoverageDistro,\
                    self.locusCoverageDistro)
                # Iterate over the LOCI of this ST
                for indexGT in range(1,nLoci+1):
                    self.appLogger.debug("Iterating over nLoci: {0}/{1}".format(indexGT,nLoci))
                    numGeneTreeDigits=len(str(nLoci))
                    filepathLoc="{0}/{1:0{2}d}/{3}_{4:0{5}d}.fasta".format(\
                        self.path,\
                        indexST,\
                        self.numSpeciesTreesDigits,\
                        self.dataprefix,\
                        indexGT,\
                        numGeneTreeDigits\
                    )
                    # dictionary. indices=variable sites, content=variable nucleotide set
                    variableSites=self.extractTrueVariantsPositions(filepathLoc)
                    nVarSites=len(variableSites.keys())
                    self.appLogger.info(\
                        "Found [{0}] variable sites.".format(\
                            nVarSites))
                    # Parse MSA files - dictionary[speciesID][tipID]
                    #                   content={description,sequence}
                    msa=parseMSAFile(filepathLoc)
                    # indices to get the sequence of the REFERENCE!
                    referenceSpeciesID=str(referenceForCurrST[1])
                    referenceTipID=str(referenceForCurrST[2])
                    referenceSeqFull=msa[referenceSpeciesID][referenceTipID]['sequence']
                    referenceFilepath=self.writeReference(\
                        indexST,indexGT,\
                        referenceSpeciesID,referenceTipID,\
                        referenceSeqFull)
                    # Coverage per SNP is the same for both true/sampled dataset
                    self.appLogger.info("Getting coverage per individual per variant")
                    DP=[\
                        self.getDepthCoveragePerIndividual(nVarSites,coverageMatrix[index][(indexGT-1)]) \
                        for index in range(0,numIndividuals)]

                    if self.settings.ploidy==1:
                        self.computeHaploid(\
                            indexST,indexGT,\
                            msa,individuals,\
                            referenceFilepath,referenceSeqFull,\
                            variableSites,DP)
                    if self.settings.ploidy==2:
                        self.computeDiploid(\
                            indexST,indexGT,\
                            msa,individuals,\
                            referenceFilepath,referenceSeqFull,\
                            variableSites,DP)


                    self.appLogger.debug("So far so good")
        except ValueError as ve:
            # If there's a wrong ploidy inserted.
            status=False
            message=ve
        except TypeError as te:
            # One of the files is not fasta.
            status=False
            message=te
        return status,message
