#!/usr/bin/env python
# -*- coding: utf-8 -*-
import argparse,csv,datetime,logging,os,subprocess,sys,threading,time,warnings
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
    if data!=0:
        return -10*np.log10(data)
    else:
        return 0
def getScoreMatrix(data):
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
        numSTs=self.settings.parser.getint("general","number_ST")

        self.refereceFilepath=self.settings.parser.get("read-count","reference")
        if self.refereceFilepath=="None":
            self.refereceFilepath=None
        else:
            self.refereceFilepath=os.path.abspath(self.refereceFilepath)

        try:
            self.seqerror=float(self.settings.parser.get("read-count","error"))
        except:
            self.seqerror=0
        self.generateFolderStructureGeneral()

    def run(self):
        status=True
        message="Read counts finished ok."

        self.appLogger.debug( "Run - read count")
        # generating folder structure for this part
        self.generateFolderStructureDetail()
        try:
            # Have ONE and ONLY ONE file with the relation of the references
            # if such file is missing then reference by default is 1_0_0
            # for all species tree replicates
            # STID, SPID, INDID
            referenceList=self.parseReferenceList(self.refereceFilepath)
            # extractTrueVariantsPositions - may raise TypeError if wrong file format
            for i in range(0,len(self.filteredST)):
                indexST=self.filteredST[i]
                self.appLogger.debug("Iterating over filteredST ({0}/{1})".format(indexST,len(self.filteredST)))
                # I have one "Individual-description relation" file per species tree
                filepath="{0}/tables/{1}.{2:0{3}d}.individuals.csv".format(\
                    self.output,\
                    self.projectName,\
                    indexST,\
                    self.numSpeciesTreesDigits
                )
                individuals=self.parseIndividualRelationFile(filepath)
                # [[row["indID"],row["indexST"],row["seqDescription"]]
                # [[row["indID"],row["indexST"],row["speciesID"], row["mateID1"],row["mateID2"]]
                numIndividuals=len(individuals)
                nLoci=self.numLociPerST[indexST-1]
                coverageMatrix=self.computeCoverageMatrix(\
                    numIndividuals,nLoci,\
                    self.experimentCoverageDistro,\
                    self.individualCoverageDistro,\
                    self.locusCoverageDistro)
                referenceForCurrST=referenceList[(indexST-1)]

                for indexGT in range(1,nLoci+1):
                    self.appLogger.debug("Iterating over nLoci: {0}/{1}".format(indexGT,nLoci))
                    numGeneTreeDigits=len(str(nLoci))
                    fastapath="{0}/{1:0{2}d}/{3}_{4:0{5}d}.fasta".format(\
                        self.path,\
                        indexST,\
                        self.numSpeciesTreesDigits,\
                        self.dataprefix,\
                        indexGT,\
                        numGeneTreeDigits\
                    )
                    # this is a dictionary - 'POS': {'Nuc': freq, ... }
                    variableSites=self.extractTrueVariantsPositions(fastapath)
                    # for item in np.sort(variableSites.keys()):
                    #     print item, variableSites[item]
                    self.appLogger.info("Found {0} variable sites.".format(len(variableSites.keys())))
                    # Parsing msa file
                    self.appLogger.debug("Parsing MSA file")
                    msa=parseMSAFile(fastapath)
                    referenceSeq=msa[str(referenceForCurrST[1])][str(referenceForCurrST[2])]['sequence']
                    # print referenceSeq[0]
                    # descriptionRefSeq=msa[str(referenceForCurrST[1])][str(referenceForCurrST[2])]['description']
                    # organizing to make the vcf table content
                    positions=np.sort([int(pos) for pos in variableSites.keys()])
                    # ref=[referenceSeq[pos] for pos in positions]
                    # altTmp=[set(variableSites[pos]) for pos in variableSites.keys()]
                    indices=np.sort([int(item) for item in variableSites.keys()])
                    alt=dict()
                    for item in indices: alt[str(item)]=[]
                    for pos in indices:
                        alt[str(pos)]=list(set(variableSites[str(pos)])-set(referenceSeq[pos]))

                    # alt=[list(set(variableSites[str(pos)])-set(referenceSeq[pos])) for pos in indices]
                    # print "ALT 0","| A=",variableSites[str(0)], "|B=",referenceSeq[0], "| A-B= ",set(variableSites[str(0)])-set(referenceSeq[0]), "| alt[0]=", alt[0]
                    self.appLogger.debug("Number of individuals: {0}".format(numIndividuals))

                    rcsTrue=[]; DP=[]; altInds=[]
                    # This coverage per SNP is the same for both true/sampled dataset
                    DP=[ self.getDepthCoveragePerIndividual(positions,coverageMatrix[index][indexGT]) for index in range(0,numIndividuals)]
                    # General info varaibles
                    HTGeneral=dict();HLGeneral=dict();GQGeneral=dict();ADGeneral=dict()
                    HTGeneralSampled=dict();HLGeneralSampled=dict();GQGeneralSampled=dict();ADGeneralSampled=dict()
                    # Initialization
                    for index in range(0,numIndividuals):
                        HTGeneral[str(index)]=[];HLGeneral[str(index)]=[];ADGeneral[str(index)]=[]
                    # first get true data
                    for index in range(0,numIndividuals):
                        self.appLogger.debug("Iterating over individuals: {0}/{1}".format(index,numIndividuals))
                        self.appLogger.debug("True")
                        # ind is a dic entry
                        # [[row["indID"],row["indexST"],row["seqDescription"]]
                        # [[row["indID"],row["indexST"],row["speciesID"], row["mateID1"],row["mateID2"]]
                        ind=individuals[index]
                        # I get single sequence if haplod/get matrix with both sequences if diploid
                        indSeq=self.getIndividualSequence(msa,ind,positions)
                        # here i get the individual sequence filtered
                        # Read Count - Per Individual - different if true/sampledInd
                        # function returns both arrays
                        ADTrue,ADSampled=self.getAllelicDepthPerIndividual(indSeq,positions,DP[index])
                        rcTrue,rcSampled=self.getReadCountPerIndividual(ADTrue,ADSampled)
                        rcsTrue+=[rcTrue]
                        newAltInd=self.getAltUpdatedPerIndividual(referenceSeq,alt,ADSampled)
                        for i in indices:
                            newAlt=set(alt[str(i)]+newAltInd[str(i)])-set(referenceSeq[i])
                            altInds+=[list(np.sort(list(newAlt)))]
                        # Get pe TRUE
                        # HT True is a dict
                        HTTrue=self.gettingHaplotype(referenceSeq,indSeq,alt)
                        HLTrue=self.haplotypeLikehood(rcTrue,ADTrue,variableSites 0)

                        # filter HLTrue / GQTrue  for the trueALT
                        # HTGeneral[HTPERINDIVIDUAL][VARIABLEPOSITION]
                        HTGeneral[str(index)]=HTTrue
                        HLGeneral[str(index)]=HLTrue
                        ADGeneral[str(index)]=ADTrue
                        self.appLogger.debug("Sampled")




                    # here corresponds to the INDEXST that is going to be written
                    self.writeVCFFile(indexST,indexGT,referenceSeq,alt,variableSites,HTGeneral,HLGeneral,ADGeneral,DP,True)
                    print "check format haploid sampled process"
                    print "check diploid process"
                    print "check format diploid sampled process"
                    print "check indels"

                    sys.exit()



        except ValueError as ve:
            # If there's a wrong ploidy inserted.
            status=False
            message=ve
        except TypeError as te:
            # One of the files is not fasta.
            status=False
            message=te
        return status,message

    def getIndividualSequence(self,msa,ind, positions):
        seqSize=len(msa[str(1)][str(1)]['sequence'])
        fullInd=None; speciesID=None; tipID=None; tmp=None
        if self.settings.ploidy==1:
            fullInd=np.chararray(shape=(1,len(positions)), itemsize=1)
            # [[row["indID"],row["indexST"],row["seqDescription"]]
            speciesID=ind[1]
            tipID=ind[2].strip().split("_")[2]
            tmp=list(msa[str(speciesID)][str(tipID)]['sequence'])
            fullInd=[tmp[pos] for pos in positions]
        if self.settings.ploidy==2:
            fullInd=np.chararray(shape=(2,len(position)), itemsize=1)
            tmp=list(msa[str(ind[2])][str(ind[3])]['sequence'])
            fullInd[0,:]=[tmp[pos] for pos in positions]
            tmp=list(msa[str(ind[2])][str(ind[4])]['sequence'])
            fullInd[1,:]=[tmp[pos] for pos in positions]
        return fullInd


    def getDepthCoveragePerIndividual(self, variableSites,startingCoverage):
        distro=ngsphydistro(0,"nb:{0},{1}".format(startingCoverage,startingCoverage))
        # since the value for SNPS is fixed to a Negative Binomial distribution
        # I iterate over the number of variants for this specific locus and
        # since everythime i hit play! (called value()) I'll sample a new value
        # from the previously set distribution, i'll then get a random value from a
        # negative bionmial distribution with mean max-coverage
        DP=[ distro.value()  for index in variableSites ]
        return DP

    def getAllelicDepthPerIndividual(self,fullInd,variableSites,DP):
        nVariants=len(variableSites)
        ADTrue=np.ndarray(shape=(4,nVariants))
        ADSampled=np.ndarray(shape=(4,nVariants))
        for indexPos in range(0,nVariants):
            # getting general coverage per position
            posCoverage=DP[indexPos]
            # information of the position whether the individual is haploid or diploid
            finalRC=[];indNucs=[]
            if self.settings.ploidy==1:
                indNucs=fullInd[indexPos]
                finalRC=[fullInd[indexPos]]*(posCoverage)
            if self.settings.ploidy==2:
                diploid=ngsphydistro(0,"b:{0},{1}".format(posCoverage,0.5))
                indNucs=fullInd[0,indexPos]+fullInd[1,indexPos]
                finalRC=[fullInd[0,indexPos]*(posCoverage/2)]+[fullInd[1,indexPos]*(posCoverage/2)]

            counter=Counter(finalRC)

            # TRUE READ count
            ADTrue[0,indexPos]=counter["A"];ADTrue[1,indexPos]=counter["C"]
            ADTrue[2,indexPos]=counter["G"];ADTrue[3,indexPos]=counter["T"]
            # SAMPLE READ COUNT
            # need to know error distribution
            errorDistro=ngsphydistro(0,"b:{0},{1}".format(posCoverage,self.seqerror))
            errorD=errorDistro.value()
            errorPositions=[]
            # need to know possible nucleotides to substitute my position with error
            possibleNucs=list(set(self.__NUCLEOTIDES)-set([indNucs]))
            # I have some positions (at least 1) that is an error
            # errorD= array with coded error nucleotides that will be modified
            if (errorD>0):
                errorChoices=np.random.choice(possibleNucs, int(errorD), replace=True)
                errorPositions=np.random.randint(posCoverage,size=int(errorD))
                for item in range(0,len(errorPositions)):
                    posToChange=errorPositions[item]
                    finalRC[posToChange]=errorChoices[item]
            counter=Counter(finalRC)
            # if any of the nucleotides does not have a counter retrieving it will be 0
            ADSampled[0,indexPos]=counter["A"];ADSampled[1,indexPos]=counter["C"]
            ADSampled[2,indexPos]=counter["G"];ADSampled[3,indexPos]=counter["T"]
        return ADTrue,ADSampled


    def getReadCountPerIndividual(self,ADTrue,ADSampled):
        nPos=ADTrue.shape[1]
        rcTrue=dict();rcSampled=dict()
        # init variantsRC
        for indexVar in range(0,nPos):
            rcTrue[str(indexVar)]=[]
            rcSampled[str(indexVar)]=[]
        # getting read count
        for indexVar in range(0,nPos):
            for indexNuc in range(0,4):
                rcTrue[str(indexVar)]+=[self.__NUCLEOTIDES[indexNuc]]*ADTrue[indexNuc,indexVar]
                rcSampled[str(indexVar)]+=[self.__NUCLEOTIDES[indexNuc]]*ADSampled[indexNuc,indexVar]
        return rcTrue, rcSampled

    def getAltUpdatedPerIndividual(self,ref,alt,AD):
        # update alt values
        # Sampled - ALT is a DICT
        # altUpdate is going to be a dict too
        altUpdated=dict()
        for item in alt.keys(): altUpdated[item]=[]
        sortedAltKeys=np.sort([int(item) for item in alt.keys()])
        for index in range(0,len(alt.keys())):
            possibleNucsAlt=[]
            if (AD[0,index]>0): possibleNucsAlt+=["A"]
            if (AD[1,index]>0): possibleNucsAlt+=["C"]
            if (AD[2,index]>0): possibleNucsAlt+=["G"]
            if (AD[3,index]>0): possibleNucsAlt+=["T"]
            newAlt=set(alt[str(sortedAltKeys[index])]+possibleNucsAlt)-set(ref[sortedAltKeys[index]])
            altUpdated[str(sortedAltKeys[index])]+=list(np.sort(list(newAlt)))
        return altUpdated


    def computeVCFColumnPerIndividualHaploid(self,msa, ind, positions,ref,alt, coverage):
        self.appLogger.debug("Starting")
        # haploid
        seqSize=len(msa[str(1)][str(1)]['sequence'])
        fullInd=np.chararray(shape=(1,len(positions)), itemsize=1)
        # [[row["indID"],row["indexST"],row["seqDescription"]]
        speciesID=ind[1]
        tipID=ind[2].strip().split("_")[2]
        tmp=list(msa[str(speciesID)][str(tipID)]['sequence'])
        fullInd=[tmp[pos] for pos in positions]
        # Need to create the distribution of the SNP based on the GLOBAL
        # coverage value
        distro=ngsphydistro(0,"nb:{0},{1}".format(coverage,coverage))
        # since the value for SNPS is fixed to a Negative Binomial distribution
        # I iterate over the number of variants for this specific locus and
        # since everythime i hit play! (called value()) I'll sample a new value
        # from the previously set distribution, i'll then get a random value from a
        # negative bionmial distribution with mean max-coverage
        DP=[ distro.value()  for index in positions ]
        # Now I have the coverage for all the variant positions
        # Initialization of AD matrices
        ADSampled=np.zeros(shape=(4,len(positions)),dtype=int)
        ADTrue=np.zeros(shape=(4,len(positions)),dtype=int)
        # iterate over variants

        for indexPos in range(0,len(positions)):
            # getting general coverage per position
            posCoverage=DP[indexPos]
            finalRC=fullInd[indexPos]*(posCoverage)
            counter=Counter(finalRC)
            # TRUE READ count
            ADTrue[0,indexPos]=counter["A"];ADTrue[1,indexPos]=counter["C"]
            ADTrue[2,indexPos]=counter["G"];ADTrue[3,indexPos]=counter["T"]
            # SAMPLE READ COUNT
            # need to know error distribution
            errorDistro=ngsphydistro(0,"b:{0},{1}".format(posCoverage,self.seqerror))
            errorD=errorDistro.value()
            errorPositions=[]
            # need to know possible nucleotides to substitute my position with error
            possibleNucs=list(set(self.__NUCLEOTIDES)-set([fullInd[indexPos]]))
            # I have some positions (at least 1) that is an error
            # errorD= array with coded error nucleotides that will be modified
            if (errorD>0):
                errorPositions=np.random.choice(possibleNucs, int(errorD), replace=True)
                finalRC=[fullInd[indexPos]]*(posCoverage-int(errorD))
                finalRC+=list(errorPositions)
            counter=Counter(finalRC)
            # if any of the nucleotides does not have a counter retrieving it will be 0
            ADSampled[0,indexPos]=counter["A"];ADSampled[1,indexPos]=counter["C"]
            ADSampled[2,indexPos]=counter["G"];ADSampled[3,indexPos]=counter["T"]

        # Get pe TRUE
        HTTrue=self.gettingHaplotype(fullInd,ref,list(alt))
        # true alt positions are on alt variable
        # Sampled - ALT can be modified
        # update alt values
        altUpdated=[]
        for index in range(0,len(alt)):
            possibleNucsAlt=[]
            if (ADSampled[0,index]>0):  possibleNucsAlt+=["A"]
            if (ADSampled[1,index]>0):  possibleNucsAlt+=["C"]
            if (ADSampled[2,index]>0):  possibleNucsAlt+=["G"]
            if (ADSampled[3,index]>0):  possibleNucsAlt+=["T"]
            newAlt=set(alt[index]+possibleNucsAlt)-set(ref[index])
            altUpdated+=[list(np.sort(list(newAlt)))]

        HTSampled=self.gettingHaplotype(fullInd,ref,altUpdated)
        # up until here I have the haplotype for the individual and the RC - AD var.
        # i need now, likelihoods
        HLSampled=self.haplotypeLikehood(ADSampled, ref,self.seqerror)
        HLTrue=self.haplotypeLikehood(ADTrue, ref, 0)
        # PLTrue=getScoreMatrix(HLTrue).astype(int)
        # PLSampled=getScoreMatrix(HLSampled).astype(int)
        # "GT:GL:AD:DP"
        trueRes=[]
        sampledRes=[]
        for index in range(0,len(positions)):
            # alt son letras necesito numeros
            trueRows=self.codifySequences(alt[index])
            sampledRows=self.codifySequences(altUpdated[index])
            trueResLine="{0}:{1}:{2}:{3}:{4}:{5}".format(\
                HTTrue[index],\
                ",".join(HLTrue[trueRows,index].astype(dtype=np.str)),\
                ",".join(ADTrue[trueRows,index].astype(dtype=np.str)),\
                DP[index])
            sampledResLine="{0}:{1}:{2}:{3}:{4}:{5}".format(\
                HTSampled[index],\
                ",".join(HLSampled[:,index].astype(dtype=np.str)),\
                ",".join(ADSampled[:,index].astype(dtype=np.str)),\
                DP[index])
            trueRes+=[trueResLine]
            sampledRes+=[sampledResLine]
        return altUpdated,trueRes,sampledRes

        #
        # for index in range(0,len(alt)):
        #     print index, ref[index],alt[index], "\t|", ref[index],altUpdated[index], "\t|| ",\
        #         fullInd[index],"HT:",HTTrue[index], " HTSampled: ", HTSampled[index],"\t|| HL:",\
        #         ",".join(HLTrue[:,index].astype(dtype=np.str)), " HLS: ",\
        #         "{0}".format(HLSampled[0,index]),\
        #         "{0}".format(HLSampled[1,index]),\
        #         "{0}".format(HLSampled[2,index]),\
        #         "{0}".format(HLSampled[3,index])
        #
        # print finalRC
        # print "end"
        # sys.exit()

    def getHaplotypeGenotypeQuality(self,HL):
        # genotype qualities:
        # https://www.biostars.org/p/115507/
        # unphred(0) / (unphred(0) + unphred(255) + unphred(255)).
        # best score/(sum of unphred all scores)
        # need to change these values to integers
        nrows,ncols=HL.shape
        HQ=[]
        for indexVar in range(0,ncols):
            nominator=min(HL[:,indexVar])
            denominator=sum(HL[:,indexVar])
            HQ+=[nominator/denominator]
        return HQ


    # returns position that coincide
    def compareSeqs(self,seqA,seqB):
        positions=[]
        if len(seqA)==len(seqB):
            for index in range(0,len(seqA)):
                if (seqA[index]==seqB[index]):
                    positions+=[index]
        return positions

    # Computes GT
    # returns a DICT[POS]=GT of that specific position
    def gettingHaplotype(self,ref,seq,alt, variableSites):
        nVariants=len(variableSites)
        indices=np.sort([int(i) for i in variableSites.keys()])
        GT=dict()
        # init variantsRC
        for indexVar in range(0,nVariants):
            GT[str(indexVar)]=0

        filteredSeq=seq[indices]
        for pos in range(0,nVariants):
            try:
                refPos=ref[indices[pos]]
                seqPos=seq[indices[pos]]
                altValues=alt[str(pos)]
                for index in range(0,len(altValues)):
                    altToCompare=altValues[index]
                    if refPos==seqPos:
                        GT[pos]=0
                    else:
                        if (index==0) and (altToCompare==seqPos):
                            GT[pos]=1
                        if (index==1) and (altToCompare==seqPos):
                            GT[pos]=2
                        if (index==2) and (altToCompare==seqPos):
                            GT[pos]=3
            except:
                pass
        return GT


    # Returs matrix 4xNPositionsVariables
    def haplotypeLikehood(self,variantsRC,observed,variableSites,error):
        # observed - read count per nuc per pos
        # matrix - rows (ACGT) cols(variantsSize)
        # alleles - possible alleles set(ref/alt[pos])
        indices=np.sort([int(i) for i in variableSites.keys()])
        nVariants=observed.shape[1]
        HL=np.ones(shape=(4,nVariants), dtype=np.float)
        error=float(error)
        for indexVar in range(0,nVariants):
            reads=variantsRC[str(indices[indexVar])]
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

    def getPossibleGenotypesPerVariableSite(self,altPos):
        d=[]; possibleGenotypes=[]
        for pos1 in altPos:
            for pos2 in altPos:
                if not (set([pos1,pos2]) in d):
                    d+=[set([pos1,pos2])]
                    possibleGenotypes+=[[pos1,pos2]]
        return possibleGenotypes

    def genotypeLikehood(self, observed, ref,alt,error):
        # observed - read count per nuc per pos
        # matrix - rows (ACGT) cols(variantsSize)
        # alleles - possible alleles set(ref/alt[pos])
        # Get Possible Genotypes
        # Get Read Counts
        nVariants=observed.shape[1]
        GL=dict() # Initialization of genotypeLikehood variable
        variantsRC=dict() # Initialization of readcount variable
        for indexVar in range(0,nVariants):
            variantsRC[str(indexVar)]=[]
            GL[str(indexVar)]={}
            possibleGenotypes=self.getPossibleGenotypesPerVariableSite(alt[indexVar])
            for pg in possibleGenotypes:
                gen="".join(pg)
                GL[str(indexVar)][gen]=1

        # getting read count
        for indexVar in range(0,nVariants):
            for indexNuc in range(0,4):
                variantsRC[str(indexVar)]+=[self.__NUCLEOTIDES]*observed[indexNuc,indexVar]

        # getting likelihoods
        for indexVar in range(0,nVariants):
            possibleGenotypes=self.getPossibleGenotypesPerVariableSite(alt[indexVar])
            for pg in possibleGenotypes:
                gen="".join(pg)
                reads=variantsRC[str(indexVar)]
                for b in range(0,len(reads)):
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
                    GL[str(indexVar)][gen]=GL[str(indexVar)][gen]*(valA1+valA2)
        return GL



    def computeVCFColumnPerIndividual(self,msa, ind, positions,ref,alt, coverage):
        # ind if haploid: 3 elemnts per cel (indexST,indID,seqDescription)
        # ind if diploid: 5 elemnts per cel (indexST,indID,speciesID,mateID1,mateID2)
        column1=None;column2=None;column3=None
        if self.settings.ploidy==1:
            column1,column2,column3=self.computeVCFColumnPerIndividualHaploid(msa, ind, positions,ref,alt, coverage)
        elif self.settings.ploidy==2:
            column1,column2,column3=self.computeVCFColumnPerIndividualDiploid(msa, ind, positions,ref,alt, coverage)
        else:
            raise ValueError("Ploidy assigned is incorrect. Please verify. Exciting.")
            # There has been a verification in Settings class, but just in case.
        return column1,column2,column3

    def codifySequences(self,ref):
        codedRef=[]
        for item in ref:
            if "A"==item:
                codedRef+=[0]
            if "C"==item:
                codedRef+=[1]
            if "G"==item:
                codedRef+=[2]
            if "T"==item:
                codedRef+=[3]
        return np.array(codedRef)


    # Parsing referenceList (STID,SPID,TIPID)
    def parseReferenceList(self, filename):
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

    # get pos
    # filename: msa fasta file from where to extract the variable positions
    def extractTrueVariantsPositions(self, filename):
        filepath=os.path.abspath(filename)
        lines=[];variants=dict();seqDescriptions=[]
        numTotalSeqs=0;lenSeq=0; matrix=None
        self.appLogger.debug("Extracting variable positions from: {0}".format(\
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
        except:
            pass
            """
            Traceback (most recent call last):
            File "<stdin>", line 1, in <module>
            ValueError: list.remove(x): x not in list
            """
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
        individuals=dict()
        if (self.settings.ploidy>0 and self.settings.ploidy<=2):
            csvfile=open(os.path.abspath(filename))
            d = csv.DictReader(csvfile)
            if (self.settings.ploidy==1):
                individuals = [[row["indID"],row["indexST"],row["seqDescription"]] for row in d]
            if (self.settings.ploidy==2):
                # indexST,indID,speciesID,mateID1,mateID2
                individuals = [[row["indID"],row["indexST"],row["speciesID"], row["mateID1"],row["mateID2"]] for row in d]
            csvfile.close()
        else:
            # There has been a verification in Settings class, but just in case.
            raise ValueError("Ploidy assigned is incorrect. Please verify. Exciting.")
        return individuals

    # PASIASDHFASGFASG!!!!
    def formatIndividualDataForVCF(self,ref,alt,variableSites,HT,HL,AD,DP):
        indices=np.sort([int(pos) for pos in variableSites.keys()])
        nVariants=len(indices)
        nInds=len(HT.keys())
        allVariants=dict()
        for indexVAR in variableSites.keys():
            allVariants[indexVAR]=[]
        # print "Init allvariants dict"
        # print len(HT.keys()),len(indices)
        for indexVAR in range(0,nVariants):
            for indexIND in range(0,nInds):
                htPerInd=HT[str(indexIND)]
                trueRows=None
                # print indexVAR, indexIND,indexVAR,"/",len(htPerInd), variableSites[str(indexVAR)]
                if htPerInd[indexVAR]==0:
                    # trueRows is the ref nucleotides
                    trueRows=self.codifySequences(ref[indices[indexVAR]])
                else:
                    trueRows=self.codifySequences(alt[str(indices[indexVAR])])
                # print "Before truerows ",DP[indexIND][indexVAR]
                ind="{0}:{1}:{2}:{3}".format(\
                    htPerInd[indexVAR],\
                    ",".join(HL[str(indexIND)][trueRows,indexVAR].astype(dtype=int).astype(dtype=np.str)),\
                    ",".join(AD[str(indexIND)][trueRows,indexVAR].astype(dtype=int).astype(dtype=np.str)),\
                    DP[indexIND][indexVAR])
                # print indexVAR,ind
                allVariants[str(indices[indexVAR])]+=[ind]
        return allVariants


    def writeVCFFile(self, indexST,indexGT,REF,alt,variableSites,HT,HL,AD,DP,flag):
        # flag is either true or sampled
        nInds=len(HT.keys())
        self.appLogger.info("Writing VCF file")
        self.refereceFilepath=""
        header="{0}\n{1}={2}\n{3}\n{4}={5}".format(\
            "##fileformat=VCFv4.0",\
            "##fileDate",\
            datetime.datetime.now(),\
            "##source=ngsphy.py",
            "##reference",\
            self.refereceFilepath
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
        nVariants=len(variableSites.keys())
        POS=np.sort([pos for pos in range(1,(nVariants+1))])
        # ID
        ID=["ST.{0:0{1}d}.GT.{2:0{3}d}.ID.{4}".format(\
            indexST,\
            self.numSpeciesTreesDigits,\
            indexGT,
            numGeneTreeDigits,\
            ID) for ID in range(1, (nVariants+1))]

        # ALT
        indices=np.sort([int(item) for item in variableSites.keys()])
        ALT=[ ",".join(alt[str(pos)]) for pos in indices ]
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
        allVariants=self.formatIndividualDataForVCF(REF,alt,variableSites,HT,HL,AD,DP)
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
        colWidths=self.getColWidhts(chromName,POS,ID,REF,alt,QUAL,FILTER,INFO,FORMAT,allVariants,variableSites)
        # (sizeChrom,sizePOS,sizeID,sizeREF,sizeALT,sizeQUAL,sizeFILTER,sizeINFO,sizeFORMAT,sizeInds)
        maxLenIndName=max([len(elem) for elem in indnames])
        maxLenIndName=max(colWidths[9], maxLenIndName)
        print maxLenIndName
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
                REF[index],colWidths[3],\
                ",".join(alt[str(indices[index])]),colWidths[4],\
                QUAL[index].encode("UTF-8"),colWidths[5],\
                FILTER[index].encode("UTF-8"),colWidths[6],\
                INFO[index].encode("UTF-8"),colWidths[7],\
                FORMAT[index].encode("UTF-8"),colWidths[8],\
                "\t".join(\
                    ["{0:{1}s}".format(indVar,maxLenIndName) for indVar in allVariants[str(indices[index])]]\
                )
            )
            filevcf.write(line)
        filevcf.close()

    def getColWidhts(self,chromName,POS,ID,REF,ALT,QUAL,FILTER,INFO,FORMAT,allVariants, variableSites):
        #CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO,FORMAT
        # allVariants is a dict
        # alt is a dict
        # all the other variables are arrays
        indices=np.sort([int(i) for i in range(0,len(variableSites.keys()))])
        sizeChrom=len(chromName)
        sizePOS=max([len(str(item)) for item in POS])
        sizeID=max([len(item) for item in ID])
        sizeREF=max([len(item) for item in REF])
        tmpALT=[",".join(ALT[var]) for var in variableSites.keys()]
        sizeALT=max([len(item) for item in tmpALT])
        sizeQUAL=max([len(item) for item in QUAL])
        sizeFILTER=max([len(item) for item in FILTER])
        sizeINFO=max([len(item) for item in INFO])
        sizeFORMAT=max([len(item) for item in FORMAT])
        tmpInds=[]
        for item in allVariants.keys():
            tmpInds+=[max([ len(elem) for elem in allVariants[item]])]
        # CHECK LEN OF STRINGS GENERATED HERE!!!
        sizeInds=max(tmpInds)
        return [sizeChrom,sizePOS,sizeID,sizeREF,sizeALT,sizeQUAL,sizeFILTER,sizeINFO,sizeFORMAT,sizeInds]

    def generateFolderStructureDetail(self):
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

    def generateFolderStructureGeneral(self):
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

    def computeCoverageMatrix(self, nInds, nLoci,expCov, indCov, locCov):
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

    def writeCoverageMatrixIntoFile(self, coverageMatrix, filename):
        filepath=os.path.abspath(filename)
        with open(filepath, 'w') as csvfile:
            writer = csv.writer(csvfile)
            [writer.writerow(r) for r in coverageMatrix]
