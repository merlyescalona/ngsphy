def getScoreSingle(data):
def getScoreMatrix(data):
    def __init__(self, filename):
    def addLine(self, line):
    def __init__(self,settings):
    def generateFolderStructureGeneral(self):
    def generateFolderStructureDetail(self):
    def computeCoverageMatrix(self, nInds, nLoci,expCov, indCov, locCov):
    def writeCoverageMatrixIntoFile(self, coverageMatrix, filename):
        If "None" inputted (file is missing) then reference by default is 1_0_0
    def parseReferenceList(self, filename):
                        "A default reference has been introduced.\n"+\
    def extractTrueVariantsPositions(self, filename):
    def parseIndividualRelationFile(self,filename):
    def getDepthCoveragePerIndividual(self, numVarSites,startingCoverage):
    def getHaploidIndividualSequence(self,msa,ind):
    def getDiploidIndividualSequence(self,msa,ind):
    def computeHaploid(self,indexST,indexGT,msa,individuals,referenceFilepath,referenceSeqFull,variableSites,DP):
    def haplotypeLikehood(self,variantsRC,variableSitesPositionIndices,error):
    def gettingHaplotype(self,ref,seq,alt, variableSitesPositionIndices):
    def getAltUpdatedPerIndividual(self,ref,alt,AD):
    def getReadCountPerIndividual(self,ADTrue,ADSampled, variableSitesPositionIndices):
    def getAllelicDepthPerHaploidIndividual(self,individualSeq,variableSitesPositionIndices,DP):
    def computeDiploid(self,indexST,indexGT,msa,individuals,referenceFilepath,referenceSeqFull,variableSites,DP):
    def genotypeLikehood(self, variantsRC,variableSitesPositionIndices,error):
        def getAllPossibleGenotypes(self):
    def getAllPossibleGenotypes(self):
    def getPossibleGenotypesPerVariableSite(self,ref,alt, variableSitesPositionIndices):
    def genotypeOrder(self,alleles):
    def getAllelicDepthPerDiploidIndividual(self,individualSeq,variableSitesPositionIndices,DP):
    def formatIndividualDataForVCF(self,ref,alt,variableSitesPositionIndices,HT,HL,AD,DP):
    def codifySequences(self,seq):
    def writeVCFFile(self, indexST,indexGT,referenceFilepath,REF,alt,variableSitesPositionIndices,HT,HL,AD,DP,flag):
    def getColWidhts(self,chromName,POS,ID,REF,ALT,QUAL,FILTER,INFO,FORMAT,allVariants, variableSitesPositionIndices):
    def writeReference(self,indexST,indexGT,referenceSpeciesID,referenceTipID,referenceSeqFull):
    def launchCommand(self, referenceForCurrST, indexST,indexGT, individuals,coverageMatrix):
    def run(self):
