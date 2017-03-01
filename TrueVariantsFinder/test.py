import argparse,datetime,logging,os, sys
from collections  import Counter
from MELoggingFormatter import MELoggingFormatter as mlf
import numpy as np

# Reading file list
filelistPath="files.txt"
fl=open(filelistPath)
filelist=fl.readlines()
fl.close()
del fl

numFiles=len(filelist)
for index in range(0,numFiles):
  filelist[index]=filelist[index][0:-1]

# Checking sequence length
filepath=os.path.abspath(filelist[np.random.randint(0,numFiles,1)])
f=open(filepath)
seq=f.readline()
seq=f.readline()
f.close()

lenSeqs=len(seq.strip())
numTotalSeqs=numFiles*2
matrix=np.chararray((numTotalSeqs,lenSeqs), itemsize=1)
variants=dict()
matrixIndex=0

individualsRelation=dict()
# Iterating over filelist
for indexFile in range(0,len(filelist)):
    filepath=os.path.abspath(filelist[indexFile])
    f=open(filepath)
    lines=f.readlines()
    f.close()
    print(filelist[indexFile])
    for index in range(0,len(lines)):
      lines[index]=lines[index].strip()
    try:
      lines.remove("")
    except:
      pass
    numLinesFile=len(lines)
    numSeqs=len(lines)/2
    index=0
    # If I got here, i can do what i want
    indexSeqs=range(1,numLinesFile,2)
    for index in range(0,numSeqs):
        try:
            val=individualsRelation[os.path.basename(filepath)]
            individualsRelation[os.path.basename(filepath)]+=[matrixIndex]
        except:
            individualsRelation[os.path.basename(filepath)]=[matrixIndex]
        matrix[matrixIndex,:]=list(lines[indexSeqs[index]])
        matrixIndex+=1



for indexCol in range(0,matrix.shape[1]):
    c=Counter(matrix[:,indexCol])
    l=np.unique(matrix[:,indexCol])
    if (len(l)>1):
        variants[str(indexCol)]=dict(c)


refFile="reference.fasta"
refF=open(os.path.abspath(refFile))
refseq=refF.readlines()
refF.close()
refseq=refseq[1][0:-1]


data=np.chararray(shape=(len(variants),9+numFiles), itemsize=50)
pos=np.sort([int(newpos) for newpos in variants.keys()])
dashes=["-"]*len(pos)
data[:,0]=[os.path.splitext(os.path.basename(refFile))[0]]*len(pos)

data[:,2]=dashes;data[:,5]=dashes;data[:,6]=dashes;data[:,7]=dashes;
data[:,8]=["GT"]*len(pos)
data[:,1]=pos
ref=[];alt=[]
for var in range(0,len(pos)):
    varPos=pos[var]
    ref+=[refseq[varPos]]
    alt+=[",".join(variants[str(varPos)].keys())]

data[:,3]=ref
data[:,4]=alt

for var in range(0,len(pos)):
  # var is pos
  for ind in range(0,len(individualsRelation.keys())):
    indSeqs=matrix[individualsRelation[individualsRelation.keys()[ind]]]
    varInds=indSeqs[:,str(pos[var])]
    data[var,ind+9]="/".join(varInds)
    print varInds


################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
import copy

truedata=copy.deepcopy(data)
moddata=copy.deepcopy(data)
moddata[:,8]=["GT:GL"]*len(moddata)
gt=list(moddata[:,4])
maxDepth=10
numVars=len(moddata)
# Coverage
M=np.random.poisson(maxDepth, 1)
matrixVariantSampling=np.chararray((numVars,M), itemsize=1)
for p in range(0,numVars):
    pair=gt[p]
    nucs=pair.split(",")
    variationNucs=len(nucs)
    samples=[ n*(int(M)/variationNucs) for n in nucs]
    samples=[ s for s in "".join(samples)]
    matrixVariantSampling











################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
# WRITING VCF
# Mandatory fields
# The header line names the 8 fixed, mandatory columns. These columns are as follows:
#CHROM POS ID REF ALT QUAL FILTER INFO
# header
header="""
##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=Integer,Description="Genotype">
"""
indnames=[os.path.splitext(os.path.basename(ff))[0] for ff in filelist]
# filename, file_extension = os.path.splitext('/path/to/somefile.ext')
headerCols=["#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT"]+indnames

import csv
# write it
with open('test_file.vcf', 'w') as filevcf:
    filevcf.write(header)

with open('test_file.vcf', 'a') as csvfile:
    writer = csv.writer(csvfile, delimiter="\t")
    writer.writerow(headerCols)
    [writer.writerow(r) for r in data]
