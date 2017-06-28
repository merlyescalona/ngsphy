import matplotlib.pyplot as plt
import argparse,datetime,logging,os,subprocess,sys
import numpy as np
import random as rnd
import Mating as mat
import Settings as sp
import NGSReads as ngs
from MELoggingFormatter import MELoggingFormatter as mlf
from select import select

# Variable for poisson distro
loc=0
class NGSReadsReadCount:
    def __init__(self,refSequence, errorRate, coverage):
        self.refSequence=refSequence
        self.sequenceSize=len(refSequence)
        self.errorRate=errorRate
        self.expCoverage=coverage
        self.coverageTable=[0]*self.sequenceSize

    def getScore(self):
        q=str(unichr(-10*np.log10(self.errorRate)+33))

    def computeCoveragePoisson(self):
        from scipy.stats import poisson
        mu = self.sequenceSize/2
        x = map(int,np.arange(poisson.ppf(0.1, mu),poisson.ppf(0.9, mu)))
        print x
        covDistro=(poisson.pmf(x, mu)*self.expCoverage)/max(poisson.pmf(x, mu))
        for index in range(0,len(x)):
            self.coverageTable[x[index]]=covDistro[index]
        for index in range(0,x[0]):
            self.coverageTable[index]=covDistro[0]
        for index in range(x[len(x)-1],self.sequenceSize):
            self.coverageTable[index]=covDistro[len(x)-1]
        fig, ax = plt.subplots(1, 1)
        print self.coverageTable
        print range(1,self.sequenceSize+1)
        ax.plot(range(1,self.sequenceSize+1),self.coverageTable, 'bo', ms=8, label='poisson pmf')
        plt.show()

    def computeCoverageBinom(self):
        from scipy.stats import binom
        mu = self.sequenceSize/2
        fig, ax = plt.subplots(1, 1)
        n, p = 200, 0.5
        x = map(np.arange(binom.ppf(0.001, n, p),binom.ppf(0.999, n, p)))
        print x
        covDistro=(binom.pmf(x, n, p)*self.expCoverage)/max(binom.pmf(x, n, p))
        for index in range(0,len(x)):
            self.coverageTable[x[index]]=covDistro[index]
        for index in range(0,x[0]):
            self.coverageTable[index]=covDistro[0]
        for index in range(x[len(x)-1],self.sequenceSize):
            self.coverageTable[index]=covDistro[len(x)-1]
        fig, ax = plt.subplots(1, 1)
        print self.coverageTable
        print range(1,self.sequenceSize+1)
        ax.plot(range(1,self.sequenceSize+1),self.coverageTable, 'bo', ms=8, label='poisson pmf')
        plt.show()

    def computeCoverageLogNormal(self):
        from scipy.stats import lognorm
        fig, ax = plt.subplots(1, 1)
        s = 0.954
        mean, var, skew, kurt = lognorm.stats(s, moments='mvsk')
        x = np.linspace(lognorm.ppf(0.01, s),lognorm.ppf(0.99, s), 100)
        mu = self.sequenceSize/2
        fig, ax = plt.subplots(1, 1)
        n, p = 200, 0.5
        x = map(np.arange(binom.ppf(0.001, n, p),binom.ppf(0.999, n, p)))
        print x
        covDistro=(binom.pmf(x, n, p)*self.expCoverage)/max(binom.pmf(x, n, p))
        for index in range(0,len(x)):
            self.coverageTable[x[index]]=covDistro[index]
        for index in range(0,x[0]):
            self.coverageTable[index]=covDistro[0]
        for index in range(x[len(x)-1],self.sequenceSize):
            self.coverageTable[index]=covDistro[len(x)-1]
        fig, ax = plt.subplots(1, 1)
        print self.coverageTable
        print range(1,self.sequenceSize+1)
        ax.plot(range(1,self.sequenceSize+1),self.coverageTable, 'bo', ms=8, label='poisson pmf')
        plt.show()

if __name__ == '__main__':
    referenceSequence=["A"]*100
    error=0.1
    coverage=10
    d=NGSReadsReadCount(referenceSequence,error,coverage)
    d.computeCoveragePoisson()



from scipy.stats import exponential
mu = self.sequenceSize/2
fig, ax = plt.subplots(1, 1)
n, p = 200, 0.5
x = map(np.arange(binom.ppf(0.001, n, p),binom.ppf(0.999, n, p)))
print x
covDistro=(binom.pmf(x, n, p)*self.expCoverage)/max(binom.pmf(x, n, p))
for index in range(0,len(x)):
    self.coverageTable[x[index]]=covDistro[index]
for index in range(0,x[0]):
    self.coverageTable[index]=covDistro[0]
for index in range(x[len(x)-1],self.sequenceSize):
    self.coverageTable[index]=covDistro[len(x)-1]
fig, ax = plt.subplots(1, 1)
print self.coverageTable
print range(1,self.sequenceSize+1)
ax.plot(range(1,self.sequenceSize+1),self.coverageTable, 'bo', ms=8, label='poisson pmf')
plt.show()


import matplotlib.pyplot as plt
from scipy.stats import nbinom
import numpy as np
mu=80
r=50
p=(r*1.0)/(r+mu)
distro=nbinom(r,p)
fig, ax = plt.subplots(1, 1)
x = np.arange(distro.ppf(0.001),distro.ppf(0.999))
ax.plot(x, distro.pmf(x), 'bo', ms=8, label='nbinom pmf')
ax.vlines(mu, 0, distro.pmf(mu), colors='b', lw=5, alpha=0.5)
ax.set_title(distro.stats(moments="m"))
plt.show()


import matplotlib.pyplot as plt
from scipy.stats import gamma
import numpy as np
mu=80
r=50
p=(r*1.0)/(r+mu)
distro=nbinom(r,p)
fig, ax = plt.subplots(1, 1)
x = np.arange(distro.ppf(0.001),distro.ppf(0.999))
ax.plot(x, distro.pmf(x), 'bo', ms=8, label='nbinom pmf')
ax.vlines(mu, 0, distro.pmf(mu), colors='b', lw=5, alpha=0.5)
ax.set_title(distro.stats(moments="m"))
plt.show()
