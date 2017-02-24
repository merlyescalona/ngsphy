#!/usr/bin/env python
# -*- coding: utf-8 -*-
import argparse,csv,datetime,logging,os,subprocess,sys,threading,time
import numpy as np
import random as rnd
from scipy.stats import  binom,expon,gamma,lognorm,norm,nbinom,poisson,uniform

class NGSPhyDistribution:

    __relationNumParams=dict({\
        "b":2,\
        "e":1,\
        "f":1,\
        "g":2,\
        "ln":2,\
        "n":2,\
        "nb":2,\
        "p":1,\
        "r":1,\
        "u":1\
    })
    __DISTRIBUTIONS=__relationNumParams.keys()
    __value=0
    __params=None
    __distro=None
    __type=None
    __dependsOnExperimentDistro=-1
    __dependsOnIndividualDistro=-1

    def __init__(self,type,distroline, flag=False):
        self.__type=type
        distroline=distroline.split(":")
        if not (len(distroline)==1):
            self.__distro=distroline[0].lower()
            self.__params=distroline[1].split(",")
        if flag:
            self.checkDistribution()

    def params(self):
        print(self.__params)

    def binom(self):
        #n=number of times tested
        #p=probability
        n=float(self.__params[0])
        p=float(self.__params[1])
        distro=binom(n,p)
        f=distro.rvs(size=1)
        self.__value=int(f)
        return self.__value

    def exponential(self):
        self.__value=int(np.random.exponential(float(self.__params[0]),1)[0])
        return self.__value

    def fixed(self):
        self.__value=int(self.__params[0])
        return self.__value

    def gamma(self):
        # The parameterization with alpha and beta is more common in Bayesian statistics,
        # where the gamma distribution is used as a conjugate prior distribution
        # for various types of inverse scale (aka rate) parameters, such as the
        #lambda of an exponential distribution or a Poisson distribution[4]  or for
        #t hat matter, the beta of the gamma distribution itself.
        #(The closely related inverse gamma distribution is used as a conjugate
        # prior for scale parameters, such as the variance of a normal distribution.)
        # shape, scale = 2., 2. # mean=4, std=2*sqrt(2)
        # s = np.random.gamma(shape, scale, 1000)
        shape=float(self.__params[0]*1.0)
        beta=float(1/self.__params[1]*1.0)
        distro=gamma(shape,beta)
        self.__value=int(distro.rvs(size=1))
        return self.__value

    def lognorm(self):
        # m=mean, s=sigma - standard deviation
        mean=float(self.__params[0]*1.0)
        sigma=float(self.__params[1]*1.0)
        distro=lognorm(mean,sigma)
        self.__value=int(distro.rvs(size=1))
        return self.__value

    def nbinom(self):
        # mu= possion mean
        # r controls the deviation from the poisson
        # This makes the negative binomial distribution suitable as a robust alternative to the Poisson,
        # which approaches the Poisson for large r, but which has larger variance than the Poisson for small r.
        r=float(self.__params[0])
        mu=float(self.__params[1])
        p=(mu*1.0)/(mu+r)
        distro=nbinom(r,p)
        self.__value=int(distro.rvs(size=1))
        return self.__value

    def normal(self):
        # mean (location) - variance (squared scale)
        locMean=float(self.__params[0]*1.0)
        scaleVariance=float(self.__params[1]*1.0)
        distro=norm(loc=locMean,scale=scaleVariance)
        self.__value=int(distro.rvs(size=1))
        return self.__value

    def poisson(self):
        l=float(self.__params[0]*1.0)
        distro=poisson(l)
        self.__value=int(distro.rvs(size=1))
        return self.__value

    def random(self):
        self.__value=int(np.random.sample(1)*float(self.__params[0]))
        return self.__value

    def uniform(self):
        # mean
        mean=float(self.params[0]*1.0)
        u=uniform(mean)
        self.__value=int(u.rvs(size=1))
        return self.__value

    def checkDistribution(self):
        if self.__type==0:
            #TYPE 0 - EXPERIMENT WISE Coverage
            return self.experimentCoverageDistroCheck()
        elif self.__type==1:
            #TYPE 1 - Individual WISE Coverage
            return self.individualCoverageDistroCheck()
        elif self.__type==2:
            #TYPE 1 - Individual WISE Coverage
            return self.lociCoverageDistroCheck()
        else:
            return False, "Error with the coverage type. Please verify. Exiting."


    def experimentCoverageDistroCheck(self):
        messageOk="Experiment Coverage distribution parsing OK. "
        distroOk=False
        if self.__distro in self.__DISTRIBUTIONS:
            distroOk=True
            # check cas distro - numparameters
            if (self.__relationNumParams[self.__distro]==len(self.__params)):
                # check parameters are numbers
                for index in range(0,len(self.__params)):
                    try:
                        self.__params[index]=float(self.__params[index])
                    except ValueError:
                        messageOk="Experiment coverage distribution: Not a correct parameter value. Please verify. Exiting"
                        distroOk=False
                        break
            else:
                messageOk="Not the correct number of parameters for the selected distribution. Please verify. Exiting"
                distroOk=False
        else:
            messageOk="Distribution selected is not correct or unavailable. Please verify. Exiting"
        return distroOk, messageOk

    def individualCoverageDistroCheck(self):
        distroOk=False
        messageOk="Individual Coverage distribution parsing OK. "
        if self.__distro in self.__DISTRIBUTIONS:
            distroOk=True
            # check cas distro - numparameters
            if (self.__relationNumParams[self.__distro]==len(self.__params)):
                # print("got the number of parameter")
                # check parameters are numbers
                for index in range(0,len(self.__params)):
                    if (self.__params[index] in ["ec","experimentcoverage","expcov"]):
                        self.__dependsOnExperimentDistro=index
                        # print("Distro based")
                    else:
                        try:
                            self.__params[index]=float(self.__params[index])
                        except ValueError:
                            messageOk="Individual coverage distribution: Not a correct parameter value. Please verify. Exiting"
                            distroOk=False
                            break
                # print(self.__distro,self.__params)
                # Finished checking all the parameters, but there's no depdency - makes no sense - quit.
                if self.__dependsOnExperimentDistro <0:
                    messageOk="Individual coverage distribution: Not a correct parameter value. Must be depdendent of *experimentCoverage*. Please verify. Exiting"
                    distroOk=False
            else:
                messageOk="Not the correct number of parameters for the selected distribution. Please verify. Exiting"
                distroOk=False
        else:
            messageOk="Distribution selected is not correct or unavailable. Please verify. Exiting"
        return distroOk, messageOk

    def lociCoverageDistroCheck(self):
        distroOk=False
        messageOk="Individual Coverage distribution parsing OK. "
        if self.__distro in self.__DISTRIBUTIONS:
            distroOk=True
            # check cas distro - numparameters
            if (self.__relationNumParams[self.__distro]==len(self.__params)):
                # print("got the number of parameter")
                # check parameters are numbers
                for index in range(0,len(self.__params)):
                    if (self.__params[index] in ["ic","individualcoverage","indcov"]):
                        self.__dependsOnIndividualDistro=index
                        # print("Distro based")
                    else:
                        try:
                            self.__params[index]=float(self.__params[index])
                        except ValueError:
                            messageOk="Loci coverage distribution: Not a correct parameter value. Please verify. Exiting"
                            distroOk=False
                            break
                # print(self.__distro,self.__params)
                # Finished checking all the parameters, but there's no depdency - makes no sense - quit.
                if self.__dependsOnExperimentDistro <0:
                    messageOk="Loci coverage distribution: Not a correct parameter value. Must be depdendent of *individualCoverage*. Please verify. Exiting"
                    distroOk=False
            else:
                messageOk="Not the correct number of parameters for the selected distribution. Please verify. Exiting"
                distroOk=False
        else:
            messageOk="Distribution selected is not correct or unavailable. Please verify. Exiting"
        return distroOk, messageOk


    def dependsOnExperimentDistro(self):
        return self.__dependsOnExperimentDistro

    def updateDependingParamValue(self,value):
        if self.__type==1:
            if self.__dependsOnExperimentDistro >= 0:
                self.__params[self.__dependsOnExperimentDistro]=value
        if self.__type==2:
            if self.__dependsOnIndividualDistro >= 0:
                self.__params[self.__dependsOnIndividualDistro]=value
        # print(self.__params)

    def value(self):
        if self.__distro=="b":
            self.binom()
        elif self.__distro=="e":
            self.exponential()
        elif self.__distro=="f":
            self.fixed()
        elif self.__distro=="g":
            self.gamma()
        elif self.__distro=="ln":
            self.lognormal()
        elif self.__distro=="n":
            self.normal()
        elif self.__distro=="nb":
            self.nbinom()
        elif self.__distro=="p":
            self.poisson()
        elif self.__distro=="u":
            self.uniform()
        else:
            # if i got to this point I have:
            # 1) checked that the distribution is correct
            # 2) that params belong to the distribution
            # 3) that params are numbers
            # Also the only distro missing is F (Fixed value)
            self.__value=0
        return self.__value

    def distroDescription(self):
        d=""
        if self.__distro=="b":
            d="binom"
        elif self.__distro=="e":
            d="exponential"
        elif self.__distro=="f":
            d="fixed value"
        elif self.__distro=="g":
            d="gamma"
        elif self.__distro=="ln":
            d="lognormal"
        elif self.__distro=="n":
            d="normal"
        elif self.__distro=="nb":
            d="negative binomial"
        elif self.__distro=="p":
            d="poisson"
        elif self.__distro=="u":
            d="uniform"
        else:
            # if i got to this point I have:
            # 1) checked that the distribution is correct
            # 2) that params belong to the distribution
            # 3) that params are numbers
            d="unknown"

        return "{} - {}".format(d,",".join(self.__params))
