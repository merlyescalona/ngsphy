#!/usr/bin/home/python
import argparse,copy,datetime,logging,os,sys, sqlite3
import numpy as np
import random as rnd
db="test_wrapper/test_wrapper.db"
con = sqlite3.connect(db)
query="select SID, Leaves, Ind_per_sp from Species_Trees WHERE Ind_per_sp % 2 = 0"
res=con.execute(query).fetchall()
con.close()



# NOTES!!!

# I'm assuming that there's a unique outgroup.
# Simphy generates either 1 or none
# have to check how to get information ong the outgroup Generation
# probably have to add a new parameter to the settings in case not using
# simphy whether if the user have an outgroup or not.
for triplet in res:
    indexST=triplet[0]
    leaves=triplet[1]
    nIndsPerSp=triplet[2]
    nInds=(leaves-1)/2
    inds=range(0,nIndsPerSp)
    species=range(1,leaves)
    mates=[]
    print "indexST: {0} / inds:{1} ".format(indexST,inds)
    # I'm always assuming there's an outgroup
    for sp in species:
        t=copy.deepcopy(inds)
        while not t==[]:
            p1=0
            p2=0
            try:
                p1=t.pop(rnd.sample(range(0,len(t)),1)[0])
                p2=t.pop(rnd.sample(range(0,len(t)),1)[0])
            except Exception as e:
                break
            pair=(indexST,sp,p1,p2)
            print pair
