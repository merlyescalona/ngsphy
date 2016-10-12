#!/usr/bin/home/python
import argparse,copy,datetime,logging,os,sys, sqlite3
import numpy as np
import random as rnd
db="test_wrapper/test_wrapper.db"
con = sqlite3.connect(db)
query="select SID from Species_Trees WHERE Ind_per_sp % 2 = 0"
res=con.execute(query).fetchall()
con.close()
for triplet in res:
    indexST=triplet[0]
    leaves=triplet[1]
    nIndsPerSp=triplet[2]
    nInds=(leaves-1)/2
# sample without replacement
# remove data from sample array
