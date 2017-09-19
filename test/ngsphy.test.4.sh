#!/bin/sh
CWD=$(pwd)
# Creating test folder
mkdir ${CWD}/ngsphy-test4
# Gathering all the data in a single folder
ln -s ${CWD}/data/settings/ngsphy.settings.4.txt ${CWD}/ngsphy-test4/
ln -s ${CWD}/data/indelible/control.4.txt ${CWD}/ngsphy-test4/
# Moving to the working directory
cd ${CWD}/ngsphy-test4
# Simulating gene/species trees with SimPhy
simphy -rs 2 -rl f:10 -sb ln:-15,1 -st u:200000,20000000 -sl f:5 -so f:1 -sp f:100000 -su f:0.00001 -si f:6 -hh ln:1.2,1 -hl ln:1.4,1 -hg f:200 -v 1 -o testwsimphy -cs 6656 -od 1 -op 1 -on 1
# Simulating DNA sequences from the gene/species trees above
perl INDELIble_wrapper.pl testwsimphy/ control.4.txt $RANDOM 1
# Running NGSPHY
ngsphy -s ngsphy.settings.4.txt
