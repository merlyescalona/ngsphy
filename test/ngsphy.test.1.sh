#!/bin/sh
CWD=$(pwd)
# Creating test folder
mkdir ${CWD}/ngsphy-test1
# Gathering all the data in a single folder
ln -s ${CWD}/data/settings/ngsphy.settings.1.txt ${CWD}/ngsphy-test1/
ln -s ${CWD}/data/indelible/control.1.txt ${CWD}/ngsphy-test1/
ln -s ${CWD}/data/trees/t1.tree ${CWD}/ngsphy-test1/
# Moving to the working directory
cd ${CWD}/ngsphy-test1
# Running NGSPHY
ngsphy -s ngsphy.settings.1.txt
