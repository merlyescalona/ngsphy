#!/bin/sh
CWD=$(pwd)
# Creating test folder
mkdir ${CWD}/ngsphy-test2
# Gathering all the data in a single folder
ln -s ${CWD}/data/settings/ngsphy.settings.2.txt ${CWD}/ngsphy-test2/
ln -s ${CWD}/data/indelible/control.2.txt ${CWD}/ngsphy-test2/
ln -s ${CWD}/data/trees/t2.tree ${CWD}/ngsphy-test2/
ln -s ${CWD}/data/sequences/my_ancestral_sequence.tar.gz ${CWD}/ngsphy-test2/
# Moving to the working directory
cd ${CWD}/ngsphy-test2
# Extracting ancestral sequence
tar -xzf my_ancestral_sequence.tar.gz
# Running NGSPHY
ngsphy -s ngsphy.settings.2.txt
