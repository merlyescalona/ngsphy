#!/bin/sh
CWD=$(pwd)
# Creating test folder
mkdir ${CWD}/ngsphy-test3
# Gathering all the data in a single folder
ln -s ${CWD}/data/settings/ngsphy.settings.3.txt ${CWD}/ngsphy-test3/
ln -s ${CWD}/data/indelible/control.3.txt ${CWD}/ngsphy-test3/
ln -s ${CWD}/data/trees/t3.tree ${CWD}/ngsphy-test3/
ln -s ${CWD}/data/sequences/my_anchor_sequence.tar.gz ${CWD}/ngsphy-test3/
ln -s ${CWD}/data/reference_alleles/my_reference_allele_file.txt ${CWD}/ngsphy-test3/
# Moving to the working directory
cd ${CWD}/ngsphy-test3
# Extracting anchor sequence
tar -xzf my_anchor_sequence.tar.gz
# Running NGSPHY
ngsphy -s ngsphy.settings.3.txt
