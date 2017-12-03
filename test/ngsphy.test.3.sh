#!/bin/sh
ngsphyPATH=""
if [[ $# -eq 0 ]]; then
    echo -e "Error: ngsphy source path is missing.\nPlease verify. Exiting."
    exit -1
elif [[ $# -eq 1 ]]; then
    ngsphyPATH=$1
else
    echo -e "Error: wrong number of parameters.\nPlease verify. Exiting."
    exit -1
fi
if [[ ! -d $ngsphyPATH ]]; then
    echo -e "Error: ngsphy source path is incorrect or does not exist.\nPlease verify. Exiting."
    exit -1
fi
echo "Creating test folder"
mkdir $(pwd)/ngsphy-test3
echo "Gathering all the data in a single folder"
cp ${ngsphyPATH}/data/settings/ngsphy.settings.3.txt $(pwd)/ngsphy-test3/
cp ${ngsphyPATH}/data/indelible/control.3.txt $(pwd)/ngsphy-test3/
cp ${ngsphyPATH}/data/trees/t3.tree $(pwd)/ngsphy-test3/
cp ${ngsphyPATH}/data/sequences/my_anchor_sequence.tar.gz $(pwd)/ngsphy-test3/
cp ${ngsphyPATH}/data/reference_alleles/my_reference_allele_file.txt $(pwd)/ngsphy-test3/
echo "Moving to the working directory"
cd $(pwd)/ngsphy-test3
echo "Extracting anchor sequence"
tar -xzf my_anchor_sequence.tar.gz
echo "Running NGSPHY"
ngsphy -s ngsphy.settings.3.txt
