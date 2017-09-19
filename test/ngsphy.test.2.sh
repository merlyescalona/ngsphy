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
mkdir $(pwd)/ngsphy-test2
echo "Gathering all the data in a single folder"
cp ${ngsphyPATH}/data/settings/ngsphy.settings.2.txt $(pwd)/ngsphy-test2/
cp ${ngsphyPATH}/data/indelible/control.2.txt $(pwd)/ngsphy-test2/
cp ${ngsphyPATH}/data/trees/t2.tree $(pwd)/ngsphy-test2/
cp ${ngsphyPATH}/data/sequences/my_ancestral_sequence.tar.gz $(pwd)/ngsphy-test2/
echo "Moving to the working directory"
cd $(pwd)/ngsphy-test2
echo "Extracting anchor sequence"
tar -xzf my_ancestral_sequence.tar.gz
echo "Running NGSPHY"
ngsphy -s ngsphy.settings.2.txt
