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
mkdir $(pwd)/ngsphy-test1
echo "Gathering all the data in a single folder"
cp ${ngsphyPATH}/data/settings/ngsphy.settings.1.txt $(pwd)/ngsphy-test1/
cp ${ngsphyPATH}/data/indelible/control.1.txt $(pwd)/ngsphy-test1/
cp ${ngsphyPATH}/data/trees/t1.tree $(pwd)/ngsphy-test1/
echo "Moving to the working directory"
cd $(pwd)/ngsphy-test1
echo "Running NGSPHY"
ngsphy -s ngsphy.settings.1.txt
