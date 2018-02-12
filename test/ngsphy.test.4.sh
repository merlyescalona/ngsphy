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
mkdir $(pwd)/ngsphy-test4
echo "Gathering all the data in a single folder"
cp ${ngsphyPATH}/data/settings/ngsphy.settings.4.txt $(pwd)/ngsphy-test4/
cp ${ngsphyPATH}/data/indelible/control.4.txt $(pwd)/ngsphy-test4/
echo "Moving to the working directory"
cd $(pwd)/ngsphy-test4
echo "Simulating gene/species trees with SimPhy"
simphy -rs 2 -rl f:10 -sb ln:-15,1 -st u:200000,20000000 -sl f:5 -so f:1 -sp f:100000 -su f:0.00001 -si f:6 -hh ln:1.2,1 -hl ln:1.4,1 -hg f:200 -v 1 -o testwsimphy -cs 6656 -od 1 -op 1 -on 1
echo -e "Downloading INDELIble_wrapper.pl file from repository"
indelibleWrapperURL="https://raw.githubusercontent.com/adamallo/SimPhy/master/scripts/INDELIble_wrapper.pl"
if [[ $(uname -s) -eq "Linux" ]]; then
    wget $indelibleWrapperURL
elif [[ $(uname -s) -eq "Darwin" ]]; then
    curl -0 $indelibleWrapperURL
fi
echo "Simulating DNA sequences from the gene/species trees above"
perl INDELIble_wrapper.pl testwsimphy/ control.4.txt $RANDOM 1

if [[ $(uname -s) -eq "Darwin" ]]; then
    echo $(pwd)
    cd $(pwd)/ngsphy-test4/testwsimphy/1
    echo $(pwd)
    indelible
    cd $(pwd)/ngsphy-test4/testwsimphy/2
    echo $(pwd)
    indelible;
    cd $(pwd)/ngsphy-test4/
    echo $(pwd)
fi
echo "Running NGSPHY"
ngsphy -s ngsphy.settings.4.txt
