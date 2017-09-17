#!/bin/sh
CWD=$(pwd)
mkdir ${CWD}/ngsphy-test1
ln -s ${CWD}/data/settings/ngsphy.settings.1.txt ${CWD}/ngsphy-test1/
ln -s ${CWD}/data/indelible/control.1.txt ${CWD}/ngsphy-test1/
ln -s ${CWD}/data/trees/t1.tree ${CWD}/ngsphy-test1/
cd ${CWD}/ngsphy-test1
${CWD}/scripts/ngsphy -s ngsphy.settings.1.txt
