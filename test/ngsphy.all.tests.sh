#!/bin/sh
ngsphyPATH=$1
echo "Running test 1"
echo "-------------------------------------------------------------------------"
bash  $ngsphyPATH/test/ngsphy.test.1.sh $ngsphyPATH
echo  $ngsphyPATH/"Running test 1"
echo "-------------------------------------------------------------------------"
bash  $ngsphyPATH/test/ngsphy.test.2.sh $ngsphyPATH
echo  $ngsphyPATH/"Running test 1"
echo "-------------------------------------------------------------------------"
bash  $ngsphyPATH/test/ngsphy.test.3.sh $ngsphyPATH
echo  $ngsphyPATH/"Running test 1"
echo "-------------------------------------------------------------------------"
bash  $ngsphyPATH/test/ngsphy.test.4.sh $ngsphyPATH
