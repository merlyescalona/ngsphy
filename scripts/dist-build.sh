#!/bin/bash
version="1.0.0"
cd $HOME/git/ngsphy
sudo rm -rf dist ngsphy.egg-info
sudo rm -rf "/usr/local/lib/python2.7/dist-packages/ngsphy-${version}-py2.7.egg/ngsphy"
sudo rm /usr/local/bin/ngsphy
python setup.py sdist
cd $HOME/git/ngsphy/dist
tar -xzvf "ngsphy-${version}.tar.gz"
cd "ngsphy-${version}"
sudo python setup.py install
