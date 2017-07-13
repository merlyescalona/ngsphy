#!/bin/bash
cd $HOME/git/ngsphy
sudo rm -rf dist ngsphy.egg-info
sudo rm -rf /usr/local/lib/python2.7/dist-packages/ngsphy-0.1-py2.7.egg/ngsphy
sudo rm /usr/local/bin/ngsphy
python setup.py sdist
cd $HOME/git/ngsphy/dist
tar -xzvf ngsphy-0.1.tar.gz
cd ngsphy-0.1
sudo python setup.py install
