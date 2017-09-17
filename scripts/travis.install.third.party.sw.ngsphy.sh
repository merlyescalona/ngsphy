#!/bin/sh
echo "Installing art"
wget https://www.niehs.nih.gov/research/resources/assets/docs/artbinmountrainier20160605linux64tgz.tgz
tar -zxvf  artbinmountrainier20160605linux64tgz.tgz
export PATH=$PATH:$PWD/art_bin_MountRainier/
git clone https://github.com/merlyescalona/indelible-ngsphy
cd indelible-ngsphy && make && popd
export PATH=$PATH:$PWD/indelible-ngsphy/bin
