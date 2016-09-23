# Data simulation for SimPhyNGSWrapper test
# Assuming "SimPhy" is fully installed.
# Also, assuming "INDELible" have been ran.

WD="/home/merly/git/simphy-ngs-wrapper"
controlFile="$WD/scripts/indelible.control.txt"
wrapper="$WD/scripts/INDELIble_wrapper.pl"
###############################################################################
pipelinesName="test_wrapper"
simphy_folder="$WD/$pipelinesName"
stReplicates=2
numDigits=${#stReplicates}
nGTs="f:10"
numGeneTrees="f:${nGTs}"
numDigitsGTs=${#nGTs}
sizeLoci="500bp"
prefix="data"
timeRange="u:200000,20000000" # 200k to 20My
nIndsTaxa="f:6" # Number of individuals per taxa (actual num of ind will be 3)
nTaxa="f:5" # Number of taxa
subsRate="e:10000000"

simphy -rs $stReplicates -rl $nGTs -sb ln:-15,1 -st $timeRange -sl $nTaxa -so f:1 -sp f:100000 -su $subsRate -si $nIndsTaxa -hh ln:1.2,1 -hl ln:1.4,1 -hg f:200 -v 1 -o $pipelinesName -cs $RANDOM -om 1 -od 1 -op 1 -oc 1 -on 1

# Move to simphyproject folder
cd $simphy_folder

perl $wrapper $simphy_folder $controlFile $RANDOM 1

for replicate in $(seq 1 $stReplicates); do
  st=$(printf "%0$numDigits""g" $replicate)
  echo "Running indelible... $simphy_folder/$st/"
  cd $simphy_folder/$st/
  indelible
done

# have a profile
python simphy.ngs.wrapper -l DEBUG
