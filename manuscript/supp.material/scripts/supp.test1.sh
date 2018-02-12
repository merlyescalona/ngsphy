################################################################################kjkj
# PATHS
################################################################################
ngsphyPATH="$HOME/git/ngsphy/"
CURRENT_DIR="$(pwd)"
CASE_NAME="test1"
MYRANDOMNUM=50426717
################################################################################
# Data organization
################################################################################
echo "Creating test folder"
mkdir -p    ${CURRENT_DIR}/${CASE_NAME}/files/ \
            ${CURRENT_DIR}/${CASE_NAME}/output/ \
            $CURRENT_DIR/${CASE_NAME}/img

echo "Gathering all the data in a single folder"
cp ${ngsphyPATH}/data/settings/ngsphy.settings.supp.test1.txt ${CURRENT_DIR}/${CASE_NAME}/files/
cp ${ngsphyPATH}/data/indelible/control.supp.test1.txt ${CURRENT_DIR}/${CASE_NAME}/files/
cp ${ngsphyPATH}/data/trees/supp.test1.tree ${CURRENT_DIR}/${CASE_NAME}/files/
cp ${ngsphyPATH}/data/sequences/anchor.supp.test1.fasta.tar.gz ${CURRENT_DIR}/${CASE_NAME}/files/
echo "Unzipping sequence file"
cd ${CURRENT_DIR}/${CASE_NAME}/files/
tar -xzf ${CURRENT_DIR}/${CASE_NAME}/files/anchor.supp.test1.fasta.tar.gz
echo "Moving to the working directory"
cd ${CURRENT_DIR}/${CASE_NAME}
################################################################################
# 1. NGSphy - input mode 3
################################################################################
for item in $(seq 1 100); do
    echo "Running NGSphy replicate $item"
    { time ngsphy -s files/ngsphy.settings.supp.test1.txt &> ${CURRENT_DIR}/${CASE_NAME}/files/ngsphy.all.output  ;} 2>> ${CURRENT_DIR}/${CASE_NAME}/files/ngsphy.timings
done
################################################################################
# 2. raxml-ng
################################################################################
for item in $(find ${CURRENT_DIR}/${CASE_NAME}/output/ -mindepth 3 -maxdepth 3 -type d | grep alignments); do
    echo "Running raxml-ng for replicate $item"
    { time raxml-ng --all --msa $item/ngsphydata_1_TRUE.fasta --model JC --seed $MYRANDOMNUM &> ${CURRENT_DIR}/${CASE_NAME}/files/raxml-ng.all.output ;} 2>> ${CURRENT_DIR}/${CASE_NAME}/files/raxml-ng.timing
done
################################################################################
# Files for R scripts
################################################################################
find ${CURRENT_DIR}/${CASE_NAME}/output/ -mindepth 3 -type f | grep bestTree | sort > ${CURRENT_DIR}/${CASE_NAME}/files/best.tree.files
find ${CURRENT_DIR}/${CASE_NAME}/output/ -mindepth 3 -type f | grep ngsphy.tree  | sort > ${CURRENT_DIR}/${CASE_NAME}/files/rerooted.trees.files

cd ${CURRENT_DIR}/${CASE_NAME}/
<<RCODE
library(ape)
library(ggplot2)
library(phangorn)
library(reshape2)
###############################################################################
bestTreeFilelistName="files/best.tree.files"
rerootedTreeFilelistName="files/rerooted.trees.files"
originalTreeFilename="files/supp.test1.tree"
treefilelist=read.table(bestTreeFilelistName, colClasses="character")
rerootedtreelist=read.table(rerootedTreeFilelistName, colClasses="character")
colnames(treefilelist)="filename"
colnames(rerootedtreelist)="filename"

original=read.tree(originalTreeFilename)
treefilelist$filename=sort(treefilelist$filename)
treefilelist$RF=rep(0,nrow(treefilelist))
treefilelist$KF=rep(0,nrow(treefilelist))
treefilelist$AnchorDistToRef=rep(0,nrow(treefilelist))
treelist=c()
for(index in 1:length(treefilelist$filename)){
    tree=read.tree(treefilelist[index,]$filename)
    treefilelist[index, ]$RF=RF.dist(original,tree)
    treefilelist[index, ]$KF=KF.dist(original,tree)
}
png("img/rf.kf.AnchorDistToRef.png", height=400, width=800)
ggplot(treefilelist, aes(x=RF)) + geom_histogram(bins=2, fill="steelblue") + theme_classic() +
    theme(
        text = element_text(size=20),
        panel.grid.major.x=element_line(color="grey20",linetype="dotted", size=.3),
        panel.grid.major.y=element_line(color="grey20",linetype="dotted", size=.3)) +
    xlab("Robinson Foulds distance") +
    ylab("Count")

m=melt(treefilelist[,c(1,3,2,4)])
ggplot(m,aes(x=variable, y=value,fill=variable)) +
  theme_classic() +
  geom_violin() +
  geom_boxplot(width=0.1) +
  theme(legend.position="none",
        panel.grid.major.x=element_line(color="grey60",linetype="dotted", size=.3),
        panel.grid.major.y=element_line(color="grey60",linetype="dotted", size=.3),
        panel.grid.minor.y=element_line(color="grey60",linetype="dotted", size=.3),
        axis.line.x = element_line(color="grey60"),
        axis.line.y = element_line(color="grey60")
        ) +
  scale_x_discrete(labels=c("Branch Score","Robinson Foulds","P-distance (Anchor given vs. simulated)")) +
  xlab("Distances") + ylab("Values")
dev.off()
RCODE
