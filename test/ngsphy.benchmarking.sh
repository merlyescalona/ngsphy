################################################################################
# Paths
################################################################################
ngsphyPATH="$HOME/git/ngsphy/"
CURRENT_DIR="$(pwd)"
CASE_NAME="benchmarking2"
referenceFile="$CURRENT_DIR/${CASE_NAME}/reference/reference.fasta"
MYRANDOMNUM=50426717
################################################################################
# Organization (Folder generation and getting proper files)
################################################################################
echo "Creating test folder"
mkdir -p    ${CURRENT_DIR}/${CASE_NAME}/settings/ ${CURRENT_DIR}/${CASE_NAME}/files/ \
            ${CURRENT_DIR}/${CASE_NAME}/output/ ${CURRENT_DIR}/${CASE_NAME}/src/ \
            $CURRENT_DIR/${CASE_NAME}/reference $CURRENT_DIR/${CASE_NAME}/img \
            $CURRENT_DIR/${CASE_NAME}/msaWreference

echo "Gathering all the data in a single folder"
cp ${ngsphyPATH}/data/settings/ngsphy.settings.benchmarking.100x.txt ${CURRENT_DIR}/${CASE_NAME}/settings/
cp ${ngsphyPATH}/data/indelible/control.benchmarking.txt ${CURRENT_DIR}/${CASE_NAME}/files/
cp ${ngsphyPATH}/data/sequences/anchor.benchmarking.fasta  ${CURRENT_DIR}/${CASE_NAME}/files/
cp ${ngsphyPATH}/data/trees/benchmarking.tree  ${CURRENT_DIR}/${CASE_NAME}/files/original.tree
cp $CURRENT_DIR/${CASE_NAME}/files/anchor.benchmarking.fasta  $CURRENT_DIR/${CASE_NAME}/reference/reference.fasta
echo "Moving to the working directory"
cd ${CURRENT_DIR}/${CASE_NAME}
################################################################################
# 1. NGSphy
# need to rename the outgrup tip within the genetrees taking into account
# that NGSphy  does not allow the 0_0_0 label on any tip.
################################################################################
for index in $(seq 1 100); do
    echo -e "$index"
    ngsphySettings="${CURRENT_DIR}/${CASE_NAME}/settings/ngsphy.settings.benchmarking.100x.txt"
    { time ngsphy -s $ngsphySettings; } 2>> ${CURRENT_DIR}/${CASE_NAME}/files/time.ngsphy.txt
done
################################################################################
# Tree inference
################################################################################
# Corre 100 bootstraps auto
for index in $(seq 1 100); do
    echo $index
    if [[ $index -eq 1 ]]; then
        repFolder="${CURRENT_DIR}/${CASE_NAME}/output/BenchmarkGT"
    else
        repFolder="${CURRENT_DIR}/${CASE_NAME}/output/BenchmarkGT_$index"
    fi
    echo -e "$repFolder"
    { time raxml-ng --redo --all --msa $repFolder/alignments/1/ngsphydata_1.fasta \
         --model JC --tree pars{10} --seed $MYRANDOMNUM; } 2>> ${CURRENT_DIR}/${CASE_NAME}/files/time.all.raxml.ng.txt
done

################################################################################
# Checking timings:
################################################################################
cat ${CURRENT_DIR}/${CASE_NAME}/files/time.all.raxml.ng.txt | grep real | tr "m" " "| tr -d "s" | awk 'BEGIN{ m=0; s=0; ms=0; total=0}{ m=$2; s=$3; total+=((m*60)+s)}END{print (total/60), "mins"}'
cat ${CURRENT_DIR}/${CASE_NAME}/files/time.all.raxml.ng.txt | grep user | tr "m" " "| tr -d "s" | awk 'BEGIN{ m=0; s=0; ms=0; total=0}{ m=$2; s=$3; total+=((m*60)+s)}END{print (total/60), "mins"}'
cat ${CURRENT_DIR}/${CASE_NAME}/files/time.all.raxml.ng.txt | grep sys  | tr "m" " "| tr -d "s" | awk 'BEGIN{ m=0; s=0; ms=0; total=0}{ m=$2; s=$3; total+=((m*60)+s)}END{print (total/60), "mins"}'

find ${CURRENT_DIR}/${CASE_NAME}/output/ -type f -name "ngsphydata_1.fasta"  > ${CURRENT_DIR}/${CASE_NAME}/files/sequences.files.txt
find ${CURRENT_DIR}/${CASE_NAME}/output/ -type f -name "ngsphydata_1.fasta.raxml.bootstraps" > ${CURRENT_DIR}/${CASE_NAME}/files/bootstraps.files.txt
find ${CURRENT_DIR}/${CASE_NAME}/output/ -name "*.bestTree" > ${CURRENT_DIR}/${CASE_NAME}/files/best.tree.files
find ${CURRENT_DIR}/${CASE_NAME}/output/ -name "ngsphy.tree" > ${CURRENT_DIR}/${CASE_NAME}/files/rerooted.trees.files
for sequence in $(cat files/sequences.files.txt); do
    echo $sequence
    filename=$(basename $(dirname $(dirname $(dirname $sequence))))
    cat ${CURRENT_DIR}/${CASE_NAME}/reference/reference.fasta $sequence >${CURRENT_DIR}/${CASE_NAME}/msaWreference/${filename}.fasta
done
find ${CURRENT_DIR}/${CASE_NAME}/msaWreference -type f -name "*.fasta"  > ${CURRENT_DIR}/${CASE_NAME}/files/sequences.wr.files.txt

<<R_RF_DIST
library(ape)
library(ggplot2)
library(stringr)
library(ggtree)
library(phangorn)
###############################################################################
originalTreeFile="/home/merly/test/benchmarking2/files/original.tree"
bestTreeFilelistName="/home/merly/test/benchmarking2/files/best.tree.files"
rerootedTreeFilelistName="/home/merly/test/benchmarking2/files/rerooted.trees.files"
treefilelist=read.table(bestTreeFilelistName, colClasses="character")
rerootedtreelist=read.table(rerootedTreeFilelistName, colClasses="character")
colnames(treefilelist)="filename"
colnames(rerootedtreelist)="filename"

original=read.tree(originalTreeFile)
treefilelist$filename=sort(treefilelist$filename)
treefilelist$RF=rep(0,nrow(treefilelist))
treefilelist$KF=rep(0,nrow(treefilelist))
treefilelist$AnchorDistToRef=rep(0,nrow(treefilelist))
treelist=c()
for(index in 1:nrow(treefilelist)){
    tree=read.tree(treefilelist[index,]$filename)
    treefilelist[index, ]$RF=RF.dist(original,tree)
    treefilelist[index, ]$KF=KF.dist(original,tree)
}


png("img/kf.branchscores.original.inferred.trees.boxplot.png", height=800, width=800)
ggplot(treefilelist, aes(1, KF)) +
    geom_boxplot(varwidth=T) +
    geom_jitter(width=0.1) +
    theme_classic() +
    theme(
        text = element_text(size=20),
        panel.grid.major.y=element_line(color="grey10",linetype="dotted", size=.3),
        panel.grid.minor.y=element_line(color="grey60",linetype="dotted", size=.3)
    ) +
    xlab("") + ylab("Branch score distance (KF)")
dev.off()


###############################################################################

png("img/kf.branchscores.original.inferred.trees.png", height=400, width=800)
ggplot(treefilelist, aes(x=1:100,y=KF)) + geom_point() +
    geom_hline(data=treefilelist,yintercept=mean(treefilelist$KF)) +
    geom_hline(data=treefilelist,yintercept=min(treefilelist$KF), colour="red", linetype="dashed" ) +
    geom_hline(data=treefilelist,yintercept=max(treefilelist$KF), colour="red", linetype="dashed" ) +
    theme(text = element_text(size=20)) +
    ylab("Branch score distance (KF)") + xlab("Replicates")
dev.off()


png("img/rf.distances.original.inferred.trees.png", height=400, width=800)
ggplot(treefilelist, aes(x=RF)) + geom_histogram(bins=2, fill="steelblue") + theme_classic() +
    theme(
        text = element_text(size=20),
        panel.grid.major.x=element_line(color="grey20",linetype="dotted", size=.3),
        panel.grid.major.y=element_line(color="grey20",linetype="dotted", size=.3)) +
    xlab("Robinson Foulds distance") +
    ylab("Count")
dev.off()

###############################################################################
fileListName="/home/merly/test/benchmarking2/files/sequences.wr.files.txt"
filelist=sort(as.list(read.table(fileListName, colClasses = "character"))[[1]])
seqsize=1000
dd=c()
for( index in 1:length(filelist)){
    baseComposition=c()
    sequence=filelist[index]
    seqStrSplit=str_split(sequence,"\\/")[[1]]
    d=str_split(seqStrSplit[7], "_")[[1]]
    dna=read.dna(sequence, format="fasta",as.character=T)
    pdist=phangorn::dist.p(as.phyDat(dna))
    treefilelist[index,]$AnchorDistToRef=pdist[1][1]
}
png("img/pdist.anchor.to.tip.png", height=400, width=800)
ggplot(treefilelist, aes(x=1:100, y=AnchorDistToRef)) + geom_line() + theme_classic() +
    theme(
        text = element_text(size=20),
        panel.grid.major.x=element_line(color="grey20",linetype="dotted", size=.3),
        panel.grid.major.y=element_line(color="grey20",linetype="dotted", size=.3)) +
    ylab("P-Distance") +
    xlab("Replicates") +
    ggtitle("Pairwise Polymorphism (P-Distances)")+
    labs(subtitle="From reference (anchor) sequence to corresponding anchor tip (simulated) sequence")
dev.off()

###############################################################################
png("img/tree.original.png", height=400, width=800)
ggtree(original) + theme_classic() + geom_tiplab(size=4) +
    geom_point(aes(shape=isTip)) +
    theme(legend.position="none") +
    theme(axis.line.x=element_line(), axis.line.y=element_line()) +
    theme(
        text = element_text(size=20),
        panel.grid.major.x=element_line(color="grey20", linetype="dotted", size=.3),
        panel.grid.major.y=element_blank()) +
    ggtitle("Original tree")
dev.off()
# ------------------------------------------------------------------------------
tree=read.tree(rerootedtreelist$filename[46])
png("img/tree.rerooted.rep1.png", height=400, width=800)
ggtree(tree) + theme_classic() + geom_tiplab(size=4) +
    geom_point(aes(shape=isTip)) +
    theme(legend.position="none") +
    theme(
        text = element_text(size=20),
        axis.line.x=element_line(),
        axis.line.y=element_line(),
        panel.grid.major.x=element_line(color="grey20", linetype="dotted", size=.3),
        panel.grid.major.y=element_blank())+
    ggtitle("Re-rooted tree (Replicate 1)")
dev.off()
R_RF_DIST
