################################################################################kjkj
# PATHS
################################################################################
ngsphyPATH="$HOME/git/ngsphy/"
CURRENT_DIR="$(pwd)"
CASE_NAME="test2"
MYRANDOMNUM=50426717
GATK="$HOME/apps/gatk/3.8-0-ge9d806836/GenomeAnalysisTK.jar"
PICARD="$HOME/apps/picard/picard.jar"
referenceFile="$CURRENT_DIR/${CASE_NAME}/reference/reference.fasta"
coverages=( "2x" "10x" "20x" "100x" "200x")
################################################################################
# Data organization
################################################################################
echo "Creating test folder"
mkdir -p    ${CURRENT_DIR}/${CASE_NAME}/files/ \
            ${CURRENT_DIR}/${CASE_NAME}/output/ ${CURRENT_DIR}/${CASE_NAME}/src/ \
            $CURRENT_DIR/${CASE_NAME}/reference $CURRENT_DIR/${CASE_NAME}/img

echo "Gathering all the data in a single folder"
cp ${ngsphyPATH}/data/settings/ngsphy.settings.supp.test2.2x.txt ${CURRENT_DIR}/${CASE_NAME}/files/
cp ${ngsphyPATH}/data/settings/ngsphy.settings.supp.test2.10x.txt ${CURRENT_DIR}/${CASE_NAME}/files/
cp ${ngsphyPATH}/data/settings/ngsphy.settings.supp.test2.20x.txt ${CURRENT_DIR}/${CASE_NAME}/files/
cp ${ngsphyPATH}/data/settings/ngsphy.settings.supp.test2.100x.txt ${CURRENT_DIR}/${CASE_NAME}/files/
cp ${ngsphyPATH}/data/settings/ngsphy.settings.supp.test2.200x.txt ${CURRENT_DIR}/${CASE_NAME}/files/
cp ${ngsphyPATH}/data/settings/ngsphy.settings.supp.test2.200x.rc.txt ${CURRENT_DIR}/${CASE_NAME}/files/
cp ${ngsphyPATH}/data/indelible/control.supp.test2.txt ${CURRENT_DIR}/${CASE_NAME}/files/
cp ${ngsphyPATH}/data/trees/supp.test2.tree ${CURRENT_DIR}/${CASE_NAME}/files/
cp ${ngsphyPATH}/data/reference_alleles/my_reference_allele_file.test2.txt ${CURRENT_DIR}/${CASE_NAME}/files/
echo "Moving to the working directory"
cd ${CURRENT_DIR}/${CASE_NAME}
################################################################################
# 1. NGSphy read counts - For true Variants
################################################################################
ngsphy -s files/ngsphy.settings.supp.test2.200x.rc.txt
################################################################################
# 2. Reference selection
################################################################################
cat ${CURRENT_DIR}/${CASE_NAME}/NGSphy_test2_200x_RC/alignments/1/ngsphydata_1_TRUE.fasta | grep -a1 "1_0_0" | tail -2 | tr -d " "  > ${CURRENT_DIR}/${CASE_NAME}/reference/reference.fasta
cp ${CURRENT_DIR}/${CASE_NAME}/NGSphy_test2_200x_RC/reads/no_error/REPLICATE_1/ngsphydata_1_1_NOERROR.vcf ${CURRENT_DIR}/${CASE_NAME}/files/true.vcf
################################################################################
# 3. Running NGSphy
################################################################################
echo "Running NGSphy - 100 replicates - Coverage 2x"
for replicate in $(seq 1 100); do { time ngsphy -s ${CURRENT_DIR}/${CASE_NAME}/files/ngsphy.settings.supp.test2.2x.txt &> ${CURRENT_DIR}/${CASE_NAME}/files/ngsphy.2x.output; } 2>> ${CURRENT_DIR}/${CASE_NAME}/files/ngsphy.2x.timings; done
echo "Running NGSphy - 100 replicates - Coverage 10x"
for replicate in $(seq 1 100); do { time ngsphy -s ${CURRENT_DIR}/${CASE_NAME}/files/ngsphy.settings.supp.test2.10x.txt &> ${CURRENT_DIR}/${CASE_NAME}/files/ngsphy.10x.output; } 2>> ${CURRENT_DIR}/${CASE_NAME}/files/ngsphy.10x.timings; done
echo "Running NGSphy - 100 replicates - Coverage 20x"
for replicate in $(seq 1 100); do { time ngsphy -s ${CURRENT_DIR}/${CASE_NAME}/files/ngsphy.settings.supp.test2.20x.txt &> ${CURRENT_DIR}/${CASE_NAME}/files/ngsphy.20x.output; } 2>> ${CURRENT_DIR}/${CASE_NAME}/files/ngsphy.20x.timings; done
echo "Running NGSphy - 100 replicates - Coverage 100x"
for replicate in $(seq 1 100); do { time ngsphy -s ${CURRENT_DIR}/${CASE_NAME}/files/ngsphy.settings.supp.test2.100x.txt &> ${CURRENT_DIR}/${CASE_NAME}/files/ngsphy.100x.output; } 2>> ${CURRENT_DIR}/${CASE_NAME}/files/ngsphy.100x.timings; done
echo "Running NGSphy - 100 replicates - Coverage 200x"
for replicate in $(seq 1 100); do { time ngsphy -s ${CURRENT_DIR}/${CASE_NAME}/files/ngsphy.settings.supp.test2.200x.txt &> ${CURRENT_DIR}/${CASE_NAME}/files/ngsphy.200x.output; } 2>> ${CURRENT_DIR}/${CASE_NAME}/files/ngsphy.200x.timings; done
################################################################################
# 4. Indexing reference
################################################################################
bwa index $referenceFile
samtools faidx $referenceFile
java -jar -Xmx4G $PICARD CreateSequenceDictionary REFERENCE=$referenceFile OUTPUT="$CURRENT_DIR/${CASE_NAME}/reference/reference.dict"
################################################################################
# 5. Mapping
################################################################################
# Organizational purposes
for coverageLevel in ${coverages[*]}; do
    mkdir -p $CURRENT_DIR/${CASE_NAME}/mappings/$coverageLevel
done
for ngsphyoutput in $(find ${CURRENT_DIR}/${CASE_NAME}/output -mindepth 1 -maxdepth 1 -type d); do
    coverageFolder=$(basename ${ngsphyoutput})
    for ngsphyreplicate in $(ls ${ngsphyoutput}| sort); do
        numInds=$(cat ${ngsphyoutput}/${ngsphyreplicate}/ind_labels/${SIMPHY_PROJECT_NAME}.1.individuals.csv | wc -l)
        let numInds=numInds-2 # This file has a header
        mkdir -p "$CURRENT_DIR/${CASE_NAME}/mappings/${coverageFolder}/${ngsphyreplicate}/"
        for ind in $(seq 0 $numInds); do
            echo "$ngsphyreplicate/$ind"
            infile="${ngsphyoutput}/$ngsphyreplicate/reads/REPLICATE_1/LOCUS_1/${SIMPHY_PROJECT_NAME}_1_1_data_${ind}_"
            outfile="$CURRENT_DIR/${CASE_NAME}/mappings/${coverageFolder}/${ngsphyreplicate}/${ngsphyreplicate}_${ind}.sam"
            RGID="${ngsphyreplicate}-I${ind}"
            machine="HiSeq2500"
            echo "bwa mem -M -t 4 -R \"@RG\tID:${RGID}\tSM:${RGID}\tPL:Illumina\tLB:${RGID}\tPU:${machine}\" ${referenceFile} ${infile}R1.fq ${infile}R2.fq > $outfile"  >> ${CURRENT_DIR}/${CASE_NAME}/src/mappings.sh
        done
    done
done
bash $CURRENT_DIR/${CASE_NAME}/src/mappings.sh
################################################################################
# 6. Sorting + bamming
################################################################################
for samFile in $(find ${CURRENT_DIR}/${CASE_NAME}/mappings -type f | grep sam$); do
    echo $samFile
    outputDIR=$(dirname $samFile)
    outputFILE="$(basename $samFile .sam).sorted.bam"
    echo "samtools view -bSh $samFile | samtools sort - -f $outputDIR/${outputFILE} -@ 4" >> $CURRENT_DIR/${CASE_NAME}/src/bamming.sh
    echo "samtools index $outputDIR/$outputFILE" >> $CURRENT_DIR/${CASE_NAME}/src/bamming.sh
    echo "rm $samFile"  >> $CURRENT_DIR/${CASE_NAME}/src/bamming.sh
done
bash $CURRENT_DIR/${CASE_NAME}/src/bamming.sh

################################################################################
# 7. Mark Duplicates
################################################################################
summaryFile="$CURRENT_DIR/${CASE_NAME}/files/duplicates.summary.txt"
for bamFile in $(find ${CURRENT_DIR}/${CASE_NAME}/mappings -type f | grep sorted.bam$); do
    coverageFolder=$(basename $(dirname $(dirname $bamFile)))
    outputDIR=$(dirname $bamFile)
    values=($(basename $bamFile | tr "_" " " | tr "." " " ))
    indID=${values[-4]}
    repID=1
    if [[ ${#values} -eq 6 ]]; then
        repID=${values[-3]}
    fi
    dedupOutput="$outputDIR/$(basename $bamFile .sorted.bam).dedup.bam"
    metricsOutput="$outputDIR/$(basename $bamFile .sorted.bam).metrics.txt"
    histogramOutput="$outputDIR/$(basename $bamFile .sorted.bam).histogram.txt"
    echo "picard MarkDuplicates I=$bamFile O=$dedupOutput M=$metricsOutput"
    java -jar -Xmx4G $HOME/apps/picard/picard.jar MarkDuplicates INPUT=$bamFile OUTPUT=$dedupOutput METRICS_FILE=$metricsOutput

    header=$(head -7 $metricsOutput | tail -n+7)
    summaryInfo=$(head -8 $metricsOutput | tail -n+8)
    if [[ ! -f $summaryFile ]]; then
        echo -e "COVERAGE\tREPLICATE\tINDIVIDUAL_ID\tNUM_MAPPED_READS_SAMTOOLS\tNUM_RECORDS_SAMTOOLS\tNUM_N_READS\t$header" > $summaryFile
    fi
    numReads=$(samtools view -c $bamFile)
    numRecords=$(samtools view $bamFile | wc -l)
    numMappedReads=$(samtools view -F 0x4 $bamFile | cut -f 1 | sort | uniq | wc -l)
    echo -e "$coverageFolder\t$repID\t$indID\t$numMappedReads\t$numRecords\t$numReads\t$summaryInfo" >> $summaryFile
    tail -n+11 $metricsOutput > $histogramOutput
    samtools index $dedupOutput
done
################################################################################
# 8. INDEL REALIGNMENT
################################################################################
for bamFile in $(find ${CURRENT_DIR}/${CASE_NAME}/mappings -type f | grep dedup.bam$| tail -n+1700); do
    echo "$bamFile"
    outputDIR=$(dirname $bamFile)
    mkdir -p $outputDIR/target/
    targetOutput="$outputDIR/target/$(basename $bamFile .dedup.bam).target.intervals.list"
    realignedBam="$outputDIR/$(basename $bamFile .dedup.bam).realigned.bam"
    java -jar -Xmx4g $GATK \
        -T RealignerTargetCreator \
        -R $referenceFile \
        -I $bamFile \
        -o $targetOutput

    java -jar -Xmx4g $GATK \
        -T IndelRealigner \
        -R $referenceFile \
        -I $bamFile \
        -targetIntervals $targetOutput \
        -o $realignedBam
    samtools index $realignedBam
done

################################################################################
# 9. GATK - single call joint genotyping
################################################################################
for bamFile in $(find ${CURRENT_DIR}/${CASE_NAME}/mappings -type f | grep realigned.bam$); do
    echo "$bamFile"
    outputDIR=$(dirname $bamFile)
    mkdir -p $outputDIR/vcf-singlevc-joint-gt/
    OUTPUTVCF="$outputDIR/vcf-singlevc-joint-gt/$(basename $bamFile .realigned.bam).g.vcf"
    { time java -jar -Xmx4g $GATK \
        -T HaplotypeCaller \
        -R $referenceFile \
        -I $bamFile \
        -ERC GVCF \
        -o $OUTPUTVCF;  } 2>> ${CURRENT_DIR}/${CASE_NAME}/files/time.gatk.HaplotypeCaller.g.vcf.txt
done

for coverageLevel in ${coverages[*]}; do
    coverageFolder="${CURRENT_DIR}/${CASE_NAME}/mappings/$coverageLevel"
    for replicate in $(ls $coverageFolder); do
        individuals=""
        replicateFolder="${CURRENT_DIR}/${CASE_NAME}/mappings/$coverageLevel/$replicate"
        for indFile in $(find ${CURRENT_DIR}/${CASE_NAME}/mappings/$coverageLevel/$replicate -type f | grep .g.vcf$); do
            individuals+=" -V $indFile"
        done
        OUTPUTVCF="$replicateFolder/vcf-singlevc-joint-gt/$replicate.vcf"
        { time java -jar -Xmx4g $GATK \
            -T GenotypeGVCFs \
            -R $referenceFile \
            $individuals \
            -o $OUTPUTVCF ;} 2>> ${CURRENT_DIR}/${CASE_NAME}/files/time.gatk.genotypeGVCF.txt
    done
done
################################################################################
# 10 - Count discovered variants
################################################################################
mkdir ${CURRENT_DIR}/${CASE_NAME}/varsites/
numVariantsSummary="${CURRENT_DIR}/${CASE_NAME}/files/numvariants.summary.txt"
echo -e "COVERAGE\tREPLICATE\tNUM_VARIANTS" >  ${CURRENT_DIR}/${CASE_NAME}/files/numvariants.summary.txt
for coverageLevel in ${coverages[*]}; do
    for vcffile in $(find ${CURRENT_DIR}/${CASE_NAME}/mappings/$coverageLevel -name "*.vcf"| grep -v g.vcf | grep vcf-singlevc-joint-gt); do
        base=$(basename $vcffile)
        repID=$(echo $base |tr "_" " " | tr "." " " | awk '{print $4}' )
        if [[ repID -eq "vcf" ]]; then
            repID=1
        fi
        numVariants=$(cat $vcffile | grep -v "^#" |wc -l)
        mkdir -p ${CURRENT_DIR}/${CASE_NAME}/varsites/$coverageLevel/
        cat $vcffile | grep -v "^#" | awk '{print $2}'  >  ${CURRENT_DIR}/${CASE_NAME}/varsites/$coverageLevel/${base}.varsites
        echo -e "$coverageLevel\t$repID\t$numVariants" >> $numVariantsSummary
    done
done
################################################################################
# 11. get information per coverage on the varibale sites
################################################################################
for coverageLevel in ${coverages[*]}; do
    find ${CURRENT_DIR}/${CASE_NAME}/varsites/$coverageLevel -name "*.varsites" > ${CURRENT_DIR}/${CASE_NAME}/files/varsites.$coverageLevel.files
done
