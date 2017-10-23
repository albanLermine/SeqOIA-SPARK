#!/bin/bash          

DATASET=$1
IP=$3
START_DATE=$(date +'%s')
NUMBER_OF_SLAVES=6
NUMBER_OF_CORES_BY_SLAVE=16
CORES=${2:-$NUMBER_OF_CORES_BY_SLAVE}

MASTER="spark://$IP:7077" # Spark Standalone
BASE=hdfs://$IP/user/ogirardot/seqoia
LOCAL_BASE=/user/ogirardot/seqoia
SPARK_CMD=spark-submit
REF_BASE=/data/seqoia/hg19/

echo "launching Spark BWA for $DATASET..."
time spark-submit \
	--master $MASTER \
	--driver-memory 30g \
	--executor-memory 20g \
	--num-executors $NUMBER_OF_SLAVES \
	--executor-cores 1 \
	--class com.github.sparkbwa.SparkBWA \
	SparkBWA-0.2.jar \
	-m -r -k -p \
	--index $REF_BASE/hg19.fasta \
	-w "-t $(($NUMBER_OF_CORES_BY_SLAVE-1)) -R @RG\tID:foo\tLB:bar\tPL:illumina\tPU:illumina\tSM:DRR002346" \
	$LOCAL_BASE/input/$DATASET/*_1.fastq \
	$LOCAL_BASE/input/$DATASET/*_2.fastq \
	$LOCAL_BASE/output/$DATASET/sam

echo "BWA MEM with Full SAM creation done."

echo "Starting SAM => Sorted BAM for $DATASET..."

time ./gatk-4.beta.6/gatk-launch SortReadFileSpark \
        --sparkMaster $MASTER \
        --input $BASE/output/$DATASET/sam/FullOutput.sam \
        --output $BASE/output/$DATASET/sorted.bam \
        -- --sparkRunner SPARK --driver-memory 10G --executor-memory 5G\
        --num-executors $NUMBER_OF_SLAVES\
        --executor-cores $CORES

echo "Finished sorted bam creation."

echo "Starting VCF creation with HaplotypeCallerSpark for $DATASET"

time ./gatk-4.beta.6/gatk-launch HaplotypeCallerSpark \
	--sparkMaster $MASTER \
	--input $BASE/output/$DATASET/sorted.bam \
	--output $BASE/output/$DATASET/output.vcf \
	--reference $REF_BASE/hg19.2bit -- --sparkRunner SPARK --driver-memory 10G --executor-memory 10G\
	--num-executors $NUMBER_OF_SLAVES\
	--executor-cores $CORES

echo "Finished HaplotypeCallerSpark."

"Whole Genome Sequencing took  $(($(date +'%s')-$START_DATE)) seconds"
