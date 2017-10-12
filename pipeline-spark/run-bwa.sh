./spark-2.2.0-bin-hadoop2.7/bin/spark-submit --master local[*] --driver-memory 10g --class com.github.sparkbwa.SparkBWA SparkBWA-0.2.jar -m -r -k -p --index /var/data/hg19/hg19.fasta -w "-R '@RG\tID:Sample1\tLB:library1\tSM:Sample1\tPL:ILLUMINA' " /var/data/data/DRR002346_1.fastq /var/data/data/DRR002346_2.fastq /var/data/output/Output_DRR002346

