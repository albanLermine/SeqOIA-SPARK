# Dockerfile du pipeline de base - Projet Seqoia-spar

## Commandes

- Création de l'image Docker (à réaliser 1 fois): 

```
docker build --rm --build-arg http_proxy=YOUR_PROXY:PORT/ --build-arg https_proxy=YOUR_PROXY:PORT/ -t seqoia-spark:1.0.0 .
```

- Lancement du Docker:

```
docker run --rm -v /your/data/path/:/data -it seqoia-spark:1.0.0 /bin/bash
```

- Création des index BWA (à réaliser 1 fois):

```
bwa index /data/hg19.fasta
```

- Création de l'index fasta (à réaliser 1 fois):

```
samtools faidx /data/hg19.fasta
```

- Création du fichier dictionary (à réaliser 1 fois):

```
./gatk-4.beta.1/gatk-launch --javaOptions "-Xmx4G" CreateSequenceDictionary -R /data/hg19.fasta -O /data/hg19.dict
```

- Lancement de BWA sur les fichiers fastq (de-zippé): 

```
spark-submit --driver-memory 60g --executor-memory 31g --num-executors 32 --class com.github.sparkbwa.SparkBWA SparkBWA-0.2.jar -m -r -k -p --index /data/seqoia/hg19/hg19.fasta -w "-R @RG\tID:foo\tLB:bar\tPL:illumina\tPU:illumina\tSM:DRR002346" data/DRR002346_1.fastq data/DRR002346_2.fastq output/DRR002346_SAM
```

- Création du fichier image bwa (à réaliser 1 fois):

```
./gatk-4.beta.5/gatk-launch BwaMemIndexImageCreator -I hg19.fasta
```

- Conversion SAM en BAM:

```
./gatk-4.beta.5/gatk-launch BwaSpark --sparkMaster spark://134.158.75.222:7077 --input output/DRR002346_SAM/FullOutput.sam --output output/DRR002346_BAM -R /data/seqoia/hg19/hg19.fasta -image /data/seqoia/hg19/hg19.fasta.img
```

END OF WIP

- Sort du fichier SAM

```
./gatk-4.beta.1/gatk-launch --javaOptions "-Xmx4G" SortSam -I /data/Sample1.bam -O /data/Sample1.sorted.bam -SO coordinate
```

- Index du fichier BAM:

```
samtools index /data/Sample1.sam
```

- Lancement de GATK 4 HaplotypeCaller non Spark:

```
./gatk-4.beta.1/gatk-launch --javaOptions "-Xmx4G" HaplotypeCaller -I /data/Sample1.sorted.bam -O /data/Sample1.vcf -R /data/hg19.fasta
```

