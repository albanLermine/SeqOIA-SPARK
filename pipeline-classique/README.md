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

- Lancement de BWA sur les fichiers fastq.gz: 

```
bwa mem -R '@RG\tID:Sample1\tLB:library1\tSM:Sample1\tPL:ILLUMINA' /data/hg19.fasta /data/Sample1_R1.fastq.gz /data/Sample1_R2.fastq.gz > /data/Sample1.sam 
```

- Conversion SAM en BAM:

```
./gatk-4.beta.1/gatk-launch --javaOptions "-Xmx4G" SamFormatConverter -I /data/Sample1.sam -O /data/Sample1.bam
```

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

