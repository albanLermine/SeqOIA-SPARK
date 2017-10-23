#################################################################
# Dockerfile
#
# Docker version:   1
# Software:         seqoia-spark : SPARK + SPARK BWA + GATK
# Software Version: Spark (2.2.0) + SparkBWA (0.2) + 4.beta.6 (GATK)
# Description:      BWA is a software package for mapping low-divergent sequences against a large reference genome, such as the human genome. 
#                   GATK: Variant Discovery in High-Throughput Sequencing Data
# Website:          https://software.broadinstitute.org/gatk/
# URL:              https://github.com/broadinstitute/gatk
# Tags:             Genomics
# Build Cmd:        docker build --rm -t seqoia-spark:1.0.0 .
# Pull Cmd:         docker pull seqoia-spark:1.0.0
# Run Cmd:          docker run --rm -it seqoia-spark:1.0.0 <options>
#################################################################

FROM openjdk:8

# File Author / Maintainer
LABEL authors="Olivier GIRARDOT <o.girardot@lateral-thoughts.com>"

################## BEGIN INSTALLATION ######################
RUN wget -q -O gatk.zip https://github.com/broadinstitute/gatk/releases/download/4.beta.6/gatk-4.beta.6.zip && \
    unzip gatk.zip && \
    rm gatk.zip

RUN wget -q -O spark.tgz http://apache.mirrors.ovh.net/ftp.apache.org/dist/spark/spark-2.2.0/spark-2.2.0-bin-hadoop2.7.tgz && \
    tar -xzvf spark.tgz && \
    rm spark.tgz

ADD SparkBWA-0.2.jar
ADD run-*.sh

# CHANGE USER
RUN useradd -ms /bin/bash docker
#USER docker

##################### INSTALLATION END #####################

