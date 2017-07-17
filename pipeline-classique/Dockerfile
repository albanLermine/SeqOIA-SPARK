#################################################################
# Dockerfile
#
# Docker version:   1
# Software:         seqoia-spark : BWA + GATK
# Software Version: 0.7.15 (BWA) + 4.beta.1 (GATK)
# Description:      BWA is a software package for mapping low-divergent sequences against a large reference genome, such as the human genome. 
#                   GATK: Variant Discovery in High-Throughput Sequencing Data
# Website:          http://bio-bwa.sourceforge.net/
#                   https://software.broadinstitute.org/gatk/
# URL:              https://sourceforge.net/projects/bio-bwa/files/bwa-0.7.15.tar.bz2/download
#                   https://github.com/broadinstitute/gatk
# Tags:             Genomics
# Build Cmd:        docker build --rm -t seqoia-spark:1.0.0 .
# Pull Cmd:         docker pull seqoia-spark:1.0.0
# Run Cmd:          docker run --rm -it seqoia-spark:1.0.0 <options>
#################################################################

FROM debian:jessie

RUN echo "deb http://httpredir.debian.org/debian jessie-backports main contrib non-free" >> /etc/apt/sources.list

RUN apt-get update && \
    apt-get install -y wget default-jre-headless unzip python samtools && \
    apt-get -t jessie-backports install -y bwa=0.7.15-2~bpo8+1 openjdk-8-jdk && \
    apt-get clean && \
    apt-get purge && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

RUN wget https://github.com/broadinstitute/gatk/releases/download/4.beta.1/gatk-4.beta.1.zip  && \
    unzip gatk-4.beta.1.zip && \
    rm gatk-4.beta.1.zip

RUN update-alternatives --set java /usr/lib/jvm/java-8-openjdk-amd64/jre/bin/java

# CHANGE USER
RUN useradd -ms /bin/bash docker
USER docker

