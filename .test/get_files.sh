#!/bin/bash

SCRIPT=`realpath $0`
SCRIPTPATH=`dirname $SCRIPT`
DATADIR="${SCRIPTPATH}/data_test"

HS37D5="http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz"
HG00096="http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG00096/exome_alignment/HG00096.chrom20.ILLUMINA.bwa.GBR.exome.20120522.bam"

wget "${HS37D5}" -O "${DATADIR}/hs37d5.fa.gz"
wget "${HS37D5}.fai" -O "${DATADIR}/hs37d5.fa.gz.fai"

wget "${HG00096}" -O "${DATADIR}/Test.bam"
wget "${HG00096}.bai" -O "${DATADIR}/Test.bam.bai"

