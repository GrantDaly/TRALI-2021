#!/bin/bash
set -e
source /opt/sentieonStartup.sh
sentieon driver --version


echo "***************Aligning " $FASTQ1 $FASTQ2
(sentieon bwa mem -M -R "@RG\tID:"$ID"\tSM:"$SAMPLE"\tPL:ILLUMINA" \
  -t $NUMBER_THREADS /mnt/data/reference/$REFERENCE  /mnt/data/fastqs/$SAMPLE"_"$ID".trimmed.R1.fastq.gz" /mnt/data/fastqs/$SAMPLE"_"$ID".trimmed.R2.fastq.gz" || echo -n 'error during alignment') \
  | sentieon util sort -r /mnt/data/reference/$REFERENCE -o /mnt/data/output/$SAMPLE"_"$ID".sorted.bam" -t $NUMBER_THREADS --sam2bam -i -
