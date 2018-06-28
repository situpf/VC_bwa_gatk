#!/bin/env bash

ID=$1
BAM="BWA_RESULTS/"$ID
DIR=$2


### get initial lines from fastq
total_lines=`zcat $DIR/$ID\_R1_001.fastq.gz | wc -l`
initial=$(($total_lines*2/4))

### trimmed reads and properly paired
pair1=`grep properly $BAM.sort.bam.flag | cut -d " " -f 1`
pair1_perc=$((($pair1*100)/$initial))
### filtered reads and properly paired
pair2=`grep properly $BAM.*.recal.bam.flag | cut -d " " -f 1`
pair2_perc=$((($pair2*100)/$initial))


echo -e $ID"\t"$initial"\t"$pair1"("$pair1_perc"%)\t"$pair2"("$pair2_perc"%)"
