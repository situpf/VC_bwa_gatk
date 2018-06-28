#!/bin/bash

REF=$1
FW=$2
RV=$3
ID=$4
known=$5
dupli=$6
if [ "$7" == "-" ]; then
    l_bed=""
else
    l_bed="-L $7"
fi

date 

### define paths
TRIM="TRIMMED/"$ID
BAM="BWA_RESULTS/"$ID
VCF="VCF/single/"$ID




### load modules
module load BWA/0.7.16a-foss-2016b
module load SAMtools/1.6-foss-2016b
module load tabix/0.2.6-foss-2016b
module load FastQC/0.11.5-Java-1.8.0_74
module load Trimmomatic/0.36-Java-1.8.0_92
module load GATK/3.4-46-Java-1.8.0_92
module load picard/2.18.6-Java-1.8.0_92
module load R/3.4.2-foss-2016b

###################################################################
#### 1. TRIM SEQUENCES (FASTQ)
###################################################################

#### Get good quality sequences from fastq files obtained in NGS: in this case, we trim paired reads with a quality of 25 and a minimum length of 35
java -Xmx2g -jar $EBROOTTRIMMOMATIC/trimmomatic-0.36.jar PE -trimlog $TRIM-trimmo.log $FW $RV -baseout $TRIM.fq.gz ILLUMINACLIP:$EBROOTTRIMMOMATIC/adapters/NexteraPE-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:25 MINLEN:35


#### Check quality of trimmed sequences
fastqc --nogroup -o TRIMMED/ $FW &
fastqc --nogroup -o TRIMMED/ $RV &
fastqc --nogroup -o TRIMMED/ $TRIM\_1P.fq.gz &
fastqc --nogroup -o TRIMMED/ $TRIM\_2P.fq.gz &
wait


###################################################################
### 2. ALIGN TRIMMED READS (FASTQ -> SAM/BAM)
##################################################################

RUN_GATK="java -Xmx4g -XX:+UseSerialGC -jar $EBROOTGATK/GenomeAnalysisTK.jar"

if [ $dupli == "Y" ]; then
    #### Align with bwa mem with -M for picard compatibility
    bwa mem -t 4 $REF $TRIM\_1P.fq.gz $TRIM\_2P.fq.gz -M > $BAM.sam
else
    #### Align with bwa mem 
    bwa mem -t 4 $REF $TRIM\_1P.fq.gz $TRIM\_2P.fq.gz > $BAM.sam
fi

###################################################################
#### 3. GATK BEST PRACTICES: READS PRE-PROCESSING
###################################################################

### sam to bam
samtools view -@ 4 -hSb $BAM.sam -o $BAM.bam

### sort bam 
samtools sort -@ 4 $BAM.bam -o $BAM.sort.bam

### Add/replace read groups
java -Xmx4g -XX:+UseSerialGC -jar $EBROOTPICARD/picard.jar AddOrReplaceReadGroups I=$BAM.sort.bam O=$BAM.sort.RG.bam ID=$ID SM=$ID LB=1 PL=illumina PU=1 VALIDATION_STRINGENCY=SILENT

### Mark duplicates (NECESSARY???)
if [ $dupli == "Y" ]; then
    java -Xmx4g -XX:+UseSerialGC -jar $EBROOTPICARD/picard.jar MarkDuplicates I=$BAM.sort.RG.bam O=$BAM.sort.RG.dup.bam M=$BAM\_dup_metrics.txt
    bam_sort_RG=$BAM.sort.RG.dup
else
    bam_sort_RG=$BAM.sort.RG
fi
    
### base quality recalibration
samtools index -@ 4 $bam_sort_RG.bam
if [ $known == "-" ]; then
    $RUN_GATK -T BaseRecalibrator -R $REF -I $bam_sort_RG.bam -o $BAM\_recal_data.table $l_bed
    $RUN_GATK -T BaseRecalibrator -R $REF -I $bam_sort_RG.bam -BQSR $BAM\_recal_data.table -o $BAM\_post_recal_data.table  $l_bed
else
    $RUN_GATK -T BaseRecalibrator -R $REF -I $bam_sort_RG.bam -knownSites $known -o $BAM\_recal_data.table $l_bed
    $RUN_GATK -T BaseRecalibrator -R $REF -I $bam_sort_RG.bam -knownSites $known -BQSR $BAM\_recal_data.table -o $BAM\_post_recal_data.table $l_bed
fi
### generate before/after-recalibration plots
$RUN_GATK -T AnalyzeCovariates -R $REF -before $BAM\_recal_data.table -after $BAM\_post_recal_data.table -plots $BAM\_recalibration_plots.pdf $l_bed
$RUN_GATK -T PrintReads -R $REF -I $bam_sort_RG.bam -BQSR $BAM\_recal_data.table -o $bam_sort_RG.recal.bam $l_bed


### generate stats
samtools flagstat -@ 4 $BAM.sort.bam > $BAM.sort.bam.flag
samtools flagstat -@ 4 $bam_sort_RG.recal.bam > $bam_sort_RG.recal.bam.flag


###################################################################
#### 4. GATK BEST PRACTICES: VARIANT CALLING (GVCF))
###################################################################

if [ $known == "-" ]; then
    $RUN_GATK -T HaplotypeCaller -R $REF -I $bam_sort_RG.recal.bam --emitRefConfidence GVCF -o $VCF.g.vcf $l_bed
else
    $RUN_GATK -T HaplotypeCaller -R $REF -I $bam_sort_RG.recal.bam --emitRefConfidence GVCF --dbsnp $known -o $VCF.g.vcf $l_bed
fi 

###################################################################
#### 5. CLEAN DIRECTORY
###################################################################


rm -f $BAM.sam
rm -f $BAM.bam*
rm -f $BAM.sort.RG.bam
rm -f $bam_sort_RG.bam*
rm -f $BAM\_post_recal_data.table
rm -f $BAM\_recal_data.table

date
