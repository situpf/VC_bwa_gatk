#!/bin/bash
#SBATCH --job-name=x0
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=8G 
# #SBATCH -o log/slurm.x0.%j.out 
# #SBATCH -e log/slurm.x0.%j.err 
# #SBATCH --mail-type=END
# #SBATCH --mail-user=

### GATK best practices pipeline:
### https://software.broadinstitute.org/gatk/best-practices/bp_3step.php?case=GermShortWGS

export PATH=$PATH:/homes/users/mtormo/opt/ANALYSIS/VC_bwa_gatk/


usage() {
    NAME=$(basename $0)
    cat <<EOF
Usage:
    RUN AS: if [ ! -d log/ ];then mkdir log/;fi ; sbatch ${0} [1/2] [PATH_TO_FASTQ] [REFERENCE] [DBSNP] [Y/N] [BEDFILE]
    [1/2] -> 1:Trim and pair, 2:Join individual VCFs and get read stats
    [PATH_TO_FASTQ] -> Path with raw fastq
    [REFERENCE] -> reference fasta (change x1 pipeline whether the reference is SMALL -bwa index is- or BIG -bwa index bwtsw-)
    [DBSNP] -> VCF with known variants or -
    [Y/N] -> Mark PCR duplicates (Yes/No)
    [BEDFILE] -> Bed file with regions to analyse or -
EOF
}


if [ "$#" -ne 6 ]; then
    usage
    exit 1
fi


analysis=$1
path=$2
ref=$3
known=$4
dupli=$5
bed=$6

date


function part1 {
    ### load modules
    module load GATK/3.4-46-Java-1.8.0_92
    module load picard/2.18.6-Java-1.8.0_92
    module load BWA/0.7.16a-foss-2016b
    module load SAMtools/1.6-foss-2016b
    module load tabix/0.2.6-foss-2016b

    ### run GATK with 4G of memory and using one thread for collecting "garbage"
    RUN_GATK="java -Xmx4g -XX:+UseSerialGC -jar $EBROOTGATK/GenomeAnalysisTK.jar"

    if [ ! -d TRIMMED ];then mkdir TRIMMED;fi
    if [ ! -d BWA_RESULTS/ ];then mkdir BWA_RESULTS/;fi
    if [ ! -d VCF/ ];then mkdir VCF/;fi
    if [ ! -d VCF/single ];then mkdir VCF/single/;fi

    ### keep basename for dict (for compressed/uncompressed - now, gatk only works with uncompressed references)    
    REF_PATH=$(dirname "${ref}")
    FILENAME=$(basename "$ref")
    if [ "${ref##*.}" == "gz" ]; then 
        FILENAME="${FILENAME%.*}"
        REF_BASENAME=$REF_PATH"/""${FILENAME%.*}"
    else
        REF_BASENAME=$REF_PATH"/""${FILENAME%.*}"
    fi
    ### Index reference sequence. Choose between "is" or "bwtsw", for short or long references
    ##bwa index -a !!!bwtsw / is!!!  $REF
    if [ ! -e $ref.amb ]; then
        bwa index -a bwtsw $ref
    fi
    if [ ! -e $ref.fai ]; then
        samtools faidx $ref
    fi
    if [ ! -e $REF_BASENAME.dict ]; then
        java -Xmx2g -XX:+UseSerialGC -jar $EBROOTPICARD/picard.jar CreateSequenceDictionary R=$ref O=$REF_BASENAME.dict
    fi


    for i in  $path/*R1_001.fastq.gz; do
        ### don't analyze samples:
        if [[ $i == *Undetermined* ]] ;then echo $i; continue;fi

        FW=$i
        RV=`echo $i | sed 's/_R1_001.fastq.gz/_R2_001.fastq.gz/g'`
        ID=`echo $(basename "$i") | sed 's/_R1_001.fastq.gz//g'`
        echo $ID

        sbatch --cpus-per-task=2 --mem-per-cpu=8G -o log/slurm.x1.%j.out -e log/slurm.x1.%j.err x1_Trim_map_call.v1.0.sh $ref $FW $RV $ID $known $dupli $bed
    done
}




function part2 {
    ### load modules
    module load GATK/3.4-46-Java-1.8.0_92
    module load picard/2.18.6-Java-1.8.0_92
    module load tabix/0.2.6-foss-2016b
    module load VCFtools/0.1.14-foss-2016b-Perl-5.22.1
    module load Python/2.7.12-foss-2016b
    
    ### list gvcf files
    ls -1 VCF/single/*g.vcf | sort > gvcf2analyse.list

    ### run GATK with 4G of memory and using one thread for collecting "garbage"
    RUN_GATK="java -Xmx4g -XX:+UseSerialGC -jar $EBROOTGATK/GenomeAnalysisTK.jar"

    ### CALL VARIANTS
    $RUN_GATK -T GenotypeGVCFs -R $ref -V gvcf2analyse.list -o VCF/All_raw_variants.vcf --dbsnp $known $l_bed

    ### FILTER VARIANTS (Hard-Filtering)
    ### Filter variants (hard-filtering) - SNPs
    $RUN_GATK -T SelectVariants -R $ref -V VCF/All_raw_variants.vcf -selectType SNP -o VCF/All_raw_snps.vcf  $l_bed

    $RUN_GATK -T VariantFiltration -R $ref -V VCF/All_raw_snps.vcf  \
        --filterExpression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
        --filterName "snp_filter" \
        -o VCF/All_filt_snps.vcf $l_bed
        
    ### Filter variants (hard-filtering) - Indels
    $RUN_GATK -T SelectVariants -R $ref -V VCF/All_raw_variants.vcf -selectType INDEL -o VCF/All_raw_indels.vcf $l_bed 
    $RUN_GATK -T VariantFiltration -R $ref -V VCF/All_raw_indels.vcf \
        --filterExpression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0" \
        --filterName "indel_filter" \
        -o VCF/All_filt_indels.vcf $l_bed

        
    ### concatenate vcfs
    bgzip VCF/All_filt_snps.vcf
    tabix -p vcf VCF/All_filt_snps.vcf.gz 
    bgzip VCF/All_filt_indels.vcf
    tabix -p vcf VCF/All_filt_indels.vcf.gz 
    

    export PERL5LIB=$EBROOTVCFTOOLS/lib/perl5/site_perl/5.22.1/
    vcf-concat  VCF/All_filt_snps.vcf.gz VCF/All_filt_indels.vcf.gz | vcf-sort   > VCF/All_filt_variants_sorted.vcf
    $RUN_GATK -T VariantsToTable -R $ref -V VCF/All_filt_variants_sorted.vcf -o VCF/All_filt_variants_sorted_ext.csv -F CHROM -F POS -F REF -F ALT -F QD -F FS -F ReadPosRankSum -GF GT -GF GQ -GF DP -GF AD -AMD $l_bed

    ### get table with all variants and AF
    Get_variants_from_multivcf.py VCF/All_filt_variants_sorted.vcf > VCF/All_filt_variants_sorted_table.csv


    ### rm intermediate files
    rm -f VCF/All_raw_indels.vcf*
    rm -f VCF/All_raw_snps.vcf*
    rm -f VCF/All_filt_indels.vcf*
    rm -f VCF/All_filt_snps.vcf*


    ### get stats after/before mapping
    echo -e "ID\tINITIAL\tPP_AFTER_MAP\tPP_FOR_VARIANT_CALL" > reads_counts.csv
    for i in  $path/*R1_001.fastq.gz; do
        ### don't analyze samples
        if [[  $i == *Undetermined* ]] ;then echo $i; continue;fi
        ID=`echo $(basename "$i") | sed 's/_R1_001.fastq.gz//g'`
        reads_count.v1.0.sh $ID $path >> reads_counts.csv
    done

}


if [ $analysis == 1 ]; then
    part1
elif [ $analysis == 2 ]; then
    if [ "$6" == "-" ]; then
        l_bed=""
    else
        l_bed="-L $6"
    fi
    part2
fi

date
