#!/bin/bash 

set -ex

SAMPLE=$1
CELL=$2 
FASTQ1=$3 
FASTQ2=$4


if test $FASTQ1 -nt temp_stage1_${SAMPLE}.sam ; then
        ~/software/bwa/bwa mem -t 8 -M \
        -R "@RG\tID:${CELL}\tPL:ILLUMINA\tSM:${SAMPLE}\tLB:${SAMPLE}_${CELL}" \
        ~/course_data/resources/reference/human_g1k_v37.fasta $FASTQ1 $FASTQ2 \
        >temp_partial_stage1_${SAMPLE}.sam 2>log_${SAMPLE}_stage1_bwa.log
        mv temp_partial_stage1_${SAMPLE}.sam temp_stage1_${SAMPLE}.sam
fi

if test temp_stage1_${SAMPLE}.sam -nt temp_stage2_${SAMPLE}.bam ; then
        java -jar ~/software/picard/picard.jar FixMateInformation I=temp_stage1_${SAMPLE}.sam \
        O=temp_partial_stage2_${SAMPLE}.bam VALIDATION_STRINGENCY=SILENT \
        COMPRESSION_LEVEL=0 CREATE_INDEX=true SORT_ORDER=coordinate TMP_DIR=. \
        2>log_${SAMPLE}_stage2_FixMateInformation.log
        mv temp_partial_stage2_${SAMPLE}.bam temp_stage2_${SAMPLE}.bam
        mv temp_partial_stage2_${SAMPLE}.bai temp_stage2_${SAMPLE}.bai
fi

if test temp_stage2_${SAMPLE}.bam -nt temp_stage3_${SAMPLE}.bam ; then 

echo "identifying reads that are exact duplicated"

java -jar ~/software/picard/picard.jar MarkDuplicates I=temp_stage2_SAMP1.bam O=temp_stage3_${SAMPLE}.bam 
METRICS_FILE=outpu
t.duplicates REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=SILENT COMPRESSION_LEVEL=0 CREATE_INDEX=true TMP_DIR=.

echo "running realigner tool" 
java -jar ~/software/gatk/gatk-3.6.0/GenomeAnalysisTK.jar -T RealignerTargetCreator -I temp_stage3_${SAMPLE}.bam -o 
output.i
ntervals -known ~/course_data/resources/filtering_annotation/Mills_and_1000G_gold_standard.indels.b37.vcf -R 
~/course_data/r
esources/reference/human_g1k_v37.fasta

echo "running indel realigner"
java -jar ~/software/gatk/gatk-3.6.0/GenomeAnalysisTK.jar -T IndelRealigner -I temp_stage3_SAMP1.bam -o 
final_${SAMPLE}.bam 
-targetIntervals output.intervals -known 
~/course_data/resources/filtering_annotation/Mills_and_1000G_gold_standard.indels.b
37.vcf -R ~/course_data/resources/reference/human_g1k_v37.fasta

echo "running BQSR process"
java -jar ~/software/gatk/gatk-3.6.0/GenomeAnalysisTK.jar -T BaseRecalibrator -I final_${SAMPLE}.bam -R 
~/course_data/resour
ces/reference/human_g1k_v37.fasta -knownSites ~/course_data/resources/filtering_annotation/ExAC.r1.sites.vep.vcf.gz -o 
test_
output.grp

echo "calling variants"
java -jar ~/software/gatk/gatk-3.6.0/GenomeAnalysisTK.jar -T HaplotypeCaller -I final_${SAMPLE}.bam --BQSR 
test_output.grp -
o output.g.vcf --emitRefConfidence GVCF -L ~/course_data/resources/intervals/v501_all_covered_bases.interval_list -R 
~/cours
e_data/resources/reference/human_g1k_v37.fasta -stand_call_conf 10.0

echo "running genotypeGVCF tool"
java -jar ~/software/gatk/gatk-3.6.0/GenomeAnalysisTK.jar -T GenotypeGVCFs -R  
~/course_data/resources/reference/human_g1k_v
37.fasta -V output.g.vcf -o gvcf.vcf

echo "filtering common variants"
java -jar ~/software/gatk/gatk-3.6.0/GenomeAnalysisTK.jar -T SelectVariants -R 
~/course_data/resources/reference/human_g1k_v
37.fasta -V gvcf.vcf -o filter.vcf --discordance 
~/course_data/resources/filtering_annotation/common_no_known_medical_impact
_20160302-edited.vcf

echo "filtering artefacts"
java -jar ~/software/gatk/gatk-3.6.0/GenomeAnalysisTK.jar -T SelectVariants -R 
~/course_data/resources/reference/human_g1k_v
37.fasta -V filter.vcf -o second_filter.vcf --discordance 
~/course_data/resources/filtering_annotation/V501_19052014_common_
artefacts.vcf

echo "running Alamut batch"
~/software/alamut/alamut-batch-standalone-1.11/alamut-batch --in second_filter.vcf --ann annotated.txt --unann 
output_unanno
tated.txt --ssIntronicRange 2 --outputVCFInfo VariantType AC AF AN DP FS MQ QD --outputVCFGenotypeData GT AD DP GQ PL 
--outp
utVCFQuality --outputVCFFilter
  
fi
echo "pipeline complete"
