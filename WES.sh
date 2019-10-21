#!/bin/bash
set -o pipefail
## step1_filter_all
SOAPnuke1.5.6 filter -n 0.1 -q 0.5 -i -l 12 -Q 2 -5 1 -E 50 -G -A 0.3 -1 Cleandata/SAMPLE_R1.fq.gz -2 Cleandata/SAMPLE_R2.fq.gz -f AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -r AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA -M 2 -o process/SiHa -C SiHa_1.fq -D SiHa_2.fq && \
fqcheck -r process/SiHa/SiHa_1.fq.gz -c process/SiHa/SiHa_1.fqcheck && \
fqcheck -r process/SiHa/SiHa_2.fq.gz -c process/SiHa/SiHa_2.fqcheck && \
perl fqcheck_distribute.pl process/SiHa/SiHa_1.fqcheck process/SiHa/SiHa_2.fqcheck -o process/SiHa/SiHa. && \
mv process/SiHa/SiHa_1.fq.gz process/SiHa/SiHa_2.fq.gz process/SiHa/*.png process/clean_data/ && \
if [-e process/SiHa/Raw_SiHa_1.fq.gz ];then rm -rf process/SiHa/Raw_SiHa_1.fq.gz;fi && \
if [-e process/SiHa/Raw_SiHa_2.fq.gz ];then rm -rf process/SiHa/Raw_SiHa_2.fq.gz;fi 

## step2_bwa_mem
bwa mem -t 30 -M -R "@RG\tID:SiHa\tPL:ILLUMINA\tPU:SiHa_DNA\tLB:SiHa_DNA\tSM:SiHa_DNA\tCN:BGI" hg19.fasta result/clean_data/SiHa_1.fq.gz result/clean_data/SiHa_2.fq.gz | samtools view -Sb -o process/SiHa/SiHa.bam - && \
java -Xmx3g -Djava.io.tmpdir=process/java_tmp -XX:MaxPermSize=512m -XX:-UseGCOverheadLimit -jar picard.jar SortSam I=process/SiHa/SiHa.bam O=process/SiHa/SiHa.sort.bam SO=coordinate MAX_RECORDS_IN_RAM=1000000 

## step3_rmdup_SiHa_DNA
java -Xmx2G -Djava.io.tmpdir=process/java_tmp -XX:MaxPermSize=512m -XX:-UseGCOverheadLimit -jar picard.jar MarkDuplicates REMOVE_DUPLICATES=true I=process/SiHa/SiHa.sort.bam O=result/result_alignment/SiHa_DNA.bam METRICS_FILE=process/SiHa_DNA.bam.mat TMP_DIR=process/tmp && \
java -Xmx2G -jar picard.jar BuildBamIndex I=result/result_alignment/SiHa_DNA.bam O=result/result_alignment/SiHa_DNA.bam.bai

## step4_GATKreallnRecal
java -Xmx4g -Djava.io.tmpdir=process/java_tmp -jar GenomeAnalysisTK.jar -T RealignerTargetCreator -R hg19.fasta -I result/result_alignment/SiHa_DNA.bam -L shell/cvr_ex_region.sort.bed -known Database/hg19/gatk/Mills_and_1000G_gold_standard.indels.hg19.vcf -known Database/hg19/gatk/1000G_phase1.indels.hg19.vcf -o process/SiHa_DNA.intervals && \
java -Xmx4g -Djava.io.tmpdir=process/java_tmp -jar GenomeAnalysisTK.jar -T IndelRealigner -R hg19.fasta -I result/result_alignment/SiHa_DNA.bam -known Database/hg19/gatk/Mills_and_1000G_gold_standard.indels.hg19.vcf -known Database/hg19/gatk/1000G_phase1.indels.hg19.vcf -targetIntervals process/SiHa_DNA.intervals -o process/SiHa_DNA.realign.bam && \
java -Xmx4g -Djava.io.tmpdir=process/java_tmp -jar GenomeAnalysisTK.jar -nct 5 -T BaseRecalibrator -R hg19.fasta -I process/SiHa_DNA.realign.bam -knownSites Database/hg19/gatk/dbsnp_138.hg19.vcf -knownSites Database/hg19/gatk/Mills_and_1000G_gold_standard.indels.hg19.vcf -knownSites Database/hg19/gatk/1000G_phase1.indels.hg19.vcf -o process/SiHa_DNA.realign.recal.table && \
java -Xmx4g -Djava.io.tmpdir=process/java_tmp -jar GenomeAnalysisTK.jar -nct 5 -T PrintReads -R hg19.fasta -I process/SiHa_DNA.realign.bam -BQSR process/SiHa_DNA.realign.recal.table -o process/SiHa_DNA.realign.recal.bam 

## step5_callGVCF_GATK
java -Xmx5G -Djava.io.tmpdir=process/java_tmp -jar GenomeAnalysisTK.jar -T HaplotypeCaller -R hg19.fasta -I process/SiHa_DNA.realign.recal.bam -L shell/cvr_ex_region.sort.bed --emitRefConfidence GVCF --variant_index_type LINEAR --variant_index_parameter 128000 -o process/callGVCF_GATK/SiHa_DNA.g.vcf.gz && \
java -Xmx2G -Djava.io.tmpdir=process/java_tmp -jar GenomeAnalysisTK.jar -T GenotypeGVCFs -R hg19.fasta --variant process/callGVCF_GATK/SiHa_DNA.g.vcf.gz -o process/callGVCF_GATK/SiHa_DNA.vcf.gz -stand_call_conf 10 -allSites 

## step6_GATK_snp
export PERL5LIB="/pipeline/DNA_Human_WES/lib:$PERL5LIB"
java -Xmx3G -Djava.io.tmpdir=process/java_tmp -jar GenomeAnalysisTK.jar  -T SelectVariants -R hg19.fasta -V process/callGVCF_GATK/SiHa_DNA.vcf.gz -selectType SNP --excludeNonVariants -o process/snp_GATK/SiHa_DNA.raw.snp.vcf.gz && \
java -Xmx3G -Djava.io.tmpdir=process/java_tmp -jar GenomeAnalysisTK.jar -T VariantFiltration -R hg19.fasta -V process/snp_GATK/SiHa_DNA.raw.snp.vcf.gz \
--filterExpression "QD < 2.0 || FS > 60.0 || MQ <40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" --filterName "filter"  -o process/snp_GATK/SiHa_DNA.filtered_snp.vcf.gz && \
bcftools view -e 'ALT=="*"' -f "PASS" -o process/snp_GATK/anno/SiHa_DNA.filtered_snp.vcf.gz -O z process/snp_GATK/SiHa_DNA.filtered_snp.vcf.gz 

## step7_GATK_indel
export PERL5LIB="/pipeline/DNA_Human_WES/lib:$PERL5LIB"
export PATH="pipeline/DNA_Human_WES/bin:$PATH"
java -Xmx3G -Djava.io.tmpdir=java_tmp -jar GenomeAnalysisTK.jar -T VariantFiltration -R hg19.fasta -V process/indel_GATK/SiHa_DNA.raw.indel.vcf.gz \
--filterExpression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0" --filterName "filter"  -o process/indel_GATK/SiHa_DNA.filtered_indel.vcf.gz && \
bcftools view -f "PASS" -o process/indel_GATK/anno/SiHa_DNA.filtered_indel.vcf.gz -O z process/indel_GATK/SiHa_DNA.filtered_indel.vcf.gz 
