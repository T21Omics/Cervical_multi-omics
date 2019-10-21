## We used pipeline developed by Kundaje’s Lab.
CODE=/software/atac_dnase_pipelines
DATA=ATAC-SEQ_DATA_PATH
WORK_ROOT=/cervical; SUFFIX=ATAC; WORK=$WORK_ROOT/$SUFFIX
bds $CODE/atac.bds \
-species hg19 -pe \
-nth 32 \
-enable_idr \
-out_dir ./ATAC_out \
-fastq1_1 $SAMPLE1_R1.fq.gz \
-fastq1_2 $SAMPLE1_R2.fq.gz \
-fastq2_1 $SAMPLE2_R1.fq.gz \
-fastq2_2 $SAMPLE2_R2.fq.gz 
