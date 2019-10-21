# Multi-omics cervical dataset
## Raw Data


## Instructions
This project is composed of for pipelines, including WES, RNA-seq, ATAC-seq and HiC-pro.

### Required Softwares
#### WES
SoapNuke1.5.6, bwa 0.7.15, Picard 2.5.0, GATK 3.6

#### RNA-seq
SOAPnuke 1.5.6, HISAT, RSEM, Bowtie2

#### ATAC-seq
ATAC-seq data were analyzed according to the pipeline developed by Kundaje’s Lab. (https://github.com/kundajelab/atac_dnase_pipelines)

Bowtie2, Picard, SAMtools, bedtools (2.16), phantompeakqualtools (v1.2), MACS2 and IDR were needed in this pipeline.

#### Hic-pro
HiC dataset were processed by standard HiC-Pro (v2.11.1) analysis pipeline. (https://github.com/nservant/HiC-Pro)

Bowtie2 and ICE were needed in this pipeline.
