# RNA-seq-DESeq-Workflow
A workflow for analysis of paired-end RNA-seq data and differential gene expression analysis.

**Pipeline: fastp → STAR → RSEM → DESeq2**

## Pipeline Overview

The pipeline uses **SLURM HPC** for alignment & quantification and **local R** for DESeq2.

1. QC and Trimming            | fastp | HPC & local |
2. Build STAR index           | STAR | HPC |
3. Read alignment             | STAR | HPC |
4. Build RSEM reference       | RSEM | HPC |
5. Quantification             | RSEM | HPC |
6. Differential Expression    | DESeq2 + tximport | Local R |

## Step 1: QC and Trimming with fastp

fastp -i {forward_read_file}_R1.fastq.gz   -I {reverse_read_file}_R2.fastq.gz   -o {forward_read_trimmed}_R1_trimmed.fq.gz   -O {reverse_read_trimmed}_R2_trimmed.fq.gz   -P 10   --thread 5   -h {fastp_report}_fastp.html
