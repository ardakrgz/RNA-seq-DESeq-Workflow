# RNA-seq-DESeq-Workflow
A workflow for analysis of paired-end RNA-seq data and differential gene expression analysis.

## Overview

**Pipeline: fastp → STAR → RSEM → DESeq2**


The pipeline uses **SLURM HPC** for alignment & quantification and **local R** for DESeq2.

1. QC and Trimming | fastp | HPC & local |
2. Build STAR index | STAR | HPC |
3. Read alignment | STAR | HPC |
4. Build RSEM reference | RSEM | HPC |
5. Quantification | RSEM | HPC |
6. Differential Expression | DESeq2 | Local R |

## Step 1: QC and Trimming with fastp

Locate to the directory containing the raw read files and run fastp with the code below: 

```
cd /path/to/raw/reads

fastp -i <forward_read_file>_R1.fastq.gz -I <reverse_read_file>_R2.fastq.gz -o <forward_read_trimmed>_R1_trimmed.fq.gz -O <reverse_read_trimmed>_R2_trimmed.fq.gz -P 10 --thread 5 -h <fastp_report>_fastp.html
```

## Step 2: Build STAR index (on HPC)

Assembly .fa file: GRCh38.primary_assembly.genome.fa

Annotation gtf file: gencode.v46.annotation.gtf

```
#!/bin/bash
#SBATCH --job-name=STAR_index
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --mem=100G
#SBATCH --partition=mid # CHECK AVAILABILITY
#SBATCH --time=10:00:00
#SBATCH --output=star_index.log
#SBATCH --mail-user=example@mail.com # CONFIGURE MAIL
#SBATCH --mail-type=ALL

module load anaconda
source activate <your_environment>

GENOME_FA=/path/to/genome/assembly/fa_file
GTF_FILE=/path/to/annotation/gtf_file
GENOME_DIR=/path/to/star/genome_dir

STAR \
  --runThreadN 16 \
  --runMode genomeGenerate \
  --genomeDir "$GENOME_DIR" \
  --genomeFastaFiles "$GENOME_FA" \
  --sjdbGTFfile "$GTF_FILE" \
  --sjdbOverhang 100

```

## Step 3: Read Alignment with STAR (on HPC)

```
#!/bin/bash

#SBATCH --job-name=star_cnt_batch
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --mem=100G
#SBATCH --partition=mid # CHECK AVAILABILITY
#SBATCH --time=20:00:00
#SBATCH --output=star_cnt_batch.log
#SBATCH --mail-user=example@mail.com # CONFIGURE MAIL
#SBATCH --mail-type=ALL

module load anaconda
source activate <your_environment>

GENOME_DIR=/path/to/star/genome_dir
GTF_FILE=/path/to/annotation/gtf_file
INPUT_DIR=/path/to/trimmed/read/files/directory
OUTPUT_DIR=/path/to/output/directory

mkdir -p "$OUTPUT_DIR"

for R1 in ${INPUT_DIR}/*_R1_trimmed.fq.gz; do
    SAMPLE=$(basename "$R1" _R1_trimmed.fq.gz)
    R2=${INPUT_DIR}/${SAMPLE}_R2_trimmed.fq.gz
    SAMPLE_OUT=${OUTPUT_DIR}/${SAMPLE}
    mkdir -p "$SAMPLE_OUT"

    STAR \
      --runThreadN 16 \
      --genomeDir "$GENOME_DIR" \
      --readFilesIn "$R1" "$R2" \
      --readFilesCommand zcat \
      --sjdbGTFfile "$GTF_FILE" \
      --quantMode TranscriptomeSAM \
      --outSAMtype BAM SortedByCoordinate \
      --outFileNamePrefix "${SAMPLE_OUT}/"
done
```


## Step 4: Build RSEM Reference (on HPC)

```
#!/bin/bash
#SBATCH --job-name=rsem_ref_build
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --mem=100G 
#SBATCH --partition=mid # CHECK AVAILABILITY
#SBATCH --time=20:00:00
#SBATCH --output=rsem_ref_build.log
#SBATCH --mail-user=example@mail.com # CONFIGURE MAIL
#SBATCH --mail-type=ALL

module load anaconda
source activate <your_environment>

GENOME_FA=/path/to/genome/assembly/fa_file
GTF_FILE=/path/to/annotation/gtf_file
RSEM_REF=/path/to/RSEM/rsem_hg38_gencode46
RSEM_BIN=/path/to/RSEM/software/RSEM-1.2.25 #if RSEM is not in $PATH

STAR_PATH=/path/to/.conda/envs/your_environment/bin

$RSEM_BIN/rsem-prepare-reference \
    --gtf "$GTF_FILE" \
    --star \
    --star-path "$STAR_PATH" \
    --num-threads 10 \
    "$GENOME_FA" \
    "$RSEM_REF"

```

## Step 5: Quantification with RSEM (on HPC)

```
#!/bin/bash
#SBATCH --job-name=rsem_MDA_OE
#SBATCH --output=rsem_calculate_MDA_OE.log
#SBATCH --partition=mid # CHECK AVAILABILITY
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --mem=100G
#SBATCH --mail-user=example@mail.com # CONFIGURE MAIL
#SBATCH --mail-type=ALL

module load anaconda
source activate <your_environment>

INPUT_DIR=/path/to/STAR/output/directory
RSEM_REF=/path/to/RSEM/rsem_hg38_gencode46
RSEM_BIN=/path/to/RSEM/software/RSEM-1.2.25 #if RSEM is not in $PATH
THREADS=16

#type sample names as they appear in the directory
SAMPLES=("<sample1>" "<sample2>" "<sample3>")

for SAMPLE in "${SAMPLES[@]}"; do
    SAMPLE_DIR="${INPUT_DIR}/${SAMPLE}"
    BAM_FILE="${SAMPLE_DIR}/Aligned.toTranscriptome.out.bam"
    OUT_PREFIX="${SAMPLE_DIR}/${SAMPLE}"
    LOG_FILE="${SAMPLE_DIR}/rsem_${SAMPLE}.log"

    "$RSEM_BIN/rsem-calculate-expression" \
        --bam \
        --paired-end \
        -p $THREADS \
        "$BAM_FILE" \
        "$RSEM_REF" \
        "$OUT_PREFIX" \
        > "$LOG_FILE" 2>&1
done
```

## Step 6: Differential Expression Analysis with DESeq2 on R (local)

