# RNA-seq-DESeq-Workflow
A workflow for analysis of paired-end RNA-seq data and differential gene expression analysis.

## Requirements

* **fastp**

    * GitHub Repo: https://github.com/OpenGene/fastp$0

    * Anaconda Download: https://anaconda.org/bioconda/fastp$0

* **STAR**

    * GitHub Repo: https://github.com/alexdobin/STAR/tree/master$0

    * Anaconda Download: https://anaconda.org/bioconda/star$0

* **RSEM**

    * GitHub Repo: https://github.com/deweylab/RSEM$0

* **DESeq2** (R)

    * https://bioconductor.org/packages/release/bioc/html/DESeq2.html$0

* **Assembly FASTA and annotation GTF files:** 

    * https://www.gencodegenes.org/human/$0

-or-

Download via wget: 
```
# The FASTA genome file (release 49 up-to-date as of Oct 2025)
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/GRCh38.primary_assembly.genome.fa.gz

# The GTF annotation file (release 49 up-to-date as of Oct 2025)
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.annotation.gtf.gz
```


## Overview


The pipeline uses **SLURM HPC** for alignment & quantification and **local R** for differential expression analysis and further visualization.

1. QC and Trimming with fastp on HPC or locally.
2. Build STAR index with STAR on HPC.
3. Read alignment with STAR on HPC.
4. Build RSEM reference with RSEM on HPC.
5. Quantification with RSEM on HPC.
6. Differential Expression with DESeq2 on local R.

## Step 1: QC and Trimming with fastp on HPC or locally

Locate to the directory containing the raw read files and run fastp with the code below: 


```
cd /path/to/raw/reads

fastp -i <forward_read>_R1.fastq.gz -I <reverse_read>_R2.fastq.gz \
    -o <fwd_out>_R1_trimmed.fq.gz -O <rev_out>_R2_trimmed.fq.gz \
    -P 10 --thread 5 -h <fastp_report>_fastp.html
```
At the end, trimmed.fq.gz files will be generated, which will be used in **step 3.** Additionally, an **html report** will be generated, which contains general information of the reads, before and after statistics, and filtering results.

## Additional: Running jobs on SLURM
Job scripts that are going to be run on HPC can be prepared and run as the following:

```
touch <example_script>.sh #create a blank text file

nano <example_script>.sh #open the blank file, copy-paste and modify the script

sbatch <example_script>.sh #submit the script to the cluster
```

If configured right, the script will be queued in the cluster, and information mail(s) will be sent to the configured mail adress.

## Step 2: Build STAR index (on HPC)

A genome index is needed to be generated (only once) so that STAR can utilize it during the alignment/mapping (next) step. A reference genome sequence file (FASTA) and an annotation file (GTF) are required installed.
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

### CONFIGURE THESE LINES ###

module load anaconda
source activate <your_environment>

GENOME_FA=/path/to/genome/assembly/fa_file
GTF_FILE=/path/to/annotation/gtf_file
GENOME_DIR=/path/to/star/genome_dir

### CONFIGURE THESE LINES ###

STAR \
  --runThreadN 16 \
  --runMode genomeGenerate \
  --genomeDir "$GENOME_DIR" \
  --genomeFastaFiles "$GENOME_FA" \
  --sjdbGTFfile "$GTF_FILE" \
  --sjdbOverhang 100 # This value must be set to ReadLength - 1
```
A star genome directory (~30 GB) will be generated, which contains necessary information of the genome and annotations. It can be used for other analyses with the same genome reference and annotation.

By following the **additional step** above, one can run this STAR step as the following:

```
touch star_index.sh #create the script

nano star_index.sh #a blank file will open, copy the script above and configure 

sbatch star_index.sh #submit the script
```


## Step 3: Read Alignment with STAR (on HPC)
In this step, the trimmed read files are aligned to the reference genome generated in the second step.

STAR genome directory, the annotation file, and the input directory containing the trimmed read files need to be specified below.
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

### CONFIGURE THESE LINES ###

module load anaconda
source activate <your_environment> 

GENOME_DIR=/path/to/star/genome_dir
GTF_FILE=/path/to/annotation/gtf_file
INPUT_DIR=/path/to/trimmed/read/files/directory
OUTPUT_DIR=/path/to/output/directory 

### CONFIGURE THESE LINES ###

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
Similar to the STAR, RSEM also needs a reference, which can be generated with the code below.

The exact assembly genome file (FASTA) and annotation file (GTF) should be used and specified below.

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

### CONFIGURE THESE LINES ###

module load anaconda
source activate <your_environment>

GENOME_FA=/path/to/genome/assembly/fa_file
GTF_FILE=/path/to/annotation/gtf_file
RSEM_REF=/path/to/RSEM/rsem_hg38_gencode46
RSEM_BIN=/path/to/RSEM/software/RSEM-1.2.25 #if RSEM is not in $PATH

STAR_PATH=/path/to/.conda/envs/your_environment/bin # useful in case RSEM can't locate STAR

### CONFIGURE THESE LINES ###

$RSEM_BIN/rsem-prepare-reference \
    --gtf "$GTF_FILE" \
    --star \
    --star-path "$STAR_PATH" \
    --num-threads 10 \
    "$GENOME_FA" \
    "$RSEM_REF"

```

## Step 5: Quantification with RSEM (on HPC)
The directory containing STAR output files, the generated RSEM reference directory, and RSEM bin directory containing the RSEM scripts (under software directory) should be specified. Sample names/types should also be specified.


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

### CONFIGURE THESE LINES ###

module load anaconda
source activate <your_environment>

INPUT_DIR=/path/to/STAR/output/directory
RSEM_REF=/path/to/RSEM/rsem_hg38_gencode46
RSEM_BIN=/path/to/RSEM/software/RSEM-1.2.25 #if RSEM is not in $PATH
THREADS=16

#type sample names as they appear in the directory
SAMPLES=("<sample1>" "<sample2>" "<sample3>")

### CONFIGURE THESE LINES ###

for SAMPLE in "${SAMPLES[@]}"; do
    SAMPLE_DIR="${INPUT_DIR}/${SAMPLE}"
    BAM_FILE="${SAMPLE_DIR}/Aligned.toTranscriptome.out.bam"
    OUT_PREFIX="${SAMPLE_DIR}/${SAMPLE}"
    LOG_FILE="${SAMPLE_DIR}/rsem_${SAMPLE}.log"

    "$RSEM_BIN/rsem-calculate-expression" \
        --bam \
        --paired-end \ ### for paired-end reads
        -p $THREADS \
        "$BAM_FILE" \
        "$RSEM_REF" \
        "$OUT_PREFIX" \
        > "$LOG_FILE" 2>&1
done
```

RSEM will generate result files with *.genes.results extension, which are needed for differential expression analysis in R with the code below.


## Step 6: Differential Expression Analysis with DESeq2 on R (local)

```
library(tximport)
library(DESeq2)
library(org.Hs.eg.db)

### -------- Importing RSEM results --------

dir <- "/path/to/rsem/results" # the directory containing .genes.results files (RSEM result files)
samples <- c("control1", "control2", "control3", ### CONFIGURE HERE
             "sample1", "sample2", "sample3")

files <- file.path(dir, paste0(samples, ".genes.results"))
names(files) <- samples

txi.rsem <- tximport(files, type = "rsem", txIn = FALSE, txOut = FALSE)

### -------- Metadata --------

sample_info <- data.frame(
  row.names = samples,
  condition = c("CNT", "CNT", "CNT", "SAMPLE", "SAMPLE", "SAMPLE") # specify conditions here
)

### -------- DESeq2 Analysis --------
dds <- DESeqDataSetFromMatrix(
  countData = round(txi.rsem$counts),
  colData = sample_info,
  design = ~ condition
)

dds <- DESeq(dds)
results <- results(dds)

### -------- Gene Annotation --------
rownames(results) <- gsub("\..*", "", rownames(results))
results$gene_symbol <- mapIds(org.Hs.eg.db,
                              keys = rownames(results),
                              column = "SYMBOL",
                              keytype = "ENSEMBL",
                              multiVals = "first")


### -------- Save All Results --------
results_cleaned <- na.omit(results)
write.csv(results_cleaned, "deseq_all_results.csv")


### -------- Save Significant Results --------
sig_results <- results_cleaned[
  !is.na(results_cleaned$padj) & results_cleaned$padj < 0.05 &
    abs(results_cleaned$log2FoldChange) > 1, ]

write.csv(sig_results, file.path(dir, "deseq_sig_results.csv"))


### -------- Visualization (Optional) --------
library(pheatmap)

vsd <- vst(dds, blind = FALSE)

# Order the (cleaned) results and get the ENSEMBL IDs of the top 30 genes based on padj values
top_genes_ids <- rownames(head(results_cleaned[order(results_cleaned$padj), ], 30))

# Subset the transformed count matrix to get the data for these top 30 genes
vsd_matrix <- assay(vsd)
rownames(vsd_matrix) <- gsub("\\..*", "", rownames(vsd_matrix))
mat <- vsd_matrix[top_genes_ids, ]

# Center the data by subtracting the row-wise means
mat <- mat - rowMeans(mat)

# Get the gene symbols for the matrix row names
# Since we are using 'results_cleaned', we are guaranteed to find a symbol for every gene.
final_labels <- results_cleaned[rownames(mat), "gene_symbol"]
rownames(mat) <- make.unique(final_labels)


# Create the heatmap PDF
pdf(file.path(dir, "Heatmap_DESeq.pdf"))
pheatmap(
  mat,
  annotation_col = NULL, #sample_info[, "Expression Levels", drop = FALSE],
  show_rownames = TRUE,
  main = "Top 30 DEGs",
  color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100)
)
dev.off()
```
