01_RNAseq_QC_and_Alignment
================
Jack Chiang, Louis Lax-Roseman
2025-04-04

1.  Introduction

This document outlines the RNA-seq data processing pipeline for
Marchantia polymorpha. The steps include:

    Quality Control (QC): using fastp, FastQC, and MultiQC.

    Alignment: of trimmed reads using STAR.

    Summarizing alignment logs.

Adjust file paths and sample names as needed for your environment.

2.  Data Preparation and Quality Control

2.1. Locating Raw FASTQ Files

Assume your raw FASTQ files (e.g., IM_Tak1_1.fastq, IM_Tak1_2.fastq,
etc.) are located in data/raw_fastq/. 2.2. Trimming and Filtering with
fastp

The following shell command demonstrates how to trim a paired-end sample
using fastp. Adjust sample names and paths accordingly.

``` r
fastp \
  -i data/raw_fastq/IM_Tak1_1.fastq \
  -I data/raw_fastq/IM_Tak1_2.fastq \
  -o data/fastp_processed/IM_Tak1_R1.filt.fastq \
  -O data/fastp_processed/IM_Tak1_R2.filt.fastq \
  --detect_adapter_for_pe \
  --thread 4 \
  --html data/fastp_processed/IM_Tak1_fastp.html
```

Repeat the above command for other samples (e.g., IM_Tak2, M_Tak1,
etc.).

2.3. Running FastQC and MultiQC

After trimming, run FastQC on the filtered FASTQ files:

``` r
fastqc data/fastp_processed/*_R1.filt.fastq data/fastp_processed/*_R2.filt.fastq -o data/fastqc_processed/
```

Then, combine the FastQC reports with MultiQC:

``` r
multiqc data/fastqc_processed/ -o data/multiqc_processed/
```

MultiQC will generate a comprehensive report (e.g.,
multiqc_report.html). 3. Alignment with STAR 3.1. Genome Index
Generation

If a STAR index has not yet been built, run:

``` r
STAR --runThreadN 8 \
     --runMode genomeGenerate \
     --genomeDir data/star_index \
     --genomeFastaFiles data/reference/MpTak1v5.1.fasta \
     --sjdbGTFfile data/reference/MpTak1v5.1_r1.gtf \
     --sjdbOverhang 99
```

3.2. Aligning Trimmed Reads

Below is an example alignment command for one sample (IM_Tak1):

``` r
STAR --runThreadN 8 \
     --genomeDir data/star_index \
     --readFilesIn data/fastp_processed/IM_Tak1_R1.filt.fastq \
                   data/fastp_processed/IM_Tak1_R2.filt.fastq \
     --outFileNamePrefix data/aligned/IM_Tak1_ \
     --outFilterMultimapNmax 10 \
     --outFilterMismatchNoverLmax 0.05 \
     --quantMode GeneCounts \
     --outSAMtype BAM SortedByCoordinate \
     --outReadsUnmapped Fastx
```

Repeat for all samples to generate sorted BAM files and alignment log
files.

3.3. Summarizing STAR Alignment Logs

The following R code reads STAR log files and extracts basic metrics.

``` r
library(dplyr)
library(readr)

# Define log files (adjust paths as needed)
log_files <- c(
  "data/aligned/IM_Tak1_Log.final.out",
  "data/aligned/IM_Tak2_Log.final.out",
  "data/aligned/M_Tak1_Log.final.out",
  "data/aligned/M_Tak2_Log.final.out",
  "data/aligned/V_Tak1_Log.final.out",
  "data/aligned/V_Tak2_Log.final.out"
)

# Function to extract metrics from a STAR log file
extract_star_metrics <- function(file) {
  lines <- readLines(file)
  uniquely_mapped <- as.numeric(gsub(".*\\|\\s*", "", lines[grep("Uniquely mapped reads number", lines)]))
  mapped_percent <- as.numeric(gsub("%", "", gsub(".*\\|\\s*", "", lines[grep("Uniquely mapped reads %", lines)])))
  
  data.frame(
    sample = gsub("_Log.final.out", "", basename(file)),
    uniquely_mapped_reads = uniquely_mapped,
    uniquely_mapped_percent = mapped_percent
  )
}

star_log_data <- do.call(rbind, lapply(log_files, extract_star_metrics))
print(star_log_data)
```

Generate a bar plot to visualize uniquely mapped reads per sample:

``` r
library(ggplot2)

ggplot(star_log_data, aes(x = sample, y = uniquely_mapped_reads, fill = sample)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "Uniquely Mapped Reads per Sample", x = "Sample", y = "Uniquely Mapped Reads")
```

4.  Conclusions

The outputs from this pipeline (BAM files, QC reports, and alignment
summaries) will be used in downstream analyses (TF annotation,
expression analysis, and advanced plots) in subsequent R Markdown
documents. How to Use This Document

# References

- Andrews, S. (2010). FastQC: a quality control tool for high throughput
  sequence data.
- Chen, S. et al. (2018). fastp: an ultra-fast all-in-one FASTQ
  preprocessor.
- Dobin, A. et al. (2013). STAR: ultrafast universal RNA-seq aligner.
