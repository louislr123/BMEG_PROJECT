02_TF_Annotation_Expression
================
Jack Chiang, Louis Lax-Roseman
2025-04-02

1.  Introduction

This document imports transcript quantifications from Salmon, annotates
genes using a reference GTF and PlantTFDB, and summarizes transcription
factor (TF) families with basic plots.

2.  Building the Expression Matrix 2.1. Importing Salmon Quantifications

``` r
library(tximport)
library(readr)
library(dplyr)
library(tidyr)

# Define sample names and file paths (adjust to match your directory structure)
samples <- c("IM_Tak1", "IM_Tak2", "M_Tak1", "M_Tak2", "V_Tak1", "V_Tak2")
files <- file.path("data/salmon_quant", samples, "quant.sf")
names(files) <- samples

# Import transcript-level quantifications
txi_tx <- tximport(files, type = "salmon", txOut = TRUE)
head(txi_tx$abundance)
```

2.2. Creating a Transcript-to-Gene Map

``` r
library(GenomicFeatures)

# Path to reference GTF file (adjust as needed)
gtf_file <- "data/reference/MpTak1v5.1_r1.gtf"
txdb <- makeTxDbFromGFF(file = gtf_file, format = "gtf")
tx2gene <- select(txdb,
                  keys = keys(txdb, keytype = "TXNAME"),
                  columns = c("TXNAME", "GENEID"),
                  keytype = "TXNAME")
head(tx2gene)
write.csv(tx2gene, "data/tx2gene.csv", row.names = FALSE)
```

If the file tx2gene.csv is available, import gene-level quantifications:

``` r
if (file.exists("data/tx2gene.csv")) {
  tx2gene <- read.csv("data/tx2gene.csv")
  txi_gene <- tximport(files, type = "salmon", tx2gene = tx2gene)
  expression_matrix <- as.data.frame(txi_gene$abundance)
  head(expression_matrix)
}
```

3.  Merging TF Annotations

``` r
# Load PlantTFDB data for M. polymorpha (assumes file "Mpo_TF_list.txt" with columns Gene_ID and Family)
marchantia_tfs <- read.delim("data/Mpo_TF_list.txt", sep = "\t", header = TRUE)
head(marchantia_tfs)

# Merge expression data with TF annotation. Adjust matching columns as needed.
if (exists("expression_matrix")) {
  expression_matrix$GeneID <- rownames(expression_matrix)
  tf_expression <- merge(expression_matrix, marchantia_tfs,
                         by.x = "GeneID", by.y = "Gene_ID", all.x = TRUE)
  head(tf_expression)
}
```

4.  Summarizing TF Families

``` r
library(dplyr)

tf_only <- tf_expression %>% filter(!is.na(Family))
non_tf <- tf_expression %>% filter(is.na(Family))

# Count TF families
tf_family_counts <- tf_only %>%
  group_by(Family) %>%
  summarize(Count = n_distinct(GeneID)) %>%
  arrange(desc(Count))
knitr::kable(tf_family_counts, caption = "TF Family Counts in Marchantia Polymorpha")
```

4.3. Basic Bar Chart of TF Family Distribution

``` r
library(ggplot2)

ggplot(tf_family_counts, aes(x = reorder(Family, -Count), y = Count)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Distribution of TF Families in Marchantia", x = "TF Family", y = "Gene Count")
```

5.  (Optional) Basic Expression Patterns

``` r
library(tidyr)

long_tf_exp <- tf_only %>%
  select(GeneID, Family, IM_Tak1, IM_Tak2, M_Tak1, M_Tak2, V_Tak1, V_Tak2) %>%
  pivot_longer(cols = c(IM_Tak1, IM_Tak2, M_Tak1, M_Tak2, V_Tak1, V_Tak2),
               names_to = "Sample",
               values_to = "TPM")

ggplot(long_tf_exp, aes(x = Family, y = TPM)) +
  geom_boxplot() +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "TF Expression Levels Across Samples", x = "TF Family", y = "TPM")
```

6.  Conclusions / Next Steps

This document covers the import of Salmon quantifications,
transcript-to-gene mapping, TF annotation, and basic TF family
summaries.

# References

- Patro, R. et al. (2017). Salmon provides fast and bias-aware
  quantification of transcript expression.
- Lawrence, M. et al. (2013). Software for Computing and Annotating
  Genomic Ranges.
- Jin, J. et al. (2017). PlantTFDB 4.0: toward a central hub for
  transcription factors and regulatory interactions in plants.
