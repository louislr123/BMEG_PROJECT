03_GO_Analysis_and_Advanced_Plots
================
Jack Chiang, Louis Lax-Roseman
2025-04-02

1.  Introduction

This document performs downstream analyses including GO enrichment,
advanced visualizations (upset plot, PCA, and correlation analysis), and
interprets TF expression and evolution in Marchantia polymorpha.

2.  GO Enrichment Analysis

``` r
library(dplyr)
library(tidyr)
library(ggplot2)
library(GO.db)
```

Load TF expression data with GO annotations (assumed to be in
tf_expression_GO.csv):

``` r
tf_go_data <- read.csv("data/tf_expression_GO.csv", stringsAsFactors = FALSE)
head(tf_go_data)
```

If necessary, merge external GO annotations here. Assume go_annotations
is a data frame with columns: Gene_ID and GO_Term.

``` r
go_annotations <- read.delim("data/go_annotation.tsv", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
```

Define a function to perform GO enrichment analysis:

``` r
perform_go_enrichment <- function(gene_set, go_annotations, background_genes) {
  go_counts <- go_annotations %>% 
    filter(Gene_ID %in% gene_set) %>% 
    group_by(GO_Term) %>% 
    summarise(Count = n())
  
  background_counts <- go_annotations %>% 
    group_by(GO_Term) %>% 
    summarise(Background_Count = n())
  
  enrichment <- merge(go_counts, background_counts, by = "GO_Term", all.x = TRUE)
  enrichment <- enrichment %>% mutate(Total_Genes = length(gene_set), Total_Background = length(background_genes))
  
  enrichment <- enrichment %>% rowwise() %>% 
    mutate(P_Value = fisher.test(matrix(c(Count, Total_Genes - Count, Background_Count, Total_Background - Background_Count), nrow = 2))$p.value) %>% 
    ungroup() %>% 
    mutate(Adjusted_P_Value = p.adjust(P_Value, method = "BH"))
  
  sig_go <- enrichment %>% filter(Adjusted_P_Value < 0.05) %>% arrange(Adjusted_P_Value)
  return(sig_go)
}
```

Example usage for a condition (e.g., IM_Tak1):

``` r
im_tak1_genes <- tf_go_data %>% filter(IM_Tak1 > 0) %>% pull(GeneID)
background_genes <- unique(tf_go_data$GeneID)
sig_go_terms <- perform_go_enrichment(im_tak1_genes, go_annotations, background_genes)
head(sig_go_terms)
```

Visualize enriched GO terms if significant terms are found:

``` r
if(nrow(sig_go_terms) > 0){
  ggplot(sig_go_terms, aes(x = reorder(GO_Term, -log10(Adjusted_P_Value)), y = -log10(Adjusted_P_Value))) +
    geom_point(size = 3, color = "darkred") +
    coord_flip() +
    labs(title = "Enriched GO Terms (IM_Tak1)", x = "GO Term", y = "-log10(Adjusted P-Value)") +
    theme_minimal()
} else {
  message("No significant GO terms found for IM_Tak1.")
}
```

3.  Upset Plot for TF Presence/Absence

``` r
library(ComplexHeatmap)
```

Create a binary presence/absence matrix for TF genes:

``` r
tf_binary <- tf_go_data %>%
  mutate(across(c(IM_Tak1, IM_Tak2, M_Tak1, M_Tak2, V_Tak1, V_Tak2), ~ ifelse(. > 0, 1, 0))) %>%
  select(GeneID, IM_Tak1, IM_Tak2, M_Tak1, M_Tak2, V_Tak1, V_Tak2)
```

Define the condition list:

``` r
condition_list <- list(
  IM_Tak1 = tf_binary %>% filter(IM_Tak1 == 1) %>% pull(GeneID),
  IM_Tak2 = tf_binary %>% filter(IM_Tak2 == 1) %>% pull(GeneID),
  M_Tak1  = tf_binary %>% filter(M_Tak1  == 1) %>% pull(GeneID),
  M_Tak2  = tf_binary %>% filter(M_Tak2  == 1) %>% pull(GeneID),
  V_Tak1  = tf_binary %>% filter(V_Tak1  == 1) %>% pull(GeneID),
  V_Tak2  = tf_binary %>% filter(V_Tak2  == 1) %>% pull(GeneID)
)
```

Create the combination matrix and generate the Upset plot:

``` r
comb_mat <- make_comb_mat(condition_list, mode = "distinct")
UpSet(comb_mat,
      set_order = c("IM_Tak1", "IM_Tak2", "M_Tak1", "M_Tak2", "V_Tak1", "V_Tak2"),
      comb_order = order(comb_size(comb_mat), decreasing = TRUE),
      top_annotation = HeatmapAnnotation("Intersection Size" = anno_barplot(comb_size(comb_mat), add_numbers = TRUE)),
      right_annotation = rowAnnotation("Set Size" = anno_barplot(set_size(comb_mat), add_numbers = TRUE)),
      column_title = "Upset Plot: TF Expression Across Conditions",
      row_title = NULL)
```

4.  Advanced Visualizations: PCA and Correlation Analysis

``` r
library(ggrepel)
```

Perform PCA on TF expression data:

``` r
species_cols <- c("IM_Tak1", "IM_Tak2", "M_Tak1", "M_Tak2", "V_Tak1", "V_Tak2")
tf_expression_matrix <- tf_go_data %>% select(all_of(species_cols)) %>% as.matrix()
pca_result <- prcomp(tf_expression_matrix, scale. = TRUE)
pca_df <- as.data.frame(pca_result$x)
pca_df$GeneID <- tf_go_data$GeneID

ggplot(pca_df, aes(x = PC1, y = PC2)) +
  geom_point(alpha = 0.7, color = "darkblue") +
  geom_text_repel(aes(label = GeneID), size = 2) +
  labs(title = "PCA of TF Expression",
       x = paste0("PC1 (", round(summary(pca_result)$importance[2,1] * 100, 1), "%)"),
       y = paste0("PC2 (", round(summary(pca_result)$importance[2,2] * 100, 1), "%)")) +
  theme_minimal()
```

Perform correlation analysis between TF family count and average TPM:

``` r
avg_tpm <- tf_go_data %>% group_by(Family) %>% 
  summarize(Avg_TPM = mean(c(IM_Tak1, IM_Tak2, M_Tak1, M_Tak2, V_Tak1, V_Tak2), na.rm = TRUE))
cor_data <- merge(tf_family_counts, avg_tpm, by = "Family")

ggplot(cor_data, aes(x = Count, y = Avg_TPM, label = Family)) +
  geom_point(color = "forestgreen", size = 3) +
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  geom_text_repel() +
  labs(title = "Correlation Between TF Family Count and Average TPM",
       x = "TF Family Count",
       y = "Average TPM") +
  theme_minimal()
```

## References

- Ashburner, M. et al. (2000). Gene Ontology: tool for the unification
  of biology.
- Fisher, R. A. (1922). On the interpretation of χ² from contingency
  tables.
- Gu, Z. et al. (2016). Complex heatmaps reveal patterns in genomic
  data.
