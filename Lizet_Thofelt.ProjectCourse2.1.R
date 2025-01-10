####################################################################################################
################################# Bioinformatics Project Course ####################################
################################# RNA-Seq Data Analysis Workflow ###################################
##################################        LIZET THOFELT          ###################################
####################################################################################################

# Set Working Directory
setwd("/Users/lizetthofelt/Desktop/Bioinformatics") 
print(getwd())

# Install and Load Required Libraries
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("edgeR", "DESeq2", "affy", "sva", "assertthat", "clusterProfiler", "org.Hs.eg.db", "pathview", "enrichplot"))
install.packages(c("tidyverse", "ggplot2", "pheatmap", "RColorBrewer", "reshape2"))

# Load Libraries
library(edgeR)
library(tidyverse)
library(ggplot2)
library(pheatmap)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(RColorBrewer)
library(reshape2)
library(ggrepel)

# Set Seed for Reproducibility
set.seed(666)

##########################################
#### Data Loading and Preparation ########
##########################################

# Load Annotation and Normalized Counts
annotation <- read.delim("annotation_edgeR_Treated-Mock.tabular", header = TRUE, sep = "\t")
norm_counts <- read.delim("edgeR_normcounts.tabular", header = TRUE, sep = "\t")

# Merge Annotation with Normalized Counts
combined_table <- merge(annotation, norm_counts, by.x = "ENTREZID", by.y = "GeneID", all = TRUE)

# Aggregate Duplicate Gene Symbols
aggregated_table <- aggregate(. ~ SYMBOL, data = combined_table, FUN = sum)

# Remove Missing Values
combined_table <- na.omit(combined_table)

# Calculate Mean Expression for Mock and Treated Samples
combined_table$Mock_Mean <- rowMeans(combined_table[, 3:5])     # Mock Samples
combined_table$Treated_Mean <- rowMeans(combined_table[, 6:8])  # Treated Samples

##########################################
#### Metadata Creation ###################
##########################################a

metadata <- data.frame(
  SampleID = colnames(norm_counts)[-1],
  Condition = c("Mock", "Mock", "Mock", "Treated", "Treated", "Treated") 
)
rownames(metadata) <- metadata$SampleID


##########################################
#### Exploratory Data Analysis ###########
##########################################

# Scatter Plot: Treated vs Mock Mean Expression
ggplot(combined_table, aes(x = Treated_Mean, y = Mock_Mean)) +
  geom_point(color = "orange", size = 2, alpha = 0.3) +
  geom_abline(intercept = 0, slope = 1, color = "blue", linetype = "dashed") +
  labs(title = "Scatter Plot of Treated vs Mock Means",
       x = "Treated Mean Expression",
       y = "Mock Mean Expression") +
  theme_classic()

##########################################
#### Principal Component Analysis (PCA) ##
##########################################

# Prepare Data for PCA
gene_expression <- norm_counts[, -1]  # Remove GeneID Column
rownames(gene_expression) <- norm_counts$GeneID
pca_result <- prcomp(t(gene_expression), scale. = TRUE)

# Plot PCA with Percentage Variance
percent_var <- round((pca_result$sdev^2 / sum(pca_result$sdev^2)) * 100, 2)
ggplot(metadata, aes(x = pca_result$x[, 1], y = pca_result$x[, 2], color = Condition)) +
  geom_point(size = 4) +
  labs(title = "PCA Plot of Normalized Gene Expression",
       x = paste0("Principal Component 1 (", percent_var[1], "%)"),
       y = paste0("Principal Component 2 (", percent_var[2], "%)")) +
  theme_classic()

#######################
##### violin plot #####
#######################


# Reshape data to long format for ggplot
gene_count_long <- melt(norm_counts, id.vars = "GeneID",
                        variable.name = "Sample", value.name = "GeneCount")

# Add Condition column from metadata
gene_count_long$Condition <- metadata$Condition[match(gene_count_long$Sample, metadata$SampleID)]

# Check if Condition column is properly added
print(head(gene_count_long))  # Verify 'Condition' is not NULL

# Generate the violin plot
ggplot(gene_count_long, aes(x = Sample, y = GeneCount, fill = Condition)) +
  geom_violin(trim = TRUE, alpha = 0.6) +  # Violin plot with transparency
  geom_boxplot(width = 0.1, outlier.shape = NA, position = position_dodge(0.9)) +
  theme_classic() +  # Clean theme
  labs(title = "Violin Plot of Gene Counts per Sample",
       x = "Samples",
       y = "Gene Count") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        text = element_text(size = 15))  # Rotate x-axis labels


##########################################
#### Differential Expression Analysis ####
##########################################

# Perform t-tests
expression_data <- norm_counts[, -1]
rownames(expression_data) <- norm_counts$GeneID
results <- apply(expression_data, 1, function(row) {
  t.test(row[metadata$Condition == "Treated"], row[metadata$Condition == "Mock"])$p.value
})

# Adjust P-values
adjusted_pvalues <- p.adjust(results, method = "fdr")

# Calculate log fold change (logFC) for each gene
logFC <- rowMeans(expression_data[, c(4, 5, 6)]) - rowMeans(expression_data[, c(1, 2, 3)])  # Subtract Treated mean from Mock mean


# Update DEG_results with the calculated logFC
DEG_results <- data.frame(
  GeneID = rownames(expression_data),   # Gene identifiers
  P.Value = results,                    # P-values from t-tests
  Adj.P.Value = adjusted_pvalues,       # Adjusted P-values
  logFC = logFC                         # Use calculated logFC
)

# View the updated DEG_results
head(DEG_results)

  

# Volcano Plot
DEG_results$diffexpressed <- "NO"
DEG_results$diffexpressed[DEG_results$logFC > 1 & DEG_results$Adj.P.Value < 0.05] <- "UP"
DEG_results$diffexpressed[DEG_results$logFC < -1 & DEG_results$Adj.P.Value < 0.05] <- "DOWN"

ggplot(DEG_results, aes(x = logFC, y = -log10(Adj.P.Value), color = diffexpressed)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_hline(yintercept = -log10(0.05), col = "red", linetype = "dashed") +
  geom_vline(xintercept = c(-1, 1), col = "red", linetype = "dashed") +
  scale_color_manual(values = c("DOWN" = "blue", "UP" = "red", "NO" = "black")) +
  labs(title = "Volcano Plot of Differential Gene Expression",
       x = "Log2 Fold Change", y = "-Log10 Adjusted P-Value") +
  theme_minimal()

##########################################
#### Heatmap of Differential Genes #######
##########################################

top_genes <- rownames(DEG_results[DEG_results$Adj.P.Value < 0.05, ])[1:20]
heatmap_data <- expression_data[top_genes, ]
pheatmap(heatmap_data, scale = "row",
         annotation_col = metadata[, "Condition", drop = FALSE],
         main = "Heatmap of Top 20 Differentially Expressed Genes")


##########################################
#### Pathway Enrichment Analysis #########
##########################################

# Define the universe (background gene set)
universe_genes <- annotation$ENTREZID  

# Extract ENTREZ IDs for differentially expressed genes
entrez_ids <- DEG_results$GeneID[DEG_results$Adj.P.Value < 0.05]  # Only DEGs

# GO Enrichment Analysis
go_enrich <- enrichGO(
  gene = entrez_ids, 
  OrgDb = org.Hs.eg.db, 
  keyType = "ENTREZID", 
  ont = "BP",  # Biological Process
  pAdjustMethod = "fdr", 
  pvalueCutoff = 0.05, 
  universe = universe_genes  # Background gene set
)

# KEGG Enrichment Analysis
kegg_enrich <- enrichKEGG(
  gene = entrez_ids, 
  organism = "hsa",  # Human
  pAdjustMethod = "fdr", 
  pvalueCutoff = 0.05, 
  universe = universe_genes  # Background gene set
)


# View KEGG enrichment results
print("KEGG Enrichment Results:")
head(kegg_enrich)

# Visualize Results
# GO Enrichment Dotplot
dotplot(go_enrich, showCategory = 10, title = "Top 10 GO Enrichment with Universe")

# KEGG Enrichment Dotplot
dotplot(kegg_enrich, showCategory = 10, title = "Top 10 KEGG Enriched Pathways with Universe")

# Notes to myself 
# is there only 5 pathways that is making the cut in the KEGG? 
# To display more pathways i would need to change the thresholds or increase the number of Differentially Expressed Genes included .
head(kegg_enrich)  # View top results
dim(kegg_enrich)  # Check the total number of enriched pathways



