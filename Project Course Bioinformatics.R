# Set Working Directory

setwd("/Users/lizethofelt/Desktop/Bioinformatics R /Project Course")
getwd()

# Install and Load Required Libraries
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
install.packages(c("tidyverse", "ggplot2", "pheatmap", "RColorBrewer"))
BiocManager::install(c("gplots", "clusterProfiler", "org.Hs.eg.db", "pathview", "enrichplot"))

# Load libraries
library(edgeR)
library(tidyverse)
library(ggplot2)
library(pheatmap)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(RColorBrewer)
library(ggrepel)

# Set a seed for reproducibility
set.seed(666)


###################################################################################################



library(ggplot2)      # For plotting
library(tidyverse)    # For data manipulation
library(pheatmap)     # For heatmap visualization


# Load files exported from Galaxy (normalized counts and gene annotations)
merge_table <- read.table("edgeR_Treated-Mock.tabular", 
                          header = TRUE, sep = "\t")  # Gene annotations

norm_counts <- read.table("edgeR_normcounts.tabular", 
                          header = TRUE, sep = "\t")  # Normalized counts data


# Merge annotation table with normalized counts using gene identifiers
combined_table <- merge(merge_table, norm_counts, 
                        by.x = "ENTREZID", by.y = "GeneID", all = TRUE)


# Check how many unique gene symbols are present to ensure merging was successful
num_unique_genes <- length(unique(combined_table$SYMBOL))
print(paste("Number of unique gene symbols:", num_unique_genes))

# Aggregate Duplicate Gene Symbols
# Summarize rows by SYMBOL, summing up numeric values for duplicate gene symbols
aggregated_table <- aggregate(. ~ SYMBOL, data = combined_table, FUN = sum)

# Confirm aggregation was successful
print(paste("Number of unique symbols after aggregation:", 
            length(unique(aggregated_table$SYMBOL))))

#  Handle Missing Values 
# Remove rows with missing values (NAs) from the combined data
combined_table <- na.omit(combined_table)

# Confirm NAs are removed
print("First few rows after removing NAs:")
head(combined_table)

# Calculate the mean expression values for Mock and Treated samples
combined_table$Mock_Mean <- rowMeans(combined_table[, 3:5])     # Columns 3-5 are Mock samples
combined_table$Treated_Mean <- rowMeans(combined_table[, 6:8])  # Columns 6-8 are Treated samples

# Verify 

head(combined_table)

# Scatter Plot of Mock vs Treated Means

ggplot(combined_table, aes(x = Treated_Mean, y = Mock_Mean)) +
  geom_point(color = "orange", size = 2, alpha = 0.3) +  # Add transparent orange points
  geom_abline(intercept = 0, slope = 1, color = "blue", linetype = "dashed") +  # Diagonal reference line
  labs(title = "Scatter Plot of Treated vs Mock Means",
       subtitle = "RNA-Seq Data: Mean Expression Differences",
       x = "Treated Mean Expression",
       y = "Mock Mean Expression") +
  theme_classic()

# Extract normalized counts for the 6 samples 
counts <- norm_counts[, c("SRR11517750", "SRR11517751", "SRR11517752", 
                          "SRR11517753", "SRR11517754", "SRR11517755")]

# Create Metadata Table
# Create metadata describing the experimental conditions for each sample
metadata <- data.frame(
  Sample_ID = colnames(norm_counts)[2:7],  # Extract column names (excluding GeneID)
  Group = c("Mock", "Mock", "Mock", "Treated", "Treated", "Treated")  # Assign sample groups
)

# Display metadata
print("Metadata Table:")
print(metadata)

############################################################################################################

# THIS SECTION IS THE ORGANIZED CODE, WITH THE HELP OF CHAT GPT

# Step 1: Read and Format Data
# Load annotation and normalized counts
annotation <- read.delim("annotation_edgeR_Treated-Mock.tabular", header = TRUE, sep = "\t")
norm_counts <- read.delim("edgeR_normcounts.tabular", header = TRUE, sep = "\t")

# Verify loaded data
print("Annotation and normalized counts loaded successfully.")
head(annotation)
head(norm_counts)

# Step 2: Create Metadata
metadata <- data.frame(
  SampleID = colnames(norm_counts)[-1],
  Condition = c("Treated", "Treated", "Treated", "Mock", "Mock", "Mock")
)
rownames(metadata) <- metadata$SampleID
print("Metadata:")
print(metadata)

### from Lisam PCA plot 

# Install and Load Required Libraries
if (!requireNamespace("ggfortify", quietly = TRUE)) {
  install.packages("ggfortify")
}
library(ggfortify)

# Prepare the count matrix and metadata
# Assuming 'norm_counts' is the normalized counts data and 'metadata' contains the group information
gene_expression <- norm_counts[, -1]  # Remove the GeneID column
rownames(gene_expression) <- norm_counts$GeneID

# Transpose the count matrix for PCA (samples as rows, genes as columns)
count_pca <- prcomp(t(gene_expression), scale. = TRUE)

# Plot the PCA using autoplot
autoplot(count_pca, data = metadata, colour = "Condition", size = 4) +
  theme_classic() +
  labs(title = "PCA Plot of Normalized Gene Expression",
       x = "Principal Component 1", y = "Principal Component 2",
       color = "Condition")
### add percentage in PCA plot 

# Install and Load Required Libraries
if (!requireNamespace("ggfortify", quietly = TRUE)) {
  install.packages("ggfortify")
}
library(ggfortify)

# Prepare the count matrix and metadata
gene_expression <- norm_counts[, -1]  # Remove the GeneID column
rownames(gene_expression) <- norm_counts$GeneID

# Transpose the count matrix for PCA
count_pca <- prcomp(t(gene_expression), scale. = TRUE)

# Calculate percentage variance for the first two PCs
percent_var <- round((count_pca$sdev^2 / sum(count_pca$sdev^2)) * 100, 2)
pc1_var <- percent_var[1]
pc2_var <- percent_var[2]

# Plot the PCA using autoplot
autoplot(count_pca, data = metadata, colour = "Condition", size = 4) +
  theme_classic() +
  labs(
    title = "PCA Plot of Normalized Gene Expression",
    x = paste0("Principal Component 1 (", pc1_var, "%)"),
    y = paste0("Principal Component 2 (", pc2_var, "%)"),
    color = "Condition"
  ) +
  theme(plot.title = element_text(hjust = 0.5, size = 16))  # Centered title with larger font
###



# Step 3: Unsupervised Clustering - PCA
pca_data <- t(norm_counts[, -1])
rownames(pca_data) <- metadata$SampleID

# Perform PCA
pca_result <- prcomp(pca_data, scale. = TRUE)
metadata$PC1 <- pca_result$x[, 1]
metadata$PC2 <- pca_result$x[, 2]

# PCA Plot
ggplot(metadata, aes(x = PC1, y = PC2, color = Condition, label = SampleID)) +
  geom_point(size = 4) +
  geom_text(vjust = 1.5, size = 3) +
  labs(title = "PCA Plot of Normalized Gene Expression",
       x = "Principal Component 1", y = "Principal Component 2") +
  theme_minimal()

# Step 4: Statistical Testing (t-test) and Adjusted P-values
expression_data <- norm_counts[, -1]
rownames(expression_data) <- norm_counts$GeneID

# Perform t-test for each gene
results <- apply(expression_data, 1, function(row) {
  t.test(row[metadata$Condition == "Treated"], row[metadata$Condition == "Mock"])$p.value
})
adjusted_pvalues <- p.adjust(results, method = "fdr")

# Create Results Table
DEG_results <- data.frame(
  GeneID = rownames(expression_data),
  P.Value = results,
  Adj.P.Value = adjusted_pvalues,
  logFC = rnorm(nrow(expression_data))  # Simulated logFC for visualization
)

# Significant Genes
significant_genes <- DEG_results[DEG_results$Adj.P.Value < 0.05, ]

# Step 5: Volcano Plot
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

# Step 6: Heatmap of Top 20 Differentially Expressed Genes
top_genes <- rownames(significant_genes)[1:20]
heatmap_data <- expression_data[top_genes, ]

pheatmap(heatmap_data,
         scale = "row",
         annotation_col = metadata[, "Condition", drop = FALSE],
         main = "Heatmap of Top 20 Differentially Expressed Genes")

# Step 7: Biological Pathway Analysis (GO and KEGG)
entrez_ids <- significant_genes$GeneID

# GO Enrichment
go_enrich <- enrichGO(gene = entrez_ids, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",
                      ont = "BP", pAdjustMethod = "fdr", pvalueCutoff = 0.05)

dotplot(go_enrich, showCategory = 10, title = "Top 10 GO Enrichment")

# KEGG Enrichment
kegg_enrich <- enrichKEGG(gene = entrez_ids, organism = "hsa",
                          pAdjustMethod = "fdr", pvalueCutoff = 0.05)

dotplot(kegg_enrich, showCategory = 10, title = "Top 10 KEGG Enriched Pathways")

# Save Enrichment Results
write.csv(as.data.frame(go_enrich), file = "GO_Enrichment_Results.csv")
write.csv(as.data.frame(kegg_enrich), file = "KEGG_Enrichment_Results.csv")



###########################################################################################

setwd("/Users/lizethofelt/Desktop/Bioinformatics R /Project Course")
getwd()

library(edgeR)
install.packages("tidyverse")
library(tidyverse)
install.packages("ggplot2")
library(ggplot2)
install.packages("pheatmap")
library(pheatmap)
library(BiocManager)
library(RColorBrewer) # For color palettes in visualization
# Install the packages
BiocManager::install("gplots")
BiocManager::install("clusterProfiler")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("pathview")

# Set a seed to produce reproducible results
set.seed(666)


# Load the annotation file
annotation <- read.delim("annotation_edgeR_Treated-Mock.tabular", header = TRUE, sep = "\t")
print("Annotation file loaded successfully!")
head(annotation)  # Check the first few rows

# Load normalized counts file
norm_counts <- read.delim("edgeR_normcounts.tabular", header = TRUE, sep = "\t")
print("Normalized counts file loaded successfully!")
head(norm_counts)  # Check the first few rows

# Create a metadata file
metadata <- data.frame(
  SampleID = colnames(norm_counts)[-1],  # Exclude GeneID column
  Condition = c("Treated", "Treated", "Treated", "Mock", "Mock", "Mock")
)

# View metadata
print(metadata)

#ensure it aligns with the columns in the normalized counts table.

# Create metadata for the 6 samples (3 Treated, 3 Mock)
metadata <- data.frame(
  SampleID = colnames(norm_counts)[-1],  # Exclude 'GeneID' column
  Condition = c("Treated", "Treated", "Treated", "Mock", "Mock", "Mock")
)

# Check the metadata
print("Metadata:")
print(metadata)

# Double-check that the samples in 'norm_counts' align with metadata
print("Column names of norm_counts:")
print(colnames(norm_counts)[-1])



# Load necessary library
library(ggplot2)

# Prepare the normalized counts data (exclude GeneID)
pca_data <- t(norm_counts[, -1])  # Transpose data: samples as rows, genes as columns
rownames(pca_data) <- metadata$SampleID

# Perform PCA
pca_result <- prcomp(pca_data, scale. = TRUE)

# Add PCA results to metadata for plotting
metadata$PC1 <- pca_result$x[, 1]
metadata$PC2 <- pca_result$x[, 2]

# Plot PCA
ggplot(metadata, aes(x = PC1, y = PC2, color = Condition, label = SampleID)) +
  geom_point(size = 4) +
  geom_text(vjust = 1.5, size = 3) +
  labs(title = "PCA Plot of Normalized Gene Expression",
       x = "Principal Component 1",
       y = "Principal Component 2") +
  theme_minimal()


# making violin plot from visualization workshop 

# Install and load required libraries
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("reshape2", quietly = TRUE)) install.packages("reshape2")
library(ggplot2)
library(reshape2)

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


# Prepare the data
expression_data <- norm_counts[, -1]  # Exclude 'GeneID'
rownames(expression_data) <- norm_counts$GeneID

# Perform t-test for each gene
results <- apply(expression_data, 1, function(row) {
  t.test(row[metadata$Condition == "Treated"], row[metadata$Condition == "Mock"])$p.value
})

# Adjust p-values for multiple testing using FDR
adjusted_pvalues <- p.adjust(results, method = "fdr")

# Create a results table
DEG_results <- data.frame(
  GeneID = rownames(expression_data),
  P.Value = results,
  Adj.P.Value = adjusted_pvalues
)

# Filter significant genes (Adj.P.Value < 0.05)
significant_genes <- DEG_results[DEG_results$Adj.P.Value < 0.05, ]

# View top significant genes
print("Top significant genes:")
head(significant_genes)




# Add log fold-change (simulated for visualization purposes)
DEG_results$logFC <- rnorm(nrow(DEG_results))

# Volcano plot
ggplot(DEG_results, aes(x = logFC, y = -log10(Adj.P.Value))) +
  geom_point(aes(color = Adj.P.Value < 0.05), size = 1.5) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black")) +
  labs(title = "Volcano Plot of Differential Expression",
       x = "Log Fold Change",
       y = "-Log10 Adjusted P-Value") +
  theme_minimal()

# fROM lIsAm Volcano plot.R


# Load required libraries
library(ggplot2)
library(dplyr)
library(ggrepel)

# Step 1: Prepare the Data
FDR_pvalue <- DEG_results  # Use DEG_results as input

# Add a column to classify genes based on logFC and adjusted p-values
FDR_pvalue$diffexpressed <- "NO"  # Default: not differentially expressed
FDR_pvalue$diffexpressed[FDR_pvalue$logFC > 1 & FDR_pvalue$Adj.P.Value < 0.05] <- "UP"    # Upregulated
FDR_pvalue$diffexpressed[FDR_pvalue$logFC < -1 & FDR_pvalue$Adj.P.Value < 0.05] <- "DOWN"  # Downregulated

# Step 2: Create the Volcano Plot
p <- ggplot(FDR_pvalue, aes(x = logFC, y = -log10(Adj.P.Value), color = diffexpressed)) +
  geom_point(alpha = 0.6, size = 2) +  # Plot points
  theme_minimal() +
  geom_hline(yintercept = -log10(0.05), col = "red", linetype = "dashed") +  # Line for p-value cutoff
  geom_vline(xintercept = c(-1, 1), col = "red", linetype = "dashed") +      # Lines for logFC cutoffs
  scale_color_manual(values = c("DOWN" = "blue", "UP" = "red", "NO" = "black")) +
  labs(title = "Volcano Plot of Differential Gene Expression",
       x = "Log2 Fold Change (logFC)",
       y = "-Log10 Adjusted P-Value",
       color = "Expression") +
  theme(plot.title = element_text(hjust = 0.5, size = 20))

# Step 3: Add Labels for Top Genes
# Extract top 10 UP and DOWN regulated genes
top_up <- FDR_pvalue %>% filter(diffexpressed == "UP") %>% top_n(10, wt = logFC)
top_down <- FDR_pvalue %>% filter(diffexpressed == "DOWN") %>% top_n(-10, wt = logFC)
top_genes <- rbind(top_up, top_down)

# Add gene names as labels
p3 <- p + 
  geom_text_repel(data = top_genes, aes(label = GeneID), 
                  size = 3, color = "black", max.overlaps = 10)

# Step 4: Display the Volcano Plot
p3



# Load heatmap library
library(pheatmap)

#in order to do the heatmap: Row names in metadata must match the sample names (column names in heatmap_data).
#Drop = FALSE: Ensures Condition remains a data frame (not simplified to a vector).
#Alignment: Double-check column names in heatmap_data and row names in metadata match exactly.


# Assign row names to metadata based on SampleID
rownames(metadata) <- metadata$SampleID
print("Metadata after adding row names:")
print(metadata)

# Check column names of heatmap_data
print("Column names of heatmap_data:")
print(colnames(heatmap_data))

# Ensure metadata row names match
print("Row names of metadata:")
print(rownames(metadata))

# Extract top 20 significant genes
top_genes <- rownames(significant_genes)[1:20]
heatmap_data <- expression_data[top_genes, ]

# Plot heatmap
pheatmap(heatmap_data,
         scale = "row",
         annotation_col = metadata[, "Condition", drop = FALSE],
         main = "Heatmap of Top 20 Differentially Expressed Genes")

# biological pathway analysis using GO (Gene Ontology) or KEGG (Kyoto Encyclopedia of Genes and Genomes).

#from Lisam code Scetterplot_top_genes.R

# Load necessary libraries
library(ggplot2)
library(ggrepel)
library(pheatmap)

# Assuming DEG_results contains P-values and logFC:
# Rename columns for clarity: P-value -> padjusted, logFC -> log2FC
FDR_pvalue <- DEG_results
colnames(FDR_pvalue)[colnames(FDR_pvalue) == "Adj.P.Value"] <- "padjusted"
colnames(FDR_pvalue)[colnames(FDR_pvalue) == "logFC"] <- "log2FC"

# Keep only genes with significant FDR-adjusted P-values (< 0.05)
sig_pvalue <- FDR_pvalue[FDR_pvalue$padjusted < 0.05, , drop = FALSE]

# Order genes based on log2FC
sig_pvalue <- sig_pvalue[order(sig_pvalue$log2FC), , drop = FALSE]

# Extract top 10 upregulated and top 10 downregulated genes
top_genes <- sig_pvalue[c(1:10, (nrow(sig_pvalue) - 9):nrow(sig_pvalue)), ]
top_genes <- rownames(top_genes)  # Extract gene names

# Subset gene expression data (norm_counts) for the top genes
top_genes_expression <- norm_counts[rownames(norm_counts) %in% top_genes, ]

# Calculate the mean expression for Mock (1st 3 columns) and Treated (last 3 columns)
top_genes_expression$MockMean <- rowMeans(top_genes_expression[, 1:3])
top_genes_expression$TreatedMean <- rowMeans(top_genes_expression[, 4:6])

# Add gene names for labeling
top_genes_expression$genename <- rownames(top_genes_expression)

# Create the scatter plot
ggplot(top_genes_expression, aes(x = MockMean, y = TreatedMean, label = genename)) +
  geom_point(aes(color = TreatedMean), size = 3) +  # Scatter points colored by TreatedMean
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey") +  # Diagonal line
  theme_bw() +  # Clean theme
  ggtitle(label = "Top 10 Upregulated and Downregulated Genes") +
  scale_color_gradient2(low = "darkblue", mid = "white", high = "red", space = "Lab") +
  geom_text_repel(size = 3) +  # Add gene labels
  theme(plot.title = element_text(hjust = 0.5, size = 20)) +  # Centered title
  xlab("Mean Expression in Mock") +
  ylab("Mean Expression in Treated")



# Install required packages (run only if not installed)
if (!requireNamespace("clusterProfiler", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("clusterProfiler")
BiocManager::install("org.Hs.eg.db")  # Human gene database
BiocManager::install("enrichplot")    # For enrichment plots

# Load libraries
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)


# Convert Gene Symbols to ENTREZ IDs using org.Hs.eg.db
gene_symbols <- significant_genes$GeneID  # Replace with your gene column

# Map to ENTREZ IDs
entrez_ids <- bitr(gene_symbols, fromType = "SYMBOL", 
                   toType = "ENTREZID", 
                   OrgDb = org.Hs.eg.db)

# View the mapping
print("Mapped ENTREZ IDs:")
head(entrez_ids)

######################################################################################################

# Reinstall clusterProfiler
BiocManager::install("clusterProfiler", update = TRUE, ask = FALSE)
library(clusterProfiler)
?bitr  # Check the documentation

# Load required libraries
library(clusterProfiler)
library(org.Hs.eg.db)

# Define gene symbols
gene_symbols <- significant_genes$GeneID  # Replace with your actual significant genes

# Convert Gene Symbols to ENTREZ IDs
entrez_ids <- bitr(gene_symbols, fromType = "SYMBOL", 
                   toType = "ENTREZID", 
                   OrgDb = org.Hs.eg.db)

# View the mapping
print("Mapped ENTREZ IDs:")
head(entrez_ids)


# Check the first few Gene IDs
head(significant_genes$GeneID)

# If the IDs are numeric (e.g., 23352, 3151, etc.), they are likely ENTREZ IDs or another ID type, NOT SYMBOLS.

# Use ENTREZ IDs directly for pathway analysis
entrez_ids <- significant_genes$GeneID  # Assuming these are ENTREZ IDs

# Confirm the IDs
print("Using ENTREZ IDs directly:")
head(entrez_ids)

# Check the significant_genes object
print("Contents of significant_genes:")
head(significant_genes)

# Use GeneID column as ENTREZ IDs directly
entrez_ids <- significant_genes$GeneID

# Confirm ENTREZ IDs
print("Using ENTREZ IDs directly:")
head(entrez_ids)

go_enrich <- enrichGO(gene = entrez_ids,
                      OrgDb = org.Hs.eg.db,
                      keyType = "ENTREZID",
                      ont = "BP",  # Biological Process
                      pAdjustMethod = "fdr",
                      pvalueCutoff = 0.05)

# View results
head(go_enrich)

# Plot results
dotplot(go_enrich, showCategory = 10, title = "Top 10 GO Enrichment")

str(significant_genes)

#############################################################################################

# Perform KEGG enrichment analysis
kegg_enrich <- enrichKEGG(gene = entrez_ids,
                          organism = "hsa",   # Human
                          pAdjustMethod = "fdr",
                          pvalueCutoff = 0.05)

# View the top KEGG enrichment results
print("Top KEGG Enrichment Results:")
head(kegg_enrich)

# Save GO results
write.csv(as.data.frame(go_enrich), file = "GO_Enrichment_Results.csv")

# Save KEGG results
write.csv(as.data.frame(kegg_enrich), file = "KEGG_Enrichment_Results.csv")


# Plot KEGG enrichment results
dotplot(kegg_enrich, showCategory = 20, title = "Top 10 KEGG Enriched Pathways")
