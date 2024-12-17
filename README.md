# Project-Course-R
In Galaxy we downloaded RNA-seq data using accession numbers and preprocessed it by performing quality control with FastQC and fastp to clean the reads. The cleaned data was aligned to the human reference genome using HISAT2, and gene expression levels were quantified with featureCounts.

The RNA-seq data was generated from A549 lung adenocarcinoma cells treated under mock conditions for 24 hours. Total polyA RNA was extracted using the RNeasy Mini Kit (Qiagen), and cDNA libraries were prepared using the TruSeq Stranded mRNA protocol. Sequencing was performed on the Illumina NextSeq 500 platform, and reads were aligned to the human genome (hg19). This information was obtained from the NCBI GEO dataset

In R, I analyzed the quantified data to identify differentially expressed genes (DEGs) between Mock and RSV-treated samples. I performed PCA to visualize sample clustering, created volcano plots to highlight up- and downregulated genes, and generated heatmaps for the top differentially expressed genes. Finally, pathway enrichment analysis using GO and KEGG provided biological insights into the affected pathways.
