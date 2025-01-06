Galaxy Analysis

RNA-Seq raw data was processed in Galaxy to ensure high-quality inputs for downstream analysis. This included quality control, adapter trimming, and alignment of reads to the reference genome with tools like FastQC, Cutadapt and HISAT2. 

Gene Quantification: After alignment, gene expression levels were quantified to generate a normalized count matrix. This step ensures comparability across samples and accounts for differences in sequencing depth.
Annotation: Gene IDs were annotated with corresponding symbols and biological information to make possible easier interpretation of the results.

R Analysis

Data Preprocessing: Annotation and normalized counts were imported into R for further analysis.
Gene annotations were merged with normalized counts to create a unified dataset for interpretation.
Duplicate entries were aggregated, and missing values were removed to ensure data quality.

PCA (Principal Component Analysis): Performed to visualize overall differences in gene expression between treated and mock samples, providing insights into sample clustering and variability.
Scatter Plot: Treated and mock means were compared to observe expression patterns.

Differential Gene Expression: Genes showing significant differences between treated and mock samples were identified using statistical tests.
Adjusted p-values were calculated to control for false discovery rates.

Pathway Enrichment Analysis: GO (Gene Ontology) and KEGG (Kyoto Encyclopedia of Genes and Genomes): Identified biological pathways enriched with significantly expressed genes.
A background gene set (universe) was used to ensure reliable pathway analysis.
Results were visualized in dot plots showing key pathways with gene counts and statistical significance.

Visualization: Volcano plots highlighted differentially expressed genes, categorizing them as upregulated, downregulated, or non-significant.
Heatmaps displayed the top 20 differentially expressed genes across samples, offering a visual summary of expression changes.
A violin plot is also included. 




