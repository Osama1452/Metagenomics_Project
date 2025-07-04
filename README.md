Metagenomics Analysis Project

Project Overview

This project focuses on the metagenomic analysis of 16S rRNA gene amplicon sequencing data to investigate microbial community diversity, the impact of antibiotic treatment on microbial populations, and comparison of different sample types. The DADA2 pipeline was utilized for raw data processing, followed by phyloseq analysis for alpha and beta diversity exploration, microbial community composition, and identification of differentially abundant taxa.

Data

The raw sequencing data (FASTQ files) used in this analysis are publicly available on Figshare (https://doi.org/10.6084/m9.figshare.29458793.v1).

Tools and Technologies Used

The analysis was performed using the R programming language, leveraging the following key packages:

DADA2: For processing high-throughput amplicon sequencing data (quality filtering, denoising, merging reads, chimera removal, ASV inference).
phyloseq: For microbiome data analysis and visualization (alpha/beta diversity, taxonomic plots, heatmaps).
ggplot2: For creating high-quality data visualizations.
vegan: For statistical ecological analyses (often used in conjunction with phyloseq).
DESeq2: For differential abundance analysis between groups.
pheatmap / ComplexHeatmap: For generating heatmaps.
corrplot: For visualizing correlation matrices.
microbiome: For core microbiota analysis.

Project Structure

This repository contains the following main files and directories:

dada2_script.R: The primary R script containing all commands used for metagenomic data processing, analysis, and visualization.
mapping.txt: A metadata file containing sample information crucial for phyloseq analyses.
Results/: This directory stores all the output from the analysis, including:
filtered/: FASTQ files after quality filtering and trimming.
dataset/: original FASTQ files are stored here.
Alpha_diversity/: Directory containing plots and tables related to alpha diversity analyses.
beta/: Directory containing plots and tables related to beta diversity analyses.
Taxonomy Composition/: Directory for plots illustrating microbial community composition at various taxonomic levels (Phylum, Class, Order, Family, Genus).
Taxonomy distribution/: Additional directory for taxonomic results.
aligned_sequences.fasta: Aligned sequence file.
error_F.pdf, error_R.pdf: Error estimation plots for forward and reverse reads from DADA2.
heatmap_sample.pdf: A sample heatmap of the data.
quality_F.pdf, quality_R.pdf: Quality profiles of raw reads.
Volcano Plot of Differential Abundance.pdf: Volcano plot summarizing differential abundance analysis results.
Generated data tables such as taxa_results.csv, dada5.asv_table.txt, dada5.asv_track.txt, dada5.asv_taxa.txt, alpha_diversity_measurements.csv, and core_taxa_taxonomy.csv.
.Rhistory, .RData, .RDataTmp4: Temporary files generated by R.

How to Run the Script:

Software Requirements: R and RStudio installed.
R Packages: Install all necessary R packages mentioned in the "Tools and Technologies Used" section (specifically dada2, phyloseq, ggplot2, vegan, DESeq2, microbiome, pheatmap, ComplexHeatmap, corrplot, RColorBrewer). These can be installed using install.packages() or BiocManager::install() as shown in your dada2_script.R.
Data Setup: Download the FASTQ files from the Figshare link provided above and place them into the meta directory within your project folder. Ensure the mapping.txt file is in the same directory. 
Path Adjustments: In the dada2_script.R file: Review all specified paths (e.g., setwd(), file.path(), read.delim2(), assignTaxonomy(), write.csv(), write.table()). 


Contact

If you have any questions or inquiries, please feel free to reach out.
Name: Osama
Email: osama2172003aldesoky@gmail.com
