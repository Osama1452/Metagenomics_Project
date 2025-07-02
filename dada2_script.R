#DADA2 Tutorial

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("dada2", version = "3.21")

library("dada2")
getwd()
list.files()
setwd("D:/Osama/biotec/Bioinformatics/metagenomeics/meta")
list.files()
#sort out the .fastq files into forward and reverse reads and separate them into
# fnFs and fnRs groups
fnFs <- sort(list.files("D:/Osama/biotec/Bioinformatics/metagenomeics/meta", pattern = "_R1_001.fastq"))
fnRs <- sort(list.files("D:/Osama/biotec/Bioinformatics/metagenomeics/meta", pattern = "_R2_001.fastq"))

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(fnFs, "_"), `[`, 1)
sample.names <- sapply(strsplit(fnRs, "_"), `[`, 1)

# Specify the full path to the fnFs and fnRs
fnFs <- file.path("D:/Osama/biotec/Bioinformatics/metagenomeics/meta", fnFs)
fnRs <- file.path("D:/Osama/biotec/Bioinformatics/metagenomeics/meta", fnRs)

#examine quality profiles
plotQualityProfile(fnFs[1:6])
plotQualityProfile(fnRs[1:6])

#filtering and trimming
#filt_path <- file.path("D:/Osama/biotec/Bioinformatics/metagenomeics/meta/", "D:/Osama/biotec/Bioinformatics/metagenomeics/meta/filtered") # Place filtered files in filtered/ subdirectory
filtFs <- file.path("D:/Osama/biotec/Bioinformatics/metagenomeics/meta/filtered", paste0(sample.names, "_F_filt.fastq"))
filtRs <- file.path("D:/Osama/biotec/Bioinformatics/metagenomeics/meta/filtered", paste0(sample.names, "_R_filt.fastq"))


out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(285,225),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=FALSE) # On Windows set multithread=FALSE
head(out)

out

print(filtFs)
print(filtRs)

#learn error
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)


#plot the errors as a sanity check 
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)



#dereplicate sequences
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)

# Name the derep-Species objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names

#sample inference using sequence variant inference New_project/algorithm on dereped data
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errF, multithread=TRUE)

#Inspecting the dada-Species object returned by dada:
dadaFs

#Merge paired reads
#Spurious sequence variants are further reduced by merging overlapping reads. The core function here is mergePairs, which depends on the forward and reverse reads being in matching Species at the FOV4 they were dereplicated.

#Merge the denoised forward and reverse reads:
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)

# Inspect the merger data.frame from the first sample
head(mergers[[1]])

#We now have a data.frame for each sample with the merged $sequence, its $abundance, and the indices of the merged $forward and  $reverse denoised sequences. Paired reads that did not exactly overlap were removed by mergePairs.

#Construct sequence table
#We can now construct a â€œsequence tableâ€ of our mouse samples, a higher-resolution version of the â€œOTU tableâ€ produced by Speciesical methods:

seqtab <- makeSequenceTable(mergers)
dim(seqtab)
seqtab


#Remove chimeras
#The core dada method removes substitution and indel errors, but chimeras remain. Fortunately, the accuracy of the sequences after denoising makes identifying chimeras simpler than it is when dealing with fuzzy OTUs: all sequences which can be exactly reconstructed as a bimera (two-parent chimera) from more abundant sequences.

#Remove chimeric sequences:

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)

dim(seqtab.nochim)

sum(seqtab.nochim)/sum(seqtab)

#The fraction of chimeras varies based on factors including experimental procedures and sample complexity, but can be substantial. Here chimeras make up about 20% of the inferred sequence variants, but those variants account for only about 4% of the total sequence reads.

#Track reads through the pipeline
#As a final check of our progress, weâ€™ll look at the number of reads that made it through each step in the pipeline:

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(mergers, getN), rowSums(seqtab), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoised", "merged", "tabled", "nonchim")
rownames(track) <- sample.names
head(track)


#Looks good, we kept the majority of our raw reads, and there is no over-large drop associated with any single step.

#Assign taxonomy
#It is common at this point, especially in 16S/18S/ITS amplicon sequencing, to Speciesify sequence variants taxonomically. The DADA2 package provides a native implementation of the RDPâ€™s naive Bayesian Speciesifier for this purpose. The assignTaxonomy function takes a set of sequences and a training set of taxonomically Speciesified sequences, and outputs the taxonomic assignments with at least minBoot bootstrap confidence.

#We maintain formatted training fastas for the RDP training set, GreenGenes clustered at 97% identity, and the Silva reference database. For fungal taxonomy, the General Fasta release files from the UNITE ITS database can be used as is. To follow along, download the  rdp_train_set_16.fa.gz file, and place it in the directory with the fastq files.

#Assign taxonomy:


taxa <- assignTaxonomy(seqtab.nochim, "D:/Osama/biotec/Bioinformatics/metagenomeics/silva_nr99_v138.1_train_set.fa.gz", multithread=TRUE, tryRC = TRUE)

#Optional: The dada2 package also implements a method to make species level assignments based on exact matching between ASVs and sequenced reference strains.

#taxa <- addSpecies(taxa, "D:/Osama/biotec/Bioinformatics/metagenomeics/meta/silva_species_assignment_v128.fa.gz")


# A quick taxonomic inspection:
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

taxa_df <- as.data.frame(taxa)

write.csv(taxa_df, "D:/Osama/biotec/Bioinformatics/metagenomeics/meta/taxa_results.csv", row.names = FALSE)


#save all your files as tables externally> write.table(seqtab.nochim, file = "/Volumes/SBPD/grazing_Lewis/lewis_fastq/uparse_otus/lewis.asv_table.txt", col.names = TRUE, row.names = TRUE)
write.table(seqtab.nochim, file = "D:/Osama/biotec/Bioinformatics/metagenomeics/meta/dada5.asv_table.txt", col.names = TRUE, row.names = TRUE, sep = "\t")
write.table(track, file = "D:/Osama/biotec/Bioinformatics/metagenomeics/meta/dada5.asv_track.txt", col.names = TRUE, row.names = TRUE, sep = "\t")
write.table(taxa, file = "D:/Osama/biotec/Bioinformatics/metagenomeics/meta/dada5.asv_taxa.txt", col.names = TRUE, row.names = TRUE, sep = "\t")


#phyloseq Tutorial
if(!requireNamespace("BiocManager")){
  install.packages("BiocManager")
}
BiocManager::install("phyloseq")

install.packages("ggplot2")
install.packages("magrittr")
install.packages("ape")
install.packages("vegan")
install.packages("plyr")
install.packages("scales")
install.packages("grid")
install.packages("reshape2")

library("dada2")
library("phyloseq")
library("ggplot2")
library("magrittr")
library("Biostrings")
library("ape")
library("vegan")
library("plyr")
library("scales")
library("grid")
library("reshape2")

mapping <- read.delim2("D:/Osama/biotec/Bioinformatics/metagenomeics/meta/mapping.txt")
mapping <- as.data.frame(mapping)
print(mapping)
rownames(mapping) <- mapping[ , 1]
Species(mapping)
head(mapping)
mapping


asv_nochim = otu_table(seqtab.nochim, taxa_are_rows = F)
tax = tax_table(taxa) # make 'tax_table' into a matrix
mapping <- sample_data(mapping)
#merge everything into one object:
ps = merge_phyloseq(asv_nochim, tax, mapping)
ps 


#Alpha_diversity:#######################
indices <- c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson", "Fisher")
plot_richness(ps, x= "Treatment", color = "Treatment", measures = indices) +
  geom_boxplot(aes(fill = Treatment), alpha = 0.4)+ labs(title = "Bacterial Alpha Diversity Measurements")

indices <- c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson", "Fisher")
plot_richness(ps, x= "Sample", color = "Sample", measures = indices) +
  geom_boxplot(aes(fill = Sample), alpha = 0.4)+ labs(title = "Bacterial Alpha Diversity Measurements")

library(phyloseq)
library(ggplot2)
library(vegan)
alpha_measures <- c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson", "Fisher")
alpha_diversity_data <- estimate_richness(ps, measures = alpha_measures)
head(alpha_diversity_data)
sample_data(ps) <- cbind(sample_data(ps), alpha_diversity_data)
head(sample_data(ps))


p1 <- plot_richness(ps, x = "Treatment", color = "Treatment", measures = alpha_measures) +
  geom_boxplot(aes(fill = Treatment), alpha = 0.4, outlier.shape = NA) + 
  geom_jitter(aes(color = Treatment), width = 0.2, alpha = 0.7, size = 2) + 
  labs(title = "Bacterial Alpha Diversity by Treatment") +
  theme_bw() 

print(p1)

p2 <- plot_richness(ps, x = "Sample", color = "Sample", measures = alpha_measures) +
  geom_boxplot(aes(fill = Sample), alpha = 0.4, outlier.shape = NA) +
  geom_jitter(aes(color = Sample), width = 0.2, alpha = 0.7, size = 2) +
  labs(title = "Bacterial Alpha Diversity by Sample") +
  theme_bw()

print(p2)
ggplot(sample_data(ps), aes(x = Shannon, fill = Treatment)) +
  geom_density(alpha = 0.5) +
  labs(title = "Density Plot of Shannon Diversity by Treatment") +
  theme_bw()

write.csv(alpha_diversity_data, file = "alpha_diversity_measurements.csv", row.names = TRUE)
print("Alpha diversity data saved to 'alpha_diversity_measurements.csv'")

ps

library(phyloseq)
library(ggplot2)
library(vegan)
library(pheatmap) # For a simple, elegant heatmap
library(ComplexHeatmap)
# Re-calculating alpha diversity data just in case, using the same measures
alpha_measures <- c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson", "Fisher")
alpha_diversity_data <- estimate_richness(ps, measures = alpha_measures)

# Ensure the data is in a format suitable for heatmaps (numeric matrix)
# It's already in the correct format, but sometimes subsetting might be needed
# If you want to include sample metadata as annotations, we need to prepare it.
sample_info <- as.data.frame(sample_data(ps))

# We need to make sure the row names of alpha_diversity_data match the sample names from ps
# which should already be the case if it's derived directly from ps.
if (!all(rownames(alpha_diversity_data) == rownames(sample_info))) {
  stop("Row names of alpha_diversity_data and sample_info do not match. Check data integrity.")
}

print("Alpha diversity data prepared for heatmaps.")
head(alpha_diversity_data)



# Convert alpha_diversity_data to a matrix, which pheatmap expects
alpha_div_matrix <- as.matrix(alpha_diversity_data)

# Extract the 'Treatment' information for annotation
# Ensure 'Treatment' is a factor if it isn't already
sample_info$Treatment <- as.factor(sample_info$Treatment)
annotation_col <- data.frame(Treatment = sample_info$Treatment)
rownames(annotation_col) <- rownames(sample_info)

# Create the heatmap using pheatmap
pheatmap(
  alpha_div_matrix,
  main = "Alpha Diversity Heatmap (pheatmap)",
  cluster_rows = TRUE,       # Cluster the diversity indices
  cluster_cols = TRUE,       # Cluster the samples
  show_rownames = TRUE,      # Show names of diversity indices
  show_colnames = TRUE,      # Show sample names
  scale = "column",          # Scale values within each column (per sample) to highlight relative differences
  # Or scale = "row" (per index) or "none"
  color = colorRampPalette(c("blue", "white", "red"))(50), # Custom color scale
  annotation_col = annotation_col, # Add column annotations for Treatment
  gaps_col = NULL,           # No gaps in columns by default
  fontsize_row = 8,          # Font size for row labels
  fontsize_col = 8           # Font size for column labels
)

print("Basic pheatmap generated.")



# Re-calculate alpha diversity data, as this is the matrix we're annotating
alpha_measures <- c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson", "Fisher")
alpha_diversity_data <- estimate_richness(ps, measures = alpha_measures)

# Convert alpha_diversity_data to a matrix, which pheatmap expects
alpha_div_matrix <- as.matrix(alpha_diversity_data)

# --- CRUCIAL FIX: Ensure sample_info is derived directly from the *current* phyloseq object ---
# This ensures that sample_info's rows correspond exactly to the samples in alpha_div_matrix.
# We'll use the original 'ps' object for sample data, as alpha_diversity_data is calculated from it.
sample_info <- as.data.frame(sample_data(ps))

# Now, filter sample_info to only include samples that are actually in alpha_div_matrix,
# and ensure the order matches. Although usually they should match if alpha_diversity_data
# was directly derived from 'ps'.
# It's good practice to explicitly ensure names match for robustness.
sample_names_in_matrix <- colnames(alpha_div_matrix)
sample_info <- sample_info[sample_names_in_matrix, , drop = FALSE] # Ensure order and presence

# Extract the 'Treatment' information for annotation
# Ensure 'Treatment' is a factor. This is good practice.
sample_info$Treatment <- as.factor(sample_info$Treatment)

# Create the annotation_col dataframe.
# It should ONLY contain the columns you want to use as annotations.
# And its rownames MUST match the colnames of your alpha_div_matrix.
annotation_col <- data.frame(Treatment = sample_info$Treatment)
rownames(annotation_col) <- rownames(sample_info) # This will be the sample names

# --- Verification steps (highly recommended for debugging) ---
print(paste("Number of columns in alpha_div_matrix:", ncol(alpha_div_matrix)))
print(paste("Number of rows in annotation_col:", nrow(annotation_col)))
print("Column names of alpha_div_matrix (first 5):")
print(head(colnames(alpha_div_matrix), 5))
print("Row names of annotation_col (first 5):")
print(head(rownames(annotation_col), 5))

# Check if the names match and are in the same order
if (all(rownames(annotation_col) == colnames(alpha_div_matrix))) {
  print("Sample names in annotation_col and alpha_div_matrix match and are in order. Proceeding with heatmap.")
} else {
  stop("Mismatch in sample names or order between annotation_col and alpha_div_matrix. Heatmap cannot be generated correctly.")
}

# --- Create the heatmap using pheatmap with the corrected annotation_col ---
pheatmap(
  alpha_div_matrix,
  main = "Alpha Diversity Heatmap (pheatmap)",
  cluster_rows = TRUE,       # Cluster the diversity indices
  cluster_cols = TRUE,       # Cluster the samples
  show_rownames = TRUE,      # Show names of diversity indices
  show_colnames = TRUE,      # Show sample names
  scale = "column",          # Scale values within each column (per sample) to highlight relative differences
  # Or scale = "row" (per index) or "none"
  color = colorRampPalette(c("blue", "white", "red"))(50), # Custom color scale
  annotation_col = annotation_col, # Add column annotations for Treatment
  gaps_col = NULL,           # No gaps in columns by default
  fontsize_row = 8,          # Font size for row labels
  fontsize_col = 8           # Font size for column labels
)
print("Basic pheatmap generated successfully after error fix.")




# (Assume alpha_diversity_data is already defined correctly from previous steps)
# (Also assume libraries like ComplexHeatmap and circlize are loaded)
library(ComplexHeatmap)
library(circlize)
library(grid) # For grid.text for a global title

# Convert alpha_diversity_data to a matrix
alpha_div_matrix_ch <- as.matrix(alpha_diversity_data)

# Create the ComplexHeatmap object
ht1 <- Heatmap(
  alpha_div_matrix_ch,
  name = "Diversity Value",
  row_title = "Alpha Diversity Index",
  # *** FIX IS HERE: Using column_title for the main title, or adjusting column_title_gp ***
  column_title = "Alpha Diversity Heatmap (ComplexHeatmap - Simple)", # This acts as the main title for columns
  column_title_gp = gpar(fontsize = 14, fontface = "bold"), # Customize title appearance
  heatmap_legend_param = list(title = "Value"),
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 8),
  col = circlize::colorRamp2(c(min(alpha_div_matrix_ch), median(alpha_div_matrix_ch), max(alpha_div_matrix_ch)), c("blue", "white", "red")) # Custom color scale
)

# Draw the heatmap
draw(ht1)

print("Simple ComplexHeatmap generated.")

if (!requireNamespace("corrplot", quietly = TRUE)) {
  install.packages("corrplot")
}
if (!requireNamespace("pheatmap", quietly = TRUE)) {
  install.packages("pheatmap")
}
if (!requireNamespace("RColorBrewer", quietly = TRUE)) {
  install.packages("RColorBrewer")
}

library(phyloseq)
library(ggplot2)
library(vegan)
library(corrplot)
library(pheatmap)
library(RColorBrewer) 



alpha_measures <- c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson", "Fisher")
alpha_diversity_data <- estimate_richness(ps, measures = alpha_measures)
correlation_matrix <- cor(alpha_diversity_data, method = "pearson")
head(correlation_matrix)

corrplot(correlation_matrix,
         method = "circle", # أو "number", "pie", "shade", "color"
         type = "upper", 
         order = "hclust", 
         col = brewer.pal(n = 8, name = "RdBu"), 
         addCoef.col = "black", 
         tl.col = "black", 
         tl.srt = 45, 
         diag = FALSE, 
         title = "Correlation Heatmap of Alpha Diversity Measures" #
)

pheatmap(correlation_matrix,
         cluster_rows = TRUE, 
         cluster_cols = TRUE, 
         show_rownames = TRUE,
         show_colnames = TRUE,
         color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100), 
         display_numbers = TRUE, 
         number_color = "black", 
         fontsize_number = 8, 
         main = "Correlation Heatmap of Alpha Diversity Measures (pheatmap)" 
)


ann <- data.frame(Treatment = sample_data(ps)$Treatment)
rownames(ann) <- rownames(sample_data(ps))
pheatmap(alpha_diversity_data, 
         annotation = ann, 
         scale = "none", 
         main = "Heatmap of Alpha Diversity Metrics with Treatment Annotations",
         fontsize = 8,
         color = colorRampPalette(c("lightyellow", "red"))(50))

###############################
## Beta diversity
M_ps.relabund <- transform_sample_counts(ps, function (OTU) OTU/sum(OTU))
M_ps.relabund

#view("M_ps.relabund")
## make Bray-Curtis PCoA plot
bray_ord_M <- ordinate(M_ps.relabund, method = "PCoA", distance = "bray")
bray_pcoa_M <- plot_ordination(M_ps.relabund, bray_ord_M, color = "Sample") + ggtitle("PCoA plot of Bacterial Communities - Bray Curtis") + theme(plot.title = element_text(hjust = 0.5, size = 10), legend.title = element_text(size = 10))
bray_pcoa_M

## make Bray-Curtis NMDS plot
bray_ord_M <- ordinate(M_ps.relabund, method = "NMDS", distance = "bray")
bray_nmds_M <- plot_ordination(M_ps.relabund, bray_ord_M, color = "Sample") + ggtitle("NMDS plot of Bacterial Communities - Bray-Curtis") + theme(plot.title = element_text(hjust = 0.5, size = 10), legend.title = element_text(size = 10))
bray_nmds_M

## make Bray-Curtis NMDS plot
bray_ord_M <- ordinate(M_ps.relabund, method = "NMDS", distance = "bray")
bray_nmds_M <- plot_ordination(M_ps.relabund, bray_ord_M, color = "Sample") + 
  ggtitle("NMDS plot of Bacterial Communities - Bray-Curtis") + 
  theme(plot.title = element_text(hjust = 0.5, size = 10), legend.title = element_text(size = 10))
bray_nmds_M

jaccard_ord <- ordinate(M_ps.relabund, method = "PCoA", distance = "jaccard")
jaccard_pcoa <- plot_ordination(M_ps.relabund, jaccard_ord, color = "Sample") + ggtitle("PCoA Plot - Jaccard Distance") + theme(plot.title = element_text(hjust = 0.5, size = 10))
jaccard_pcoa




################################
#Heatmap
ps
plot_heatmap(ps, sample.label="Sample", sample.Species = "Sample")


top_N_taxa <- names(sort(taxa_sums(ps), TRUE)[1:500]) 
ps_topN <- prune_taxa(top_N_taxa, ps)
plot_heatmap(ps_topN, sample.label = "Sample", sample.Species = "Sample")

#Taxonomy bar plot################

#Species taxplots
dr_fam <- tax_glom(ps, taxrank = "Phylum", NArm = F)
ps0<-names(sort(taxa_sums(dr_fam), TRUE)[1:600]) #get most abundant ones
ps1<-prune_taxa(ps0, dr_fam)

asv_Phylum <- ps1 %>% 
  tax_glom(taxrank = "Phylum") %>%                     # agglomerate at Species level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  #filter(Abundance > 0.005) %>%                         # Filter out low abundance taxa
  arrange(Phylum)                                    # Sort data frame alphabetically by Species

Phylum_colors <- c(
  "purple", "green", "blue", "red","orange", "yellow","#CBD588", "#5F7FC7", "#DA5724", "#508578", "#CD9BCD",
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248", 
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861", "goldenrod", 
  "mediumpurple1", "lightcoral", "brown")



ggplot(asv_Phylum, aes(x =Sample, y = Abundance, fill = Phylum)) + 
  #facet_grid(Diet~.) +
  geom_bar(stat = "identity", position = "Fill") +
  
  # Remove x axis title
  theme(axis.title.x = element_blank()) + 
  #
  guides(fill = guide_legend(reverse = TRUE, keywidth = 1.5, keyheight = 1.5)) +
  ylab("Phylum Relative Abundance \n") +
  ggtitle("Phylum Composition of Bacterial Communities Across Different Sample ")

ggplot(asv_Phylum, aes(x =Treatment, y = Abundance, fill = Phylum)) + 
  #facet_grid(Diet~.) +
  geom_bar(stat = "identity", position = "Fill") +
  
  # Remove x axis title
  theme(axis.title.x = element_blank()) + 
  #
  guides(fill = guide_legend(reverse = TRUE, keywidth = 1.5, keyheight = 1.5)) +
  ylab("Phylum Relative Abundance \n") +
  ggtitle("Phylum Composition of Bacterial Communities Across Different Treatment ")



dr_fam <- tax_glom(ps, taxrank = "Class", NArm = F)
ps0<-names(sort(taxa_sums(dr_fam), TRUE)[1:670]) #get most abundant ones
ps1<-prune_taxa(ps0, dr_fam)

asv_Class <- ps1 %>% 
  tax_glom(taxrank = "Class") %>%                     # agglomerate at Species level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  #filter(Abundance > 0.005) %>%                         # Filter out low abundance taxa
  arrange(Class)                                    # Sort data frame alphabetically by Species

Class_colors <- c(
  "purple", "green", "blue", "red","orange", "yellow","#CBD588", "#5F7FC7", "#DA5724", "#508578", "#CD9BCD",
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248", 
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861", "goldenrod", 
  "mediumpurple1", "lightcoral", "brown")



ggplot(asv_Class, aes(x =Sample, y = Abundance, fill = Class)) + 
  #facet_grid(Diet~.) +
  geom_bar(stat = "identity", position = "Fill") +
  
  # Remove x axis title
  theme(axis.title.x = element_blank()) + 
  #
  guides(fill = guide_legend(reverse = TRUE, keywidth = 1.5, keyheight = 1.5)) +
  ylab("Class Relative Abundance \n") +
  ggtitle("Class Composition of Bacterial Communities Across Different Sample ")

ggplot(asv_Class, aes(x =Treatment, y = Abundance, fill = Class)) + 
  #facet_grid(Diet~.) +
  geom_bar(stat = "identity", position = "Fill") +
  
  # Remove x axis title
  theme(axis.title.x = element_blank()) + 
  #
  guides(fill = guide_legend(reverse = TRUE, keywidth = 1.5, keyheight = 1.5)) +
  ylab("Class Relative Abundance \n") +
  ggtitle("Class Composition of Bacterial Communities Across Different Treatment ")





dr_fam <- tax_glom(ps, taxrank = "Order", NArm = F)
ps0<-names(sort(taxa_sums(dr_fam), TRUE)[1:670]) #get most abundant ones
ps1<-prune_taxa(ps0, dr_fam)

asv_Order <- ps1 %>% 
  tax_glom(taxrank = "Order") %>%                     # agglomerate at Species level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  #filter(Abundance > 0.005) %>%                         # Filter out low abundance taxa
  arrange(Order)                                    # Sort data frame alphabetically by Species

Order_colors <- c(
  "purple", "green", "blue", "red","orange", "yellow","#CBD588", "#5F7FC7", "#DA5724", "#508578", "#CD9BCD",
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248", 
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861", "goldenrod", 
  "mediumpurple1", "lightcoral", "brown")



ggplot(asv_Order, aes(x =Sample, y = Abundance, fill = Order)) + 
  #facet_grid(Diet~.) +
  geom_bar(stat = "identity", position = "Fill") +
  
  # Remove x axis title
  theme(axis.title.x = element_blank()) + 
  #
  guides(fill = guide_legend(reverse = TRUE, keywidth = 1.5, keyheight = 1.5)) +
  ylab("Order Relative Abundance \n") +
  ggtitle("Order Composition of Bacterial Communities Across Different Sample ")

ggplot(asv_Order, aes(x =Treatment, y = Abundance, fill = Order)) + 
  #facet_grid(Diet~.) +
  geom_bar(stat = "identity", position = "Fill") +
  
  # Remove x axis title
  theme(axis.title.x = element_blank()) + 
  #
  guides(fill = guide_legend(reverse = TRUE, keywidth = 1.5, keyheight = 1.5)) +
  ylab("Order Relative Abundance \n") +
  ggtitle("Order Composition of Bacterial Communities Across Different Treatment ")





dr_fam <- tax_glom(ps, taxrank = "Family", NArm = F)
ps0<-names(sort(taxa_sums(dr_fam), TRUE)[1:670]) #get most abundant ones
ps1<-prune_taxa(ps0, dr_fam)

asv_Family <- ps1 %>% 
  tax_glom(taxrank = "Family") %>%                     # agglomerate at Species level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  #filter(Abundance > 0.005) %>%                         # Filter out low abundance taxa
  arrange(Family)                                    # Sort data frame alphabetically by Species

Family_colors <- c(
  "purple", "green", "blue", "red","orange", "yellow","#CBD588", "#5F7FC7", "#DA5724", "#508578", "#CD9BCD",
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248", 
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861", "goldenrod", 
  "mediumpurple1", "lightcoral", "brown")



ggplot(asv_Family, aes(x =Sample, y = Abundance, fill = Family)) + 
  #facet_grid(Diet~.) +
  geom_bar(stat = "identity", position = "Fill") +
  
  # Remove x axis title
  theme(axis.title.x = element_blank()) + 
  #
  guides(fill = guide_legend(reverse = TRUE, keywidth = 1.5, keyheight = 1.5)) +
  ylab("Family Relative Abundance \n") +
  ggtitle("Family Composition of Bacterial Communities Across Different Sample ")

ggplot(asv_Family, aes(x =Treatment, y = Abundance, fill = Family)) + 
  #facet_grid(Diet~.) +
  geom_bar(stat = "identity", position = "Fill") +
  
  # Remove x axis title
  theme(axis.title.x = element_blank()) + 
  #
  guides(fill = guide_legend(reverse = TRUE, keywidth = 1.5, keyheight = 1.5)) +
  ylab("Family Relative Abundance \n") +
  ggtitle("Family Composition of Bacterial Communities Across Different Treatment ")





dr_fam <- tax_glom(ps, taxrank = "Genus", NArm = F)
ps0<-names(sort(taxa_sums(dr_fam), TRUE)[1:670]) #get most abundant ones
ps1<-prune_taxa(ps0, dr_fam)

asv_Genus <- ps1 %>% 
  tax_glom(taxrank = "Genus") %>%                     # agglomerate at Species level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  #filter(Abundance > 0.005) %>%                         # Filter out low abundance taxa
  arrange(Genus)                                    # Sort data frame alphabetically by Species

Genus_colors <- c(
  "purple", "green", "blue", "red","orange", "yellow","#CBD588", "#5F7FC7", "#DA5724", "#508578", "#CD9BCD",
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248", 
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861", "goldenrod", 
  "mediumpurple1", "lightcoral", "brown")



ggplot(asv_Genus, aes(x =Sample, y = Abundance, fill = Genus)) + 
  #facet_grid(Diet~.) +
  geom_bar(stat = "identity", position = "Fill") +
  
  # Remove x axis title
  theme(axis.title.x = element_blank()) + 
  #
  guides(fill = guide_legend(reverse = TRUE, keywidth = 1.5, keyheight = 1.5)) +
  ylab("Genus Relative Abundance \n") +
  ggtitle("Genus Composition of Bacterial Communities Across Different Sample ")

ggplot(asv_Genus, aes(x =Treatment, y = Abundance, fill = Genus)) + 
  #facet_grid(Diet~.) +
  geom_bar(stat = "identity", position = "Fill") +
  
  # Remove x axis title
  theme(axis.title.x = element_blank()) + 
  #
  guides(fill = guide_legend(reverse = TRUE, keywidth = 1.5, keyheight = 1.5)) +
  ylab("Genus Relative Abundance \n") +
  ggtitle("Genus Composition of Bacterial Communities Across Different Treatment ")





library(microbiome)
core_taxa <- core(ps, detection = 0.001, prevalence = 0.8)
# Print a summary of the core_taxa object
print(core_taxa)

# Check the number of ASVs in the core microbiota
ntaxa(core_taxa)

# View the taxa names in the core microbiota
taxa_names(core_taxa)
# Extract the taxonomy table from core_taxa
core_taxa_table <- tax_table(core_taxa)

# Convert to a data frame for easier viewing
core_taxa_df <- as.data.frame(core_taxa_table)

# View the first few rows
head(core_taxa_df)

# Save the taxonomy table to a file
write.csv(core_taxa_df, "D:/Osama/biotec/Bioinformatics/metagenomeics/meta/core_taxa_taxonomy.csv", row.names = TRUE)




library(phyloseq)
library(ggplot2)

# Transform core_taxa to relative abundance
core_taxa_rel <- transform_sample_counts(core_taxa, function(x) x / sum(x))

# Melt the data for plotting
core_taxa_melt <- psmelt(core_taxa_rel)

# Plot a barplot for core taxa at the Genus level
ggplot(core_taxa_melt, aes(x = Sample, y = Abundance, fill = Genus)) +
  geom_bar(stat = "identity", position = "fill") +
  theme(axis.title.x = element_blank()) +
  guides(fill = guide_legend(reverse = TRUE, keywidth = 1.5, keyheight = 1.5)) +
  ylab("Relative Abundance") +
  ggtitle("Core Microbiota Composition (Genus Level)") +
  theme_bw()

ggplot(core_taxa_melt, aes(x = Treatment, y = Abundance, fill = Genus)) +
  geom_bar(stat = "identity", position = "fill") +
  theme(axis.title.x = element_blank()) +
  guides(fill = guide_legend(reverse = TRUE, keywidth = 1.5, keyheight = 1.5)) +
  ylab("Relative Abundance") +
  ggtitle("Core Microbiota Composition (Genus Level)") +
  theme_bw()





ggplot(core_taxa_melt, aes(x = Sample, y = Abundance, fill = Phylum)) +
  geom_bar(stat = "identity", position = "fill") +
  theme(axis.title.x = element_blank()) +
  guides(fill = guide_legend(reverse = TRUE, keywidth = 1.5, keyheight = 1.5)) +
  ylab("Relative Abundance") +
  ggtitle("Core Microbiota Composition (Phylum Level) in Samples") +
  theme_bw()

ggplot(core_taxa_melt, aes(x = Treatment, y = Abundance, fill = Phylum)) +
  geom_bar(stat = "identity", position = "fill") +
  theme(axis.title.x = element_blank()) +
  guides(fill = guide_legend(reverse = TRUE, keywidth = 1.5, keyheight = 1.5)) +
  ylab("Relative Abundance") +
  ggtitle("Core Microbiota Composition (Phylum Level) in Trearment") +
  theme_bw()





ggplot(core_taxa_melt, aes(x = Sample, y = Abundance, fill = Class)) +
  geom_bar(stat = "identity", position = "fill") +
  theme(axis.title.x = element_blank()) +
  guides(fill = guide_legend(reverse = TRUE, keywidth = 1.5, keyheight = 1.5)) +
  ylab("Relative Abundance") +
  ggtitle("Core Microbiota Composition (Class Level) in Samples") +
  theme_bw()

ggplot(core_taxa_melt, aes(x = Treatment, y = Abundance, fill = Class)) +
  geom_bar(stat = "identity", position = "fill") +
  theme(axis.title.x = element_blank()) +
  guides(fill = guide_legend(reverse = TRUE, keywidth = 1.5, keyheight = 1.5)) +
  ylab("Relative Abundance") +
  ggtitle("Core Microbiota Composition (Class Level) in Trearment") +
  theme_bw()





# Plot core taxa by Treatment
ggplot(core_taxa_melt, aes(x = Treatment, y = Abundance, fill = Genus)) +
  geom_bar(stat = "identity", position = "fill") +
  theme(axis.title.x = element_blank()) +
  guides(fill = guide_legend(reverse = TRUE, keywidth = 1.5, keyheight = 1.5)) +
  ylab("Relative Abundance") +
  ggtitle("Core Microbiota by Treatment") +
  theme_bw()


library(DESeq2)
ps_deseq <- phyloseq_to_deseq2(ps, ~ Treatment)
deseq_result <- DESeq(ps_deseq, test = "Wald")
res <- results(deseq_result)
res <- res[!is.na(res$padj), ]  
res_ordered <- res[order(res$padj), ]
sig_res <- res_ordered[res_ordered$padj < 0.01, ]
if (nrow(sig_res) > 0) {
  print(sig_res)
} else {
  print("No significant results found with padj < 0.01")
}
res_df <- as.data.frame(res)
res_df$diffexpressed <- "NO"
res_df$diffexpressed[res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1] <- "YES"
ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = diffexpressed)) +
  geom_point(size = 2) +
  scale_color_manual(values = c("NO" = "gray", "YES" = "red")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  labs(title = "Volcano Plot of Differential Abundance", x = "log2 Fold Change", y = "-log10(padj)") +
  theme_bw()


ps_deseq <- phyloseq_to_deseq2(ps, ~ Treatment)
print(dim(ps_deseq))
print(colData(ps_deseq))
dds <- DESeq(ps_deseq, test = "Wald", fitType = "parametric", betaPrior = FALSE)
plotDispEsts(dds, main = "Dispersion Plot")

ps_deseq <- phyloseq_to_deseq2(ps, ~ Sample)
print(dim(ps_deseq))
print(colData(ps_deseq))
dds <- DESeq(ps_deseq, test = "Wald", fitType = "parametric", betaPrior = FALSE)
plotDispEsts(dds, main = "Dispersion Plot")
