# Analyzing CRP CFG results and quantifying cis sequence features
# Authors: Ruby Froom
# Date: June 11, 2024

#### Libraries -----------------------------------------------------------------

library(tidyverse)
library(stringr)
library(stringi)
library(Biostrings)
library(seqinr)
library(ggseqlogo)
library(DESeq2)
library(EnhancedVolcano) 

## DESeq2 disp ests ------------------------------------------------------------

design_matrix_string = '5enrich_CRP/selectThreshold/DESeq2/noCRP_CRP_Wald_design_matrix.csv'
count_data_string = '5enrich_CRP/selectThreshold/DESeq2/noCRP_CRP_countdata_ends_20counts.csv'
formula_string= '~condition'

# Read in design matrix
design_matrix = as.data.frame(read.table(design_matrix_string, sep=",",
                                         header=TRUE,colClasses = "factor"))
conditions <- levels(design_matrix[['condition']])

# Read in count data
count_data_all <- read.csv(count_data_string)

# Remove first extraneous column from count data
counts_only <- count_data_all[,c(2:ncol(count_data_all))]

# Add in the coordinate ID (location and strand) as row names
rownames(counts_only) <- count_data_all[,1]

# Generate DESeqDataSet
deseq_dataset <- DESeqDataSetFromMatrix(countData = counts_only,
                                        colData = design_matrix,
                                        design = formula(formula_string))

# Estimate size factors based on spike counts from downsampling ----------------
# Note: the total amounts of RNA didn't happen to differ by very much,
# so the size factor differences aren't dramatic
num_spike_coords <- nrow(counts_only %>%
                           filter(grepl("^Eco_",rownames(counts_only))))
spike_row_indices <- ((nrow(counts_only)-num_spike_coords)+1):(nrow(counts_only))
dds_sizeFactors <- estimateSizeFactors(deseq_dataset,
                                       controlGenes = spike_row_indices)
# Run LRT
dds <- DESeq(dds_sizeFactors, parallel=TRUE, fitType = 'local')
vsd <- varianceStabilizingTransformation(dds, fitType = 'local')

# Visually inspect model fits 
# Visually inspect model fits 
png("5enrich_CRP/initiation_figures/extData_fig5e_CRP_disp_ests.png",
    res = 300,
    width = 5,
    height = 5,
    units = 'in')
plotDispEsts(dds)
dev.off()

CRP_DE <- results(dds,
                  contrast=c('condition',
                             'CRP',
                             'noCRP'),
                  test = 'Wald')

CRP_results <- data.frame(CRP_DE)
CRP_results$coordID <- rownames(CRP_DE)
CRP_noEco <- CRP_results[!grepl("Eco",rownames(CRP_results)),]

#### Fig 1D: volcano plot of CRP-DE hits with motif coloring -------------------

CRP_DF <- read.csv('5enrich_CRP/selectThreshold/allMotifResults/noCRP_CRP_min100plus30_20counts_DESeq2.csv')
CRP_DF_toKeep <- CRP_DF[,c(5:ncol(CRP_DF))]
CRP_annotated <- left_join(CRP_noEco, CRP_DF_toKeep,
                           by = c("coordID" = "coordID"))
CRP_sig <- CRP_annotated %>%
  filter(padj <= 0.05)

CRP_DF <- CRP_annotated

strong_motif_color = '#0075B6'
weak_motif_color = '#54B3E2'
no_motif_color = '#D3D3D3'
black = '#000000'

# Assign colors based on motif strength
keyvals <- ifelse(
  (CRP_DF$padj > 0.05), black,
  ifelse(
    (CRP_DF$score >= median(na.omit(CRP_DF$score))), strong_motif_color,
    ifelse(
      (CRP_DF$score < median(na.omit(CRP_DF$score))), weak_motif_color,
        no_motif_color)))
keyvals[is.na(keyvals)] <- 'No motif'
names(keyvals)[keyvals == weak_motif_color] <- 'Weak'
names(keyvals)[keyvals == no_motif_color] <- 'No motif'
names(keyvals)[keyvals == strong_motif_color] <- 'Strong'
names(keyvals)[keyvals == black] <- 'Not sig'

CRP_volcano <- EnhancedVolcano(CRP_DF,
                lab = CRP_DF$coordID,
                x = 'log2FoldChange',
                y = 'padj',
                pCutoff = 0.05,
                colAlpha = 1,
                labSize = 0,
                title = NULL,
                subtitle = NULL,
                xlim = c(-6,6),
                xlab = NULL,
                ylab = NULL,
                legendPosition = "none",
                gridlines.minor = FALSE,
                colCustom = keyvals,
                pointSize =  2) +
  theme_minimal() +
  theme(legend.position = "none")

ggsave(filename = '5enrich_CRP/initiation_figures/fig3c_crp_volcano.png',
       plot = CRP_volcano,
       width = 101.6,
       height = 80,
       units = 'mm',
       dpi = 300)

CRP_activated <- CRP_DF[CRP_DF$padj <= 0.05 & CRP_DF$log2FoldChange > 0,]
CRP_repressed <- CRP_DF[CRP_DF$padj <= 0.05 & CRP_DF$log2FoldChange < 0,]

#nrow(CRP_activated[CRP_activated$score != 'No motif',])
nrow(CRP_activated[is.na(CRP_activated$score),])
nrow(CRP_activated)

#nrow(CRP_repressed[CRP_repressed$score != 'No motif',])
nrow(CRP_repressed[is.na(CRP_repressed$score),])
nrow(CRP_repressed)

# Plot of correlation between motif strength and CRP L2FC ----------------------

CRP_DF_scatterplot <- data.frame(CRP_DF) %>%
  mutate(CRP_DF,
         motif_strength = case_when(score > median(na.omit(CRP_DF$score)) ~ "strong",
                                    score < median(na.omit(CRP_DF$score)) ~ "weak",
                                    is.na(score) ~ "none"))

CRP_DF_scatterplot_noNA <- filter(CRP_DF_scatterplot,
                                  motif_strength != "none")

scatterplot <- ggplot(CRP_DF_scatterplot_noNA, aes(x = score, y = abs(log2FoldChange))) +
  geom_point(data = subset(CRP_DF_scatterplot_noNA, motif_strength == 'strong'),
             size = 1.5,
             color = strong_motif_color) +
  geom_point(data = subset(CRP_DF_scatterplot_noNA, motif_strength == 'weak'),
             size = 1.5,
             color = weak_motif_color) +
  theme_minimal() +
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE, color = 'black',
              size = 0.6)

cor.test(CRP_DF_scatterplot_noNA$score, abs(CRP_DF_scatterplot_noNA$log2FoldChange))

ggsave(filename = '5enrich_CRP/initiation_figures/fig3d_motifSim_vs_effectSize_scatterplot.png',
       plot = scatterplot,
       width = 80,
       height = 80,
       units = 'mm',
       dpi = 300)

#### Distributions of motif spacing in CRP-DE hits ---------------------

# Add coordinate and direction by splitting up coordID
CRP_DF$coordinate <- as.integer(stri_sub(CRP_DF$coordID, from = 1, to = -2))
CRP_DF$direction <- sapply(strsplit(CRP_DF$coordID,""), tail, 1)

# Filter for promoters differentially-regulated by CRP 
CRP_sig <- CRP_DF[CRP_DF['padj'] <= 0.05,]

CRP_sig$CRP_effect <- "activated"
CRP_sig[CRP_sig['log2FoldChange'] < 0,"CRP_effect"] <- "repressed"

CRP_repressed_col <- '#CC6677'
CRP_activated_col <- '#117733'

motif_hist <- ggplot(CRP_sig, aes(x=motif_center)) +
  geom_histogram(data = subset(CRP_sig, CRP_effect=='activated'),
                 bins = 21,
                 fill = CRP_activated_col,
                 alpha = 0.8) +
  geom_histogram(data = subset(CRP_sig, CRP_effect=='repressed'),
                 bins = 21,
                 fill = CRP_repressed_col,
                 alpha = 0.8) +
  scale_x_continuous(breaks= seq(-85,20,5),
                     lim = c(-85,20)) +
  theme_minimal() +
  ylim(0,8)+
  theme(panel.grid.minor = element_blank())

motif_density <- ggplot(CRP_sig,
       aes(x = motif_center)) +
  geom_density(data = subset(CRP_sig, CRP_effect=='activated'),
               alpha = 0.5, fill = CRP_activated_col, size = 0.3) +
  geom_density(data = subset(CRP_sig, CRP_effect=='repressed'),
               alpha = 0.5, fill = CRP_repressed_col, size = 0.3) +
  theme_minimal() +
  scale_x_continuous(limits = c(-100,40), breaks = seq(-100, 40, 20)) +
  theme(panel.grid.minor = element_blank())

ggsave(filename = '5enrich_CRP/initiation_figures/extData_fig5f_motif_hist.png',
       plot = motif_hist,
       width = 200,
       height = 40,
       units = 'mm',
       dpi = 300)

ggsave(filename = '5enrich_CRP/initiation_figures/fig3g_motif_density.png',
       plot = motif_density,
       width = 83,
       height = 40,
       units = 'mm',
       dpi = 300)

#### Raw data for top CRP-activated and -repressed hits ------------------------
library(Gviz)
library(rtracklayer)

options(ucscChromosomeNames = FALSE)

genome <- readDNAStringSet('5enrich_CRP/initiation_figures/genomic_info/Eco_Mtb_genome.fasta')
sTrack <- SequenceTrack(genome)
sTrack_comp <- SequenceTrack(complement(genome))

# Identify top activated hit
most_sig_activated <- data.frame(CRP_sig) %>%
  filter(CRP_effect == 'activated') %>%
  filter(padj == min(padj))

# If the most activated hit is in the + orientation (yes, in this case)
if (most_sig_activated$direction == '+') {
  
  # Import positive-stranded raw data (bigWigs)
  noCRP1_plus_bw <- import.bw("5enrich_CRP/identifyEnrichedEnds/bigWigs_vis/noCRP1_plus_vis.bw")
  noCRP2_plus_bw <- import.bw("5enrich_CRP/identifyEnrichedEnds/bigWigs_vis/noCRP2_plus_vis.bw")
  noCRP3_plus_bw <- import.bw("5enrich_CRP/identifyEnrichedEnds/bigWigs_vis/noCRP3_plus_vis.bw")
  CRP1_plus_bw <- import.bw("5enrich_CRP/identifyEnrichedEnds/bigWigs_vis/CRP1_plus_vis.bw")
  CRP2_plus_bw <- import.bw("5enrich_CRP/identifyEnrichedEnds/bigWigs_vis/CRP2_plus_vis.bw")
  CRP3_plus_bw <- import.bw("5enrich_CRP/identifyEnrichedEnds/bigWigs_vis/CRP3_plus_vis.bw")
  
  # For activated hits, expect CRP motif to be upstream,
  # so include a long upstream region
  coord = most_sig_activated$coordinate
  start_coord = coord - 73
  end_coord = coord + 7
  alpha_val = 0.7
  noCRP_col = "#ED7121"
  CRP_col = "#0075B8"
  
  # Generate DataTracks from each individual replicate
  noCRP1_DT <- DataTrack(noCRP1_plus_bw, start = start_coord, end = end_coord,
                          name = "noCRP", groups = factor("noCRP",
                                                           levels = c("noCRP","CRP")),
                          alpha = alpha_val, col = noCRP_col)
  noCRP2_DT <- DataTrack(noCRP2_plus_bw, start = start_coord, end = end_coord,
                          name = "noCRP", groups = factor("noCRP",
                                                           levels = c("noCRP","CRP")),
                          alpha = alpha_val, col = noCRP_col)
  noCRP3_DT <- DataTrack(noCRP3_plus_bw, start = start_coord, end = end_coord,
                          name = "noCRP", groups = factor("noCRP",
                                                           levels = c("noCRP","CRP")),
                          alpha = alpha_val, col = noCRP_col)
  CRP1_DT <- DataTrack(CRP1_plus_bw, start = start_coord, end = end_coord,
                       name = "CRP", groups = factor("CRP",
                                                     levels = c("noCRP","CRP")),
                       alpha = alpha_val, col = CRP_col)
  CRP2_DT <- DataTrack(CRP2_plus_bw, start = start_coord, end = end_coord,
                       name = "CRP", groups = factor("CRP",
                                                     levels = c("noCRP","CRP")),
                       alpha = alpha_val, col = CRP_col)
  CRP3_DT <- DataTrack(CRP3_plus_bw, start = start_coord, end = end_coord,
                       name = "CRP", groups = factor("CRP",
                                                     levels = c("noCRP","CRP")),
                       alpha = alpha_val, col = CRP_col)
  
  ## Overlay each replicate, but continue to separate by condition
  noCRP_overlay <- OverlayTrack(trackList = list(noCRP3_DT, noCRP1_DT,
                                                  noCRP2_DT))
  CRP_overlay <- OverlayTrack(trackList = list(CRP1_DT, CRP2_DT, CRP3_DT))
  
  colForLetters <- c(A = "black", C = "black",
                     T = "black", G = "black")
  
  # Note: for now, I cannot find a way to automatically scale the y-axis.
  # Need to manually adjust to observe the max height for each replicate.
  # Note that the y-axis values are 1000-fold higher than observed, because
  # I needed to incorporate DESeq2 size factors when scaling the bigWigs-  
  # so the relative scaling between samples is appropriate (based on spike),
  # but all values are 1000-fold less than the axes here show.
  # Good rule of thumb is to start around 10^6 for max_y
  max_y <- 1500000
  
  # In figure, final y axis values are calculated by dividing by 1000
  # then converting from counts to counts per million 
  # (dividing by 7.18mill reads)
  # Resize in Illustrator to align with sequence in 5pt Courier font:
  # W = 98.565mm
  # H = 57.496mm
  png("5enrich_CRP/initiation_figures/fig3f_CFG_activated_rawData.png",
      res = 300,
      width = 12,
      height = 12,
      units = 'in')
  plotTracks(list(noCRP_overlay, CRP_overlay, sTrack),
             from = start_coord,
             to = end_coord,
             type = c("hist"),
             featureAnnotation = "feature",
             legend = FALSE,
             col.axis = "black",
             background.title = "transparent",
             fontcolor.item = "black",
             fontcolor = colForLetters,
             ylim = c(0,max_y))
  dev.off()
}

# Identify top repressed hit
most_sig_repressed <- data.frame(CRP_sig) %>%
  filter(CRP_effect == 'repressed') %>%
  filter(padj == min(padj))

# If the most activated hit is in the - orientation (yes, in this case)
if (most_sig_repressed$direction == '-') {
  
  # Import positive-stranded raw data (bigWigs)
  noCRP1_minus_bw <- import.bw("5enrich_CRP/identifyEnrichedEnds/bigWigs_vis/noCRP1_minus_vis.bw")
  noCRP2_minus_bw <- import.bw("5enrich_CRP/identifyEnrichedEnds/bigWigs_vis/noCRP2_minus_vis.bw")
  noCRP3_minus_bw <- import.bw("5enrich_CRP/identifyEnrichedEnds/bigWigs_vis/noCRP3_minus_vis.bw")
  CRP1_minus_bw <- import.bw("5enrich_CRP/identifyEnrichedEnds/bigWigs_vis/CRP1_minus_vis.bw")
  CRP2_minus_bw <- import.bw("5enrich_CRP/identifyEnrichedEnds/bigWigs_vis/CRP2_minus_vis.bw")
  CRP3_minus_bw <- import.bw("5enrich_CRP/identifyEnrichedEnds/bigWigs_vis/CRP3_minus_vis.bw")
  
  # For repressed hits, expect CRP motif to be within the promoter,
  # so don't need as long an upstream section, but may need more downstream
  coord = most_sig_repressed$coordinate
  start_coord = coord - 25
  end_coord = coord + 55
  alpha_val = 0.7
  noCRP_col = "#ED7121"
  CRP_col = "#0075B8"
  
  # Generate DataTracks from each individual replicate
  noCRP1_DT <- DataTrack(noCRP1_minus_bw, start = start_coord, end = end_coord,
                          name = "noCRP", groups = factor("noCRP",
                                                           levels = c("noCRP","CRP")),
                          alpha = alpha_val, col = noCRP_col)
  noCRP2_DT <- DataTrack(noCRP2_minus_bw, start = start_coord, end = end_coord,
                          name = "noCRP", groups = factor("noCRP",
                                                           levels = c("noCRP","CRP")),
                          alpha = alpha_val, col = noCRP_col)
  noCRP3_DT <- DataTrack(noCRP3_minus_bw, start = start_coord, end = end_coord,
                          name = "noCRP", groups = factor("noCRP",
                                                           levels = c("noCRP","CRP")),
                          alpha = alpha_val, col = noCRP_col)
  CRP1_DT <- DataTrack(CRP1_minus_bw, start = start_coord, end = end_coord,
                       name = "CRP", groups = factor("CRP",
                                                     levels = c("noCRP","CRP")),
                       alpha = alpha_val, col = CRP_col)
  CRP2_DT <- DataTrack(CRP2_minus_bw, start = start_coord, end = end_coord,
                       name = "CRP", groups = factor("CRP",
                                                     levels = c("noCRP","CRP")),
                       alpha = alpha_val, col = CRP_col)
  CRP3_DT <- DataTrack(CRP3_minus_bw, start = start_coord, end = end_coord,
                       name = "CRP", groups = factor("CRP",
                                                     levels = c("noCRP","CRP")),
                       alpha = alpha_val, col = CRP_col)
  
  ## Overlay each replicate, but continue to separate by condition
  noCRP_overlay <- OverlayTrack(trackList = list(noCRP3_DT, noCRP1_DT,
                                                  noCRP2_DT))
  CRP_overlay <- OverlayTrack(trackList = list(CRP1_DT, CRP2_DT, CRP3_DT))
  
  colForLetters <- c(A = "black", C = "black",
                     T = "black", G = "black")
  
  # Note: for now, I cannot find a way to automatically scale the y-axis.
  # Need to manually adjust to observe the max height for each replicate.
  # Note that the y-axis values are 1000-fold higher than observed, because
  # I needed to incorporate DESeq2 size factors when scaling the bigWigs-  
  # so the relative scaling between samples is appropriate (based on spike),
  # but all values are 1000-fold less than the axes here show.
  # Good rule of thumb is to start around 10^6 for max_y
  max_y <- 1000000
  
  # In figure, final y axis values are calculated by dividing by 1000
  # then converting from counts to counts per million 
  # (dividing by 7.18mill reads)
  # Resize in Illustrator to align with sequence in 5pt Courier font:
  # W = 75.382mm
  # H = 75.382mm
  png("5enrich_CRP/initiation_figures/fig3e_CFG_repressed_rawData.png",
      res = 300,
      width = 12,
      height = 12,
      units = 'in')
  plotTracks(list(noCRP_overlay, CRP_overlay, sTrack_comp),
             from = start_coord,
             to = end_coord,
             type = c("hist"),
             featureAnnotation = "feature",
             legend = FALSE,
             col.axis = "black",
             background.title = "transparent",
             fontcolor.item = "black",
             fontcolor = colForLetters,
             ylim = c(0,max_y))
  dev.off()
}

## End -------------------------------------------------------------------------

