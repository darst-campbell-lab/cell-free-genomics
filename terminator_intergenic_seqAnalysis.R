#### Analyzing sequence features of intergenic terminators from cell-free genomics
#### Authors: Ruby Froom
#### Date: June 3, 2024

#### Libraries ----------------------------------------------------------------

library(dplyr)
library(DESeq2)
library(viridis)
library(RColorBrewer)
library(stringr)
library(stringi)
library(Biostrings)
library(seqinr)
library(universalmotif)
library(ggseqlogo)
library(GenomicRanges)
library(seqinr)
library(Biostrings)
library(DNAshapeR)
library(ggplot2)

## Import intergenic terminators -----------------------------------------------

intergenic_allCols <- read.csv('3enrich_NusAG/termination_figures/terminator_annotations/Mtb_annotated.csv')

## dsDNA, U tract and hairpin features of intergenic Nus-stim. candidates ------
## from each cluster -----------------------------------------------------------

# Clean up DF
colsToKeep <- c('Nontemplate_strand','coordID','cluster','coordinate',
                'strand','RNAfold_hairpin_region_44_deltaG',
                'RNAfold_hairpin_region_44_spacerLength',
                'RNAfold_hairpin_region_44_stemLength',
                'RNAfold_hairpin_region_44_loopLength',
                'U_count_tractRegion','NusA_log2FC','NusG_log2FC','NusAG_log2FC',
                'NusA_padj','NusG_padj','NusAG_padj',
                'NusG_stim_cluster','NusA_stim_cluster',
                'AT_count_downstream')

intergenic <- intergenic_allCols[,colsToKeep]

colnames(intergenic) <- c('Nontemplate_strand','coordID','cluster','coordinate',
                          'strand','deltaG',
                          'spacerLength',
                          'stemLength',
                          'loopLength',
                          'Utract_count','NusA_log2FC','NusG_log2FC','NusAG_log2FC',
                          'NusA_padj','NusG_padj','NusAG_padj',
                          'NusG_stim','NusA_stim',
                          'AT_dsDNA')

intergenic_noRandom <- intergenic %>%
  filter(!cluster == 'random') %>%
  mutate(cluster = factor(cluster,
                          levels = c("NE",'A',"B","C","D")))

random <- intergenic %>%
  filter(cluster == 'random')
NE <- intergenic %>%
  filter(cluster == 'NE')
A <- intergenic %>%
  filter(cluster == 'A')
B <- intergenic %>%
  filter(cluster == 'B')
C <- intergenic %>%
  filter(cluster == 'C')
D <- intergenic %>%
  filter(cluster == 'D')

random_col <- 'white'
nonstim_col <- 'grey'
clusterA_col <- "#FFA904"
NusA_col <- '#5F9AEC'
NusG_col <- '#F26B7A'
NusAG_col <- '#BD67F3'

intergenic_effectOnly <- intergenic_noRandom %>%
  filter(cluster != "NE")

summary_DF <- data.frame(cluster = c("NE",'A','B',
                                     "C",'D'),
                         spacer_medians = c(median(NE$spacerLength),
                                            median(A$spacerLength),
                                            median(B$spacerLength),
                                            median(C$spacerLength),
                                            median(D$spacerLength)),
                         deltaG_medians = c(median(NE$deltaG),
                                            median(A$deltaG),
                                            median(B$deltaG),
                                            median(C$deltaG),
                                            median(D$deltaG)),
                         AT_medians = c(median(NE$AT_dsDNA),
                                        median(A$AT_dsDNA),
                                        median(B$AT_dsDNA),
                                        median(C$AT_dsDNA),
                                        median(D$AT_dsDNA)),
                         spacer_means = c(mean(NE$spacerLength),
                                          mean(A$spacerLength),
                                          mean(B$spacerLength),
                                          mean(C$spacerLength),
                                          mean(D$spacerLength)),
                         spacer_SD = c(sd(NE$spacerLength),
                                       sd(A$spacerLength),
                                       sd(B$spacerLength),
                                       sd(C$spacerLength),
                                       sd(D$spacerLength)),
                         deltaG_means = c(mean(NE$deltaG),
                                          mean(A$deltaG),
                                          mean(B$deltaG),
                                          mean(C$deltaG),
                                          mean(D$deltaG)),
                         deltaG_SD = c(sd(NE$deltaG),
                                       sd(A$deltaG),
                                       sd(B$deltaG),
                                       sd(C$deltaG),
                                       sd(D$deltaG)),
                         AT_means = c(mean(NE$AT_dsDNA),
                                      mean(A$AT_dsDNA),
                                      mean(B$AT_dsDNA),
                                      mean(C$AT_dsDNA),
                                      mean(D$AT_dsDNA)),
                         AT_SD = c(sd(NE$AT_dsDNA),
                                       sd(A$AT_dsDNA),
                                       sd(B$AT_dsDNA),
                                       sd(C$AT_dsDNA),
                                       sd(D$AT_dsDNA)),
                         num_terms = c(nrow(NE),nrow(A),
                                       nrow(B),
                                       nrow(C),
                                       nrow(D)))

summary_DF$deltaG_mean_minSD <- summary_DF$deltaG_means - summary_DF$deltaG_SD
summary_DF$deltaG_mean_plusSD <- summary_DF$deltaG_means + summary_DF$deltaG_SD
summary_DF$spacer_mean_minSD <- summary_DF$spacer_means - summary_DF$spacer_SD
summary_DF$spacer_mean_plusSD <- summary_DF$spacer_means + summary_DF$spacer_SD
summary_DF$AT_mean_minSD <- summary_DF$AT_means - summary_DF$AT_SD
summary_DF$AT_mean_plusSD <- summary_DF$AT_means + summary_DF$AT_SD

write.csv(summary_DF, '3enrich_NusAG/termination_figures/intergenic_summary_DF.csv')

## Distributions (NusG stim vs. not, NusA stim vs not) -------------------------

NusA_stim <- intergenic %>%
  filter(NusA_stim == 'Yes')
NusA_nonstim <- intergenic %>%
  filter(NusA_stim == 'No')
NusG_stim <- intergenic %>%
  filter(NusG_stim == 'Yes')
NusG_nonstim <- intergenic %>%
  filter(NusG_stim == 'No')
random <- intergenic %>%
  filter(NusG_stim == 'random')

## Distribution plots ----------------------------------------------------------

# Hairpin deltaG
NusA_deltaG <- ggplot(intergenic_noRandom, aes(x = NusA_stim, y = deltaG,
                                               fill = NusA_stim)) + 
  geom_violin(linewidth = 2) +
  geom_boxplot(width = 0.2,
               linewidth = 2) +
  scale_y_continuous(limits = c(-40,0),
                     breaks = seq(-40,0,10)) +
  scale_fill_manual(values = c(nonstim_col,NusA_col)) +
  theme_minimal() +
  theme(legend.position = 'none')
ggsave(filename = '3enrich_NusAG/termination_figures/fig5b_deltaG_NusA_boxplot.png',
       plot = NusA_deltaG,
       width = 200,
       height = 200,
       units = 'mm',
       dpi = 300)
NusA_deltaG_pval <- wilcox.test(NusA_nonstim$deltaG,
                                NusA_stim$deltaG)

NusG_deltaG <- ggplot(intergenic_noRandom, aes(x = NusG_stim, y = deltaG,
                                               fill = NusG_stim)) + 
  geom_violin(linewidth = 2) +
  geom_boxplot(width = 0.2,
               linewidth = 2) +
  scale_y_continuous(limits = c(-40,0),
                     breaks = seq(-40,0,10)) +
  scale_fill_manual(values = c(nonstim_col,NusG_col)) +
  theme_minimal() +
  theme(legend.position = 'none')
ggsave(filename = '3enrich_NusAG/termination_figures/fig5b_deltaG_NusG_boxplot.png',
       plot = NusG_deltaG,
       width = 200,
       height = 200,
       units = 'mm',
       dpi = 300)
NusG_deltaG_pval <- wilcox.test(NusG_nonstim$deltaG,
                                NusG_stim$deltaG)

# Downstream AT counts
NusA_downstreamATs <- ggplot(intergenic_noRandom, aes(x = NusA_stim,
                                                      y = AT_dsDNA,
                                                      fill = NusA_stim)) + 
  geom_violin(linewidth = 2) +
  geom_boxplot(width = 0.2,
               linewidth = 2) +
  scale_y_continuous(limits = c(0,10),
                     breaks = seq(0,10,2)) +
  scale_fill_manual(values = c(nonstim_col,NusA_col)) +
  theme_minimal() +
  theme(legend.position = 'none')
ggsave(filename = '3enrich_NusAG/termination_figures/fig5d_downstreamATs_NusA_boxplot.png',
       plot = NusA_downstreamATs,
       width = 200,
       height = 200,
       units = 'mm',
       dpi = 300)
NusA_dsDNA_pval <- wilcox.test(NusA_nonstim$AT_dsDNA,
                               NusA_stim$AT_dsDNA)

NusG_downstreamATs <- ggplot(intergenic_noRandom, aes(x = NusG_stim,
                                                      y = AT_dsDNA,
                                                      fill = NusG_stim)) + 
  geom_violin(linewidth = 2) +
  geom_boxplot(width = 0.2,
               linewidth = 2) +
  scale_y_continuous(limits = c(0,10),
                     breaks = seq(0,10,2)) +
  scale_fill_manual(values = c(nonstim_col, NusG_col)) +
  theme_minimal() +
  theme(legend.position = 'none')
ggsave(filename = '3enrich_NusAG/termination_figures/fig5d_downstreamATs_NusG_boxplot.png',
       plot = NusG_downstreamATs,
       width = 200,
       height = 200,
       units = 'mm',
       dpi = 300)
NusG_dsDNA_pval <- wilcox.test(NusG_nonstim$AT_dsDNA,
                               NusG_stim$AT_dsDNA)

# Hairpin spacer lengths
NusA_spacerLength_pval <- wilcox.test(NusA_nonstim$spacerLength,
                                      NusA_stim$spacerLength)

NusG_spacerLength_pval <- wilcox.test(NusG_nonstim$spacerLength,
                                      NusG_stim$spacerLength)

spacerLength_percentage_DF <- data.frame(spacer_length = seq(0,10,1),
                              NusG_stim_spacer_percentage = c(nrow(intergenic_noRandom[(intergenic_noRandom$spacerLength == 0) &
                                                                                         (intergenic_noRandom$NusG_stim == 'Yes'),]) /
                                                                nrow(intergenic_noRandom[(intergenic_noRandom$NusG_stim == 'Yes'),]),
                                                              nrow(intergenic_noRandom[(intergenic_noRandom$spacerLength == 1) &
                                                                                         (intergenic_noRandom$NusG_stim == 'Yes'),]) /
                                                                nrow(intergenic_noRandom[(intergenic_noRandom$NusG_stim == 'Yes'),]),
                                                              nrow(intergenic_noRandom[(intergenic_noRandom$spacerLength == 2) &
                                                                                         (intergenic_noRandom$NusG_stim == 'Yes'),]) /
                                                                nrow(intergenic_noRandom[(intergenic_noRandom$NusG_stim == 'Yes'),]),
                                                              nrow(intergenic_noRandom[(intergenic_noRandom$spacerLength == 3) &
                                                                                         (intergenic_noRandom$NusG_stim == 'Yes'),]) /
                                                                nrow(intergenic_noRandom[(intergenic_noRandom$NusG_stim == 'Yes'),]),
                                                              nrow(intergenic_noRandom[(intergenic_noRandom$spacerLength == 4) &
                                                                                         (intergenic_noRandom$NusG_stim == 'Yes'),]) /
                                                                nrow(intergenic_noRandom[(intergenic_noRandom$NusG_stim == 'Yes'),]),
                                                              nrow(intergenic_noRandom[(intergenic_noRandom$spacerLength == 5) &
                                                                                         (intergenic_noRandom$NusG_stim == 'Yes'),]) /
                                                                nrow(intergenic_noRandom[(intergenic_noRandom$NusG_stim == 'Yes'),]),
                                                              nrow(intergenic_noRandom[(intergenic_noRandom$spacerLength == 6) &
                                                                                         (intergenic_noRandom$NusG_stim == 'Yes'),]) /
                                                                nrow(intergenic_noRandom[(intergenic_noRandom$NusG_stim == 'Yes'),]),
                                                              nrow(intergenic_noRandom[(intergenic_noRandom$spacerLength == 7) &
                                                                                         (intergenic_noRandom$NusG_stim == 'Yes'),]) /
                                                                nrow(intergenic_noRandom[(intergenic_noRandom$NusG_stim == 'Yes'),]),
                                                              nrow(intergenic_noRandom[(intergenic_noRandom$spacerLength == 8) &
                                                                                         (intergenic_noRandom$NusG_stim == 'Yes'),]) /
                                                                nrow(intergenic_noRandom[(intergenic_noRandom$NusG_stim == 'Yes'),]),
                                                              nrow(intergenic_noRandom[(intergenic_noRandom$spacerLength == 9) &
                                                                                         (intergenic_noRandom$NusG_stim == 'Yes'),]) /
                                                                nrow(intergenic_noRandom[(intergenic_noRandom$NusG_stim == 'Yes'),]),
                                                              nrow(intergenic_noRandom[(intergenic_noRandom$spacerLength == 10) &
                                                                                         (intergenic_noRandom$NusG_stim == 'Yes'),]) /
                                                                nrow(intergenic_noRandom[(intergenic_noRandom$NusG_stim == 'Yes'),])),
                              NusA_stim_spacer_percentage = c(nrow(intergenic_noRandom[(intergenic_noRandom$spacerLength == 0) &
                                                                                         (intergenic_noRandom$NusA_stim == 'Yes'),]) /
                                                                nrow(intergenic_noRandom[(intergenic_noRandom$NusA_stim == 'Yes'),]),
                                                              nrow(intergenic_noRandom[(intergenic_noRandom$spacerLength == 1) &
                                                                                         (intergenic_noRandom$NusA_stim == 'Yes'),]) /
                                                                nrow(intergenic_noRandom[(intergenic_noRandom$NusA_stim == 'Yes'),]),
                                                              nrow(intergenic_noRandom[(intergenic_noRandom$spacerLength == 2) &
                                                                                         (intergenic_noRandom$NusA_stim == 'Yes'),]) /
                                                                nrow(intergenic_noRandom[(intergenic_noRandom$NusA_stim == 'Yes'),]),
                                                              nrow(intergenic_noRandom[(intergenic_noRandom$spacerLength == 3) &
                                                                                         (intergenic_noRandom$NusA_stim == 'Yes'),]) /
                                                                nrow(intergenic_noRandom[(intergenic_noRandom$NusA_stim == 'Yes'),]),
                                                              nrow(intergenic_noRandom[(intergenic_noRandom$spacerLength == 4) &
                                                                                         (intergenic_noRandom$NusA_stim == 'Yes'),]) /
                                                                nrow(intergenic_noRandom[(intergenic_noRandom$NusA_stim == 'Yes'),]),
                                                              nrow(intergenic_noRandom[(intergenic_noRandom$spacerLength == 5) &
                                                                                         (intergenic_noRandom$NusA_stim == 'Yes'),]) /
                                                                nrow(intergenic_noRandom[(intergenic_noRandom$NusA_stim == 'Yes'),]),
                                                              nrow(intergenic_noRandom[(intergenic_noRandom$spacerLength == 6) &
                                                                                         (intergenic_noRandom$NusA_stim == 'Yes'),]) /
                                                                nrow(intergenic_noRandom[(intergenic_noRandom$NusA_stim == 'Yes'),]),
                                                              nrow(intergenic_noRandom[(intergenic_noRandom$spacerLength == 7) &
                                                                                         (intergenic_noRandom$NusA_stim == 'Yes'),]) /
                                                                nrow(intergenic_noRandom[(intergenic_noRandom$NusA_stim == 'Yes'),]),
                                                              nrow(intergenic_noRandom[(intergenic_noRandom$spacerLength == 8) &
                                                                                         (intergenic_noRandom$NusA_stim == 'Yes'),]) /
                                                                nrow(intergenic_noRandom[(intergenic_noRandom$NusA_stim == 'Yes'),]),
                                                              nrow(intergenic_noRandom[(intergenic_noRandom$spacerLength == 9) &
                                                                                         (intergenic_noRandom$NusA_stim == 'Yes'),]) /
                                                                nrow(intergenic_noRandom[(intergenic_noRandom$NusA_stim == 'Yes'),]),
                                                              nrow(intergenic_noRandom[(intergenic_noRandom$spacerLength == 10) &
                                                                                         (intergenic_noRandom$NusA_stim == 'Yes'),]) /
                                                                nrow(intergenic_noRandom[(intergenic_noRandom$NusA_stim == 'Yes'),])),
                              nonNusG_stim_spacer_percentage = c(nrow(intergenic_noRandom[(intergenic_noRandom$spacerLength == 0) &
                                                                                            (intergenic_noRandom$NusG_stim == 'No'),]) /
                                                                   nrow(intergenic_noRandom[(intergenic_noRandom$NusG_stim == 'No'),]),
                                                                 nrow(intergenic_noRandom[(intergenic_noRandom$spacerLength == 1) &
                                                                                            (intergenic_noRandom$NusG_stim == 'No'),]) /
                                                                   nrow(intergenic_noRandom[(intergenic_noRandom$NusG_stim == 'No'),]),
                                                                 nrow(intergenic_noRandom[(intergenic_noRandom$spacerLength == 2) &
                                                                                            (intergenic_noRandom$NusG_stim == 'No'),]) /
                                                                   nrow(intergenic_noRandom[(intergenic_noRandom$NusG_stim == 'No'),]),
                                                                 nrow(intergenic_noRandom[(intergenic_noRandom$spacerLength == 3) &
                                                                                            (intergenic_noRandom$NusG_stim == 'No'),]) /
                                                                   nrow(intergenic_noRandom[(intergenic_noRandom$NusG_stim == 'No'),]),
                                                                 nrow(intergenic_noRandom[(intergenic_noRandom$spacerLength == 4) &
                                                                                            (intergenic_noRandom$NusG_stim == 'No'),]) /
                                                                   nrow(intergenic_noRandom[(intergenic_noRandom$NusG_stim == 'No'),]),
                                                                 nrow(intergenic_noRandom[(intergenic_noRandom$spacerLength == 5) &
                                                                                            (intergenic_noRandom$NusG_stim == 'No'),]) /
                                                                   nrow(intergenic_noRandom[(intergenic_noRandom$NusG_stim == 'No'),]),
                                                                 nrow(intergenic_noRandom[(intergenic_noRandom$spacerLength == 6) &
                                                                                            (intergenic_noRandom$NusG_stim == 'No'),]) /
                                                                   nrow(intergenic_noRandom[(intergenic_noRandom$NusG_stim == 'No'),]),
                                                                 nrow(intergenic_noRandom[(intergenic_noRandom$spacerLength == 7) &
                                                                                            (intergenic_noRandom$NusG_stim == 'No'),]) /
                                                                   nrow(intergenic_noRandom[(intergenic_noRandom$NusG_stim == 'No'),]),
                                                                 nrow(intergenic_noRandom[(intergenic_noRandom$spacerLength == 8) &
                                                                                            (intergenic_noRandom$NusG_stim == 'No'),]) /
                                                                   nrow(intergenic_noRandom[(intergenic_noRandom$NusG_stim == 'No'),]),
                                                                 nrow(intergenic_noRandom[(intergenic_noRandom$spacerLength == 9) &
                                                                                            (intergenic_noRandom$NusG_stim == 'No'),]) /
                                                                   nrow(intergenic_noRandom[(intergenic_noRandom$NusG_stim == 'No'),]),
                                                                 nrow(intergenic_noRandom[(intergenic_noRandom$spacerLength == 10) &
                                                                                            (intergenic_noRandom$NusG_stim == 'No'),]) /
                                                                   nrow(intergenic_noRandom[(intergenic_noRandom$NusG_stim == 'No'),])),
                              nonNusA_stim_spacer_percentage = c(nrow(intergenic_noRandom[(intergenic_noRandom$spacerLength == 0) &
                                                                                            (intergenic_noRandom$NusA_stim == 'No'),]) /
                                                                   nrow(intergenic_noRandom[(intergenic_noRandom$NusA_stim == 'No'),]),
                                                                 nrow(intergenic_noRandom[(intergenic_noRandom$spacerLength == 1) &
                                                                                            (intergenic_noRandom$NusA_stim == 'No'),]) /
                                                                   nrow(intergenic_noRandom[(intergenic_noRandom$NusA_stim == 'No'),]),
                                                                 nrow(intergenic_noRandom[(intergenic_noRandom$spacerLength == 2) &
                                                                                            (intergenic_noRandom$NusA_stim == 'No'),]) /
                                                                   nrow(intergenic_noRandom[(intergenic_noRandom$NusA_stim == 'No'),]),
                                                                 nrow(intergenic_noRandom[(intergenic_noRandom$spacerLength == 3) &
                                                                                            (intergenic_noRandom$NusA_stim == 'No'),]) /
                                                                   nrow(intergenic_noRandom[(intergenic_noRandom$NusA_stim == 'No'),]),
                                                                 nrow(intergenic_noRandom[(intergenic_noRandom$spacerLength == 4) &
                                                                                            (intergenic_noRandom$NusA_stim == 'No'),]) /
                                                                   nrow(intergenic_noRandom[(intergenic_noRandom$NusA_stim == 'No'),]),
                                                                 nrow(intergenic_noRandom[(intergenic_noRandom$spacerLength == 5) &
                                                                                            (intergenic_noRandom$NusA_stim == 'No'),]) /
                                                                   nrow(intergenic_noRandom[(intergenic_noRandom$NusA_stim == 'No'),]),
                                                                 nrow(intergenic_noRandom[(intergenic_noRandom$spacerLength == 6) &
                                                                                            (intergenic_noRandom$NusA_stim == 'No'),]) /
                                                                   nrow(intergenic_noRandom[(intergenic_noRandom$NusA_stim == 'No'),]),
                                                                 nrow(intergenic_noRandom[(intergenic_noRandom$spacerLength == 7) &
                                                                                            (intergenic_noRandom$NusA_stim == 'No'),]) /
                                                                   nrow(intergenic_noRandom[(intergenic_noRandom$NusA_stim == 'No'),]),
                                                                 nrow(intergenic_noRandom[(intergenic_noRandom$spacerLength == 8) &
                                                                                            (intergenic_noRandom$NusA_stim == 'No'),]) /
                                                                   nrow(intergenic_noRandom[(intergenic_noRandom$NusA_stim == 'No'),]),
                                                                 nrow(intergenic_noRandom[(intergenic_noRandom$spacerLength == 9) &
                                                                                            (intergenic_noRandom$NusA_stim == 'No'),]) /
                                                                   nrow(intergenic_noRandom[(intergenic_noRandom$NusA_stim == 'No'),]),
                                                                 nrow(intergenic_noRandom[(intergenic_noRandom$spacerLength == 10) &
                                                                                            (intergenic_noRandom$NusA_stim == 'No'),]) /
                                                                   nrow(intergenic_noRandom[(intergenic_noRandom$NusA_stim == 'No'),])))

spacerLength_percentage_DF2 <- data.frame(spacerLength = rep(spacerLength_percentage_DF$spacer_length, 4),
                                          condition = c(rep('NusG_stim',11),
                                                        rep('nonNusG_stim',11),
                                                        rep('NusA_stim',11),
                                                        rep('nonNusA_stim',11)),
                                          spacer_percentage = c(spacerLength_percentage_DF$NusG_stim_spacer_percentage,
                                                                spacerLength_percentage_DF$nonNusG_stim_spacer_percentage,
                                                                spacerLength_percentage_DF$NusA_stim_spacer_percentage,
                                                                spacerLength_percentage_DF$nonNusA_stim_spacer_percentage))

spacerLength_NusG_only <- spacerLength_percentage_DF2 %>%
  filter(grepl("NusG", condition))

spacerLength_NusG_only$spacerLength_min <- -(spacerLength_NusG_only$spacerLength)

hairpin_spacer_hist_NusG <- ggplot(spacerLength_NusG_only,
                                   aes(x = spacerLength_min,
                                       y = spacer_percentage,
                                       fill = condition)) +
  geom_col(position = 'dodge',
           alpha = 0.8, width = 0.8,
           color = 'black') +
  scale_x_continuous(limits = c(-11,1),
                     breaks = seq(-10, 0, 1)) +
  scale_y_continuous(limits = c(0, 0.8),
                     breaks = seq(0, 0.8, 0.2)) +
  scale_fill_manual(values = c(nonstim_col, NusG_col))+
  theme(legend.position = "none")+
  theme_minimal()
ggsave(filename = '3enrich_NusAG/termination_figures/fig5c_spacerLength_hist_NusG.png',
       plot = hairpin_spacer_hist_NusG,
       width = 400,
       height = 300,
       units = 'mm',
       dpi = 300)


spacerLength_NusA_only <- spacerLength_percentage_DF2 %>%
  filter(grepl("NusA", condition))

spacerLength_NusA_only$spacerLength_min <- -(spacerLength_NusA_only$spacerLength)

hairpin_spacer_hist_NusA <- ggplot(spacerLength_NusA_only,
                                   aes(x = spacerLength_min,
                                       y = spacer_percentage,
                                       fill = condition)) +
  geom_col(position = 'dodge',
           alpha = 0.8, width = 0.8,
           color = 'black') +
  scale_x_continuous(limits = c(-11,1),
                     breaks = seq(-10, 0, 1)) +
  scale_y_continuous(limits = c(0, 0.8),
                     breaks = seq(0, 0.8, 0.2)) +
#  geom_area(aes(color=condition), alpha=0.3, position="identity") +
  scale_fill_manual(values = c(nonstim_col, NusA_col))+
#  scale_color_manual(values = c(nonstim_col, NusA_col))+
  theme(legend.position = "none")+
  theme_minimal()
ggsave(filename = '3enrich_NusAG/termination_figures/fig5c_spacerLength_hist_NusA.png',
       plot = hairpin_spacer_hist_NusA,
       width = 400,
       height = 300,
       units = 'mm',
       dpi = 300)

## Not significant indicators of Nus responses ---------------------------------

# U tract U-count
Utract_count_percentage_DF <- data.frame(Utract_count_length = seq(0,7,1),
                                         NusG_stim_Utract_count_percentage = c(nrow(intergenic_noRandom[(intergenic_noRandom$Utract_count == 0) &
                                                                                                    (intergenic_noRandom$NusG_stim == 'Yes'),]) /
                                                                           nrow(intergenic_noRandom[(intergenic_noRandom$NusG_stim == 'Yes'),]),
                                                                         nrow(intergenic_noRandom[(intergenic_noRandom$Utract_count == 1) &
                                                                                                    (intergenic_noRandom$NusG_stim == 'Yes'),]) /
                                                                           nrow(intergenic_noRandom[(intergenic_noRandom$NusG_stim == 'Yes'),]),
                                                                         nrow(intergenic_noRandom[(intergenic_noRandom$Utract_count == 2) &
                                                                                                    (intergenic_noRandom$NusG_stim == 'Yes'),]) /
                                                                           nrow(intergenic_noRandom[(intergenic_noRandom$NusG_stim == 'Yes'),]),
                                                                         nrow(intergenic_noRandom[(intergenic_noRandom$Utract_count == 3) &
                                                                                                    (intergenic_noRandom$NusG_stim == 'Yes'),]) /
                                                                           nrow(intergenic_noRandom[(intergenic_noRandom$NusG_stim == 'Yes'),]),
                                                                         nrow(intergenic_noRandom[(intergenic_noRandom$Utract_count == 4) &
                                                                                                    (intergenic_noRandom$NusG_stim == 'Yes'),]) /
                                                                           nrow(intergenic_noRandom[(intergenic_noRandom$NusG_stim == 'Yes'),]),
                                                                         nrow(intergenic_noRandom[(intergenic_noRandom$Utract_count == 5) &
                                                                                                    (intergenic_noRandom$NusG_stim == 'Yes'),]) /
                                                                           nrow(intergenic_noRandom[(intergenic_noRandom$NusG_stim == 'Yes'),]),
                                                                         nrow(intergenic_noRandom[(intergenic_noRandom$Utract_count == 6) &
                                                                                                    (intergenic_noRandom$NusG_stim == 'Yes'),]) /
                                                                           nrow(intergenic_noRandom[(intergenic_noRandom$NusG_stim == 'Yes'),]),
                                                                         nrow(intergenic_noRandom[(intergenic_noRandom$Utract_count == 7) &
                                                                                                    (intergenic_noRandom$NusG_stim == 'Yes'),]) /
                                                                           nrow(intergenic_noRandom[(intergenic_noRandom$NusG_stim == 'Yes'),])),
                                         NusA_stim_Utract_count_percentage = c(nrow(intergenic_noRandom[(intergenic_noRandom$Utract_count == 0) &
                                                                                                    (intergenic_noRandom$NusA_stim == 'Yes'),]) /
                                                                           nrow(intergenic_noRandom[(intergenic_noRandom$NusA_stim == 'Yes'),]),
                                                                         nrow(intergenic_noRandom[(intergenic_noRandom$Utract_count == 1) &
                                                                                                    (intergenic_noRandom$NusA_stim == 'Yes'),]) /
                                                                           nrow(intergenic_noRandom[(intergenic_noRandom$NusA_stim == 'Yes'),]),
                                                                         nrow(intergenic_noRandom[(intergenic_noRandom$Utract_count == 2) &
                                                                                                    (intergenic_noRandom$NusA_stim == 'Yes'),]) /
                                                                           nrow(intergenic_noRandom[(intergenic_noRandom$NusA_stim == 'Yes'),]),
                                                                         nrow(intergenic_noRandom[(intergenic_noRandom$Utract_count == 3) &
                                                                                                    (intergenic_noRandom$NusA_stim == 'Yes'),]) /
                                                                           nrow(intergenic_noRandom[(intergenic_noRandom$NusA_stim == 'Yes'),]),
                                                                         nrow(intergenic_noRandom[(intergenic_noRandom$Utract_count == 4) &
                                                                                                    (intergenic_noRandom$NusA_stim == 'Yes'),]) /
                                                                           nrow(intergenic_noRandom[(intergenic_noRandom$NusA_stim == 'Yes'),]),
                                                                         nrow(intergenic_noRandom[(intergenic_noRandom$Utract_count == 5) &
                                                                                                    (intergenic_noRandom$NusA_stim == 'Yes'),]) /
                                                                           nrow(intergenic_noRandom[(intergenic_noRandom$NusA_stim == 'Yes'),]),
                                                                         nrow(intergenic_noRandom[(intergenic_noRandom$Utract_count == 6) &
                                                                                                    (intergenic_noRandom$NusA_stim == 'Yes'),]) /
                                                                           nrow(intergenic_noRandom[(intergenic_noRandom$NusA_stim == 'Yes'),]),
                                                                         nrow(intergenic_noRandom[(intergenic_noRandom$Utract_count == 7) &
                                                                                                    (intergenic_noRandom$NusA_stim == 'Yes'),]) /
                                                                           nrow(intergenic_noRandom[(intergenic_noRandom$NusA_stim == 'Yes'),])),
                                         nonNusG_stim_Utract_count_percentage = c(nrow(intergenic_noRandom[(intergenic_noRandom$Utract_count == 0) &
                                                                                                       (intergenic_noRandom$NusG_stim == 'No'),]) /
                                                                              nrow(intergenic_noRandom[(intergenic_noRandom$NusG_stim == 'No'),]),
                                                                            nrow(intergenic_noRandom[(intergenic_noRandom$Utract_count == 1) &
                                                                                                       (intergenic_noRandom$NusG_stim == 'No'),]) /
                                                                              nrow(intergenic_noRandom[(intergenic_noRandom$NusG_stim == 'No'),]),
                                                                            nrow(intergenic_noRandom[(intergenic_noRandom$Utract_count == 2) &
                                                                                                       (intergenic_noRandom$NusG_stim == 'No'),]) /
                                                                              nrow(intergenic_noRandom[(intergenic_noRandom$NusG_stim == 'No'),]),
                                                                            nrow(intergenic_noRandom[(intergenic_noRandom$Utract_count == 3) &
                                                                                                       (intergenic_noRandom$NusG_stim == 'No'),]) /
                                                                              nrow(intergenic_noRandom[(intergenic_noRandom$NusG_stim == 'No'),]),
                                                                            nrow(intergenic_noRandom[(intergenic_noRandom$Utract_count == 4) &
                                                                                                       (intergenic_noRandom$NusG_stim == 'No'),]) /
                                                                              nrow(intergenic_noRandom[(intergenic_noRandom$NusG_stim == 'No'),]),
                                                                            nrow(intergenic_noRandom[(intergenic_noRandom$Utract_count == 5) &
                                                                                                       (intergenic_noRandom$NusG_stim == 'No'),]) /
                                                                              nrow(intergenic_noRandom[(intergenic_noRandom$NusG_stim == 'No'),]),
                                                                            nrow(intergenic_noRandom[(intergenic_noRandom$Utract_count == 6) &
                                                                                                       (intergenic_noRandom$NusG_stim == 'No'),]) /
                                                                              nrow(intergenic_noRandom[(intergenic_noRandom$NusG_stim == 'No'),]),
                                                                            nrow(intergenic_noRandom[(intergenic_noRandom$Utract_count == 7) &
                                                                                                       (intergenic_noRandom$NusG_stim == 'No'),]) /
                                                                              nrow(intergenic_noRandom[(intergenic_noRandom$NusG_stim == 'No'),])),
                                         nonNusA_stim_Utract_count_percentage = c(nrow(intergenic_noRandom[(intergenic_noRandom$Utract_count == 0) &
                                                                                                       (intergenic_noRandom$NusA_stim == 'No'),]) /
                                                                              nrow(intergenic_noRandom[(intergenic_noRandom$NusA_stim == 'No'),]),
                                                                            nrow(intergenic_noRandom[(intergenic_noRandom$Utract_count == 1) &
                                                                                                       (intergenic_noRandom$NusA_stim == 'No'),]) /
                                                                              nrow(intergenic_noRandom[(intergenic_noRandom$NusA_stim == 'No'),]),
                                                                            nrow(intergenic_noRandom[(intergenic_noRandom$Utract_count == 2) &
                                                                                                       (intergenic_noRandom$NusA_stim == 'No'),]) /
                                                                              nrow(intergenic_noRandom[(intergenic_noRandom$NusA_stim == 'No'),]),
                                                                            nrow(intergenic_noRandom[(intergenic_noRandom$Utract_count == 3) &
                                                                                                       (intergenic_noRandom$NusA_stim == 'No'),]) /
                                                                              nrow(intergenic_noRandom[(intergenic_noRandom$NusA_stim == 'No'),]),
                                                                            nrow(intergenic_noRandom[(intergenic_noRandom$Utract_count == 4) &
                                                                                                       (intergenic_noRandom$NusA_stim == 'No'),]) /
                                                                              nrow(intergenic_noRandom[(intergenic_noRandom$NusA_stim == 'No'),]),
                                                                            nrow(intergenic_noRandom[(intergenic_noRandom$Utract_count == 5) &
                                                                                                       (intergenic_noRandom$NusA_stim == 'No'),]) /
                                                                              nrow(intergenic_noRandom[(intergenic_noRandom$NusA_stim == 'No'),]),
                                                                            nrow(intergenic_noRandom[(intergenic_noRandom$Utract_count == 6) &
                                                                                                       (intergenic_noRandom$NusA_stim == 'No'),]) /
                                                                              nrow(intergenic_noRandom[(intergenic_noRandom$NusA_stim == 'No'),]),
                                                                            nrow(intergenic_noRandom[(intergenic_noRandom$Utract_count == 7) &
                                                                                                       (intergenic_noRandom$NusA_stim == 'No'),]) /
                                                                              nrow(intergenic_noRandom[(intergenic_noRandom$NusA_stim == 'No'),])))

Utract_count_percentage_DF2 <- data.frame(Utract_count = rep(Utract_count_percentage_DF$Utract_count_length, 4),
                                          condition = c(rep('NusG_stim',8),
                                                        rep('nonNusG_stim',8),
                                                        rep('NusA_stim',8),
                                                        rep('nonNusA_stim',8)),
                                          Utract_count_percentage = c(Utract_count_percentage_DF$NusG_stim_Utract_count_percentage,
                                                                Utract_count_percentage_DF$nonNusG_stim_Utract_count_percentage,
                                                                Utract_count_percentage_DF$NusA_stim_Utract_count_percentage,
                                                                Utract_count_percentage_DF$nonNusA_stim_Utract_count_percentage))

Utract_count_NusG_only <- Utract_count_percentage_DF2 %>%
  filter(grepl("NusG", condition))

hairpin_Utract_hist_NusG <- ggplot(Utract_count_NusG_only,
                                   aes(x = Utract_count,
                                       y = Utract_count_percentage,
                                       fill = condition)) +
  geom_col(position = 'dodge',
           alpha = 0.8, width = 0.8,
           color = 'black') +
  scale_x_continuous(limits = c(-1,8),
                     breaks = seq(0, 7, 1)) +
  scale_y_continuous(limits = c(0, 0.5),
                     breaks = seq(0, 0.5, 0.1)) +
  scale_fill_manual(values = c(nonstim_col, NusG_col))+
  theme(legend.position = "none")+
  theme_minimal()
ggsave(filename = '3enrich_NusAG/termination_figures/extData_fig8e_Utract_hist_NusG.png',
       plot = hairpin_Utract_hist_NusG,
       width = 400,
       height = 300,
       units = 'mm',
       dpi = 300)
NusG_Utract_pval <- wilcox.test(NusG_nonstim$Utract_count,
                                NusG_stim$Utract_count)


Utract_count_NusA_only <- Utract_count_percentage_DF2 %>%
  filter(grepl("NusA", condition))

hairpin_Utract_hist_NusA <- ggplot(Utract_count_NusA_only,
                                   aes(x = Utract_count,
                                       y = Utract_count_percentage,
                                       fill = condition)) +
  geom_col(position = 'dodge',
           alpha = 0.8, width = 0.8,
           color = 'black') +
  scale_x_continuous(limits = c(-1,8),
                     breaks = seq(0, 7, 1)) +
  scale_y_continuous(limits = c(0, 0.5),
                     breaks = seq(0, 0.5, 0.1)) +
  scale_fill_manual(values = c(nonstim_col, NusA_col))+
  theme(legend.position = "none")+
  theme_minimal()
ggsave(filename = '3enrich_NusAG/termination_figures/extData_fig8e_Utract_hist_NusA.png',
       plot = hairpin_Utract_hist_NusA,
       width = 400,
       height = 300,
       units = 'mm',
       dpi = 300)
NusA_Utract_pval <- wilcox.test(NusA_nonstim$Utract_count,
                                NusA_stim$Utract_count)

# Hairpin loop length differences
loopLength_percentage_DF <- data.frame(loopLength_length = seq(0,12,1),
                                         NusG_stim_loopLength_percentage = c(nrow(intergenic_noRandom[(intergenic_noRandom$loopLength == 0) &
                                                                                                          (intergenic_noRandom$NusG_stim == 'Yes'),]) /
                                                                                 nrow(intergenic_noRandom[(intergenic_noRandom$NusG_stim == 'Yes'),]),
                                                                               nrow(intergenic_noRandom[(intergenic_noRandom$loopLength == 1) &
                                                                                                          (intergenic_noRandom$NusG_stim == 'Yes'),]) /
                                                                                 nrow(intergenic_noRandom[(intergenic_noRandom$NusG_stim == 'Yes'),]),
                                                                               nrow(intergenic_noRandom[(intergenic_noRandom$loopLength == 2) &
                                                                                                          (intergenic_noRandom$NusG_stim == 'Yes'),]) /
                                                                                 nrow(intergenic_noRandom[(intergenic_noRandom$NusG_stim == 'Yes'),]),
                                                                               nrow(intergenic_noRandom[(intergenic_noRandom$loopLength == 3) &
                                                                                                          (intergenic_noRandom$NusG_stim == 'Yes'),]) /
                                                                                 nrow(intergenic_noRandom[(intergenic_noRandom$NusG_stim == 'Yes'),]),
                                                                               nrow(intergenic_noRandom[(intergenic_noRandom$loopLength == 4) &
                                                                                                          (intergenic_noRandom$NusG_stim == 'Yes'),]) /
                                                                                 nrow(intergenic_noRandom[(intergenic_noRandom$NusG_stim == 'Yes'),]),
                                                                               nrow(intergenic_noRandom[(intergenic_noRandom$loopLength == 5) &
                                                                                                          (intergenic_noRandom$NusG_stim == 'Yes'),]) /
                                                                                 nrow(intergenic_noRandom[(intergenic_noRandom$NusG_stim == 'Yes'),]),
                                                                               nrow(intergenic_noRandom[(intergenic_noRandom$loopLength == 6) &
                                                                                                          (intergenic_noRandom$NusG_stim == 'Yes'),]) /
                                                                                 nrow(intergenic_noRandom[(intergenic_noRandom$NusG_stim == 'Yes'),]),
                                                                               nrow(intergenic_noRandom[(intergenic_noRandom$loopLength == 7) &
                                                                                                          (intergenic_noRandom$NusG_stim == 'Yes'),]) /
                                                                                 nrow(intergenic_noRandom[(intergenic_noRandom$NusG_stim == 'Yes'),]),
                                                                             nrow(intergenic_noRandom[(intergenic_noRandom$loopLength == 8) &
                                                                                                        (intergenic_noRandom$NusG_stim == 'Yes'),]) /
                                                                               nrow(intergenic_noRandom[(intergenic_noRandom$NusG_stim == 'Yes'),]),
                                                                             nrow(intergenic_noRandom[(intergenic_noRandom$loopLength == 9) &
                                                                                                        (intergenic_noRandom$NusG_stim == 'Yes'),]) /
                                                                               nrow(intergenic_noRandom[(intergenic_noRandom$NusG_stim == 'Yes'),]),
                                                                             nrow(intergenic_noRandom[(intergenic_noRandom$loopLength == 10) &
                                                                                                        (intergenic_noRandom$NusG_stim == 'Yes'),]) /
                                                                               nrow(intergenic_noRandom[(intergenic_noRandom$NusG_stim == 'Yes'),]),
                                                                             nrow(intergenic_noRandom[(intergenic_noRandom$loopLength == 11) &
                                                                                                        (intergenic_noRandom$NusG_stim == 'Yes'),]) /
                                                                               nrow(intergenic_noRandom[(intergenic_noRandom$NusG_stim == 'Yes'),]),
                                                                             nrow(intergenic_noRandom[(intergenic_noRandom$loopLength == 12) &
                                                                                                        (intergenic_noRandom$NusG_stim == 'Yes'),]) /
                                                                               nrow(intergenic_noRandom[(intergenic_noRandom$NusG_stim == 'Yes'),])),
                                         NusA_stim_loopLength_percentage = c(nrow(intergenic_noRandom[(intergenic_noRandom$loopLength == 0) &
                                                                                                          (intergenic_noRandom$NusA_stim == 'Yes'),]) /
                                                                                 nrow(intergenic_noRandom[(intergenic_noRandom$NusA_stim == 'Yes'),]),
                                                                               nrow(intergenic_noRandom[(intergenic_noRandom$loopLength == 1) &
                                                                                                          (intergenic_noRandom$NusA_stim == 'Yes'),]) /
                                                                                 nrow(intergenic_noRandom[(intergenic_noRandom$NusA_stim == 'Yes'),]),
                                                                               nrow(intergenic_noRandom[(intergenic_noRandom$loopLength == 2) &
                                                                                                          (intergenic_noRandom$NusA_stim == 'Yes'),]) /
                                                                                 nrow(intergenic_noRandom[(intergenic_noRandom$NusA_stim == 'Yes'),]),
                                                                               nrow(intergenic_noRandom[(intergenic_noRandom$loopLength == 3) &
                                                                                                          (intergenic_noRandom$NusA_stim == 'Yes'),]) /
                                                                                 nrow(intergenic_noRandom[(intergenic_noRandom$NusA_stim == 'Yes'),]),
                                                                               nrow(intergenic_noRandom[(intergenic_noRandom$loopLength == 4) &
                                                                                                          (intergenic_noRandom$NusA_stim == 'Yes'),]) /
                                                                                 nrow(intergenic_noRandom[(intergenic_noRandom$NusA_stim == 'Yes'),]),
                                                                               nrow(intergenic_noRandom[(intergenic_noRandom$loopLength == 5) &
                                                                                                          (intergenic_noRandom$NusA_stim == 'Yes'),]) /
                                                                                 nrow(intergenic_noRandom[(intergenic_noRandom$NusA_stim == 'Yes'),]),
                                                                               nrow(intergenic_noRandom[(intergenic_noRandom$loopLength == 6) &
                                                                                                          (intergenic_noRandom$NusA_stim == 'Yes'),]) /
                                                                                 nrow(intergenic_noRandom[(intergenic_noRandom$NusA_stim == 'Yes'),]),
                                                                               nrow(intergenic_noRandom[(intergenic_noRandom$loopLength == 7) &
                                                                                                          (intergenic_noRandom$NusA_stim == 'Yes'),]) /
                                                                                 nrow(intergenic_noRandom[(intergenic_noRandom$NusA_stim == 'Yes'),]),
                                                                             nrow(intergenic_noRandom[(intergenic_noRandom$loopLength == 8) &
                                                                                                        (intergenic_noRandom$NusA_stim == 'Yes'),]) /
                                                                               nrow(intergenic_noRandom[(intergenic_noRandom$NusA_stim == 'Yes'),]),
                                                                             nrow(intergenic_noRandom[(intergenic_noRandom$loopLength == 9) &
                                                                                                        (intergenic_noRandom$NusA_stim == 'Yes'),]) /
                                                                               nrow(intergenic_noRandom[(intergenic_noRandom$NusA_stim == 'Yes'),]),
                                                                             nrow(intergenic_noRandom[(intergenic_noRandom$loopLength == 10) &
                                                                                                        (intergenic_noRandom$NusA_stim == 'Yes'),]) /
                                                                               nrow(intergenic_noRandom[(intergenic_noRandom$NusA_stim == 'Yes'),]),
                                                                             nrow(intergenic_noRandom[(intergenic_noRandom$loopLength == 11) &
                                                                                                        (intergenic_noRandom$NusA_stim == 'Yes'),]) /
                                                                               nrow(intergenic_noRandom[(intergenic_noRandom$NusA_stim == 'Yes'),]),
                                                                             nrow(intergenic_noRandom[(intergenic_noRandom$loopLength == 12) &
                                                                                                        (intergenic_noRandom$NusA_stim == 'Yes'),]) /
                                                                               nrow(intergenic_noRandom[(intergenic_noRandom$NusA_stim == 'Yes'),])),
                                         nonNusG_stim_loopLength_percentage = c(nrow(intergenic_noRandom[(intergenic_noRandom$loopLength == 0) &
                                                                                                             (intergenic_noRandom$NusG_stim == 'No'),]) /
                                                                                    nrow(intergenic_noRandom[(intergenic_noRandom$NusG_stim == 'No'),]),
                                                                                  nrow(intergenic_noRandom[(intergenic_noRandom$loopLength == 1) &
                                                                                                             (intergenic_noRandom$NusG_stim == 'No'),]) /
                                                                                    nrow(intergenic_noRandom[(intergenic_noRandom$NusG_stim == 'No'),]),
                                                                                  nrow(intergenic_noRandom[(intergenic_noRandom$loopLength == 2) &
                                                                                                             (intergenic_noRandom$NusG_stim == 'No'),]) /
                                                                                    nrow(intergenic_noRandom[(intergenic_noRandom$NusG_stim == 'No'),]),
                                                                                  nrow(intergenic_noRandom[(intergenic_noRandom$loopLength == 3) &
                                                                                                             (intergenic_noRandom$NusG_stim == 'No'),]) /
                                                                                    nrow(intergenic_noRandom[(intergenic_noRandom$NusG_stim == 'No'),]),
                                                                                  nrow(intergenic_noRandom[(intergenic_noRandom$loopLength == 4) &
                                                                                                             (intergenic_noRandom$NusG_stim == 'No'),]) /
                                                                                    nrow(intergenic_noRandom[(intergenic_noRandom$NusG_stim == 'No'),]),
                                                                                  nrow(intergenic_noRandom[(intergenic_noRandom$loopLength == 5) &
                                                                                                             (intergenic_noRandom$NusG_stim == 'No'),]) /
                                                                                    nrow(intergenic_noRandom[(intergenic_noRandom$NusG_stim == 'No'),]),
                                                                                  nrow(intergenic_noRandom[(intergenic_noRandom$loopLength == 6) &
                                                                                                             (intergenic_noRandom$NusG_stim == 'No'),]) /
                                                                                    nrow(intergenic_noRandom[(intergenic_noRandom$NusG_stim == 'No'),]),
                                                                                  nrow(intergenic_noRandom[(intergenic_noRandom$loopLength == 7) &
                                                                                                             (intergenic_noRandom$NusG_stim == 'No'),]) /
                                                                                    nrow(intergenic_noRandom[(intergenic_noRandom$NusG_stim == 'No'),]),
                                                                                nrow(intergenic_noRandom[(intergenic_noRandom$loopLength == 8) &
                                                                                                           (intergenic_noRandom$NusG_stim == 'No'),]) /
                                                                                  nrow(intergenic_noRandom[(intergenic_noRandom$NusG_stim == 'No'),]),
                                                                                nrow(intergenic_noRandom[(intergenic_noRandom$loopLength == 9) &
                                                                                                           (intergenic_noRandom$NusG_stim == 'No'),]) /
                                                                                  nrow(intergenic_noRandom[(intergenic_noRandom$NusG_stim == 'No'),]),
                                                                                nrow(intergenic_noRandom[(intergenic_noRandom$loopLength == 10) &
                                                                                                           (intergenic_noRandom$NusG_stim == 'No'),]) /
                                                                                  nrow(intergenic_noRandom[(intergenic_noRandom$NusG_stim == 'No'),]),
                                                                                nrow(intergenic_noRandom[(intergenic_noRandom$loopLength == 11) &
                                                                                                           (intergenic_noRandom$NusG_stim == 'No'),]) /
                                                                                  nrow(intergenic_noRandom[(intergenic_noRandom$NusG_stim == 'No'),]),
                                                                                nrow(intergenic_noRandom[(intergenic_noRandom$loopLength == 12) &
                                                                                                           (intergenic_noRandom$NusG_stim == 'No'),]) /
                                                                                  nrow(intergenic_noRandom[(intergenic_noRandom$NusG_stim == 'No'),])),
                                         nonNusA_stim_loopLength_percentage = c(nrow(intergenic_noRandom[(intergenic_noRandom$loopLength == 0) &
                                                                                                             (intergenic_noRandom$NusA_stim == 'No'),]) /
                                                                                    nrow(intergenic_noRandom[(intergenic_noRandom$NusA_stim == 'No'),]),
                                                                                  nrow(intergenic_noRandom[(intergenic_noRandom$loopLength == 1) &
                                                                                                             (intergenic_noRandom$NusA_stim == 'No'),]) /
                                                                                    nrow(intergenic_noRandom[(intergenic_noRandom$NusA_stim == 'No'),]),
                                                                                  nrow(intergenic_noRandom[(intergenic_noRandom$loopLength == 2) &
                                                                                                             (intergenic_noRandom$NusA_stim == 'No'),]) /
                                                                                    nrow(intergenic_noRandom[(intergenic_noRandom$NusA_stim == 'No'),]),
                                                                                  nrow(intergenic_noRandom[(intergenic_noRandom$loopLength == 3) &
                                                                                                             (intergenic_noRandom$NusA_stim == 'No'),]) /
                                                                                    nrow(intergenic_noRandom[(intergenic_noRandom$NusA_stim == 'No'),]),
                                                                                  nrow(intergenic_noRandom[(intergenic_noRandom$loopLength == 4) &
                                                                                                             (intergenic_noRandom$NusA_stim == 'No'),]) /
                                                                                    nrow(intergenic_noRandom[(intergenic_noRandom$NusA_stim == 'No'),]),
                                                                                  nrow(intergenic_noRandom[(intergenic_noRandom$loopLength == 5) &
                                                                                                             (intergenic_noRandom$NusA_stim == 'No'),]) /
                                                                                    nrow(intergenic_noRandom[(intergenic_noRandom$NusA_stim == 'No'),]),
                                                                                  nrow(intergenic_noRandom[(intergenic_noRandom$loopLength == 6) &
                                                                                                             (intergenic_noRandom$NusA_stim == 'No'),]) /
                                                                                    nrow(intergenic_noRandom[(intergenic_noRandom$NusA_stim == 'No'),]),
                                                                                  nrow(intergenic_noRandom[(intergenic_noRandom$loopLength == 7) &
                                                                                                             (intergenic_noRandom$NusA_stim == 'No'),]) /
                                                                                    nrow(intergenic_noRandom[(intergenic_noRandom$NusA_stim == 'No'),]),
                                                                                nrow(intergenic_noRandom[(intergenic_noRandom$loopLength == 8) &
                                                                                                           (intergenic_noRandom$NusA_stim == 'No'),]) /
                                                                                  nrow(intergenic_noRandom[(intergenic_noRandom$NusA_stim == 'No'),]),
                                                                                nrow(intergenic_noRandom[(intergenic_noRandom$loopLength == 9) &
                                                                                                           (intergenic_noRandom$NusA_stim == 'No'),]) /
                                                                                  nrow(intergenic_noRandom[(intergenic_noRandom$NusA_stim == 'No'),]),
                                                                                nrow(intergenic_noRandom[(intergenic_noRandom$loopLength == 10) &
                                                                                                           (intergenic_noRandom$NusA_stim == 'No'),]) /
                                                                                  nrow(intergenic_noRandom[(intergenic_noRandom$NusA_stim == 'No'),]),
                                                                                nrow(intergenic_noRandom[(intergenic_noRandom$loopLength == 11) &
                                                                                                           (intergenic_noRandom$NusA_stim == 'No'),]) /
                                                                                  nrow(intergenic_noRandom[(intergenic_noRandom$NusA_stim == 'No'),]),
                                                                                nrow(intergenic_noRandom[(intergenic_noRandom$loopLength == 12) &
                                                                                                           (intergenic_noRandom$NusA_stim == 'No'),]) /
                                                                                  nrow(intergenic_noRandom[(intergenic_noRandom$NusA_stim == 'No'),])))

loopLength_percentage_DF2 <- data.frame(loopLength = rep(loopLength_percentage_DF$loopLength_length, 4),
                                          condition = c(rep('NusG_stim',13),
                                                        rep('nonNusG_stim',13),
                                                        rep('NusA_stim',13),
                                                        rep('nonNusA_stim',13)),
                                          loopLength_percentage = c(loopLength_percentage_DF$NusG_stim_loopLength_percentage,
                                                                      loopLength_percentage_DF$nonNusG_stim_loopLength_percentage,
                                                                      loopLength_percentage_DF$NusA_stim_loopLength_percentage,
                                                                      loopLength_percentage_DF$nonNusA_stim_loopLength_percentage))

loopLength_NusG_only <- loopLength_percentage_DF2 %>%
  filter(grepl("NusG", condition))

hairpin_loopLength_hist_NusG <- ggplot(loopLength_NusG_only,
                                   aes(x = loopLength,
                                       y = loopLength_percentage,
                                       fill = condition)) +
  geom_col(position = 'dodge',
           alpha = 0.8, width = 0.8,
           color = 'black') +
  scale_x_continuous(limits = c(-1,13),
                     breaks = seq(0, 12, 1)) +
  scale_y_continuous(limits = c(0, 0.5),
                     breaks = seq(0, 0.5, 0.1)) +
  scale_fill_manual(values = c(nonstim_col, NusG_col))+
  theme(legend.position = "none")+
  theme_minimal()
ggsave(filename = '3enrich_NusAG/termination_figures/extData_fig8g_loopLength_hist_NusG.png',
       plot = hairpin_loopLength_hist_NusG,
       width = 400,
       height = 300,
       units = 'mm',
       dpi = 300)
NusG_loopLength_pval <- wilcox.test(NusG_nonstim$loopLength,
                                NusG_stim$loopLength)


loopLength_NusA_only <- loopLength_percentage_DF2 %>%
  filter(grepl("NusA", condition))

hairpin_loopLength_hist_NusA <- ggplot(loopLength_NusA_only,
                                       aes(x = loopLength,
                                           y = loopLength_percentage,
                                           fill = condition)) +
  geom_col(position = 'dodge',
           alpha = 0.8, width = 0.8,
           color = 'black') +
  scale_x_continuous(limits = c(-1,13),
                     breaks = seq(0, 12, 1)) +
  scale_y_continuous(limits = c(0, 0.5),
                     breaks = seq(0, 0.5, 0.1)) +
  scale_fill_manual(values = c(nonstim_col, NusA_col))+
  theme(legend.position = "none")+
  theme_minimal()
ggsave(filename = '3enrich_NusAG/termination_figures/extData_fig8g_loopLength_hist_NusA.png',
       plot = hairpin_loopLength_hist_NusA,
       width = 400,
       height = 300,
       units = 'mm',
       dpi = 300)
NusA_loopLength_pval <- wilcox.test(NusA_nonstim$loopLength,
                                    NusA_stim$loopLength)


## Hairpin stem length differences
NusA_stemLength <- ggplot(intergenic_noRandom, aes(x = NusA_stim,
                                                   y = stemLength,
                                                   fill = NusA_stim)) + 
  geom_violin() +
  geom_boxplot(width = 0.2) +
  scale_y_continuous(limits = c(0,20),
                     breaks = seq(0,20,5)) +
  scale_fill_manual(values = c(nonstim_col, NusA_col)) +
  theme_minimal()
ggsave(filename = '3enrich_NusAG/termination_figures/extData_fig8f_stemLength_NusA.png',
       plot = NusA_stemLength,
       width = 144,
       height = 80,
       units = 'mm',
       dpi = 300)
NusA_stemLength_pval <- wilcox.test(NusA_nonstim$stemLength,
                                    NusA_stim$stemLength)

NusG_stemLength <- ggplot(intergenic_noRandom, aes(x = NusG_stim,
                                                   y = stemLength,
                                                   fill = NusG_stim)) + 
  geom_violin() +
  geom_boxplot(width = 0.2) +
  scale_y_continuous(limits = c(0,20),
                     breaks = seq(0,20,5)) +
  scale_fill_manual(values = c(nonstim_col, NusG_col)) +
  theme_minimal()
ggsave(filename = '3enrich_NusAG/termination_figures/extData_fig8f_stemLength_NusG.png',
       plot = NusG_stemLength,
       width = 144,
       height = 80,
       units = 'mm',
       dpi = 300)
NusG_stemLength_pval <- wilcox.test(NusG_nonstim$stemLength,
                                    NusG_stim$stemLength)

pval_DF <- data.frame(deltaG_NusA = NusA_deltaG_pval$p.value,
                      deltaG_NusA_sig = NusA_deltaG_pval$p.value <= 0.05,
                      deltaG_NusG = NusG_deltaG_pval$p.value,
                      deltaG_NusG_sig = NusG_deltaG_pval$p.value <= 0.05,
                      dsDNA_NusA = NusA_dsDNA_pval$p.value,
                      dsDNA_NusA_sig = NusA_dsDNA_pval$p.value <= 0.05,
                      dsDNA_NusG = NusG_dsDNA_pval$p.value,
                      dsDNA_NusG_sig = NusG_dsDNA_pval$p.value <= 0.05,
                      spacerLength_NusA = NusA_spacerLength_pval$p.value,
                      spacerLength_NusA_sig = NusA_spacerLength_pval$p.value <= 0.05,
                      spacerLength_NusG = NusG_spacerLength_pval$p.value,
                      spacerLength_NusG_sig = NusG_spacerLength_pval$p.value <= 0.05,
                      Utract_NusA = NusA_Utract_pval$p.value,
                      Utract_NusA_sig = NusA_Utract_pval$p.value <= 0.05,
                      Utract_NusG = NusG_Utract_pval$p.value,
                      Utract_NusG_sig = NusG_Utract_pval$p.value <= 0.05,
                      loopLength_NusA = NusA_loopLength_pval$p.value,
                      loopLength_NusA_sig = NusA_loopLength_pval$p.value <= 0.05,
                      loopLength_NusG = NusG_loopLength_pval$p.value,
                      loopLength_NusG_sig = NusG_loopLength_pval$p.value <= 0.05,
                      stemLength_NusA = NusA_stemLength_pval$p.value,
                      stemLength_NusA_sig = NusA_stemLength_pval$p.value <= 0.05,
                      stemLength_NusG = NusG_stemLength_pval$p.value,
                      stemLength_NusG_sig = NusG_stemLength_pval$p.value <= 0.05)

write.csv(pval_DF, '3enrich_NusAG/termination_figures/pval_DF.csv')

## DNAshape predictions --------------------------------------------------------

stim_list <- c("Yes","No",'random')

# Write out dsDNA as fasta files for DNAshapeR
for (i in 1:length(stim_list)) {
  
  NusG_stim_DF <- data.frame(intergenic) %>%
    filter(NusG_stim == stim_list[i])
  
  write.fasta(sequences = strsplit(NusG_stim_DF$Nontemplate_strand, split=''),
              name = NusG_stim_DF$coordID, open = "w",
              as.string = FALSE,
              file.out = paste("3enrich_NusAG/termination_figures/extData_fig8h_shape_NusGstim_",stim_list[i],".fasta",sep=''))
  
  NusA_stim_DF <- data.frame(intergenic) %>%
    filter(NusA_stim == stim_list[i])
  
  write.fasta(sequences = strsplit(NusA_stim_DF$Nontemplate_strand, split=''),
              name = NusA_stim_DF$coordID, open = "w",
              as.string = FALSE,
              file.out = paste("3enrich_NusAG/termination_figures/extData_fig8h_shape_NusAstim_",stim_list[i],".fasta",sep=''))
}

NusG_DF <- data.frame(stim = 0,
                      position = 0,
                      MGW = 0,
                      HelT = 0,
                      ProT = 0,
                      Roll = 0,
                      EP = 0)

NusA_DF <- data.frame(stim = 0,
                      position = 0,
                      MGW = 0,
                      HelT = 0,
                      ProT = 0,
                      Roll = 0,
                      EP = 0)

# Generate NusG dataframe with DNA shape parameters
for (i in 1:length(stim_list)) {
  
  NusG_shape_features <- getShape(paste('3enrich_NusAG/termination_figures/extData_fig8h_shape_NusGstim_',stim_list[i],'.fasta',
                                        sep = ''))
  
  MGW <- NusG_shape_features$MGW[,1:ncol(NusG_shape_features$MGW) - 1]
  HelT <- NusG_shape_features$HelT
  ProT <- NusG_shape_features$ProT[,1:ncol(NusG_shape_features$ProT) - 1]
  Roll <- NusG_shape_features$Roll
  EP <- NusG_shape_features$EP[,1:ncol(NusG_shape_features$EP) - 1]
  
  MGW_vals <- MGW[,1]
  HelT_vals <- HelT[,1]
  ProT_vals <- ProT[,1]
  Roll_vals <- Roll[,1]
  EP_vals <- EP[,1]
  
  for (j in 2:ncol(MGW)) {
    
    MGW_vals <- c(MGW_vals, MGW[,j])
    HelT_vals <- c(HelT_vals, HelT[,j])
    ProT_vals <- c(ProT_vals, ProT[,j])
    Roll_vals <- c(Roll_vals, Roll[,j])
    EP_vals <- c(EP_vals, EP[,j])
  }
  
  NusG_DF_new <- data.frame(stim = rep(stim_list[i],
                                       length(MGW_vals)),
                            position = rep(-51:(-51+88),
                                           each = length(MGW[,1])),
                            MGW = MGW_vals,
                            HelT = HelT_vals,
                            ProT = ProT_vals,
                            Roll = Roll_vals,
                            EP = EP_vals)
  
  NusG_DF <- rbind(NusG_DF, NusG_DF_new)
  
}

for (i in 1:length(stim_list)) {
  
  NusA_shape_features <- getShape(paste('3enrich_NusAG/termination_figures/extData_fig8h_shape_NusAstim_',
                                        stim_list[i],'.fasta',
                                        sep = ''))
  
  MGW <- NusA_shape_features$MGW[,1:ncol(NusA_shape_features$MGW) - 1]
  HelT <- NusA_shape_features$HelT
  ProT <- NusA_shape_features$ProT[,1:ncol(NusA_shape_features$ProT) - 1]
  Roll <- NusA_shape_features$Roll
  EP <- NusA_shape_features$EP[,1:ncol(NusA_shape_features$EP) - 1]
  
  MGW_vals <- MGW[,1]
  HelT_vals <- HelT[,1]
  ProT_vals <- ProT[,1]
  Roll_vals <- Roll[,1]
  EP_vals <- EP[,1]
  
  for (j in 2:ncol(MGW)) {
    
    MGW_vals <- c(MGW_vals, MGW[,j])
    HelT_vals <- c(HelT_vals, HelT[,j])
    ProT_vals <- c(ProT_vals, ProT[,j])
    Roll_vals <- c(Roll_vals, Roll[,j])
    EP_vals <- c(EP_vals, EP[,j])
  }
  
  NusA_DF_new <- data.frame(stim = rep(stim_list[i],
                                       length(MGW_vals)),
                            position = rep(-51:(-51+88),
                                           each = length(MGW[,1])),
                            MGW = MGW_vals,
                            HelT = HelT_vals,
                            ProT = ProT_vals,
                            Roll = Roll_vals,
                            EP = EP_vals)
  
  NusA_DF <- rbind(NusA_DF, NusA_DF_new)
  
}

NusG_DF_toPlot <- NusG_DF[2:nrow(NusG_DF),]
NusA_DF_toPlot <- NusA_DF[2:nrow(NusA_DF),]

NusG_DF_subset <- NusG_DF_toPlot
NusA_DF_subset <- NusA_DF_toPlot

NusG_DF_subset <- NusG_DF_toPlot[NusG_DF_toPlot$stim != 'random',]
NusA_DF_subset <- NusA_DF_toPlot[NusA_DF_toPlot$stim != 'random',]

library(dplyr)

# Propeller twist
NusA_ProT_means <- NusA_DF_subset %>%
  group_by(position, stim) %>%
  summarise_at(vars(ProT), funs(mean(., na.rm=TRUE)))

NusA_ProT_means <- NusA_ProT_means %>%
  mutate(stim = factor(stim,
                       levels = c("No","Yes")))


NusA_ProT_lines <- ggplot(NusA_ProT_means, aes(x = position, y = ProT, color = stim,
                                               fill = stim)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "loess", se = TRUE,
              alpha = 0.3) +
  scale_color_manual(values = c(nonstim_col, NusA_col)) +
  scale_fill_manual(values = c(nonstim_col, NusA_col)) +
  theme_minimal() +
  #  ylim(4.5, 5.4) +
  theme(legend.position = 'none') +
  scale_x_continuous(limits = c(-2,12),
                     breaks = seq(-2,12,1)) +
  scale_y_continuous(limits = c(-10, -4),
                     breaks = seq(-10, -4, 2))
ggsave(filename = paste('3enrich_NusAG/termination_figures/extData_fig8h_shapeProT_NusA.png',sep=''),
       plot = NusA_ProT_lines,
       width = 144*1.74,
       height = 85,
       units = 'mm',
       dpi = 300)

NusG_ProT_means <- NusG_DF_subset %>%
  group_by(position, stim) %>%
  summarise_at(vars(ProT), funs(mean(., na.rm=TRUE)))

NusG_ProT_means <- NusG_ProT_means %>%
  mutate(stim = factor(stim,
                       levels = c("No","Yes")))

NusG_ProT_lines <- ggplot(NusG_ProT_means, aes(x = position, y = ProT, color = stim,
                                               fill = stim)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "loess", se = TRUE,
              alpha = 0.3) +
  scale_color_manual(values = c(nonstim_col, NusG_col)) +
  scale_fill_manual(values = c(nonstim_col, NusG_col)) +
  theme_minimal() +
  #  ylim(4.5, 5.4) +
  theme(legend.position = 'none') +
  scale_x_continuous(limits = c(-2,12),
                     breaks = seq(-2,12,1)) +
  scale_y_continuous(limits = c(-10, -4),
                     breaks = seq(-10, -4, 2))
ggsave(filename = paste('3enrich_NusAG/termination_figures/extData_fig8h_shapeProT_NusG.png',sep=''),
       plot = NusG_ProT_lines,
       width = 144*1.74,
       height = 85,
       units = 'mm',
       dpi = 300)

## End -------------------------------------------------------------------------

