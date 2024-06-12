#!/usr/bin/env Rscript
library(optparse)

option_list = list(
  make_option(c("-g", "--gene_model"), type="character", default=NULL, 
              help="dataset file name", metavar="character"),
  make_option(c("-s", "--summed_counts_input"), type="character", default=NULL, 
              help="dataset file name", metavar="character"),
  make_option(c("-f", "--output_5end"), type="character", default=NULL, 
              help="dataset file name", metavar="character"),
  make_option(c("-e", "--output_3end"), type="character", default=NULL, 
              help="dataset file name", metavar="character"),
); 

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

library(ggplot2)
library(dplyr)

## Colors ----------------------------------------------------------------------

inCell_col <- '#EDD3B4'
inCell_col_dark <- '#A89074'
gDNA_col <- '#666666'
core_col <- 'grey'
RC_CRP_col <- '#EC7121'
CRP_col <- '#0075B7'

noTF_col = "grey"
NusG_col = "#E26B7A"
NusA_col = "#5F8ADC"
NusA_NusG_col = "#BD67F3"

## In-cell vs. out of cell coverage --------------------------------------------

counts_all <- read.csv(opt$summed_counts_input)
geneModel_all <- read.csv(opt$gene_model)
counts_noEco <- counts_all[3:nrow(counts_all),]
geneModel <- geneModel_all[3:nrow(geneModel_all),]

gDNA_cov <- counts_noEco[,c("sorted_gDNA1_trimmed.bam",
                      "sorted_gDNA2_trimmed.bam",
                      "sorted_gDNA3_trimmed.bam")]
gDNA_cov$strand <- geneModel$strand
gDNA_cov$start <- geneModel$start_norm

genes <- gDNA_cov %>%
  filter(!grepl("_as",geneID))

as_regions <- gDNA_cov %>%
  filter(grepl("_as",geneID))

gDNA_plus <- genes[,c(1:3)] + as_regions[,c(1:3)]
gDNA_minus <- gDNA_plus
gDNA_plus$strand <- genes$strand
gDNA_plus$geneID <- genes$geneID
gDNA_plus$start <- genes$start
gDNA_minus$strand <- as_regions$strand
gDNA_minus$geneID <- as_regions$geneID
gDNA_minus$start <- as_regions$start
gDNA_all <- rbind(gDNA_plus, gDNA_minus)
gDNA_all <- data.frame(gDNA_all) %>%
  arrange(start)

all_5end_samples <- cbind(gDNA_all[,c(1:3)], counts_noEco[,c("sorted_noCRP1_trimmed.bam",
                                                             "sorted_noCRP2_trimmed.bam",
                                                             "sorted_noCRP3_trimmed.bam",
                                                             "sorted_CRP1_trimmed.bam",
                                                             "sorted_CRP2_trimmed.bam",
                                                             "sorted_CRP3_trimmed.bam",
                                                             "sorted_Exp_R1_TSS_trimmed.bam",
                                                             "sorted_Exp_R2_TSS_trimmed.bam",
                                                             "sorted_Exp_R3_TSS_trimmed.bam")])
all_3end_samples <- cbind(gDNA_all[,c(1:3)], counts_noEco[,c("sorted_noTF1_trimmed.bam",
                                                             "sorted_noTF2_trimmed.bam",
                                                             "sorted_noTF3_trimmed.bam",
                                                             "sorted_NusA1_trimmed.bam",
                                                             "sorted_NusA2_trimmed.bam",
                                                             "sorted_NusA3_trimmed.bam",
                                                             "sorted_NusG1_trimmed.bam",
                                                             "sorted_NusG2_trimmed.bam",
                                                             "sorted_NusG3_trimmed.bam",
                                                             "sorted_NusA_NusG1_trimmed.bam",
                                                             "sorted_NusA_NusG2_trimmed.bam",
                                                             "sorted_NusA_NusG3_trimmed.bam",
                                                             "sorted_termseq_expo_r1_trimmed.bam",
                                                             "sorted_termseq_expo_r2_trimmed.bam",
                                                             "sorted_termseq_expo_r3_trimmed.bam")])

#### 5' end direct count comparison across genomic windows ------------------------------

CFG_hist_color <- '#CBADC1'
gDNA_col_light <- '#A0A0A0'

mean_gDNA_samples <- rowMeans(all_5end_samples[,c("sorted_gDNA1_trimmed.bam",
                                                  "sorted_gDNA2_trimmed.bam",
                                                  "sorted_gDNA3_trimmed.bam")])
mean_inCell_5end_samples <- rowMeans(all_5end_samples[,c("sorted_Exp_R1_TSS_trimmed.bam",
                                                         "sorted_Exp_R2_TSS_trimmed.bam",
                                                         "sorted_Exp_R3_TSS_trimmed.bam")])
mean_CFG_5end_samples <- rowMeans(all_5end_samples[,c("sorted_noCRP1_trimmed.bam",
                                                      "sorted_noCRP2_trimmed.bam",
                                                      "sorted_noCRP3_trimmed.bam",
                                                      "sorted_CRP1_trimmed.bam",
                                                      "sorted_CRP2_trimmed.bam",
                                                      "sorted_CRP3_trimmed.bam")])

count_comparison <- data.frame(inCell = mean_inCell_5end_samples,
                               gDNA = mean_gDNA_samples,
                               CFG = mean_CFG_5end_samples)

hist_DF <- data.frame(values = c(count_comparison$inCell,
                                 count_comparison$gDNA,
                                 count_comparison$CFG),
                      condition = c(rep("inCell", nrow(count_comparison)),
                                    rep("gDNA", nrow(count_comparison)),
                                    rep("CFG", nrow(count_comparison))))

padj_method = 'bonf'
kruskal.test(log2(values+1) ~ condition,
             hist_DF)
pairwise.wilcox.test(log2(hist_DF$values+1), hist_DF$condition,
                     p.adjust.method = padj_method)

toPlot_5end <- hist_DF %>%
  mutate(condition = factor(condition,
                            levels = c('gDNA','CFG','inCell')))

boxplot_5end <- ggplot(toPlot_5end,
                       aes(y = log2(values+1),
                           x = condition,
                           fill = condition))+
  geom_violin(linewidth = 2) +
  geom_boxplot(width = 0.2,
               linewidth = 2) +
  scale_fill_manual(values = c(gDNA_col_light, CFG_hist_color,
                               inCell_col)) +
  theme_minimal() +
  geom_hline(aes(yintercept = log2(10)),
             linetype = 'longdash',
             linewidth = 2) +
  scale_y_continuous(limits = c(0, 17),
                     breaks = seq(0, 15, 5)) +
  theme(legend.position = 'none')

ggsave(filename = opt$output_5end,
       plot = boxplot_5end,
       width = 300,
       height = 300,
       units = 'mm',
       dpi = 300)

#### 3' end direct count comparison across genomic windows ------------------------------

mean_gDNA_samples <- rowMeans(all_3end_samples[,c("sorted_gDNA1_trimmed.bam",
                                                  "sorted_gDNA2_trimmed.bam",
                                                  "sorted_gDNA3_trimmed.bam")])
mean_inCell_3end_samples <- rowMeans(all_3end_samples[,c("sorted_termseq_expo_r1_trimmed.bam",
                                                         "sorted_termseq_expo_r2_trimmed.bam",
                                                         "sorted_termseq_expo_r3_trimmed.bam")])
mean_CFG_3end_samples <- rowMeans(all_3end_samples[,c("sorted_noTF1_trimmed.bam",
                                                      "sorted_noTF2_trimmed.bam",
                                                      "sorted_noTF3_trimmed.bam",
                                                      "sorted_NusA1_trimmed.bam",
                                                      "sorted_NusA2_trimmed.bam",
                                                      "sorted_NusA3_trimmed.bam",
                                                      "sorted_NusG1_trimmed.bam",
                                                      "sorted_NusG2_trimmed.bam",
                                                      "sorted_NusG3_trimmed.bam",
                                                      "sorted_NusA_NusG1_trimmed.bam",
                                                      "sorted_NusA_NusG2_trimmed.bam",
                                                      "sorted_NusA_NusG3_trimmed.bam"))])

count_comparison <- data.frame(inCell = mean_inCell_3end_samples,
                               gDNA = mean_gDNA_samples,
                               CFG = mean_CFG_3end_samples)

hist_DF <- data.frame(values = c(count_comparison$inCell,
                                 count_comparison$gDNA,
                                 count_comparison$CFG),
                      condition = c(rep("inCell", nrow(count_comparison)),
                                    rep("gDNA", nrow(count_comparison)),
                                    rep("CFG", nrow(count_comparison))))

padj_method = 'bonf'
kruskal.test(log2(values+1) ~ condition,
             hist_DF)
pairwise.wilcox.test(log2(hist_DF$values+1), hist_DF$condition,
                     p.adjust.method = padj_method)

toPlot_3end <- hist_DF %>%
  mutate(condition = factor(condition,
                            levels = c('gDNA','CFG','inCell')))

boxplot_3end <- ggplot(toPlot_3end,
                       aes(y = log2(values+1),
                           x = condition,
                           fill = condition))+
  geom_violin(linewidth = 2) +
  geom_boxplot(width = 0.2,
               linewidth = 2) +
  scale_fill_manual(values = c(gDNA_col_light, CFG_hist_color,
                               inCell_col)) +
  theme_minimal() +
  geom_hline(aes(yintercept = log2(10)),
             linetype = 'longdash',
             linewidth = 2) +
  scale_y_continuous(limits = c(0, 17),
                     breaks = seq(0, 15, 5)) +
  theme(legend.position = 'none')

ggsave(filename = opt$output_3end,
       plot = boxplot_5end,
       width = 300,
       height = 300,
       units = 'mm',
       dpi = 300)

