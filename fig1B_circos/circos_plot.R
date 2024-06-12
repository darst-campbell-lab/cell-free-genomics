#!/usr/bin/env Rscript
library(optparse)

option_list = list(
  make_option(c("-f", "--fiveEnd_input"), type="character", default=NULL, 
              help="dataset file name", metavar="character"),
  make_option(c("-e", "--threeEnd_input"), type="character", default=NULL, 
              help="dataset file name", metavar="character"),
  make_option(c("-s", "--TSS_input"), type="character", default=NULL, 
              help="dataset file name", metavar="character"),
  make_option(c("-t", "--TTS_input"), type="character", default=NULL, 
              help="dataset file name", metavar="character"),
  make_option(c("-g", "--gene_model"), type="character", default=NULL, 
              help="dataset file name", metavar="character"),
  make_option(c("-s", "--len_spike_genome"), type="character", default=NULL, 
              help="dataset file name", metavar="character"),
  make_option(c("-o", "--TSS_output"), type="character", default=NULL,
              help="dataset file name", metavar="character"),  
  make_option(c("-x", "--TTS_output"), type="character", default=NULL,
              help="dataset file name", metavar="character"),  
); 

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

library(Biostrings)
library(seqinr)
#library(ShortRead)
library(BiocParallel)
# Use to activate conda environment
library(Rsubread)
library(GenomicAlignments)
library(Rsamtools)
library(SummarizedExperiment)
library(ggplot2)
library(dplyr)
library(sqldf)
library(GenomicFeatures)
library(Rsamtools)
library(ggplot2)
# For vectorized string operations
library(stringr)
# For dataframe manipulation
library(ggplot2)
library(circlize)

## Read in summarizeOverlaps count matrices for all samples ---------------

counts_all <- cbind(read.csv(opt$fiveEnd_input,
                             opt$threeEnd_input))
geneModel <- read.csv(opt$gene_model)
geneModel$start_norm <- geneModel$start - as.integer(opt$len_spike_genome)
geneModel$end_norm <- geneModel$end - as.integer(opt$len_spike_genome)


counts <- counts_all[3:nrow(counts_all),]
geneModel <- geneModel[3:nrow(geneModel),]
counts$unique_row_ID <- geneModel$unique_row_ID
counts$start <- geneModel$start_norm
counts$strand <- geneModel$strand

# Process gDNA samples
gDNA_cov <- counts[,c("gDNA1_downsample.bam",
                      "gDNA2_downsample.bam",
                      "gDNA3_downsample.bam",
                      "strand","start",
                      "unique_row_ID")]

plus_regions <- gDNA_cov %>%
  filter(strand == '+')
minus_regions <- gDNA_cov %>%
  filter(strand == '-')

plus_regions$gDNA1_all <- plus_regions$gDNA1_downsample.bam + minus_regions$gDNA1_downsample.bam
plus_regions$gDNA2_all <- plus_regions$gDNA2_downsample.bam + minus_regions$gDNA2_downsample.bam
plus_regions$gDNA3_all <- plus_regions$gDNA3_downsample.bam + minus_regions$gDNA3_downsample.bam

minus_regions <- plus_regions[c("gDNA1_all",
                                "gDNA2_all",
                                "gDNA3_all",
                                "strand",
                                "unique_row_ID",
                                "start")]
minus_regions$strand <- '-'

plus_regions_toJoin <- plus_regions[,c("gDNA1_all",
                                       "gDNA2_all",
                                       "gDNA3_all",
                                       "strand",
                                       "unique_row_ID",
                                       "start")]

gDNA_all_notOrdered <- rbind(plus_regions_toJoin,
                             minus_regions)

gDNA_all <- data.frame(gDNA_all_notOrdered) %>%
  arrange(start)

all_5end_samples <- cbind(gDNA_all[,c(1:3)], counts[,c(8:13)])
all_3end_samples <- cbind(gDNA_all[,c(1:3)], counts[,c(18:29)])

log2_5end <- log2(all_5end_samples + 1)
log2_3end <- log2(all_3end_samples + 1)

coverage_5end_DF <- data.frame(gDNA = rowMeans(log2_5end[,c(1:3)]),
                               RNA = rowMeans(log2_5end[,c(4:ncol(log2_5end))]))
coverage_5end_DF <- data.frame(gDNA = rowMeans(log2_5end[,c(1:3)]),
                               RNA = rowMeans(log2_5end[,c(4:6)]))
coverage_5end_DF$strand <- gDNA_all$strand

coverage_3end_DF <- data.frame(gDNA = rowMeans(log2_3end[,c(1:3)]),
                               RNA = rowMeans(log2_3end[,c(4:ncol(log2_3end))]))
coverage_3end_DF <- data.frame(gDNA = rowMeans(log2_3end[,c(1:3)]),
                               RNA = rowMeans(log2_3end[,c(4:6)]))
coverage_3end_DF$strand <- gDNA_all$strand

## Read in high-conf calls -----------------------------------------------------

TSS_all <- read.csv(opt$TSS_input)
TSS_all[is.na(TSS_all)] <- 0
TSS_all$start <- TSS_all$start - 4641652
TSS_all$end <- TSS_all$end - 4641652
TSS_all$TSS_union <- TSS_all$X5end_enriched_end_count_20
TSS_all$log2_TSS_counts <- log2(TSS_all$TSS_union + 1)

TTS_all <- read.csv(opt$TTS_input)
TTS_all$start <- TTS_all$start - 4641652
TTS_all$end <- TTS_all$end - 4641652
TTS_all$TTS_union <- TTS_all$X3end_enriched_end_count_70
TTS_all$log2_TTS_counts <- log2(TTS_all$TTS_union + 1)

TSS_colsToKeep <- c("unique_row_ID","start","end","strand",
                    "log2_TSS_counts")
TTS_colsToKeep <- c("unique_row_ID","start","end","strand",
                    "log2_TTS_counts")

TSS <- TSS_all[,TSS_colsToKeep]
TTS <- TTS_all[,TTS_colsToKeep]

TSS_minus <- TSS[TSS$strand == '-',]
TSS_plus <- TSS[TSS$strand == '+',]
TTS_minus <- TTS[TTS$strand == '-',]
TTS_plus <- TTS[TTS$strand == '+',]

## Join coverage and TSS calls -------------------------------------------------

plus_coverage_5end_DF <- coverage_5end_DF[coverage_5end_DF$strand == '+',]
minus_coverage_5end_DF <- coverage_5end_DF[coverage_5end_DF$strand == '-',]

plus_coverage_3end_DF <- coverage_3end_DF[coverage_3end_DF$strand == '+',]
minus_coverage_3end_DF <- coverage_3end_DF[coverage_3end_DF$strand == '-',]

colnames(plus_coverage_5end_DF) <- c("gDNA_cov","5end_RNA_cov_plus", "strand")
colnames(minus_coverage_5end_DF) <- c("gDNA_cov_plus","5end_RNA_cov_minus", "strand")
colnames(plus_coverage_3end_DF) <- c("gDNA_cov_plus","3end_RNA_cov_plus", "strand")
colnames(minus_coverage_3end_DF) <- c("gDNA_cov_plus","3end_RNA_cov_minus", "strand")
colnames(TSS_plus) <- c("unique_row_ID","start","end","strand","TSS_plus")
colnames(TSS_minus) <- c("unique_row_ID","start","end","strand","TSS_minus")
colnames(TTS_plus) <- c("unique_row_ID","start","end","strand","TTS_plus")
colnames(TTS_minus) <- c("unique_row_ID","start","end","strand","TTS_minus")

circos_DF_all <- cbind(plus_coverage_5end_DF,
                       minus_coverage_5end_DF,
                       plus_coverage_3end_DF,
                       minus_coverage_3end_DF,
                       TSS_plus,
                       TSS_minus,
                       TTS_plus,
                       TTS_minus)

cols_to_keep <- c("start","end","gDNA_cov",
                  "5end_RNA_cov_plus","5end_RNA_cov_minus",
                  "TSS_plus","TSS_minus",
                  "3end_RNA_cov_plus","3end_RNA_cov_minus",
                  "TTS_plus","TTS_minus")

circos_DF <- circos_DF_all[,cols_to_keep]

# Generate circos plot for TSS data --------------------------------------------

# Colors
igv_pink <- '#E7ADAC'
igv_purple <- '#AFADD5'
track_height <- 0.165

df = data.frame(name = 'H37Rv',
                start = circos_DF$start,
                end = circos_DF$end)

png(opt$TSS_output,
    res = 300,
    width = 12,
    height = 12,
    units = 'in')

circos.genomicInitialize(df,
                         major.by = 1000000)


# Union TSS coverage
circos.genomicTrack(ylim = c(0,
                             max(circos_DF$TSS_plus)),
                    track.height = track_height)
circos.genomicLines(region = data.frame(start = circos_DF$start,
                                        end = circos_DF$end),
                    value = circos_DF$TSS_plus,
                    col = igv_pink,
                    lwd = 2,
                    area = FALSE)

circos.genomicTrack(ylim = c(0,
                             max(circos_DF$TSS_minus)),
                    track.height = track_height)
circos.genomicLines(region = data.frame(start = circos_DF$start,
                                        end = circos_DF$end),
                    value = circos_DF$TSS_minus,
                    col = igv_purple,
                    lwd = 2,
                    area = FALSE)

# Mean RNA coverage
circos.genomicTrack(ylim = c(0,
                             max(circos_DF$`5end_RNA_cov_plus`)),
                    track.height = track_height)
circos.genomicLines(region = data.frame(start = circos_DF$start,
                                        end = circos_DF$end),
                    value = circos_DF$`5end_RNA_cov_plus`,
                    col = igv_pink,
                    area = TRUE,
                    border = "NA")

circos.genomicTrack(ylim = c(0,
                             max(circos_DF$`5end_RNA_cov_minus`)),
                    track.height = track_height)
circos.genomicLines(region = data.frame(start = circos_DF$start,
                                        end = circos_DF$end),
                    value = circos_DF$`5end_RNA_cov_minus`,
                    col = igv_purple,
                    area = TRUE,
                    border = "NA")


# Mean gDNA coverage
circos.genomicTrack(ylim = c(0, max(circos_DF$gDNA_cov)),
                    track.height = track_height)
circos.genomicLines(region = data.frame(start = circos_DF$start,
                                        end = circos_DF$end),
                    value = circos_DF$gDNA_cov,
                    col = "grey",
                    area = TRUE,
                    border = "NA")
dev.off()



# Generate circos plot for TTS data --------------------------------------------

png(opt$TTS_output,
    res = 300,
    width = 12,
    height = 12,
    units = 'in')

circos.genomicInitialize(df,
                         major.by = 1000000)

# Union TTS coverage
circos.genomicTrack(ylim = c(0,
                             max(circos_DF$TTS_plus)),
                    track.height = track_height)
circos.genomicLines(region = data.frame(start = circos_DF$start,
                                        end = circos_DF$end),
                    value = circos_DF$TTS_plus,
                    col = igv_pink,
                    lwd = 2,
                    area = FALSE)

circos.genomicTrack(ylim = c(0,
                             max(circos_DF$TTS_minus)),
                    track.height = track_height)
circos.genomicLines(region = data.frame(start = circos_DF$start,
                                        end = circos_DF$end),
                    value = circos_DF$TTS_minus,
                    col = igv_purple,
                    lwd = 2,
                    area = FALSE)


# Mean RNA coverage
circos.genomicTrack(ylim = c(0,
                             max(circos_DF$`3end_RNA_cov_plus`)),
                    track.height = track_height)
circos.genomicLines(region = data.frame(start = circos_DF$start,
                                        end = circos_DF$end),
                    value = circos_DF$`3end_RNA_cov_plus`,
                    col = igv_pink,
                    area = TRUE,
                    border = "NA")

circos.genomicTrack(ylim = c(0,
                             max(circos_DF$`3end_RNA_cov_minus`)),
                    track.height = track_height)
circos.genomicLines(region = data.frame(start = circos_DF$start,
                                        end = circos_DF$end),
                    value = circos_DF$`3end_RNA_cov_minus`,
                    col = igv_purple,
                    area = TRUE,
                    border = "NA")

# Mean gDNA coverage
circos.genomicTrack(ylim = c(0, max(circos_DF$gDNA_cov)),
                    track.height = track_height)
circos.genomicLines(region = data.frame(start = circos_DF$start,
                                        end = circos_DF$end),
                    value = circos_DF$gDNA_cov,
                    col = "grey",
                    area = TRUE,
                    border = "NA")
dev.off()
