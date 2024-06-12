#!/usr/bin/env Rscript
library(optparse)

option_list = list(
  make_option(c("-g", "--gene_model"), type="character", default=NULL, 
              help="dataset file name", metavar="character"),
  make_option(c("-f", "--fiveEnd_input"), type="character", default=NULL, 
              help="dataset file name", metavar="character"),
  make_option(c("-e", "--threeEnd_input"), type="character", default=NULL, 
              help="dataset file name", metavar="character"),
  make_option(c("-d", "--gDNA_input"), type="character", default=NULL, 
              help="dataset file name", metavar="character"),
  make_option(c("-s", "--TSS_input"), type="character", default=NULL, 
              help="dataset file name", metavar="character"),
  make_option(c("-t", "--TTS_input"), type="character", default=NULL, 
              help="dataset file name", metavar="character")
  make_option(c("-o", "--fiveEnd_output"), type="character", default=NULL, 
              help="dataset file name", metavar="character")
  make_option(c("-x", "--threeEnd_output"), type="character", default=NULL, 
              help="dataset file name", metavar="character")
  make_option(c("-y", "--TSS_output"), type="character", default=NULL, 
              help="dataset file name", metavar="character")
  make_option(c("-z", "--TTS_output"), type="character", default=NULL, 
              help="dataset file name", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

library(reshape2)
library(ggplot2)
library(dplyr)

## Colors ----------------------------------------------------------------------

gDNA_col <- '#666666'
core_col <- 'grey'
noCRP_col <- '#E87121'
CRP_col <- '#0075B6'

noTF_col = "grey"
NusG_col = "#E26B7A"
NusA_col = "#5F8ADC"
NusA_NusG_col = "#BD67F3"

## Gene-level: Read in summarizeOverlaps inputs --------------------------------

TSS_sumOver_withEco <- read.csv(opt$fiveEnd_input)
TTS_sumOver_withEco <- read.csv(opt$threeEnd_input)
geneModel <- read.csv(opt$gene_model)

TSS_sumOver_withEco$strand <- geneModel$strand
TSS_sumOver_withEco$start <- geneModel$start_norm
TTS_sumOver_withEco$strand <- geneModel$strand
TTS_sumOver_withEco$start <- geneModel$start_norm

TSS_sumOver <- TSS_sumOver_withEco[c(3:nrow(TSS_sumOver_withEco)),c(5:ncol(TSS_sumOver_withEco))]
TTS_sumOver <- TTS_sumOver_withEco[c(3:nrow(TTS_sumOver_withEco)),c(5:ncol(TTS_sumOver_withEco))]

# Process gDNA
gDNA_counts <- read.csv(opt$gDNA_input)
gDNA_counts$strand <- geneModel$strand
gDNA_counts$start <- geneModel$start

gDNA_counts <- gDNA_counts[3:nrow(gDNA_counts),]

genes <- gDNA_counts %>%
  filter(!grepl("_as",geneID))

as_regions <- gDNA_counts %>%
  filter(grepl("_as",geneID))

gDNA_plus <- genes[,c(2:4)] + as_regions[,c(2:4)]
gDNA_minus <- gDNA_plus
gDNA_plus$strand <- genes$strand
gDNA_plus$geneID <- genes$geneID
gDNA_plus$start <- genes$start
gDNA_minus$strand <- as_regions$strand
gDNA_minus$geneID <- as_regions$geneID
gDNA_minus$start <- as_regions$start
gDNA_all <- rbind(gDNA_plus, gDNA_minus)
gDNA_all_stringent <- data.frame(gDNA_all) %>%
  arrange(start)

all_5end_samples <- cbind(gDNA_all_stringent[,c(1:3)],
                          TSS_sumOver[,c(1:18)])
all_3end_samples <- cbind(gDNA_all_stringent[,c(1:3)],
                          TTS_sumOver[,c(1:12)])

## Log2 transform and correlation analysis -------------------------------------

# Add a constant to prevent 0 --> -Inf upon log2 transform
log2_5end <- log2(all_5end_samples + 1)
log2_5end <- log2_5end[,c("core1_downsample.bam","core2_downsample.bam",'core3_downsample.bam',
                          "noCRP1_downsample.bam","noCRP2_downsample.bam",'noCRP3_downsample.bam',
                          "CRP1_downsample.bam","CRP2_downsample.bam",'CRP3_downsample.bam')]

#cormat <- round(cor(numeric_RNA),2)
cormat <- cor(log2_5end)
melted_cormat <- melt(cormat)

ggplot(data = melted_cormat, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile()

# Get lower triangle of the correlation matrix
get_lower_tri<-function(cormat){
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}

# Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}

upper_tri <- get_upper_tri(cormat)

melted_cormat <- melt(upper_tri, na.rm = TRUE)

log2_5endCov_corr_heatmap_noGDNA <- ggplot(data = melted_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = (0.985+1)/2, limit = c(0.985,1),
                       name="Pearson\nCorrelation") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()

ggsave(filename = opt$fiveEnd_output,
       plot = log2_5endCov_corr_heatmap,
       width = 300,
       height = 300,
       units = 'mm',
       dpi = 300)


# Add a constant to prevent 0 --> -Inf upon log2 transform
log2_3end <- log2(all_3end_samples + 1)

#cormat <- round(cor(numeric_RNA),2)
cormat <- cor(log2_3end)
melted_cormat <- melt(cormat)

ggplot(data = melted_cormat, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile()

upper_tri <- get_upper_tri(cormat)

melted_cormat <- melt(upper_tri, na.rm = TRUE)

log2_3endCov_corr_heatmap <- ggplot(data = melted_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0.985, limit = c(0.97,1),
                       name="Pearson\nCorrelation") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()

ggsave(filename = opt$threeEnd_output,
       plot = log2_3endCov_corr_heatmap,
       width = 300,
       height = 300,
       units = 'mm',
       dpi = 300)

## TSS/TTS correlation analysis ------------------------------------------------

TSS_data <- read.csv(opt$TSS_input)
TSS_counts_all <- TSS_data[,c(3:ncol(TSS_data))]
TSS_counts <- TSS_counts_all[,c("core1_count","core2_count","core3_count",
                                "noCRP1_count","noCRP2_count","noCRP3_count",
                                "CRP1_count","CRP2_count","CRP3_count")]

# Add a constant to prevent 0 --> -Inf upon log2 transform
log2_5end <- log2(TSS_counts + 1)

#cormat <- round(cor(numeric_RNA),2)
cormat <- cor(log2_5end)
melted_cormat <- melt(cormat)

ggplot(data = melted_cormat, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile()

# Get lower triangle of the correlation matrix
get_lower_tri<-function(cormat){
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}

# Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}

upper_tri <- get_upper_tri(cormat)

melted_cormat <- melt(upper_tri, na.rm = TRUE)

log2_TSS_corr_heatmap <- ggplot(data = melted_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0.75, limit = c(0.5,1),
                       name="Pearson\nCorrelation") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()

ggsave(filename = opt$TSS_output,
       plot = log2_TSS_corr_heatmap,
       width = 300,
       height = 300,
       units = 'mm',
       dpi = 300)

## TTS
TTS_data <- read.csv(opt$TTS_input)
TTS_counts <- TTS_data[,c(3:ncol(TTS_data))]

# Add a constant to prevent 0 --> -Inf upon log2 transform
log2_3end <- log2(TTS_counts + 1)

#cormat <- round(cor(numeric_RNA),2)
cormat <- cor(log2_3end)
melted_cormat <- melt(cormat)

ggplot(data = melted_cormat, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile()

upper_tri <- get_upper_tri(cormat)

melted_cormat <- melt(upper_tri, na.rm = TRUE)

melted_cormat_noGDNA <- melted_cormat %>%
  filter(!grepl("gDNA", Var1) | !grepl("gDNA", Var2))

log2_TTS_corr_heatmap <- ggplot(data = melted_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1),
                       name="Pearson\nCorrelation") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()

ggsave(filename = opt$TTS_output,
       plot = log2_TTS_corr_heatmap,
       width = 300,
       height = 300,
       units = 'mm',
       dpi = 300)

