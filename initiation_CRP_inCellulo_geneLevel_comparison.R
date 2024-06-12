# Integrating gene-level CRP CFG results with Mtb RNA-seq and ChIP-seq
# Authors: Ruby Froom, Michael DeJesus (GFF parsing)
# Date: May 29, 2024

#### Libraries -----------------------------------------------------------------

library(tidyverse)
library(stringr)
library(stringi)
library(Biostrings)
library(seqinr)
library(sqldf)
library(ggseqlogo)

# Read in CRP CFG coordinates with sequences annotated -------------------------

all_TSSs_withAS <- read.csv('fig2_inCellulo_vs_CFG/all_TSSs_needMotifs.csv')

# Since RNA-seq and ChIP-seq don't distinguish between sense and antisense,
# remove antisense TSSs from comparisons
all_TSSs <- all_TSSs_withAS[!grepl("_as",all_TSSs_withAS$geneID),]

# Add in CFG differential expression info
CFG_diff_exp <- read.csv('5enrich_CRP/selectThreshold/DESeq2/noCRP_CRP_20counts_ends_results_Wald_local.csv')
colnames(CFG_diff_exp) <- c("coordID",colnames(CFG_diff_exp[2:ncol(CFG_diff_exp)]))

CFG_diff_exp$CRP_effect <- "No effect"
CFG_diff_exp[(CFG_diff_exp$padj <= 0.05) & (CFG_diff_exp$log2FoldChange > 0),
       'CRP_effect'] <- 'Activated'
CFG_diff_exp[(CFG_diff_exp$padj <= 0.05) & (CFG_diff_exp$log2FoldChange < 0),
       'CRP_effect'] <- 'Repressed'

DE_toKeep <- CFG_diff_exp[,c("coordID","log2FoldChange","padj",'CRP_effect')]
colnames(DE_toKeep) <- paste("CFG_",
                             colnames(DE_toKeep),
                             sep='')

add_DE <- left_join(all_TSSs,
                    DE_toKeep,
                    by = c("CFG_TSS_ID_higher" = "CFG_coordID"))

# Remove columns relating to TSS overlap, just keep TSS location info
colsToKeep <- c("geneID","geneName","geneDesc","CFG_coord","CFG_direction",
                "CFG_TSS_ID_higher","CFG_log2FoldChange","CFG_padj","CFG_CRP_effect",
                "CFG_NT_strand","inCell_NT_strand",
                "inCell_TSS_ID",'TSS_distance')

CFG_DF <- add_DE[,colsToKeep]

# Number of CRP-regulated sense TSSs from CFG included in gene-level analysis
CFG_reg <- CFG_DF[CFG_DF$CFG_CRP_effect != 'No effect',]
length(unique(CFG_reg$CFG_TSS_ID_higher))

## Read in Kahramanoglou RNA-seq hits -------------------------------------------

RNAseq <- read.csv('5enrich_CRP/initiation_figures/CRP_inCellulo_hits/CRP_RNAseq.csv')
RNAseq$RNAseq_CRP_KO_effect <- 'No effect'
RNAseq[(RNAseq$padj <= 0.05) & (RNAseq$log2FoldChange > 0),
              'RNAseq_CRP_KO_effect'] <- 'Upregulated'
RNAseq[(RNAseq$padj <= 0.05) & (RNAseq$log2FoldChange < 0),
              'RNAseq_CRP_KO_effect'] <- 'Downregulated'

# Change gene name to RVBD
gene_ID_Rv = unlist(strsplit(RNAseq$id,'Rv'))
gene_suffix_Rv <- gene_ID_Rv[gene_ID_Rv != '']
RVBDGenes <- paste("RVBD",gene_suffix_Rv,sep='')
RNAseq$geneID <- RVBDGenes

RNAseq_toKeep <- RNAseq[,c(5:8)]

colnames(RNAseq_toKeep) <- c("RNAseq_log2FoldChange",'RNAseq_padj',
                             'RNAseq_CRP_KO_effect','geneID')

## Read in Kahramanoglou ChIP-seq ----------------------------------------------

CRP_ChIP <- read.csv('5enrich_CRP/initiation_figures/CRP_inCellulo_hits/CRP_ChIPseq.csv')
ChIP_peaks <- unique(c(CRP_ChIP$upstream, CRP_ChIP$internal))
ChIP_peaks <- ChIP_peaks[ChIP_peaks != '']
ChIP_peaks <- ChIP_peaks[ChIP_peaks != 'intergenic']
# Filter out tRNA genes: unclear which tRNA genes the ChIP peaks correspond to
ChIP_peaks <- ChIP_peaks[!grepl("^tRNA",ChIP_peaks)]

# Change gene name to RVBD
gene_ID_Rv = unlist(strsplit(ChIP_peaks,'Rv'))
gene_suffix_Rv <- gene_ID_Rv[gene_ID_Rv != '']
RVBDGenes <- paste("RVBD",gene_suffix_Rv,sep='')
ChIP_gene_list <- RVBDGenes[RVBDGenes != 'RVBDnc']

## Identify overlapping sets of genes between the 3 assays ---------------------

full_DF <- full_join(CFG_DF, RNAseq_toKeep,
                      by = c("geneID" = "geneID"))

# Get rid of NA values
full_DF[is.na(full_DF$CFG_log2FoldChange),
        grepl("CFG",colnames(full_DF))] <- 'Not transcribed in CFG'
full_DF[is.na(full_DF$RNAseq_log2FoldChange),
        grepl("RNAseq",colnames(full_DF))] <- 'Not measured in RNA-seq'

full_DF$ChIPseq_hit <- 'No'
full_DF[full_DF$geneID %in% ChIP_gene_list,'ChIPseq_hit'] <- 'Yes'

full_DF$Hit_type <- 'None'

# Hit in all 3: CRP-activated genes
full_DF[((full_DF$CFG_CRP_effect == 'Activated') &
           (full_DF$RNAseq_CRP_KO_effect == 'Downregulated') &
           (full_DF$ChIPseq_hit == 'Yes')),
        'Hit_type'] <- 'CFG, ChIP-seq and RNA-seq'

# Hit in all 3: CRP-repressed genes
full_DF[((full_DF$CFG_CRP_effect == 'Repressed') &
           (full_DF$RNAseq_CRP_KO_effect == 'Upregulated') &
           (full_DF$ChIPseq_hit == 'Yes')),
        'Hit_type'] <- 'CFG, ChIP-seq and RNA-seq'

# Hit in CFG and ChIP-seq only
full_DF[(((full_DF$CFG_CRP_effect == 'Activated') | (full_DF$CFG_CRP_effect == 'Repressed')) &
           ((full_DF$RNAseq_CRP_KO_effect == 'No effect') | (full_DF$RNAseq_CRP_KO_effect == 'Not measured in RNA-seq')) &
           (full_DF$ChIPseq_hit == 'Yes')),
        'Hit_type'] <- 'CFG and ChIP-seq'

# Hit in RNA-seq and ChIP-seq only
full_DF[(((full_DF$CFG_CRP_effect == 'No effect') | (full_DF$CFG_CRP_effect == 'Not transcribed in CFG')) &
           ((full_DF$RNAseq_CRP_KO_effect == 'Downregulated') | (full_DF$RNAseq_CRP_KO_effect == 'Upregulated')) &
           (full_DF$ChIPseq_hit == 'Yes')),
        'Hit_type'] <- 'RNA-seq and ChIP-seq'

# Hit in CFG and RNA-seq only: CRP-activated
full_DF[((full_DF$CFG_CRP_effect == 'Activated') &
           (full_DF$RNAseq_CRP_KO_effect == 'Downregulated') &
           (full_DF$ChIPseq_hit == 'No')),
        'Hit_type'] <- 'CFG and RNA-seq'

# Hit in CFG and RNA-seq only: CRP-repressed
full_DF[((full_DF$CFG_CRP_effect == 'Repressed') &
           (full_DF$RNAseq_CRP_KO_effect == 'Upregulated') &
           (full_DF$ChIPseq_hit == 'No')),
        'Hit_type'] <- 'CFG and RNA-seq'

# Hit in CFG only
full_DF[(((full_DF$CFG_CRP_effect == 'Activated') | (full_DF$CFG_CRP_effect == 'Repressed')) &
           ((full_DF$RNAseq_CRP_KO_effect == 'No effect') | (full_DF$RNAseq_CRP_KO_effect == 'Not measured in RNA-seq')) &
           (full_DF$ChIPseq_hit == 'No')),
        'Hit_type'] <- 'CFG only'

# Hit in ChIP-seq only
full_DF[(((full_DF$CFG_CRP_effect == 'No effect') | (full_DF$CFG_CRP_effect == 'Not transcribed in CFG')) &
           ((full_DF$RNAseq_CRP_KO_effect == 'No effect') | (full_DF$RNAseq_CRP_KO_effect == 'Not measured in RNA-seq')) &
           (full_DF$ChIPseq_hit == 'Yes')),
        'Hit_type'] <- 'ChIP-seq only'

# Hit in RNA-seq only
full_DF[(((full_DF$CFG_CRP_effect == 'No effect') | (full_DF$CFG_CRP_effect == 'Not transcribed in CFG')) &
           ((full_DF$RNAseq_CRP_KO_effect == 'Downregulated') | (full_DF$RNAseq_CRP_KO_effect == 'Upregulated')) &
           (full_DF$ChIPseq_hit == 'No')),
        'Hit_type'] <- 'RNA-seq only'

# Hit in ChIP-seq, conflicting hit in CFG and RNA-seq only (option 1)
full_DF[((full_DF$CFG_CRP_effect == 'Repressed') &
           (full_DF$RNAseq_CRP_KO_effect == 'Downregulated') &
           (full_DF$ChIPseq_hit == 'Yes')),
        'Hit_type'] <- 'ChIP-seq hit, conflict CFG vs. RNA-seq'

# Hit in ChIP-seq, conflicting hit in CFG and RNA-seq only (option 2)
full_DF[((full_DF$CFG_CRP_effect == 'Activated') &
           (full_DF$RNAseq_CRP_KO_effect == 'Upregulated') &
           (full_DF$ChIPseq_hit == 'Yes')),
        'Hit_type'] <- 'ChIP-seq hit, conflict CFG vs. RNA-seq'

# Not hit in ChIP-seq, conflicting hit in CFG and RNA-seq only (option 1)
full_DF[((full_DF$CFG_CRP_effect == 'Repressed') &
           (full_DF$RNAseq_CRP_KO_effect == 'Downregulated') &
           (full_DF$ChIPseq_hit == 'No')),
        'Hit_type'] <- 'conflict CFG vs. RNA-seq'

# Not hit in ChIP-seq, conflicting hit in CFG and RNA-seq only (option 2)
full_DF[((full_DF$CFG_CRP_effect == 'Activated') &
           (full_DF$RNAseq_CRP_KO_effect == 'Upregulated') &
           (full_DF$ChIPseq_hit == 'No')),
        'Hit_type'] <- 'conflict CFG vs. RNA-seq'

## Calculate enrichment of hit gene overlap ------------------------------------

# Make a list of all genes
no_conflicting_hits <- full_DF[!grepl("conflict", full_DF$Hit_type),]

all_genes <- unique(no_conflicting_hits$geneID)

fishers_overlap <- function(input_genes1, input_genes2, num_tests) {
  
  inBoth <- length(input_genes1[input_genes1 %in% input_genes2])
  genes1_only <- length(input_genes1[!(input_genes1 %in% input_genes2)])
  genes2_only <- length(input_genes2[!(input_genes2 %in% inBoth)])
  inAtleast1 <- unique(c(input_genes2, input_genes1))
  neither <- length(all_genes[!(all_genes %in% inAtleast1)])
  
  sig_matrix <- matrix(c(inBoth, genes1_only,
                         genes2_only, neither),
                       nrow = 2, ncol= 2)
  
  fisher_test <- fisher.test(sig_matrix, conf.level = (100-(5/num_tests))/100)
  return(fisher_test)
}

# Initialize empty list
p_values <- c()
odds_ratios <- c()
CI_lower <- c()
CI_upper <- c()


RNAseq_hits <- no_conflicting_hits[grepl("RNA-seq", no_conflicting_hits$Hit_type),
                                   'geneID']
ChIPseq_hits <- no_conflicting_hits[grepl("ChIP-seq", no_conflicting_hits$Hit_type),
                                    'geneID']
RNAseq_ChIPseq_hits <- no_conflicting_hits[(grepl("RNA-seq", no_conflicting_hits$Hit_type) & 
                                              grepl("ChIP-seq", no_conflicting_hits$Hit_type)),
                                           'geneID']

CFG_hits <- no_conflicting_hits[grepl("CFG", no_conflicting_hits$Hit_type),
                                'geneID']

gene_lists <- list(unique(RNAseq_hits),
                   unique(ChIPseq_hits),
                   unique(RNAseq_ChIPseq_hits))

for (i in 1:length(gene_lists)) {
  fishers_DF <- fishers_overlap(unique(CFG_hits), gene_lists[[i]],
                                length(gene_lists))
  p_values <- c(p_values, fishers_DF$p.value)
  odds_ratios <- c(odds_ratios, fishers_DF$estimate)
  CI_lower <- c(CI_lower, fishers_DF$conf.int[1])
  CI_upper <- c(CI_upper, fishers_DF$conf.int[2])
}

oddsRatio_DF <- data.frame(assay = c("RNA-seq","ChIP-seq","Both"),
                           pvalues = p_values,
                           odds_ratios = odds_ratios,
                           CI_lower = CI_lower,
                           CI_upper = CI_upper)

oddsRatio_DF$Bonferroni_sig_cutoff <- 0.05/length(gene_lists)

oddsRatio_DF$Sig <- oddsRatio_DF$pvalues < oddsRatio_DF$Bonferroni_sig_cutoff

oddsRatio_DF$log10_OR <- log10(oddsRatio_DF$odds_ratios)
oddsRatio_DF$log10_CI_upper <- log10(oddsRatio_DF$CI_upper)
oddsRatio_DF$log10_CI_lower <- log10(oddsRatio_DF$CI_lower)

oddsRatio_toPlot <- oddsRatio_DF %>%
  mutate(assay = factor(assay,
                        levels = c("RNA-seq","ChIP-seq","Both")))

gene_overlap_enrichment <- ggplot(oddsRatio_toPlot, aes(x = pvalues, y = assay)) +
  geom_point(data = subset(oddsRatio_toPlot, assay == 'RNA-seq'),
             size = 4, color = "#FEC102") +
  geom_point(data = subset(oddsRatio_toPlot, assay == 'ChIP-seq'),
             size = 4*((oddsRatio_toPlot$odds_ratios[2]/oddsRatio_toPlot$odds_ratios[1])/2), color = "#FE0200") +
  geom_point(data = subset(oddsRatio_toPlot, assay == 'Both'),
             size = 4*((oddsRatio_toPlot$odds_ratios[3]/oddsRatio_toPlot$odds_ratios[1])/2), color = "#FE7401") +
  theme_bw() +
  scale_x_continuous(limits = c(-0.08, 0.3))

ggsave(filename = '5enrich_CRP/initiation_figures/fig3j_gene_overlap_enrichment.png',
       plot = gene_overlap_enrichment,
       width = 75,
       height = 59,
       units = 'mm',
       dpi = 300)

## Add KEGG and Patric terms ---------------------------------------------------

kegg_all <- read.csv('genome_files_misc/KEGG_Mtb_annotations.csv')
kegg_toKeep <- kegg_all[,c("gene","Category")]
colnames(kegg_toKeep) <- c("geneID", 'kegg_category')

add_kegg <- left_join(full_DF, kegg_toKeep,
                      by = c("geneID" = "geneID"))

patric_all <- read.csv('genome_files_misc/PATRIC_Mtb_annotations.csv')
patric_toKeep <- patric_all[,c("gene","Superclass_x","Class_x","Subclass_x")]
colnames(patric_toKeep) <- c("geneID","patric_superclass",
                             'patric_class','patric_subclass')
all_terms <- left_join(add_kegg, patric_toKeep,
                       by = c("geneID" = "geneID"))

all_terms[is.na(all_terms)] <- 'None found'

## Add in CRISPRi and Tn-Seq essentiality calls --------------------------------

vulnerability <- read.csv("genome_files_misc/vulnerability_calls_CRISPRi.csv")

add_vuln <- left_join(all_terms, vulnerability[,c("locus_tag",
                                                "tnseq_ess",
                                                "crispr_ess")],
                      by = c('geneID' = 'locus_tag'))

## Collapse KEGG and PATRIC terms into one row ---------------------------------

add_vuln$KEGG_collapse <- 0
add_vuln$PATRIC_superCollapse <- 0
add_vuln$PATRIC_classCollapse <- 0
add_vuln$PATRIC_subCollapse <- 0

for (i in 1:nrow(add_vuln)) {

  coord_DF <- add_vuln[add_vuln$CFG_TSS_ID_higher == add_vuln$CFG_TSS_ID_higher[i],]
  
  kegg_categories <- unique(coord_DF$kegg_category)
  kegg_categories <- kegg_categories[!kegg_categories == 'None found']
  
  patric_superclasses <- unique(coord_DF$patric_superclass)
  patric_superclasses <- patric_superclasses[!patric_superclasses == 'None found']
  
  patric_classes <- unique(coord_DF$patric_class)
  patric_classes <- patric_classes[!patric_classes == 'None found']
  
  patric_subclasses <- unique(coord_DF$patric_subclass)
  patric_subclasses <- patric_subclasses[!patric_subclasses == 'None found']
  
  add_vuln$KEGG_collapse[i] <- paste0(kegg_categories, collapse = ', ')
  add_vuln$PATRIC_superCollapse[i] <- paste0(patric_superclasses, collapse = ', ')
  add_vuln$PATRIC_classCollapse[i] <- paste0(patric_classes, collapse = ', ')
  add_vuln$PATRIC_subCollapse[i] <- paste0(patric_subclasses, collapse = ', ')
  
}

old_KEGG_PATRIC_cols <- c("kegg_category", "patric_superclass",
                          "patric_class", "patric_subclass",
                          "CFG_geneID", "CFG_strand")

no_old_cols <- add_vuln[,!(colnames(add_vuln) %in% old_KEGG_PATRIC_cols)]

collapsed_DF <- unique(no_old_cols)

## Add in CFG motif scanning info ----------------------------------------------

CFG_results <- read.csv('5enrich_CRP/selectThreshold/allMotifResults/noCRP_CRP_min100plus30_20counts_DESeq2.csv')

motif_pvals <- unique(CFG_results$p.value)
# Identify p-value threshold based on discovered CFG motif
motif_pvals_final <- as.double(motif_pvals[!is.na(motif_pvals)])
motif_pvals_threshold <- max(motif_pvals_final)

fimo_all_CFG <- read.table('5enrich_CRP/initiation_figures/CFG_motif_scan/fimo.tsv',
                           header = TRUE)
fimo_all_CFG$top_scoring_motif <- "No"
coords <- unique(fimo_all_CFG$sequence_name)

for (i in 1:length(coords)) {
  
  fimo_subDF <- fimo_all_CFG[fimo_all_CFG$sequence_name == coords[i],]
  scores <- fimo_subDF$score
  fimo_all_CFG[(fimo_all_CFG$sequence_name == coords[i]) &
                 fimo_all_CFG$score == max(scores),"top_scoring_motif"] <- "Yes"
  
}

# For coordinates with a significant CRP motif:
fimo_pval_CFG <- fimo_all_CFG[fimo_all_CFG$p.value <= motif_pvals_threshold,]
collapsed_DF$Motif_in_CFG_TSS <- 'No'
collapsed_DF[collapsed_DF$CFG_TSS_ID_higher %in% fimo_pval_CFG$sequence_name,
             'Motif_in_CFG_TSS'] <- 'Yes'

# Only keep the top-scoring motif for each coordinate
fimo_top_motif_CFG <- unique(fimo_pval_CFG[fimo_pval_CFG$top_scoring_motif == 'Yes',])
fimo_top_motif_CFG$motif_center <- rowMeans(fimo_top_motif_CFG[,c("start","stop")]) - 101

colsToKeep <- c("sequence_name","motif_id","score","p.value",
                "matched_sequence","motif_center")
fimo_top_motif_CFG_clean <- fimo_top_motif_CFG[,colsToKeep]
colnames(fimo_top_motif_CFG_clean) <- paste("CFG_",
                                            colnames(fimo_top_motif_CFG_clean),
                                            sep='')

add_CFG_motif_info <- left_join(collapsed_DF,
                                fimo_top_motif_CFG_clean,
                                by = c("CFG_TSS_ID_higher" = "CFG_sequence_name"))
add_CFG_motif_info[is.na(add_CFG_motif_info)] <- ''
add_CFG_motif_info[add_CFG_motif_info$CFG_motif_id == 'KGUGAYBHASVUCAC',
                   'CFG_motif_id'] <- 'KGTGAYBHASVTCAC'


## Add in cellulo motif scanning info ----------------------------------------------

fimo_all_inCell <- read.table('5enrich_CRP/initiation_figures/inCellulo_motif_scan/fimo.tsv',
                           header = TRUE)
fimo_all_inCell$top_scoring_motif <- "No"
coords <- unique(fimo_all_inCell$sequence_name)

for (i in 1:length(coords)) {
  
  fimo_subDF <- fimo_all_inCell[fimo_all_inCell$sequence_name == coords[i],]
  scores <- fimo_subDF$score
  fimo_all_inCell[(fimo_all_inCell$sequence_name == coords[i]) &
                 fimo_all_inCell$score == max(scores),"top_scoring_motif"] <- "Yes"
  
}

# For coordinates with a significant CRP motif:
fimo_pval_inCell <- fimo_all_inCell[fimo_all_inCell$p.value <= motif_pvals_threshold,]
add_CFG_motif_info$Motif_in_inCell_TSS <- 'No'
add_CFG_motif_info[add_CFG_motif_info$inCell_TSS_ID %in% fimo_pval_inCell$sequence_name,
             'Motif_in_inCell_TSS'] <- 'Yes'

# Only keep the top-scoring motif for each coordinate
fimo_top_motif_inCell <- unique(fimo_pval_inCell[fimo_pval_inCell$top_scoring_motif == 'Yes',])
fimo_top_motif_inCell$motif_center <- rowMeans(fimo_top_motif_inCell[,c("start","stop")]) - 101

colsToKeep <- c("sequence_name","motif_id","score","p.value",
                "matched_sequence","motif_center")
fimo_top_motif_inCell_clean <- fimo_top_motif_inCell[,colsToKeep]
colnames(fimo_top_motif_inCell_clean) <- paste("inCell_",
                                            colnames(fimo_top_motif_inCell_clean),
                                            sep='')

add_inCell_motif_info <- left_join(add_CFG_motif_info,
                                fimo_top_motif_inCell_clean,
                                by = c("inCell_TSS_ID" = "inCell_sequence_name"))
add_inCell_motif_info[is.na(add_inCell_motif_info)] <- ''
add_inCell_motif_info[add_inCell_motif_info$inCell_motif_id == 'KGUGAYBHASVUCAC',
                   'inCell_motif_id'] <- 'KGTGAYBHASVTCAC'

## Add in gene names -----------------------------------------------------------

## Functions for GFF processing (by Michael DeJesus)
get_list_of_features = function(X, col=1){
  gene_ids = c()
  for (x in X){
    gid = get_gene_feature(x, col=col)
    gene_ids = c(gene_ids, gid )
  }
  return (gene_ids)
}
# Extract arbitrary feature based on column of split (";")
get_gene_feature = function(rawline, col=1){
  feature_section = strsplit(as.character(rawline), ";")[[1]][col]
  feature = strsplit(feature_section, "=")[[1]][2]
  return(feature)
}
# Extract gene ID
get_gene_id = function(rawline){
  gid_total = get_gene_feature(rawline, col=1)
  gid = strsplit(gid_total,'gene-')[[1]][2]
  return (gid)
}
# Extract gene name
get_gene_name = function(rawline){
  name = get_gene_feature(rawline, col=3)
  return (name)
}

GFF <- read.table('genome_files_misc/H37Rv_BD.gff', sep='\t', quote = "")

# Boolean method to get rid of mRNA, exon and CDS rows in gff file
ii_gene = GFF$V3 == "gene"
gene_DF = GFF[ii_gene,]
colnames(gene_DF) <- c("accessionNumber","source","type","start","end",
                       "dot","geneStrand","dot2","metadata")

IDs = get_list_of_features(gene_DF$metadata, col=1)
geneID_split = unlist(strsplit(IDs,'gene-'))
gene_DF$geneID <- geneID_split[geneID_split != '']
gene_DF$geneNameAll <- get_list_of_features(gene_DF$metadata, col = 2)
gene_DF$geneDescAll <- get_list_of_features(gene_DF$metadata, col = 4)

gene_names_only <- gene_DF[,c("geneID","geneNameAll","geneDescAll")]

add_gene_names <- left_join(add_inCell_motif_info,
                            gene_names_only,
                            by = c("geneID" = "geneID"))

## Organize CRP CFG supplementary table ------------------------------------------

coordinate_info <- c("geneID","geneNameAll",
                     "CFG_coord", "CFG_direction",
                     "CFG_TSS_ID_higher",'inCell_TSS_ID','TSS_distance')
CRP_effect_info <- c("CFG_log2FoldChange","RNAseq_log2FoldChange",
                     "CFG_padj","RNAseq_padj",
                     "CFG_CRP_effect","RNAseq_CRP_KO_effect",
                     "Hit_type")
gene_class_info <- c("geneDescAll","crispr_ess","tnseq_ess","KEGG_collapse",
                     "PATRIC_superCollapse","PATRIC_classCollapse","PATRIC_subCollapse")
CFG_sequence_info <- c('Motif_in_CFG_TSS',"CFG_motif_id","CFG_score","CFG_p.value",
                   "CFG_matched_sequence","CFG_motif_center",
                   'CFG_NT_strand')
inCell_sequence_info <- c("Motif_in_inCell_TSS","inCell_motif_id","inCell_score","inCell_p.value",
                       "inCell_matched_sequence","inCell_motif_center",
                       'inCell_NT_strand')

col_order <- c(coordinate_info, CRP_effect_info,
               gene_class_info, CFG_sequence_info, inCell_sequence_info)
cleaned_DF <- add_gene_names[,col_order]

sorted_DF <- cleaned_DF %>%
  arrange(Hit_type)

new_colnames <- c("Gene ID","Gene name",
                  "CFG TSS coordinate","CFG TSS strand",
                  "CFG TSS ID","In cellulo TSS ID","CFG vs. in cellulo TSS distance",
                  "CFG log2FC","CRP KO log2FC",
                  "CFG p-adj.","CRP KO p-adj.",
                  "CRP effect (CFG)","CRP effect (in cellulo RNA-seq)","Hit type",
                  "Gene description","CRISPRi call",
                  "Tn-seq call", "KEGG class", "PATRIC superclass",
                  "PATRIC class","PATRIC subclass","CRP motif found in CFG promoter",
                  "CFG de novo motif","CFG motif score","CFG motif p-val.",
                  "CFG promoter matched seq.","CFG motif center",
                  "CFG Non-template strand (-100 to +20)","CRP motif found in in cellulo promoter",
                  "CFG de novo motif","In cellulo motif score","In cellulo motif p-val.",
                  "In cellulo promoter matched seq.","In cellulo motif center",
                  "In cellulo Non-template strand (-100 to +20)")

colnames(sorted_DF) <- new_colnames
unique_export_DF <- unique(sorted_DF)

write.csv(unique_export_DF,
          '5enrich_CRP/initiation_figures/table_S5_CRP_geneLevel_supp_table.csv',
          row.names = FALSE)

## Motif distributions in RNA-seq ----------------------------------------------

CRP_repressed_col <- '#CC6677'
CRP_activated_col <- '#117733'

# Cortes motif
CRP_inCell <- cleaned_DF[cleaned_DF$RNAseq_padj != 'Not measured in RNA-seq',]
CRP_inCell$RNAseq_padj <- as.double(CRP_inCell$RNAseq_padj)
CRP_inCell$motif_center <- as.double(CRP_inCell$inCell_motif_center)
CRP_sig_inCell <- CRP_inCell[CRP_inCell$RNAseq_padj <= 0.05,]
CRP_inCell_toPlot <- unique(CRP_sig_inCell[,c('inCell_TSS_ID','motif_center',
                                              'RNAseq_CRP_KO_effect')])
inCell_motif_only <- CRP_inCell_toPlot[!is.na(CRP_inCell_toPlot$motif_center),]

inCell_motif_density <- ggplot(inCell_motif_only,
                               aes(x = motif_center)) +
  geom_density(data = subset(inCell_motif_only, RNAseq_CRP_KO_effect == 'Downregulated'),
               alpha = 0.5, fill = CRP_activated_col, size = 0.3) +
  geom_density(data = subset(inCell_motif_only, RNAseq_CRP_KO_effect == 'Upregulated'),
               alpha = 0.5, fill = CRP_repressed_col, size = 0.3) +
  theme_minimal() +
  scale_x_continuous(limits = c(-100,40), breaks = seq(-100, 40, 20)) +
  scale_y_continuous(limits = c(0,0.02), breaks = seq(0, 0.02, 0.005)) +
  theme(panel.grid.minor = element_blank())

ggsave(filename = '5enrich_CRP/initiation_figures/fig3h_inCellulo_motif_density.png',
       plot = inCell_motif_density,
       width = 83,
       height = 40,
       units = 'mm',
       dpi = 300)

nrow(inCell_motif_only[inCell_motif_only$RNAseq_CRP_KO_effect == 'Downregulated',])
nrow(inCell_motif_only[inCell_motif_only$RNAseq_CRP_KO_effect == 'Upregulated',])


## Upsetplot + bar chart to show CFG, RNA-seq and ChIP-seq overlap -------------

library(UpSetR)

CFG_only <- length(unique(cleaned_DF[(cleaned_DF$Hit_type == 'CFG only'),
                                       'geneID']))
CFG_only_withMotif <- length(unique(cleaned_DF[((cleaned_DF$Hit_type == 'CFG only') &
                                                      (cleaned_DF$Motif_in_CFG_TSS == 'Yes')),
                                                    'geneID']))
CFG_only_noMotif <- length(unique(cleaned_DF[((cleaned_DF$Hit_type == 'CFG only') &
                                                       (cleaned_DF$Motif_in_CFG_TSS == 'No')),
                                                    'geneID']))

ChIPseq_only <- length(unique(cleaned_DF[(cleaned_DF$Hit_type == 'ChIP-seq only'),
                                           'geneID']))
ChIPseq_only_withMotif <- length(unique(cleaned_DF[((cleaned_DF$Hit_type == 'ChIP-seq only') &
                                                       (cleaned_DF$Motif_in_inCell_TSS == 'Yes')),
                                                    'geneID']))
ChIPseq_only_noMotif <- length(unique(cleaned_DF[((cleaned_DF$Hit_type == 'ChIP-seq only') &
                                                     (cleaned_DF$Motif_in_inCell_TSS == 'No')),
                                                  'geneID']))

RNAseq_only <- length(unique(cleaned_DF[(cleaned_DF$Hit_type == 'RNA-seq only'),
                                           'geneID']))
RNAseq_only_withMotif <- length(unique(cleaned_DF[((cleaned_DF$Hit_type == 'RNA-seq only') &
                                                           (cleaned_DF$Motif_in_inCell_TSS == 'Yes')),
                                                        'geneID']))
RNAseq_only_noMotif <- length(unique(cleaned_DF[((cleaned_DF$Hit_type == 'RNA-seq only') &
                                                         (cleaned_DF$Motif_in_inCell_TSS == 'No')),
                                                      'geneID']))


CFG_and_ChIPseq <- length(unique(cleaned_DF[(cleaned_DF$Hit_type == 'CFG and ChIP-seq'),
                                           'geneID']))
CFG_and_ChIPseq_with_inCell_motif <- length(unique(cleaned_DF[((cleaned_DF$Hit_type == 'CFG and ChIP-seq') &
                                                          (cleaned_DF$Motif_in_inCell_TSS == 'Yes')),
                                                       'geneID']))
CFG_and_ChIPseq_no_inCell_motif <- length(unique(cleaned_DF[((cleaned_DF$Hit_type == 'CFG and ChIP-seq') &
                                                        (cleaned_DF$Motif_in_inCell_TSS == 'No')),
                                                     'geneID']))
CFG_and_ChIPseq_with_CFG_motif <- length(unique(cleaned_DF[((cleaned_DF$Hit_type == 'CFG and ChIP-seq') &
                                                                      (cleaned_DF$Motif_in_CFG_TSS == 'Yes')),
                                                                   'geneID']))
CFG_and_ChIPseq_no_CFG_motif <- length(unique(cleaned_DF[((cleaned_DF$Hit_type == 'CFG and ChIP-seq') &
                                                                    (cleaned_DF$Motif_in_CFG_TSS == 'No')),
                                                                 'geneID']))

CFG_and_RNAseq <- length(unique(cleaned_DF[(cleaned_DF$Hit_type == 'CFG and RNA-seq'),
                                           'geneID']))
CFG_and_RNAseq_with_inCell_motif <- length(unique(cleaned_DF[((cleaned_DF$Hit_type == 'CFG and RNA-seq') &
                                                                      (cleaned_DF$Motif_in_inCell_TSS == 'Yes')),
                                                                   'geneID']))
CFG_and_RNAseq_no_inCell_motif <- length(unique(cleaned_DF[((cleaned_DF$Hit_type == 'CFG and RNA-seq') &
                                                                    (cleaned_DF$Motif_in_inCell_TSS == 'No')),
                                                                 'geneID']))
CFG_and_RNAseq_with_CFG_motif <- length(unique(cleaned_DF[((cleaned_DF$Hit_type == 'CFG and RNA-seq') &
                                                                   (cleaned_DF$Motif_in_CFG_TSS == 'Yes')),
                                                                'geneID']))
CFG_and_RNAseq_no_CFG_motif <- length(unique(cleaned_DF[((cleaned_DF$Hit_type == 'CFG and RNA-seq') &
                                                                 (cleaned_DF$Motif_in_CFG_TSS == 'No')),
                                                              'geneID']))

ChIP_and_RNAseq <- length(unique(cleaned_DF[(cleaned_DF$Hit_type == 'RNA-seq and ChIP-seq'),
                                             'geneID']))
ChIP_and_RNAseq_withMotif <- length(unique(cleaned_DF[((cleaned_DF$Hit_type == 'RNA-seq and ChIP-seq') &
                                                                     (cleaned_DF$Motif_in_inCell_TSS == 'Yes')),
                                                                  'geneID']))
ChIP_and_RNAseq_noMotif <- length(unique(cleaned_DF[((cleaned_DF$Hit_type == 'RNA-seq and ChIP-seq') &
                                                                   (cleaned_DF$Motif_in_inCell_TSS == 'No')),
                                                                'geneID']))

all_three <- length(unique(cleaned_DF[(cleaned_DF$Hit_type == 'CFG, ChIP-seq and RNA-seq'),
                                             'geneID']))
all_three_with_inCell_motif <- length(unique(cleaned_DF[((cleaned_DF$Hit_type == 'CFG, ChIP-seq and RNA-seq') &
                                                                      (cleaned_DF$Motif_in_inCell_TSS == 'Yes')),
                                                                   'geneID']))
all_three_no_inCell_motif <- length(unique(cleaned_DF[((cleaned_DF$Hit_type == 'CFG, ChIP-seq and RNA-seq') &
                                                                    (cleaned_DF$Motif_in_inCell_TSS == 'No')),
                                                                 'geneID']))
all_three_with_CFG_motif <- length(unique(cleaned_DF[((cleaned_DF$Hit_type == 'CFG, ChIP-seq and RNA-seq') &
                                                                   (cleaned_DF$Motif_in_CFG_TSS == 'Yes')),
                                                                'geneID']))
all_three_no_CFG_motif <- length(unique(cleaned_DF[((cleaned_DF$Hit_type == 'CFG, ChIP-seq and RNA-seq') &
                                                                 (cleaned_DF$Motif_in_CFG_TSS == 'No')),
                                                              'geneID']))

input <- c(
  CFG = CFG_only,
  ChIPseq = ChIPseq_only,
  RNAseq = RNAseq_only,
  "CFG&ChIPseq" = CFG_and_ChIPseq,
  "CFG&RNAseq" = CFG_and_RNAseq,
  "ChIPseq&RNAseq" = ChIP_and_RNAseq,
  "CFG&ChIPseq&RNAseq" = all_three
)

png(paste("5enrich_CRP/initiation_figures/fig3i_CRP_overlap_upset_noConflicts.png",sep=''),
    res = 300,
    width = 12,
    height = 12,
    units = 'in')
upset(fromExpression(input), 
      nintersects = length(input), 
      nsets = length(input), 
      order.by = "freq", 
      decreasing = T, 
      mb.ratio = c(0.6, 0.4),
      number.angles = 0, 
      text.scale = 1.1, 
      point.size = 2.8, 
      line.size = 1
)
dev.off()

library(ggbreak)
library(ggpattern)

## in-cell only DF--------------------------------------------------------------

# Figs color
yellow = '#FEC102'
red = '#FE0200'
orange = '#FE7401'
green = '#53B04C'
blue = '#5067EF'
brown = '#A66A2E'
purple = '#C58BF7'

inCell_only_DF <- data.frame(hit_condition = c("RNAseq","RNAseq",
                                                    "ChIPseq","ChIPseq",
                                               'ChIP&RNAseq','ChIP&RNAseq'),
                                  motif = c("Yes","No","Yes","No","Yes","No"),
                             num_hits = c(RNAseq_only_withMotif,
                                          RNAseq_only_noMotif,
                                          ChIPseq_only_withMotif,
                                          ChIPseq_only_noMotif,
                                          ChIP_and_RNAseq_withMotif,
                                          ChIP_and_RNAseq_noMotif))

inCell_only_DF <- inCell_only_DF %>%
  mutate(hit_condition = factor(hit_condition,
                                levels = c("RNAseq","ChIPseq",
                                           "ChIP&RNAseq"))) %>%
  mutate(motif = factor(motif,
                        levels = c("Yes","No")))

inCell_only_bar <- ggplot(inCell_only_DF, aes(x = hit_condition,
                                      y = num_hits, fill = motif)) +
  geom_col_pattern(aes(fill = hit_condition,
                       pattern_density = motif),colour='black',
                   pattern='stripe',
                   width = 0.5) +
  scale_fill_manual(values=c('lightgrey', 'lightgrey', 'lightgrey')) +
  theme_minimal() +
  theme(legend.position = 'none') +
  scale_y_break(c(200,800), scales = 2, space = 1) +
  scale_pattern_density_manual(values = c(0, 0.1,
                                          0, 0.1,
                                          0, 0.1)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

ggsave(filename = '5enrich_CRP/initiation_figures/fig3i_inCell_only_bar.png',
       plot = inCell_only_bar,
       width = 188,
       height = 100,
       units = 'mm',
       dpi = 300)

write.csv(inCell_only_DF, '5enrich_CRP/initiation_figures/fig3i_inCell_only_geneCounts.csv')

## CFG overlap DF --------------------------------------------------------------

CFG_overlap_DF <- data.frame(hit_condition = c("CFG","CFG",
                                               "CFG&RNAseq","CFG&RNAseq",
                                               'CFG&ChIP','CFG&ChIP',
                                               'CFG&ChIP&RNAseq','CFG&ChIP&RNAseq'),
                             motif = c("Yes","No","Yes","No","Yes","No",
                                       "Yes","No"),
                             num_hits = c(CFG_only_withMotif,
                                          CFG_only_noMotif,
                                          CFG_and_RNAseq_with_CFG_motif,
                                          CFG_and_RNAseq_no_CFG_motif,
                                          CFG_and_ChIPseq_with_CFG_motif,
                                          CFG_and_ChIPseq_no_CFG_motif,
                                          all_three_with_CFG_motif,
                                          all_three_no_CFG_motif))

CFG_overlap_DF <- CFG_overlap_DF %>%
  mutate(hit_condition = factor(hit_condition,
                                levels = c("CFG","CFG&RNAseq",
                                           "CFG&ChIP",'CFG&ChIP&RNAseq'))) %>%
  mutate(motif = factor(motif,
                        levels = c("Yes","No")))

CFG_overlap_bar <- ggplot(CFG_overlap_DF, aes(x = hit_condition,
                                              y = num_hits, fill = motif)) +
  geom_col_pattern(aes(fill = hit_condition,
                       pattern_density = motif),colour='black',
                   pattern='stripe',
                   width = 0.5) +
  scale_fill_manual(values=c(blue, yellow,
                             red,orange)) +
  theme_minimal() +
  theme(legend.position = 'none') +
  scale_pattern_density_manual(values = c(0, 0.5,
                                          0, 0.5,
                                          0, 0.5,
                                          0, 0.5)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_y_continuous(limits = c(0,20), breaks = seq(0,20,5))

ggsave(filename = '5enrich_CRP/initiation_figures/fig3i_CFG_overlap_bar.png',
       plot = CFG_overlap_bar,
       width = 188,
       height = 100,
       units = 'mm',
       dpi = 300)

write.csv(CFG_overlap_DF, '5enrich_CRP/initiation_figures/fig3i_CFG_overlap_geneCounts.csv')

## End -------------------------------------------------------------------------

