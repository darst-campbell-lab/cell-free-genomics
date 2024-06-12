#### Candidate terminator validation
#### Authors: Ruby Froom
#### Date: February 23, 2024

#### Libraries ----------------------------------------------------------------

library(dplyr)
library(stringr)
library(stringi)
library(Biostrings)
library(seqinr)
library(GenomicRanges)
library(seqinr)
library(Biostrings)
library(Gviz)
library(rtracklayer)


### Gviz -----------------------------------------------------------------------
## Functions for GFF processing
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
  gid = get_gene_feature(rawline, col=1)
  return (gid)
}
# Extract gene name
get_gene_name = function(rawline){
  name = get_gene_feature(rawline, col=2)
  return (name)
}
options(ucscChromosomeNames = FALSE)
genomeAxis <- GenomeAxisTrack(name = "Eco_Mtb")

## Generate gff for Gviz
Mtb_GFF <- read.table('genome_files_misc/H37Rv_BD.gff', sep='\t')
Eco_GFF <- read.table('genome_files_misc/Eco.gff', sep='\t', quote = "")
Mtb_second_GFF <- Mtb_GFF
Mtb_second_GFF$V4 <- Mtb_GFF$V4 + 4641652
Mtb_second_GFF$V5 <- Mtb_GFF$V5 + 4641652
ii_gene = Mtb_second_GFF$V3 == "gene"
geneX = Mtb_second_GFF[ii_gene,]
IDs = get_list_of_features(geneX$V9, col=1)
Names = get_list_of_features(geneX$V9, col=2)
Description = get_list_of_features(geneX$V9, col=4)
geneX$gene_id <- IDs
geneX$gene_name <- Names
geneX$desc <- Description
column_names <- c('V1', 'V4', 'V5', 'V7', 'gene_id', 'gene_name', 'desc')
Mtb_col_curated <- geneX[,column_names]
Mtb_col_curated$V1 <- 'Eco_Mtb'
Mtb_col_curated[Mtb_col_curated$gene_name == '-','gene_name'] <- Mtb_col_curated[Mtb_col_curated$gene_name == '-','gene_id'] 
eco_feature <- data.frame(V1 = 'Eco_Mtb',
                          V4 = '1',
                          V5 = '4641652',
                          V7 = '+',
                          gene_id = 'Eco_genome',
                          gene_name = 'Eco_genome',
                          desc = 'Eco_genome')
final_DF <- rbind(eco_feature, Mtb_col_curated)
Mtb_GRanges <- makeGRangesFromDataFrame(final_DF,
                                        seqnames.field = 'V1',
                                        start.field = 'V4',
                                        end.field = 'V5',
                                        strand.field = 'V7',
                                        keep.extra.columns = TRUE)

GFF_Annot <- AnnotationTrack(Mtb_GRanges, labelAllFeatures = TRUE)
feature(GFF_Annot) <- Mtb_GRanges$gene_name

bs <- readDNAStringSet('genome_files_misc/Eco_Mtb_genome.fasta')
sTrack <- SequenceTrack(bs)
sTrack_comp <- SequenceTrack(complement(bs))

noTF1_bw <- import.bw("3enrich_NusAG/identifyEnrichedEnds/bigWigs_vis/noTF1_plus_vis.bw")
noTF2_bw <- import.bw("3enrich_NusAG/identifyEnrichedEnds/bigWigs_vis/noTF2_plus_vis.bw")
noTF3_bw <- import.bw("3enrich_NusAG/identifyEnrichedEnds/bigWigs_vis/noTF3_plus_vis.bw")
NusG1_bw <- import.bw("3enrich_NusAG/identifyEnrichedEnds/bigWigs_vis/NusG1_plus_vis.bw")
NusG2_bw <- import.bw("3enrich_NusAG/identifyEnrichedEnds/bigWigs_vis/NusG2_plus_vis.bw")
NusG3_bw <- import.bw("3enrich_NusAG/identifyEnrichedEnds/bigWigs_vis/NusG3_plus_vis.bw")
NusA1_bw <- import.bw("3enrich_NusAG/identifyEnrichedEnds/bigWigs_vis/NusA1_plus_vis.bw")
NusA2_bw <- import.bw("3enrich_NusAG/identifyEnrichedEnds/bigWigs_vis/NusA2_plus_vis.bw")
NusA3_bw <- import.bw("3enrich_NusAG/identifyEnrichedEnds/bigWigs_vis/NusA3_plus_vis.bw")
NusA_NusG1_bw <- import.bw("3enrich_NusAG/identifyEnrichedEnds/bigWigs_vis/NusA_NusG1_plus_vis.bw")
NusA_NusG2_bw <- import.bw("3enrich_NusAG/identifyEnrichedEnds/bigWigs_vis/NusA_NusG2_plus_vis.bw")
NusA_NusG3_bw <- import.bw("3enrich_NusAG/identifyEnrichedEnds/bigWigs_vis/NusA_NusG3_plus_vis.bw")

## pitB ------------------------------------------------------------------------

coord = 7196579
start_coord = coord - 50
end_coord = coord + 10
alpha_val = 0.7
noTF_col = "grey"
NusG_col = "#ED2939"
NusA_col = "#006ee6"
NusA_NusG_col = "purple"

noTF1_DT <- DataTrack(noTF1_bw, start = start_coord, end = end_coord,
                      name = "noTF", groups = factor("noTF",
                                                     levels = c("noTF","NusG","NusA","NusA + NusG")),
                      alpha = alpha_val, col = noTF_col)
noTF2_DT <- DataTrack(noTF2_bw, start = start_coord, end = end_coord,
                      name = "noTF", groups = factor("noTF",
                                                     levels = c("noTF","NusG","NusA","NusA + NusG")),
                      alpha = alpha_val, col = noTF_col)
noTF3_DT <- DataTrack(noTF3_bw, start = start_coord, end = end_coord,
                      name = "noTF", groups = factor("noTF",
                                                     levels = c("noTF","NusG","NusA","NusA + NusG")),
                      alpha = alpha_val, col = noTF_col)
NusG1_DT <- DataTrack(NusG1_bw, start = start_coord, end = end_coord,
                      name = "NusG", groups = factor("NusG",
                                                     levels = c("noTF","NusG","NusA","NusA + NusG")),
                      alpha = alpha_val, col = NusG_col)
NusG2_DT <- DataTrack(NusG2_bw, start = start_coord, end = end_coord,
                      name = "NusG", groups = factor("NusG",
                                                     levels = c("noTF","NusG","NusA","NusA + NusG")),
                      alpha = alpha_val, col = NusG_col)
NusG3_DT <- DataTrack(NusG3_bw, start = start_coord, end = end_coord,
                      name = "NusG", groups = factor("NusG",
                                                     levels = c("noTF","NusG","NusA","NusA + NusG")),
                      alpha = alpha_val, col = NusG_col)
NusA1_DT <- DataTrack(NusA1_bw, start = start_coord, end = end_coord,
                      name = "NusA", groups = factor("NusA",
                                                     levels = c("noTF","NusG","NusA","NusA + NusG")),
                      alpha = alpha_val, col = NusA_col)
NusA2_DT <- DataTrack(NusA2_bw, start = start_coord, end = end_coord,
                      name = "NusA", groups = factor("NusA",
                                                     levels = c("noTF","NusG","NusA","NusA + NusG")),
                      alpha = alpha_val, col = NusA_col)
NusA3_DT <- DataTrack(NusA3_bw, start = start_coord, end = end_coord,
                      name = "NusA", groups = factor("NusA",
                                                     levels = c("noTF","NusG","NusA","NusA + NusG")),
                      alpha = alpha_val, col = NusA_col)
NusA_NusG1_DT <- DataTrack(NusA_NusG1_bw, start = start_coord, end = end_coord,
                           name = "NusA + NusG", groups = factor("NusA + NusG",
                                                                 levels = c("noTF","NusG","NusA","NusA + NusG")),
                           alpha = alpha_val, col = NusA_NusG_col)
NusA_NusG2_DT <- DataTrack(NusA_NusG2_bw, start = start_coord, end = end_coord,
                           name = "NusA + NusG", groups = factor("NusA + NusG",
                                                                 levels = c("noTF","NusG","NusA","NusA + NusG")),
                           alpha = alpha_val, col = NusA_NusG_col)
NusA_NusG3_DT <- DataTrack(NusA_NusG3_bw, start = start_coord, end = end_coord,
                           name = "NusA + NusG", groups = factor("NusA + NusG",
                                                                 levels = c("noTF","NusG","NusA","NusA + NusG")),
                           alpha = alpha_val, col = NusA_NusG_col)

## Separate tracks by sample
noTF_overlay <- OverlayTrack(trackList = list(noTF1_DT, noTF2_DT, noTF3_DT))
NusG_overlay <- OverlayTrack(trackList = list(NusG1_DT, NusG2_DT, NusG3_DT))
NusA_overlay <- OverlayTrack(trackList = list(NusA1_DT, NusA2_DT, NusA3_DT))
NusA_NusG_overlay <- OverlayTrack(trackList = list(NusA_NusG1_DT, NusA_NusG2_DT, NusA_NusG3_DT))

colForLetters <- c(A = "black", C = "black",
                   T = "red", G = "black")

png("3enrich_NusAG/termination_figures/fig4d_pitB_term.png",
    res = 300,
    width = 12,
    height = 3,
    units = 'in')
plotTracks(list(noTF_overlay, NusA_overlay, NusG_overlay, NusA_NusG_overlay,
                sTrack),
           from = start_coord,
           to = end_coord,
           type = c("hist"),
           add53 = FALSE,
           add35 = FALSE,
           ylim = range(0,150000),
           legend = FALSE,
           col.axis = "black",
           background.title = "transparent",
           fontcolor.item = "black",
           fontcolor = colForLetters)
dev.off()

## icd1 --------------------------------------------------------------------

coord = 8366386
start_coord = coord - 10
end_coord = coord + 50
alpha_val = 0.7
noTF_col = "grey"
NusG_col = "#ED2939"
NusA_col = "#006ee6"
NusA_NusG_col = "purple"

noTF1_DT <- DataTrack(noTF1_bw, start = start_coord, end = end_coord,
                      name = "noTF", groups = factor("noTF",
                                                     levels = c("noTF","NusG","NusA","NusA + NusG")),
                      alpha = alpha_val, col = noTF_col)
noTF2_DT <- DataTrack(noTF2_bw, start = start_coord, end = end_coord,
                      name = "noTF", groups = factor("noTF",
                                                     levels = c("noTF","NusG","NusA","NusA + NusG")),
                      alpha = alpha_val, col = noTF_col)
noTF3_DT <- DataTrack(noTF3_bw, start = start_coord, end = end_coord,
                      name = "noTF", groups = factor("noTF",
                                                     levels = c("noTF","NusG","NusA","NusA + NusG")),
                      alpha = alpha_val, col = noTF_col)
NusG1_DT <- DataTrack(NusG1_bw, start = start_coord, end = end_coord,
                      name = "NusG", groups = factor("NusG",
                                                     levels = c("noTF","NusG","NusA","NusA + NusG")),
                      alpha = alpha_val, col = NusG_col)
NusG2_DT <- DataTrack(NusG2_bw, start = start_coord, end = end_coord,
                      name = "NusG", groups = factor("NusG",
                                                     levels = c("noTF","NusG","NusA","NusA + NusG")),
                      alpha = alpha_val, col = NusG_col)
NusG3_DT <- DataTrack(NusG3_bw, start = start_coord, end = end_coord,
                      name = "NusG", groups = factor("NusG",
                                                     levels = c("noTF","NusG","NusA","NusA + NusG")),
                      alpha = alpha_val, col = NusG_col)
NusA1_DT <- DataTrack(NusA1_bw, start = start_coord, end = end_coord,
                      name = "NusA", groups = factor("NusA",
                                                     levels = c("noTF","NusG","NusA","NusA + NusG")),
                      alpha = alpha_val, col = NusA_col)
NusA2_DT <- DataTrack(NusA2_bw, start = start_coord, end = end_coord,
                      name = "NusA", groups = factor("NusA",
                                                     levels = c("noTF","NusG","NusA","NusA + NusG")),
                      alpha = alpha_val, col = NusA_col)
NusA3_DT <- DataTrack(NusA3_bw, start = start_coord, end = end_coord,
                      name = "NusA", groups = factor("NusA",
                                                     levels = c("noTF","NusG","NusA","NusA + NusG")),
                      alpha = alpha_val, col = NusA_col)
NusA_NusG1_DT <- DataTrack(NusA_NusG1_bw, start = start_coord, end = end_coord,
                           name = "NusA + NusG", groups = factor("NusA + NusG",
                                                                 levels = c("noTF","NusG","NusA","NusA + NusG")),
                           alpha = alpha_val, col = NusA_NusG_col)
NusA_NusG2_DT <- DataTrack(NusA_NusG2_bw, start = start_coord, end = end_coord,
                           name = "NusA + NusG", groups = factor("NusA + NusG",
                                                                 levels = c("noTF","NusG","NusA","NusA + NusG")),
                           alpha = alpha_val, col = NusA_NusG_col)
NusA_NusG3_DT <- DataTrack(NusA_NusG3_bw, start = start_coord, end = end_coord,
                           name = "NusA + NusG", groups = factor("NusA + NusG",
                                                                 levels = c("noTF","NusG","NusA","NusA + NusG")),
                           alpha = alpha_val, col = NusA_NusG_col)

## Separate tracks by sample
noTF_overlay <- OverlayTrack(trackList = list(noTF1_DT, noTF2_DT, noTF3_DT))
NusG_overlay <- OverlayTrack(trackList = list(NusG1_DT, NusG2_DT, NusG3_DT))
NusA_overlay <- OverlayTrack(trackList = list(NusA1_DT, NusA2_DT, NusA3_DT))
NusA_NusG_overlay <- OverlayTrack(trackList = list(NusA_NusG1_DT, NusA_NusG2_DT, NusA_NusG3_DT))

colForLetters <- c(A = "black", C = "black",
                   T = "red", G = "black")

png("3enrich_NusAG/termination_figures/fig4e_icd1_term.png",
    res = 300,
    width = 12,
    height = 6,
    units = 'in')
plotTracks(list(noTF_overlay, NusA_overlay, NusG_overlay, NusA_NusG_overlay,
                sTrack_comp),
           from = start_coord,
           to = end_coord,
           type = c("hist"),
           add53 = FALSE,
           add35 = FALSE,
           ylim = range(0,150000),
           legend = FALSE,
           col.axis = "black",
           background.title = "transparent",
           fontcolor.item = "black",
           fontcolor = colForLetters)
dev.off()

## rrf -------------------------------------------------------------------------

coord = 6118708
start_coord = coord - 50
end_coord = coord + 10
alpha_val = 0.7
noTF_col = "grey"
NusG_col = "#ED2939"
NusA_col = "#006ee6"
NusA_NusG_col = "purple"

noTF1_bw <- import.bw("3enrich_NusAG/identifyEnrichedEnds/bigWigs_vis/noTF1_minus_vis.bw")
noTF2_bw <- import.bw("3enrich_NusAG/identifyEnrichedEnds/bigWigs_vis/noTF2_minus_vis.bw")
noTF3_bw <- import.bw("3enrich_NusAG/identifyEnrichedEnds/bigWigs_vis/noTF3_minus_vis.bw")
NusG1_bw <- import.bw("3enrich_NusAG/identifyEnrichedEnds/bigWigs_vis/NusG1_minus_vis.bw")
NusG2_bw <- import.bw("3enrich_NusAG/identifyEnrichedEnds/bigWigs_vis/NusG2_minus_vis.bw")
NusG3_bw <- import.bw("3enrich_NusAG/identifyEnrichedEnds/bigWigs_vis/NusG3_minus_vis.bw")
NusA1_bw <- import.bw("3enrich_NusAG/identifyEnrichedEnds/bigWigs_vis/NusA1_minus_vis.bw")
NusA2_bw <- import.bw("3enrich_NusAG/identifyEnrichedEnds/bigWigs_vis/NusA2_minus_vis.bw")
NusA3_bw <- import.bw("3enrich_NusAG/identifyEnrichedEnds/bigWigs_vis/NusA3_minus_vis.bw")
NusA_NusG1_bw <- import.bw("3enrich_NusAG/identifyEnrichedEnds/bigWigs_vis/NusA_NusG1_minus_vis.bw")
NusA_NusG2_bw <- import.bw("3enrich_NusAG/identifyEnrichedEnds/bigWigs_vis/NusA_NusG2_minus_vis.bw")
NusA_NusG3_bw <- import.bw("3enrich_NusAG/identifyEnrichedEnds/bigWigs_vis/NusA_NusG3_minus_vis.bw")

noTF1_DT <- DataTrack(noTF1_bw, start = start_coord, end = end_coord,
                      name = "noTF", groups = factor("noTF",
                                                     levels = c("noTF","NusG","NusA","NusA + NusG")),
                      alpha = alpha_val, col = noTF_col)
noTF2_DT <- DataTrack(noTF2_bw, start = start_coord, end = end_coord,
                      name = "noTF", groups = factor("noTF",
                                                     levels = c("noTF","NusG","NusA","NusA + NusG")),
                      alpha = alpha_val, col = noTF_col)
noTF3_DT <- DataTrack(noTF3_bw, start = start_coord, end = end_coord,
                      name = "noTF", groups = factor("noTF",
                                                     levels = c("noTF","NusG","NusA","NusA + NusG")),
                      alpha = alpha_val, col = noTF_col)
NusG1_DT <- DataTrack(NusG1_bw, start = start_coord, end = end_coord,
                      name = "NusG", groups = factor("NusG",
                                                     levels = c("noTF","NusG","NusA","NusA + NusG")),
                      alpha = alpha_val, col = NusG_col)
NusG2_DT <- DataTrack(NusG2_bw, start = start_coord, end = end_coord,
                      name = "NusG", groups = factor("NusG",
                                                     levels = c("noTF","NusG","NusA","NusA + NusG")),
                      alpha = alpha_val, col = NusG_col)
NusG3_DT <- DataTrack(NusG3_bw, start = start_coord, end = end_coord,
                      name = "NusG", groups = factor("NusG",
                                                     levels = c("noTF","NusG","NusA","NusA + NusG")),
                      alpha = alpha_val, col = NusG_col)
NusA1_DT <- DataTrack(NusA1_bw, start = start_coord, end = end_coord,
                      name = "NusA", groups = factor("NusA",
                                                     levels = c("noTF","NusG","NusA","NusA + NusG")),
                      alpha = alpha_val, col = NusA_col)
NusA2_DT <- DataTrack(NusA2_bw, start = start_coord, end = end_coord,
                      name = "NusA", groups = factor("NusA",
                                                     levels = c("noTF","NusG","NusA","NusA + NusG")),
                      alpha = alpha_val, col = NusA_col)
NusA3_DT <- DataTrack(NusA3_bw, start = start_coord, end = end_coord,
                      name = "NusA", groups = factor("NusA",
                                                     levels = c("noTF","NusG","NusA","NusA + NusG")),
                      alpha = alpha_val, col = NusA_col)
NusA_NusG1_DT <- DataTrack(NusA_NusG1_bw, start = start_coord, end = end_coord,
                           name = "NusA + NusG", groups = factor("NusA + NusG",
                                                                 levels = c("noTF","NusG","NusA","NusA + NusG")),
                           alpha = alpha_val, col = NusA_NusG_col)
NusA_NusG2_DT <- DataTrack(NusA_NusG2_bw, start = start_coord, end = end_coord,
                           name = "NusA + NusG", groups = factor("NusA + NusG",
                                                                 levels = c("noTF","NusG","NusA","NusA + NusG")),
                           alpha = alpha_val, col = NusA_NusG_col)
NusA_NusG3_DT <- DataTrack(NusA_NusG3_bw, start = start_coord, end = end_coord,
                           name = "NusA + NusG", groups = factor("NusA + NusG",
                                                                 levels = c("noTF","NusG","NusA","NusA + NusG")),
                           alpha = alpha_val, col = NusA_NusG_col)

## Separate tracks by sample
noTF_overlay <- OverlayTrack(trackList = list(noTF1_DT, noTF2_DT, noTF3_DT))
NusG_overlay <- OverlayTrack(trackList = list(NusG1_DT, NusG2_DT, NusG3_DT))
NusA_overlay <- OverlayTrack(trackList = list(NusA1_DT, NusA2_DT, NusA3_DT))
NusA_NusG_overlay <- OverlayTrack(trackList = list(NusA_NusG1_DT, NusA_NusG2_DT, NusA_NusG3_DT))

colForLetters <- c(A = "black", C = "black",
                   T = "red", G = "black")

png("3enrich_NusAG/termination_figures/fig4f_rrf_term.png",
    res = 300,
    width = 12,
    height = 6,
    units = 'in')
plotTracks(list(noTF_overlay, NusA_overlay, NusG_overlay, NusA_NusG_overlay,
                sTrack),
           from = start_coord,
           to = end_coord,
           type = c("hist"),
           add53 = FALSE,
           add35 = FALSE,
           ylim = range(0,800000),
           legend = FALSE,
           col.axis = "black",
           background.title = "transparent",
           fontcolor.item = "black",
           fontcolor = colForLetters)
dev.off()

### Term assay results ---------------------------------------------------------

term_values <- read.csv('3enrich_NusAG/termination_figures/singleDNA_term_values.csv')

term_values <- term_values[term_values$TF != 'noTF',]

term_values <- term_values %>%
  mutate(TF = factor(TF,
                     levels = c('NusA','NusG','NusAG')))

term_values$TE_norm_each <- term_values$TE_FC - 1
term_values$avg_TE_norm <- term_values$avg_TE_FC - 1

pitB_values <- term_values[term_values$terminator == 'pitB',]
icd1_values <- term_values[term_values$terminator == 'icd1',]
rrf_values <- term_values[term_values$terminator == 'rrf',]

noTF_col = "grey"
NusG_col = "#E26B7A"
NusA_col = "#5F8ADC"
NusA_NusG_col = "#BD67F3"

noTF_col_dark = "darkgrey"
NusG_col_dark = "#ED2939"
NusA_col_dark = "#006ee6"
NusA_NusG_col_dark = "purple"

pitB_fig <- ggplot(pitB_values, aes(x = TF,
                                    y = avg_TE_norm,
                                    fill = TF)) +
  geom_bar(stat = "identity", position = "dodge",
           width = 0.5, linewidth = 1, color = 'black') +
  theme_minimal() +
  geom_point(aes(x = TF, y = TE_norm_each, color = TF),
             pitB_values, size = 7) +
  scale_fill_manual(values = c(NusA_col,
                               NusG_col, NusA_NusG_col)) +
  scale_color_manual(values = c('black','black','black')) +
  geom_errorbar(aes(ymin = avg_TE_norm - SD_TE_FC,
                    ymax = avg_TE_norm + SD_TE_FC),
                width = 0.4, linewidth = 2,
                position = position_dodge(0.5)) +
  theme(legend.position = "none")

ggsave(filename = '3enrich_NusAG/termination_figures/extData_fig7e_pitB_bargraph_zeroed.png',
       plot = pitB_fig,
       width = 300,
       height = 300,
       units = 'mm',
       dpi = 300)

icd1_fig <- ggplot(icd1_values, aes(x = TF,
                                    y = avg_TE_norm,
                                    fill = TF)) +
  geom_bar(stat = "identity", position = "dodge",
           width = 0.5, linewidth = 1, color = 'black') +
  theme_minimal() +
  geom_point(aes(x = TF, y = TE_norm_each, color = TF),
             icd1_values, size = 7) +
  scale_fill_manual(values = c(NusA_col,
                               NusG_col, NusA_NusG_col)) +
  scale_color_manual(values = c('black','black','black')) +
  geom_errorbar(aes(ymin = avg_TE_norm - SD_TE_FC,
                    ymax = avg_TE_norm + SD_TE_FC),
                width = 0.2, linewidth = 2,
                position = position_dodge(0.5)) +
  theme(legend.position = "none")

ggsave(filename = '3enrich_NusAG/termination_figures/extData_fig7e_icd1_bargraph_zeroed.png',
       plot = icd1_fig,
       width = 300,
       height = 300,
       units = 'mm',
       dpi = 300)

rrf_fig <- ggplot(rrf_values, aes(x = TF,
                                  y = avg_TE_norm,
                                  fill = TF)) +
  geom_bar(stat = "identity", position = "dodge",
           width = 0.5, linewidth = 1, color = 'black') +
  theme_minimal() +
  geom_point(aes(x = TF, y = TE_norm_each, color = TF),
             rrf_values, size = 7) +
  scale_fill_manual(values = c(NusA_col,
                               NusG_col, NusA_NusG_col)) +
  scale_color_manual(values = c('black','black','black')) +
  geom_errorbar(aes(ymin = avg_TE_norm - SD_TE_FC,
                    ymax = avg_TE_norm + SD_TE_FC),
                width = 0.4, linewidth = 2,
                position = position_dodge(0.5)) +
  theme(legend.position = "none")

ggsave(filename = '3enrich_NusAG/termination_figures/extData_fig7e_rrf_bargraph_zeroed.png',
       plot = rrf_fig,
       width = 300,
       height = 300,
       units = 'mm',
       dpi = 300)

## Pred. l2fc of 10bp-res terms ------------------------------------------------

term_DF <- data.frame(l2fc = c(rep(0.142894,3), rep(1.023766,3), rep(0.821631,3),
                               rep(1.130831,3), rep(0.639374,3), rep(1.190423,3),
                               rep(1.613103,3), rep(1.304096,3), rep(2.225975,3)),
                      term = c(rep("pitB",9),
                               rep("icd1",9),
                               rep("rrf",9)),
                      lfcSE = c(rep(0.244870,3), rep(0.240234,3), rep(0.240348,3),
                                rep(0.248250,3), rep(0.252742,3), rep(0.248484,3),
                                rep(0.263339,3), rep(0.264172,3), rep(0.262870,3)),
                      TF = c(rep(c(rep("NusA",3),rep("NusG",3),rep("NusAG",3)),3)))

pitB_DF <- term_DF[term_DF$term == 'pitB',]
icd1_DF <- term_DF[term_DF$term == 'icd1',]
rrf_DF <- term_DF[term_DF$term == 'rrf',]

pitB_DF <- pitB_DF %>%
  mutate(TF = factor(TF,
                     levels = c("NusA","NusG","NusAG")))

pitB_fig <- ggplot(pitB_DF, aes(x = TF,
                                y = l2fc,
                                fill = TF)) +
  geom_bar(stat = 'identity', position = "dodge",
           width = 0.5, linewidth = 1) +
  theme_minimal() +
  scale_fill_manual(values = c(NusA_col,
                               NusG_col, NusA_NusG_col)) +
  geom_errorbar(aes(ymin = l2fc - lfcSE,
                    ymax = l2fc + lfcSE),
                width = 0.2, linewidth = 2,
                position = position_dodge(0.5),
                color = 'black') +
  theme(legend.position = "none")

ggsave(filename = '3enrich_NusAG/termination_figures/extData_fig7e_pitB_bargraph_pred.png',
       plot = pitB_fig,
       width = 300,
       height = 300,
       units = 'mm',
       dpi = 300)

icd1_DF <- icd1_DF %>%
  mutate(TF = factor(TF,
                     levels = c("NusA","NusG","NusAG")))

icd1_fig <- ggplot(icd1_DF, aes(x = TF,
                                y = l2fc,
                                fill = TF)) +
  geom_bar(stat = "identity", position = "dodge",
           width = 0.5, linewidth = 1) +
  theme_minimal() +
  scale_fill_manual(values = c(NusA_col,
                               NusG_col, NusA_NusG_col)) +
  geom_errorbar(aes(ymin = l2fc - lfcSE,
                    ymax = l2fc + lfcSE),
                width = 0.2, linewidth = 2,
                position = position_dodge(0.5),
                color = 'black') +
  theme(legend.position = "none")

ggsave(filename = '3enrich_NusAG/termination_figures/extData_fig7e_icd1_bargraph_pred.png',
       plot = icd1_fig,
       width = 300,
       height = 300,
       units = 'mm',
       dpi = 300)

rrf_DF <- rrf_DF %>%
  mutate(TF = factor(TF,
                     levels = c("NusA","NusG","NusAG")))

rrf_fig <- ggplot(rrf_DF, aes(x = TF,
                              y = l2fc,
                              fill = TF)) +
  geom_bar(stat = "identity", position = "dodge",
           width = 0.5, linewidth = 1) +
  theme_minimal() +
  scale_fill_manual(values = c(NusA_col,
                               NusG_col, NusA_NusG_col)) +
  geom_errorbar(aes(ymin = l2fc - lfcSE,
                    ymax = l2fc + lfcSE),
                width = 0.2, linewidth = 2,
                position = position_dodge(0.5),
                color = 'black') +
  theme(legend.position = "none")

ggsave(filename = '3enrich_NusAG/termination_figures/extData_fig7e_rrf_bargraph_pred.png',
       plot = rrf_fig,
       width = 300,
       height = 300,
       units = 'mm',
       dpi = 300)

## Terminator t-tests ----------------------------------------------------------

term_stats$log2_CFG <- log2(term_stats$CFG_counts)
pitB_stats <- term_stats[term_stats$terminator == 'pitB',]
icd1_stats <- term_stats[term_stats$terminator == 'icd1',]
rrf_stats <- term_stats[term_stats$terminator == 'rrf',]

## pitB: rel to noTF

pitB_CFG_NusA_noTF <- t.test(pitB_stats[pitB_stats$TF == 'noTF','log2_CFG'],
                             pitB_stats[pitB_stats$TF == 'NusA','log2_CFG'])

pitB_CFG_NusG_noTF <- t.test(pitB_stats[pitB_stats$TF == 'noTF','log2_CFG'],
                             pitB_stats[pitB_stats$TF == 'NusG','log2_CFG'], alternative = 't')

pitB_CFG_NusAG_noTF <- t.test(pitB_stats[pitB_stats$TF == 'noTF','log2_CFG'],
                              pitB_stats[pitB_stats$TF == 'NusAG','log2_CFG'], alternative = 't')

pitB_CFG_NusA_NusG <- t.test(pitB_stats[pitB_stats$TF == 'NusA','log2_CFG'],
                             pitB_stats[pitB_stats$TF == 'NusG','log2_CFG'])

pitB_CFG_NusAG_NusG <- t.test(pitB_stats[pitB_stats$TF == 'NusG','log2_CFG'],
                              pitB_stats[pitB_stats$TF == 'NusAG','log2_CFG'])

pitB_CFG_NusAG_NusA <- t.test(pitB_stats[pitB_stats$TF == 'NusA','log2_CFG'],
                              pitB_stats[pitB_stats$TF == 'NusAG','log2_CFG'], alternative = 'l')

pitB_TE_NusA_noTF <- t.test(pitB_stats[pitB_stats$TF == 'noTF','TE_FC'],
                            pitB_stats[pitB_stats$TF == 'NusA','TE_FC'])

pitB_TE_NusG_noTF <- t.test(pitB_stats[pitB_stats$TF == 'noTF','TE_FC'],
                            pitB_stats[pitB_stats$TF == 'NusG','TE_FC'], alternative = 't')

pitB_TE_NusAG_noTF <- t.test(pitB_stats[pitB_stats$TF == 'noTF','TE_FC'],
                             pitB_stats[pitB_stats$TF == 'NusAG','TE_FC'], alternative = 't')

pitB_TE_NusA_NusG <- t.test(pitB_stats[pitB_stats$TF == 'NusA','TE_FC'],
                            pitB_stats[pitB_stats$TF == 'NusG','TE_FC'])

pitB_TE_NusAG_NusG <- t.test(pitB_stats[pitB_stats$TF == 'NusG','TE_FC'],
                             pitB_stats[pitB_stats$TF == 'NusAG','TE_FC'])

pitB_TE_NusAG_NusA <- t.test(pitB_stats[pitB_stats$TF == 'NusA','TE_FC'],
                             pitB_stats[pitB_stats$TF == 'NusAG','TE_FC'], alternative = 't')

## icd1: rel to noTF

icd1_CFG_NusA_noTF <- t.test(icd1_stats[icd1_stats$TF == 'noTF','log2_CFG'],
                             icd1_stats[icd1_stats$TF == 'NusA','log2_CFG'], alternative = 't')

icd1_CFG_NusG_noTF <- t.test(icd1_stats[icd1_stats$TF == 'noTF','log2_CFG'],
                             icd1_stats[icd1_stats$TF == 'NusG','log2_CFG'], alternative = 't')

icd1_CFG_NusAG_noTF <- t.test(icd1_stats[icd1_stats$TF == 'noTF','log2_CFG'],
                              icd1_stats[icd1_stats$TF == 'NusAG','log2_CFG'], alternative = 't')

icd1_CFG_NusA_NusG <- t.test(icd1_stats[icd1_stats$TF == 'NusA','log2_CFG'],
                             icd1_stats[icd1_stats$TF == 'NusG','log2_CFG'], alternative = 't')

icd1_CFG_NusAG_NusG <- t.test(icd1_stats[icd1_stats$TF == 'NusG','log2_CFG'],
                              icd1_stats[icd1_stats$TF == 'NusAG','log2_CFG'], alternative = 't')

icd1_CFG_NusAG_NusA <- t.test(icd1_stats[icd1_stats$TF == 'NusA','log2_CFG'],
                              icd1_stats[icd1_stats$TF == 'NusAG','log2_CFG'])

icd1_TE_NusA_noTF <- t.test(icd1_stats[icd1_stats$TF == 'noTF','TE_FC'],
                            icd1_stats[icd1_stats$TF == 'NusA','TE_FC'], alternative = 't')

icd1_TE_NusG_noTF <- t.test(icd1_stats[icd1_stats$TF == 'noTF','TE_FC'],
                            icd1_stats[icd1_stats$TF == 'NusG','TE_FC'], alternative = 't')

icd1_TE_NusAG_noTF <- t.test(icd1_stats[icd1_stats$TF == 'noTF','TE_FC'],
                             icd1_stats[icd1_stats$TF == 'NusAG','TE_FC'], alternative = 't')

icd1_TE_NusA_NusG <- t.test(icd1_stats[icd1_stats$TF == 'NusA','TE_FC'],
                            icd1_stats[icd1_stats$TF == 'NusG','TE_FC'], alternative = 't')

icd1_TE_NusAG_NusG <- t.test(icd1_stats[icd1_stats$TF == 'NusG','TE_FC'],
                             icd1_stats[icd1_stats$TF == 'NusAG','TE_FC'], alternative = 't')

icd1_TE_NusAG_NusA <- t.test(icd1_stats[icd1_stats$TF == 'NusA','TE_FC'],
                             icd1_stats[icd1_stats$TF == 'NusAG','TE_FC'])

## rrf: rel to noTF
rrf_CFG_NusA_noTF <- t.test(rrf_stats[rrf_stats$TF == 'noTF','log2_CFG'],
                            rrf_stats[rrf_stats$TF == 'NusA','log2_CFG'], alternative = 't')

rrf_CFG_NusG_noTF <- t.test(rrf_stats[rrf_stats$TF == 'noTF','log2_CFG'],
                            rrf_stats[rrf_stats$TF == 'NusG','log2_CFG'], alternative = 't')

rrf_CFG_NusAG_noTF <- t.test(rrf_stats[rrf_stats$TF == 'noTF','log2_CFG'],
                             rrf_stats[rrf_stats$TF == 'NusAG','log2_CFG'], alternative = 't')

rrf_CFG_NusA_NusG <- t.test(rrf_stats[rrf_stats$TF == 'NusA','log2_CFG'],
                            rrf_stats[rrf_stats$TF == 'NusG','log2_CFG'])

rrf_CFG_NusAG_NusG <- t.test(rrf_stats[rrf_stats$TF == 'NusG','log2_CFG'],
                             rrf_stats[rrf_stats$TF == 'NusAG','log2_CFG'], alternative = 't')

rrf_CFG_NusAG_NusA <- t.test(rrf_stats[rrf_stats$TF == 'NusA','log2_CFG'],
                             rrf_stats[rrf_stats$TF == 'NusAG','log2_CFG'], alternative = 't')

rrf_TE_NusA_noTF <- t.test(rrf_stats[rrf_stats$TF == 'noTF','TE_FC'],
                           rrf_stats[rrf_stats$TF == 'NusA','TE_FC'], alternative = 't')

rrf_TE_NusG_noTF <- t.test(rrf_stats[rrf_stats$TF == 'noTF','TE_FC'],
                           rrf_stats[rrf_stats$TF == 'NusG','TE_FC'], alternative = 't')

rrf_TE_NusAG_noTF <- t.test(rrf_stats[rrf_stats$TF == 'noTF','TE_FC'],
                            rrf_stats[rrf_stats$TF == 'NusAG','TE_FC'], alternative = 't')

rrf_TE_NusA_NusG <- t.test(rrf_stats[rrf_stats$TF == 'NusA','TE_FC'],
                           rrf_stats[rrf_stats$TF == 'NusG','TE_FC'])

rrf_TE_NusAG_NusG <- t.test(rrf_stats[rrf_stats$TF == 'NusG','TE_FC'],
                            rrf_stats[rrf_stats$TF == 'NusAG','TE_FC'], alternative = 't')

rrf_TE_NusAG_NusA <- t.test(rrf_stats[rrf_stats$TF == 'NusA','TE_FC'],
                            rrf_stats[rrf_stats$TF == 'NusAG','TE_FC'], alternative = 't')

ttest_DF <- data.frame(term = c(rep("pitB",12),
                                rep("icd1",12),
                                rep("rrf",12)),
                       comparison = c(rep(c("NusA_noTF","NusG_noTF","NusAG_noTF",
                                            "NusG_NusA","NusAG_NusA","NusAG_NusG"),6)),
                       assay = c(rep(c(rep("CFG",6),
                                       rep("inVitro",6)),3)),
                       pvalues = c(pitB_CFG_NusA_noTF$p.value,
                                   pitB_CFG_NusG_noTF$p.value,
                                   pitB_CFG_NusAG_noTF$p.value,
                                   pitB_CFG_NusA_NusG$p.value,
                                   pitB_CFG_NusAG_NusA$p.value,
                                   pitB_CFG_NusAG_NusG$p.value,
                                   pitB_TE_NusA_noTF$p.value,
                                   pitB_TE_NusG_noTF$p.value,
                                   pitB_TE_NusAG_noTF$p.value,
                                   pitB_TE_NusA_NusG$p.value,
                                   pitB_TE_NusAG_NusA$p.value,
                                   pitB_TE_NusAG_NusG$p.value,
                                   icd1_CFG_NusA_noTF$p.value,
                                   icd1_CFG_NusG_noTF$p.value,
                                   icd1_CFG_NusAG_noTF$p.value,
                                   icd1_CFG_NusA_NusG$p.value,
                                   icd1_CFG_NusAG_NusA$p.value,
                                   icd1_CFG_NusAG_NusG$p.value,
                                   icd1_TE_NusA_noTF$p.value,
                                   icd1_TE_NusG_noTF$p.value,
                                   icd1_TE_NusAG_noTF$p.value,
                                   icd1_TE_NusA_NusG$p.value,
                                   icd1_TE_NusAG_NusA$p.value,
                                   icd1_TE_NusAG_NusG$p.value,
                                   rrf_CFG_NusA_noTF$p.value,
                                   rrf_CFG_NusG_noTF$p.value,
                                   rrf_CFG_NusAG_noTF$p.value,
                                   rrf_CFG_NusA_NusG$p.value,
                                   rrf_CFG_NusAG_NusA$p.value,
                                   rrf_CFG_NusAG_NusG$p.value,
                                   rrf_TE_NusA_noTF$p.value,
                                   rrf_TE_NusG_noTF$p.value,
                                   rrf_TE_NusAG_noTF$p.value,
                                   rrf_TE_NusA_NusG$p.value,
                                   rrf_TE_NusAG_NusA$p.value,
                                   rrf_TE_NusAG_NusG$p.value))

ttest_DF$bonf_adjusted_pval <- p.adjust(ttest_DF$pvalues, method = 'bonferroni')
ttest_DF$BH_adjusted_pval <- p.adjust(ttest_DF$pvalues, method = 'BH')
ttest_DF$fdr_adjusted_pval <- p.adjust(ttest_DF$pvalues, method = 'fdr')
ttest_DF$holm_adjusted_pval <- p.adjust(ttest_DF$pvalues, method = 'holm')
ttest_DF$hochberg_adjusted_pval <- p.adjust(ttest_DF$pvalues, method = 'hochberg')
ttest_DF$hommel_adjusted_pval <- p.adjust(ttest_DF$pvalues, method = 'hommel')
ttest_DF$BY_adjusted_pval <- p.adjust(ttest_DF$pvalues, method = 'BY')

write.csv(ttest_DF, '3enrich_NusAG/termination_figures/extData_fig7e_ttest_termDF_2sided.csv')

## End -------------------------------------------------------------------------

