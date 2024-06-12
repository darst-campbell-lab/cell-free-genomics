# Generating supplementary table combining TSS overlap & annotations
# Authors: Ruby Froom
# Date: May 28, 2024

## Import and format data ------------------------------------------------------

TSS_overlaps <- read.csv('all_TSSs_needMotifs.csv')
TSS_overlaps[is.na(TSS_overlaps$CFG_TSS_ID_higher),'CFG_TSS_ID_higher'] <- ''
TSS_overlaps[is.na(TSS_overlaps$CFG_coord),'CFG_coord'] <- ''
TSS_overlaps$CFG_TSS_ID <- paste(TSS_overlaps$CFG_TSS_ID_higher,
                                 TSS_overlaps$CFG_direction,
                                 sep = '')

CFG_annotations <- read.csv('CFG_annotated_promoters.csv')
colnames(CFG_annotations) <- paste("CFG_",
                                   colnames(CFG_annotations),
                                   sep='')
CFG_annotations[CFG_annotations == 0] <- ''

inCell_annotations <- read.csv('inCell_Cortes_annotated_promoters.csv')
colnames(inCell_annotations) <- paste("inCell_",
                                   colnames(inCell_annotations),
                                   sep='')
inCell_annotations[inCell_annotations == 0] <- ''


## Merge data ------------------------------------------------------------------

add_CFG <- left_join(TSS_overlaps,
                     CFG_annotations,
                     by = c("CFG_TSS_ID" = "CFG_coordID"))

add_inCell <- left_join(add_CFG,
                        inCell_annotations,
                        by = c("inCell_TSS_ID" = "inCell_coordID"))

add_inCell[is.na(add_inCell)] <- ''

## Clean data for supp table ---------------------------------------------------

colsToKeep <- c("geneID","geneName",
                "CFG_TSS_local_region","inCell_TSS_local_region","TSS_distance",
                "CFG_coord","CFG_direction","CFG_TSS_orientation","CFG_TSS_location",
                "CFG_TSS_ID","CFG_UPregion","CFG_Min35","CFG_Spacer","CFG_SpacerLength",
                "CFG_ExtMin10","CFG_CombinedMin10_Min11A","CFG_Disc","CFG_DiscLength",
                "CFG_TSS_nt","CFG_NT_strand","inCell_coord","inCell_direction",
                "inCell_TSS_orientation","inCell_TSS_location","inCell_TSS_ID",
                "inCell_UPregion","inCell_Min35","inCell_Spacer","inCell_SpacerLength",
                "inCell_ExtMin10","inCell_CombinedMin10_Min11A","inCell_Disc","inCell_DiscLength",
                "inCell_TSS_nt","inCell_NT_strand")

cleaned_DF <- add_inCell[,colsToKeep]

colnames_for_TSS_combined <- c("GeneID","Gene name",
                               "CFG TSS region (-9 to +10)","In cellulo TSS region (-9 to +10)",
                               "CFG vs. in cellulo TSS distance","CFG TSS coordinate","CFG TSS strand",
                               "CFG TSS orientation","CFG TSS location","CFG TSS ID",
                               "CFG UP element region (-59 to -37)",	"CFG -35",
                               "CFG spacer sequence",	"CFG spacer length",
                               "CFG ext. -10",	"CFG -10",	"CFG disc. sequence",
                               "CFG disc. length",	"CFG TSS (+1)", "CFG non-template strand (-100 to +20)",
                               "In cellulo TSS coordinate","In cellulo TSS strand",
                               "In cellulo TSS orientation","In cellulo TSS location","In cellulo TSS ID",
                               "In cellulo UP element region (-59 to -37)",	"In cellulo -35",
                               "In cellulo spacer sequence",	"In cellulo spacer length",
                               "In cellulo ext. -10",	"In cellulo -10",	"In cellulo disc. sequence",
                               "In cellulo disc. length",	"In cellulo TSS (+1)",
                               "In cellulo non-template strand (-100 to +20)")
colnames(cleaned_DF) <- colnames_for_TSS_combined

# Export 
write.csv(cleaned_DF, 'table_S2_all_TSSs_CFG_inCellulo.csv',
          row.names = FALSE)
