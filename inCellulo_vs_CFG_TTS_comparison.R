# Integrating cell-free and in cellulo profiling at single-NT resolution for TTSs
# Authors: Ruby Froom, Michael DeJesus (GFF parsing)
# Date: May 27, 2024

#### Libraries -----------------------------------------------------------------

library(tidyverse)
library(stringr)
library(stringi)
library(Biostrings)
library(seqinr)
library(sqldf)
library(ggseqlogo)

#### Read in TTS files ---------------------------------------------------------

CFG_TTS_all <- read.csv('3enrich_NusAG/selectThreshold/DESeq2/multifactor_countdata_ends_70counts.csv')
# For CFG data, split coordinate and strand
no_Eco_TTS <- CFG_TTS_all[!grepl("Eco", CFG_TTS_all$coordID),]
no_Eco_TTS$coordinate <- as.integer(stri_sub(no_Eco_TTS$coordID,
                                              from = 1, to = -2))

# For CFG data, swap direction
no_Eco_TTS$direction2 <- sapply(strsplit(no_Eco_TTS$coordID,""), tail, 1)
no_Eco_TTS$direction <- no_Eco_TTS$direction2
no_Eco_TTS[no_Eco_TTS$direction2 == '+','direction'] <- '-'
no_Eco_TTS[no_Eco_TTS$direction2 == '-','direction'] <- '+'
CFG_TTS <- no_Eco_TTS[,c("coordinate","direction")]

CFG_TTS <- CFG_TTS[order(CFG_TTS$coordinate),]

inCell_TTS <- read.csv('genome_files_misc/TTS_arnvig_S3.csv')

#### Read in reference genomes -------------------------------------------------

CFG_genome <- read.fasta('genome_files_misc/Eco_Mtb_genome.fasta', seqtype="DNA")[[1]]
inCell_TTS_genome <- read.fasta('genome_files_misc/TTS_arnvig_refGenome.fasta', seqtype="DNA")[[1]]

#### Extract flanking sequences -------------------------------------------------

upstream_end <- 50
downstream_end <- 20
TTS <- upstream_end + 1

# A function that takes in a set of genomic coordinates (e.g. transcription start sites) and the direction of transcription to extract the surrounding genome sequence
# genome_seq = reference genome, pre-loaded and processed so it is indexable (see below function)
# coordinate_DF = coordinate dataframe: 
# needs a "coordinate" column with numbers corresponding to positions in the reference genome
# and a "direction" column indicating which direction transcription is going from that coordinate. If +, then transcription is going left to right
# upstream_end = how many bases upstream of the coordinate to extract
# downstream_end = how many bases downstream of the coordinate to extract
extract_coordinate_environment = function(genome_seq, coordinate_DF, upstream_end, downstream_end){
  
  # Indicates the names of the columns for the new output dataframe
  # Includes all of the column names in the original file, as well as new columns to add 
  column_names <- c(colnames(coordinate_DF),'Non.template.strand',
                    'Template.strand',"IndexingCoordinate","coord_local_region")
  
  # Initializes the output dataframe with the correct dimensions
  new_DF <- data.frame(matrix(ncol=length(column_names), nrow=nrow(coordinate_DF)))
  
  # Fills in the new dataframe with the coordinates and direction of transcription
  new_DF[,1:ncol(coordinate_DF)] <- coordinate_DF
  
  # Names the output dataframe columns using the vector column_names
  colnames(new_DF) <- column_names
  
  new_DF$IndexingCoordinate <- new_DF$coordinate
  new_DF[(new_DF$coordinate < upstream_end) & (new_DF$direction == "+"),
         "IndexingCoordinate"] <- new_DF[(new_DF$coordinate < upstream_end) & (new_DF$direction == "+"),
                                         "coordinate"] + length(genome_seq)
  new_DF[(new_DF$coordinate < downstream_end) & (new_DF$direction == "-"),
         "IndexingCoordinate"] <- new_DF[(new_DF$coordinate < downstream_end) & (new_DF$direction == "-"),
                                         "coordinate"] + length(genome_seq)
  
  double_genome <- paste0(c(genome_seq, genome_seq),collapse='')
  
  sense_DF <- new_DF[new_DF$direction == "+",]
  antisense_DF <- new_DF[new_DF$direction == "-",]
  
  # Initializes the range of the genome to extract
  sense_DF$start <- sense_DF$IndexingCoordinate - upstream_end
  sense_DF$end <- sense_DF$IndexingCoordinate + downstream_end - 1
  sense_DF$IndexCoordMin1 <- sense_DF$IndexingCoordinate - 1
  sense_DF$IndexCoordPlus1 <- sense_DF$IndexingCoordinate + 1
  
  # Breaks up the extracted region to capitalize the TSS
  # Promoter is the area directly upstream of the coordinate
  sense_promoters <- str_sub(double_genome, sense_DF$start, sense_DF$IndexCoordMin1)
  sense_TSS <- toupper(str_sub(double_genome, sense_DF$IndexingCoordinate, sense_DF$IndexingCoordinate))
  sense_downstream <- str_sub(double_genome, sense_DF$IndexCoordPlus1, sense_DF$end)
  sense_DF$Non.template.strand <- paste(sense_promoters, sense_TSS, sense_downstream,sep='')
  
  sense_TSS_local_start <- sense_DF$IndexingCoordinate - 10
  sense_TSS_local_end <- sense_DF$IndexingCoordinate + 10
  sense_DF$coord_local_region <- tolower(paste(str_sub(double_genome, sense_TSS_local_start, sense_DF$IndexCoordMin1),
                                               sense_TSS,
                                               str_sub(double_genome, sense_DF$IndexCoordPlus1, sense_TSS_local_end),sep=''))
  
  ## Generate the template strand
  # Turn the promoter DNA into a DNAString object
  sense_promoter_DNA <- DNAStringSet(sense_promoters)
  # Generates the template strand sequence by taking the complement of the NT strand
  sense_promoter_T <- complement(sense_promoter_DNA)
  
  # Repeat for the TSS and the downstream DNA section
  sense_TSS_DNA <- DNAStringSet(sense_TSS)
  sense_TSS_T <- complement(sense_TSS_DNA)
  sense_downstream_DNA <- DNAStringSet(sense_downstream)
  sense_downstream_T <- complement(sense_downstream_DNA)
  
  # Paste the template strand back together, with uppercase TSS and the rest lowercase
  sense_extracted_regions_T <- paste(tolower(as.character(sense_promoter_T)),toupper(as.character(sense_TSS_T)),tolower(as.character(sense_downstream_T)),sep='')
  sense_DF$Template.strand <- sense_extracted_regions_T
  
  # Initially orients and extracts based on the + strand
  antisense_DF$start <- antisense_DF$IndexingCoordinate - downstream_end
  antisense_DF$end <- antisense_DF$IndexingCoordinate + upstream_end
  antisense_DF$IndexCoordMin1 <- antisense_DF$IndexingCoordinate - 1
  antisense_DF$IndexCoordPlus1 <- antisense_DF$IndexingCoordinate + 1
  
  # For antisense promoters: turns the extracted region into a DNA sequence and takes the reverse complement
  # (converts the + strand to the - strand) to get the non-template strand
  antisense_extract <- str_sub(double_genome, antisense_DF$start, antisense_DF$end)
  antisense_DNAseq <- DNAStringSet(antisense_extract)
  antisense_DNAseq_NT <- as.character(reverseComplement(antisense_DNAseq))
  
  # Breaks up the reverse-complemented promoter to capitalize the TSS letter
  antisense_promoter <- tolower(str_sub(antisense_DNAseq_NT, 1, upstream_end))
  antisense_TSS <- toupper(str_sub(antisense_DNAseq_NT, upstream_end+1, upstream_end+1))
  antisense_downstream <- tolower(str_sub(antisense_DNAseq_NT, upstream_end + 2, upstream_end + downstream_end))
  antisense_DF$Non.template.strand <- paste(antisense_promoter,antisense_TSS,antisense_downstream,sep='')
  
  antisense_extract <- str_sub(double_genome, antisense_DF$IndexingCoordinate - 10, antisense_DF$IndexingCoordinate + 10)
  antisense_DNAseq <- DNAStringSet(antisense_extract)
  antisense_DF$coord_local_region <- tolower(as.character(reverseComplement(antisense_DNAseq)))
  
  # Gets the template strand sequence (the + strand in this case)
  antisense_promoter_DNA <- DNAStringSet(antisense_promoter)
  antisense_promoter_T <- as.character(complement(antisense_promoter_DNA))
  antisense_TSS_DNA <- DNAStringSet(antisense_TSS)
  antisense_TSS_T <- as.character(complement(antisense_TSS_DNA))
  antisense_downstream_DNA <- DNAStringSet(antisense_downstream)
  antisense_downstream_T <- as.character(complement(antisense_downstream_DNA))
  
  # Paste the template strand back together, with uppercase TSS and the rest lowercase
  antisense_extracted_regions_T <- paste(tolower(antisense_promoter_T),toupper(antisense_TSS_T),tolower(antisense_downstream_T),sep='')
  antisense_DF$Template.strand <- antisense_extracted_regions_T
  
  complete_DF <- rbind(sense_DF, antisense_DF)
  complete_DF_ordered <- complete_DF[order(complete_DF$coordinate),]
  final_DF <- complete_DF_ordered[,1:(ncol(coordinate_DF)+4)]
  
  return (final_DF)
}

CFG_TTS_sequences <- data.frame(extract_coordinate_environment(CFG_genome,
                                                               CFG_TTS,
                                                               upstream_end,
                                                               downstream_end))
CFG_TTS_sequences$coordinate <- CFG_TTS_sequences$coordinate - 4641652
CFG_TTS_sequences$coordID <- paste(CFG_TTS_sequences$coordinate,
                                   CFG_TTS_sequences$direction,
                                   sep='')

inCell_TTS_sequences <- data.frame(extract_coordinate_environment(inCell_TTS_genome,
                                                                  inCell_TTS,
                                                                  upstream_end,
                                                                  downstream_end))
inCell_TTS_sequences$coordID <- paste(inCell_TTS_sequences$coordinate,
                                      inCell_TTS_sequences$direction,
                                      sep='')

## Import gene annotations -----------------------------------------------------

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
gene_DF$geneName <- get_list_of_features(gene_DF$metadata, col = 2)
gene_DF$geneDesc <- get_list_of_features(gene_DF$metadata, col = 4)

# keep start, end, strand, gene_id, gene_name, and description columns
keep_cols <- c("start","end","geneStrand","geneID","geneName","geneDesc")
gene_DF_curated <- gene_DF[,keep_cols]
gene_DF_curated$start <- gene_DF_curated$start
gene_DF_curated$end <- gene_DF_curated$end

# Generate columns needed for the "inverse" dataframe
gene_DF_curated$startMin1 <- lead(gene_DF_curated$start, n=1)
gene_DF_curated$namesMin1 <- lead(gene_DF_curated$geneID, n=1)

# Generate the "inverse" dataframe- all of the intergenic spaces
inverseDF <- data.frame(geneID = paste(gene_DF_curated$geneID, "-",
                                       gene_DF_curated$namesMin1," intergenic region",
                                       sep=''),
                        start = gene_DF_curated$end,
                        end = gene_DF_curated$startMin1,
                        geneStrand = "inter",
                        geneName = '',
                        geneDesc = '')
fullDF <- na.omit(rbind(gene_DF_curated[,c(1:6)],inverseDF))
fullDF$start <- as.integer(fullDF$start)
fullDF$end <- as.integer(fullDF$end)
# Only include rows where the end is greater than the start 
#(because of overlapping genes, sometimes some weird rows get generated)
geneDF1 <- fullDF[(fullDF$end > (fullDF$start - 1)),]
geneDF <- geneDF1[(geneDF1$end != geneDF1$start),]
ordered_geneDF <- geneDF[order(geneDF$start),]

### Add gene annotations (TTS-specific) ----------------------------------------

# Adds in annotations for upstream and downstream genes for terminators
# Very important that the gene DF is ordered by genomic start location
# Removes duplicates generated during the annotation process (overlapping ORFs)
add_gene_annotations_TTS = function(coordDF, geneDF) {
  
  # For intergenic coordinates, add one ahead
  geneDF$geneID_next <- lead(geneDF$geneID, n = 1)
  geneDF[nrow(geneDF),'geneID_next'] <- geneDF[1,'geneID']
  
  geneDF$geneID_previous <- lag(geneDF$geneID, n = 1)
  geneDF[1,'geneID_previous'] <- geneDF[nrow(geneDF),'geneID']
  
  geneDF$geneName_next <- lead(geneDF$geneName, n = 1)
  geneDF[nrow(geneDF),'geneName_next'] <- geneDF[1,'geneName']
  
  geneDF$geneName_previous <- lag(geneDF$geneName, n = 1)
  geneDF[1,'geneName_previous'] <- geneDF[nrow(geneDF),'geneName']
  
  geneDF$geneDesc_next <- lead(geneDF$geneDesc, n = 1)
  geneDF[nrow(geneDF),'geneDesc_next'] <- geneDF[1,'geneDesc']
  
  geneDF$geneDesc_previous <- lag(geneDF$geneDesc, n = 1)
  geneDF[1,'geneDesc_previous'] <- geneDF[nrow(geneDF),'geneDesc']
  
  geneDF$geneStrand_next <- lead(geneDF$geneStrand, n = 1)
  geneDF[nrow(geneDF),'geneStrand_next'] <- geneDF[1,'geneStrand']
  
  geneDF$geneStrand_previous <- lag(geneDF$geneStrand, n = 1)
  geneDF[1,'geneStrand_previous'] <- geneDF[nrow(geneDF),'geneStrand']
  
  geneDF$start_next <- lead(geneDF$start, n = 1)
  geneDF[nrow(geneDF),'start_next'] <- geneDF[1,'start']
  
  geneDF$start_previous <- lag(geneDF$start, n = 1)
  geneDF[1,'start_previous'] <- geneDF[nrow(geneDF),'start']
  
  geneDF$end_next <- lead(geneDF$end, n = 1)
  geneDF[nrow(geneDF),'end_next'] <- geneDF[1,'end']
  
  geneDF$end_previous <- lag(geneDF$end, n = 1)
  geneDF[1,'end_previous'] <- geneDF[nrow(geneDF),'end']
  
  # For intragenic coordinates, need 2 ahead
  geneDF$geneID_next2 <- lead(geneDF$geneID, n = 2)
  geneDF[nrow(geneDF),'geneID_next2'] <- geneDF[2,'geneID']
  geneDF[nrow(geneDF)-1,'geneID_next2'] <- geneDF[1,'geneID']
  
  geneDF$geneID_previous2 <- lag(geneDF$geneID, n = 2)
  geneDF[1,'geneID_previous2'] <- geneDF[nrow(geneDF)-1,'geneID']
  geneDF[2,'geneID_previous2'] <- geneDF[nrow(geneDF),'geneID']
  
  geneDF$geneName_next2 <- lead(geneDF$geneName, n = 2)
  geneDF[nrow(geneDF),'geneName_next2'] <- geneDF[2,'geneName']
  geneDF[nrow(geneDF)-1,'geneName_next2'] <- geneDF[1,'geneName']
  
  geneDF$geneName_previous2 <- lag(geneDF$geneName, n = 2)
  geneDF[1,'geneName_previous2'] <- geneDF[nrow(geneDF)-1,'geneName']
  geneDF[2,'geneName_previous2'] <- geneDF[nrow(geneDF),'geneName']
  
  geneDF$geneDesc_next2 <- lead(geneDF$geneDesc, n = 2)
  geneDF[nrow(geneDF),'geneDesc_next2'] <- geneDF[2,'geneDesc']
  geneDF[nrow(geneDF)-1,'geneDesc_next2'] <- geneDF[1,'geneDesc']
  
  geneDF$geneDesc_previous2 <- lag(geneDF$geneDesc, n = 2)
  geneDF[1,'geneDesc_previous2'] <- geneDF[nrow(geneDF)-1,'geneDesc']
  geneDF[2,'geneDesc_previous2'] <- geneDF[nrow(geneDF),'geneDesc']
  
  geneDF$geneStrand_next2 <- lead(geneDF$geneStrand, n = 2)
  geneDF[nrow(geneDF),'geneStrand_next2'] <- geneDF[2,'geneStrand']
  geneDF[nrow(geneDF)-1,'geneStrand_next2'] <- geneDF[1,'geneStrand']
  
  geneDF$geneStrand_previous2 <- lag(geneDF$geneStrand, n = 2)
  geneDF[1,'geneStrand_previous2'] <- geneDF[nrow(geneDF)-1,'geneStrand']
  geneDF[2,'geneStrand_previous2'] <- geneDF[nrow(geneDF),'geneStrand']
  
  geneDF$start_next2 <- lead(geneDF$start, n = 2)
  geneDF[nrow(geneDF),'start_next2'] <- geneDF[2,'start']
  geneDF[nrow(geneDF)-1,'start_next2'] <- geneDF[1,'start']
  
  geneDF$start_previous2 <- lag(geneDF$start, n = 2)
  geneDF[1,'start_previous2'] <- geneDF[nrow(geneDF)-1,'start']
  geneDF[2,'start_previous2'] <- geneDF[nrow(geneDF),'start']
  
  geneDF$end_next2 <- lead(geneDF$end, n = 2)
  geneDF[nrow(geneDF),'end_next2'] <- geneDF[2,'end']
  geneDF[nrow(geneDF)-1,'end_next2'] <- geneDF[1,'end']
  
  geneDF$end_previous2 <- lag(geneDF$end, n = 2)
  geneDF[1,'end_previous2'] <- geneDF[nrow(geneDF)-1,'end']
  geneDF[2,'end_previous2'] <- geneDF[nrow(geneDF),'end']
  
  addGenes <- sqldf("select *
      from coordDF inner join geneDF
        on (coordDF.coordinate between geneDF.start and geneDF.end)")
  
  addGenes$upstream_geneID <- ''
  addGenes$upstream_geneName <- ''
  addGenes$upstream_geneDesc <- ''
  addGenes$upstream_geneDistance <- ''
  
  addGenes$downstream_geneID <- ''
  addGenes$downstream_geneName <- ''
  addGenes$downstream_geneDesc <- ''
  addGenes$downstream_geneDistance <- ''
  
  # Annotate up- and downstream genes in a strand-sensitive manner
  # Positive strand intergenic terminators:
  # If coordinate and upstream gene are in same orientation:
  addGenes[(addGenes$geneStrand == 'inter') & (addGenes$direction == '+') & (addGenes$geneStrand_previous == '+'),
           'upstream_geneID'] <- addGenes[(addGenes$geneStrand == 'inter') & (addGenes$direction == '+') & (addGenes$geneStrand_previous == '+'),
                                          'geneID_previous']
  addGenes[(addGenes$geneStrand == 'inter') & (addGenes$direction == '+') & (addGenes$geneStrand_previous == '+'),
           'upstream_geneDistance'] <- addGenes[(addGenes$geneStrand == 'inter') & (addGenes$direction == '+') & (addGenes$geneStrand_previous == '+'),
                                                'coordinate'] - addGenes[(addGenes$geneStrand == 'inter') & (addGenes$direction == '+') & (addGenes$geneStrand_previous == '+'),
                                                                         'end_previous']
  # If coordinate and upstream gene are in opposite orientations:
  addGenes[(addGenes$geneStrand == 'inter') & (addGenes$direction == '+') & (addGenes$geneStrand_previous == '-'),
           'upstream_geneID'] <- paste(addGenes[(addGenes$geneStrand == 'inter') & (addGenes$direction == '+') & (addGenes$geneStrand_previous == '-'),
                                                'geneID_previous'], "_as", sep='')
  addGenes[(addGenes$geneStrand == 'inter') & (addGenes$direction == '+') & (addGenes$geneStrand_previous == '-'),
           'upstream_geneDistance'] <- addGenes[(addGenes$geneStrand == 'inter') & (addGenes$direction == '+') & (addGenes$geneStrand_previous == '-'),
                                                'coordinate'] - addGenes[(addGenes$geneStrand == 'inter') & (addGenes$direction == '+') & (addGenes$geneStrand_previous == '-'),
                                                                         'end_previous']
  
  addGenes[(addGenes$geneStrand == 'inter') & (addGenes$direction == '+') & (addGenes$geneStrand_previous == '+'),
           'upstream_geneName'] <- addGenes[(addGenes$geneStrand == 'inter') & (addGenes$direction == '+') & (addGenes$geneStrand_previous == '+'),
                                            'geneName_previous']
  addGenes[(addGenes$geneStrand == 'inter') & (addGenes$direction == '+') & (addGenes$geneStrand_previous == '-'),
           'upstream_geneName'] <- paste(addGenes[(addGenes$geneStrand == 'inter') & (addGenes$direction == '+') & (addGenes$geneStrand_previous == '-'),
                                                  'geneName_previous'], " (antisense)", sep='')
  
  addGenes[(addGenes$geneStrand == 'inter') & (addGenes$direction == '+') & (addGenes$geneStrand_previous == '+'),
           'upstream_geneDesc'] <- addGenes[(addGenes$geneStrand == 'inter') & (addGenes$direction == '+') & (addGenes$geneStrand_previous == '+'),
                                            'geneDesc_previous']
  addGenes[(addGenes$geneStrand == 'inter') & (addGenes$direction == '+') & (addGenes$geneStrand_previous == '-'),
           'upstream_geneDesc'] <- paste(addGenes[(addGenes$geneStrand == 'inter') & (addGenes$direction == '+') & (addGenes$geneStrand_previous == '-'),
                                                  'geneDesc_previous'], " (antisense)", sep='')
  
  addGenes[(addGenes$geneStrand == 'inter') & (addGenes$direction == '+') & (addGenes$geneStrand_next == '+'),
           'downstream_geneID'] <- addGenes[(addGenes$geneStrand == 'inter') & (addGenes$direction == '+') & (addGenes$geneStrand_next == '+'),
                                            'geneID_next']
  addGenes[(addGenes$geneStrand == 'inter') & (addGenes$direction == '+') & (addGenes$geneStrand_next == '+'),
           'downstream_geneDistance'] <- addGenes[(addGenes$geneStrand == 'inter') & (addGenes$direction == '+') & (addGenes$geneStrand_next == '+'),
                                                  'start_next'] - addGenes[(addGenes$geneStrand == 'inter') & (addGenes$direction == '+') & (addGenes$geneStrand_next == '+'),
                                                                           'coordinate']
  addGenes[(addGenes$geneStrand == 'inter') & (addGenes$direction == '+') & (addGenes$geneStrand_next == '-'),
           'downstream_geneID'] <- paste(addGenes[(addGenes$geneStrand == 'inter') & (addGenes$direction == '+') & (addGenes$geneStrand_next == '-'),
                                                  'geneID_next'], "_as", sep='')
  addGenes[(addGenes$geneStrand == 'inter') & (addGenes$direction == '+') & (addGenes$geneStrand_next == '-'),
           'downstream_geneDistance'] <- addGenes[(addGenes$geneStrand == 'inter') & (addGenes$direction == '+') & (addGenes$geneStrand_next == '-'),
                                                  'start_next'] - addGenes[(addGenes$geneStrand == 'inter') & (addGenes$direction == '+') & (addGenes$geneStrand_next == '-'),
                                                                           'coordinate']
  
  addGenes[(addGenes$geneStrand == 'inter') & (addGenes$direction == '+') & (addGenes$geneStrand_next == '+'),
           'downstream_geneName'] <- addGenes[(addGenes$geneStrand == 'inter') & (addGenes$direction == '+') & (addGenes$geneStrand_next == '+'),
                                              'geneName_next']
  addGenes[(addGenes$geneStrand == 'inter') & (addGenes$direction == '+') & (addGenes$geneStrand_next == '-'),
           'downstream_geneName'] <- paste(addGenes[(addGenes$geneStrand == 'inter') & (addGenes$direction == '+') & (addGenes$geneStrand_next == '-'),
                                                    'geneName_next'], " (antisense)", sep='')
  
  addGenes[(addGenes$geneStrand == 'inter') & (addGenes$direction == '+') & (addGenes$geneStrand_next == '+'),
           'downstream_geneDesc'] <- addGenes[(addGenes$geneStrand == 'inter') & (addGenes$direction == '+') & (addGenes$geneStrand_next == '+'),
                                              'geneDesc_next']
  addGenes[(addGenes$geneStrand == 'inter') & (addGenes$direction == '+') & (addGenes$geneStrand_next == '-'),
           'downstream_geneDesc'] <- paste(addGenes[(addGenes$geneStrand == 'inter') & (addGenes$direction == '+') & (addGenes$geneStrand_next == '-'),
                                                    'geneDesc_next'], " (antisense)", sep='')
  
  # Negative strand intergenic terminators:
  addGenes[(addGenes$geneStrand == 'inter') & (addGenes$direction == '-') & (addGenes$geneStrand_next == '-'),
           'upstream_geneID'] <- addGenes[(addGenes$geneStrand == 'inter') & (addGenes$direction == '-') & (addGenes$geneStrand_next == '-'),
                                          'geneID_next']
  addGenes[(addGenes$geneStrand == 'inter') & (addGenes$direction == '-') & (addGenes$geneStrand_next == '-'),
           'upstream_geneDistance'] <- addGenes[(addGenes$geneStrand == 'inter') & (addGenes$direction == '-') & (addGenes$geneStrand_next == '-'),
                                                'start_next'] - addGenes[(addGenes$geneStrand == 'inter') & (addGenes$direction == '-') & (addGenes$geneStrand_next == '-'),
                                                                         'coordinate']
  addGenes[(addGenes$geneStrand == 'inter') & (addGenes$direction == '-') & (addGenes$geneStrand_next == '+'),
           'upstream_geneID'] <- paste(addGenes[(addGenes$geneStrand == 'inter') & (addGenes$direction == '-') & (addGenes$geneStrand_next == '+'),
                                                'geneID_next'], "_as", sep='')
  addGenes[(addGenes$geneStrand == 'inter') & (addGenes$direction == '-') & (addGenes$geneStrand_next == '+'),
           'upstream_geneDistance'] <- addGenes[(addGenes$geneStrand == 'inter') & (addGenes$direction == '-') & (addGenes$geneStrand_next == '+'),
                                                'start_next'] - addGenes[(addGenes$geneStrand == 'inter') & (addGenes$direction == '-') & (addGenes$geneStrand_next == '+'),
                                                                         'coordinate']
  
  addGenes[(addGenes$geneStrand == 'inter') & (addGenes$direction == '-') & (addGenes$geneStrand_next == '-'),
           'upstream_geneName'] <- addGenes[(addGenes$geneStrand == 'inter') & (addGenes$direction == '-') & (addGenes$geneStrand_next == '-'),
                                            'geneName_next']
  addGenes[(addGenes$geneStrand == 'inter') & (addGenes$direction == '-') & (addGenes$geneStrand_next == '+'),
           'upstream_geneName'] <- paste(addGenes[(addGenes$geneStrand == 'inter') & (addGenes$direction == '-') & (addGenes$geneStrand_next == '+'),
                                                  'geneName_next'], " (antisense)", sep='')
  
  addGenes[(addGenes$geneStrand == 'inter') & (addGenes$direction == '-') & (addGenes$geneStrand_next == '-'),
           'upstream_geneDesc'] <- addGenes[(addGenes$geneStrand == 'inter') & (addGenes$direction == '-') & (addGenes$geneStrand_next == '-'),
                                            'geneDesc_next']
  addGenes[(addGenes$geneStrand == 'inter') & (addGenes$direction == '-') & (addGenes$geneStrand_next == '+'),
           'upstream_geneDesc'] <- paste(addGenes[(addGenes$geneStrand == 'inter') & (addGenes$direction == '-') & (addGenes$geneStrand_next == '+'),
                                                  'geneDesc_next'], " (antisense)", sep='')
  
  addGenes[(addGenes$geneStrand == 'inter') & (addGenes$direction == '-') & (addGenes$geneStrand_previous == '-'),
           'downstream_geneID'] <- addGenes[(addGenes$geneStrand == 'inter') & (addGenes$direction == '-') & (addGenes$geneStrand_previous == '-'),
                                            'geneID_previous']
  addGenes[(addGenes$geneStrand == 'inter') & (addGenes$direction == '-') & (addGenes$geneStrand_previous == '-'),
           'downstream_geneDistance'] <- addGenes[(addGenes$geneStrand == 'inter') & (addGenes$direction == '-') & (addGenes$geneStrand_previous == '-'),
                                                  'coordinate'] - addGenes[(addGenes$geneStrand == 'inter') & (addGenes$direction == '-') & (addGenes$geneStrand_previous == '-'),
                                                                           'end_previous']
  addGenes[(addGenes$geneStrand == 'inter') & (addGenes$direction == '-') & (addGenes$geneStrand_previous == '+'),
           'downstream_geneID'] <- paste(addGenes[(addGenes$geneStrand == 'inter') & (addGenes$direction == '-') & (addGenes$geneStrand_previous == '+'),
                                                  'geneID_previous'], "_as", sep='')
  addGenes[(addGenes$geneStrand == 'inter') & (addGenes$direction == '-') & (addGenes$geneStrand_previous == '+'),
           'downstream_geneDistance'] <- addGenes[(addGenes$geneStrand == 'inter') & (addGenes$direction == '-') & (addGenes$geneStrand_previous == '+'),
                                                  'coordinate'] - addGenes[(addGenes$geneStrand == 'inter') & (addGenes$direction == '-') & (addGenes$geneStrand_previous == '+'),
                                                                           'end_previous']
  addGenes[(addGenes$geneStrand == 'inter') & (addGenes$direction == '-') & (addGenes$geneStrand_previous == '-'),
           'downstream_geneName'] <- addGenes[(addGenes$geneStrand == 'inter') & (addGenes$direction == '-') & (addGenes$geneStrand_previous == '-'),
                                              'geneName_previous']
  addGenes[(addGenes$geneStrand == 'inter') & (addGenes$direction == '-') & (addGenes$geneStrand_previous == '+'),
           'downstream_geneName'] <- paste(addGenes[(addGenes$geneStrand == 'inter') & (addGenes$direction == '-') & (addGenes$geneStrand_previous == '+'),
                                                    'geneName_previous'], " (antisense)", sep='')
  
  addGenes[(addGenes$geneStrand == 'inter') & (addGenes$direction == '-') & (addGenes$geneStrand_previous == '-'),
           'downstream_geneDesc'] <- addGenes[(addGenes$geneStrand == 'inter') & (addGenes$direction == '-') & (addGenes$geneStrand_previous == '-'),
                                              'geneDesc_previous']
  addGenes[(addGenes$geneStrand == 'inter') & (addGenes$direction == '-') & (addGenes$geneStrand_previous == '+'),
           'downstream_geneDesc'] <- paste(addGenes[(addGenes$geneStrand == 'inter') & (addGenes$direction == '-') & (addGenes$geneStrand_previous == '+'),
                                                    'geneDesc_previous'], " (antisense)", sep='')
  
  # Positive strand intragenic terminators:
  addGenes[(addGenes$geneStrand == '+') & (addGenes$direction == '+'),
           'upstream_geneID'] <- addGenes[(addGenes$geneStrand == '+') & (addGenes$direction == '+'),
                                          'geneID']
  addGenes[(addGenes$geneStrand == '+') & (addGenes$direction == '+'),
           'upstream_geneDistance'] <- 0
  addGenes[(addGenes$geneStrand == '-') & (addGenes$direction == '+'),
           'upstream_geneID'] <- paste(addGenes[(addGenes$geneStrand == '-') & (addGenes$direction == '+'),
                                                'geneID'], "_as", sep='')
  addGenes[(addGenes$geneStrand == '-') & (addGenes$direction == '+'),
           'upstream_geneDistance'] <- 0
  
  addGenes[(addGenes$geneStrand == '+') & (addGenes$direction == '+'),
           'upstream_geneName'] <- addGenes[(addGenes$geneStrand == '+') & (addGenes$direction == '+'),
                                            'geneName']
  addGenes[(addGenes$geneStrand == '-') & (addGenes$direction == '+'),
           'upstream_geneName'] <- paste(addGenes[(addGenes$geneStrand == '-') & (addGenes$direction == '+'),
                                                  'geneName'], " (antisense)", sep='')
  
  addGenes[(addGenes$geneStrand == '+') & (addGenes$direction == '+'),
           'upstream_geneDesc'] <- addGenes[(addGenes$geneStrand == '+') & (addGenes$direction == '+'),
                                            'geneDesc']
  addGenes[(addGenes$geneStrand == '-') & (addGenes$direction == '+'),
           'upstream_geneDesc'] <- paste(addGenes[(addGenes$geneStrand == '-') & (addGenes$direction == '+'),
                                                  'geneDesc'], " (antisense)", sep='')
  
  addGenes[(addGenes$geneStrand != 'inter') & (addGenes$direction == '+') & (addGenes$geneStrand_next2 == '+'),
           'downstream_geneID'] <- addGenes[(addGenes$geneStrand != 'inter') & (addGenes$direction == '+') & (addGenes$geneStrand_next2 == '+'),
                                            'geneID_next2']
  addGenes[(addGenes$geneStrand != 'inter') & (addGenes$direction == '+') & (addGenes$geneStrand_next2 == '+'),
           'downstream_geneDistance'] <- addGenes[(addGenes$geneStrand != 'inter') & (addGenes$direction == '+') & (addGenes$geneStrand_next2 == '+'),
                                                  'start_next2'] - addGenes[(addGenes$geneStrand != 'inter') & (addGenes$direction == '+') & (addGenes$geneStrand_next2 == '+'),
                                                                            'coordinate']
  addGenes[(addGenes$geneStrand != 'inter') & (addGenes$direction == '+') & (addGenes$geneStrand_next2 == '-'),
           'downstream_geneID'] <- paste(addGenes[(addGenes$geneStrand != 'inter') & (addGenes$direction == '+') & (addGenes$geneStrand_next2 == '-'),
                                                  'geneID_next2'], " (antisense)", sep='')
  addGenes[(addGenes$geneStrand != 'inter') & (addGenes$direction == '+') & (addGenes$geneStrand_next2 == '-'),
           'downstream_geneDistance'] <- addGenes[(addGenes$geneStrand != 'inter') & (addGenes$direction == '+') & (addGenes$geneStrand_next2 == '-'),
                                                  'start_next2'] - addGenes[(addGenes$geneStrand != 'inter') & (addGenes$direction == '+') & (addGenes$geneStrand_next2 == '-'),
                                                                            'coordinate']
  
  addGenes[(addGenes$geneStrand != 'inter') & (addGenes$direction == '+') & (addGenes$geneStrand_next2 == '+'),
           'downstream_geneName'] <- addGenes[(addGenes$geneStrand != 'inter') & (addGenes$direction == '+') & (addGenes$geneStrand_next2 == '+'),
                                              'geneName_next2']
  addGenes[(addGenes$geneStrand != 'inter') & (addGenes$direction == '+') & (addGenes$geneStrand_next2 == '-'),
           'downstream_geneName'] <- paste(addGenes[(addGenes$geneStrand != 'inter') & (addGenes$direction == '+') & (addGenes$geneStrand_next2 == '-'),
                                                    'geneName_next2'], " (antisense)", sep='')
  
  addGenes[(addGenes$geneStrand != 'inter') & (addGenes$direction == '+') & (addGenes$geneStrand_next2 == '+'),
           'downstream_geneDesc'] <- addGenes[(addGenes$geneStrand != 'inter') & (addGenes$direction == '+') & (addGenes$geneStrand_next2 == '+'),
                                              'geneDesc_next2']
  addGenes[(addGenes$geneStrand != 'inter') & (addGenes$direction == '+') & (addGenes$geneStrand_next2 == '-'),
           'downstream_geneDesc'] <- paste(addGenes[(addGenes$geneStrand != 'inter') & (addGenes$direction == '+') & (addGenes$geneStrand_next2 == '-'),
                                                    'geneDesc_next2'], " (antisense)", sep='')
  
  # Negative strand intragenic terminators:
  addGenes[(addGenes$geneStrand == '-') & (addGenes$direction == '-'),
           'upstream_geneID'] <- addGenes[(addGenes$geneStrand == '-') & (addGenes$direction == '-'),
                                          'geneID']
  addGenes[(addGenes$geneStrand == '-') & (addGenes$direction == '-'),
           'upstream_geneDistance'] <- 0
  addGenes[(addGenes$geneStrand == '+') & (addGenes$direction == '-'),
           'upstream_geneID'] <- paste(addGenes[(addGenes$geneStrand == '+') & (addGenes$direction == '-'),
                                                'geneID'], "_as", sep='')
  addGenes[(addGenes$geneStrand == '+') & (addGenes$direction == '-'),
           'upstream_geneDistance'] <- 0
  
  addGenes[(addGenes$geneStrand == '-') & (addGenes$direction == '-'),
           'upstream_geneName'] <- addGenes[(addGenes$geneStrand == '-') & (addGenes$direction == '-'),
                                            'geneName']
  addGenes[(addGenes$geneStrand == '+') & (addGenes$direction == '-'),
           'upstream_geneName'] <- paste(addGenes[(addGenes$geneStrand == '+') & (addGenes$direction == '-'),
                                                  'geneName'], " (antisense)", sep='')
  
  addGenes[(addGenes$geneStrand == '-') & (addGenes$direction == '-'),
           'upstream_geneDesc'] <- addGenes[(addGenes$geneStrand == '-') & (addGenes$direction == '-'),
                                            'geneDesc']
  addGenes[(addGenes$geneStrand == '+') & (addGenes$direction == '-'),
           'upstream_geneDesc'] <- paste(addGenes[(addGenes$geneStrand == '+') & (addGenes$direction == '-'),
                                                  'geneDesc'], " (antisense)", sep='')
  
  addGenes[(addGenes$geneStrand != 'inter') & (addGenes$direction == '-') & (addGenes$geneStrand_previous2 == '-'),
           'downstream_geneID'] <- addGenes[(addGenes$geneStrand != 'inter') & (addGenes$direction == '-') & (addGenes$geneStrand_previous2 == '-'),
                                            'geneID_previous2']
  addGenes[(addGenes$geneStrand != 'inter') & (addGenes$direction == '-') & (addGenes$geneStrand_next2 == '-'),
           'downstream_geneDistance'] <- addGenes[(addGenes$geneStrand != 'inter') & (addGenes$direction == '-') & (addGenes$geneStrand_next2 == '-'),
                                                  'coordinate'] - addGenes[(addGenes$geneStrand != 'inter') & (addGenes$direction == '-') & (addGenes$geneStrand_next2 == '-'),
                                                                           'end_previous2']
  addGenes[(addGenes$geneStrand != 'inter') & (addGenes$direction == '-') & (addGenes$geneStrand_previous2 == '+'),
           'downstream_geneID'] <- paste(addGenes[(addGenes$geneStrand != 'inter') & (addGenes$direction == '-') & (addGenes$geneStrand_previous2 == '+'),
                                                  'geneID_previous2'], "_as", sep='')
  addGenes[(addGenes$geneStrand != 'inter') & (addGenes$direction == '-') & (addGenes$geneStrand_next2 == '+'),
           'downstream_geneDistance'] <- addGenes[(addGenes$geneStrand != 'inter') & (addGenes$direction == '-') & (addGenes$geneStrand_next2 == '+'),
                                                  'coordinate'] - addGenes[(addGenes$geneStrand != 'inter') & (addGenes$direction == '-') & (addGenes$geneStrand_next2 == '+'),
                                                                           'end_previous2']
  
  addGenes[(addGenes$geneStrand != 'inter') & (addGenes$direction == '-') & (addGenes$geneStrand_previous2 == '-'),
           'downstream_geneName'] <- addGenes[(addGenes$geneStrand != 'inter') & (addGenes$direction == '-') & (addGenes$geneStrand_previous2 == '-'),
                                              'geneName_previous2']
  addGenes[(addGenes$geneStrand != 'inter') & (addGenes$direction == '-') & (addGenes$geneStrand_previous2 == '+'),
           'downstream_geneName'] <- paste(addGenes[(addGenes$geneStrand != 'inter') & (addGenes$direction == '-') & (addGenes$geneStrand_previous2 == '+'),
                                                    'geneName_previous2'], " (antisense)", sep='')
  
  addGenes[(addGenes$geneStrand != 'inter') & (addGenes$direction == '-') & (addGenes$geneStrand_previous2 == '-'),
           'downstream_geneDesc'] <- addGenes[(addGenes$geneStrand != 'inter') & (addGenes$direction == '-') & (addGenes$geneStrand_previous2 == '-'),
                                              'geneDesc_previous2']
  addGenes[(addGenes$geneStrand != 'inter') & (addGenes$direction == '-') & (addGenes$geneStrand_previous2 == '+'),
           'downstream_geneDesc'] <- paste(addGenes[(addGenes$geneStrand != 'inter') & (addGenes$direction == '-') & (addGenes$geneStrand_previous2 == '+'),
                                                    'geneDesc_previous2'], " (antisense)", sep='')
  
  # Specify whether a coordinate is inter- or intragenic
  addGenes$coord_location <- 'intragenic'
  addGenes[addGenes$geneStrand == 'inter',
           'coord_location'] <- 'intergenic'
  
  withDups <- unique(addGenes)
  duplicated_coords <- unique(withDups[duplicated(withDups$coordID),'coordID'])
  withDups$Remove <- 'No'

  for (i in 1:length(duplicated_coords)) {
    
    coord_DF <- withDups[withDups$coordID == duplicated_coords[i],]
    
    # If both intergenic and intragenic annotations are possible,
    # choose intergenic
    coord_locations <- unique(coord_DF$coord_location)
    if (length(coord_locations) > 1) {
      withDups[((withDups$coordID == duplicated_coords[i]) &
                  (withDups$coord_location == 'intragenic')),
               'Remove'] <- 'Yes'
    }
    
    # Otherwise, choose row with both upstream and downstream geneIDs annotated
    else {
      withDups[(withDups$coordID == duplicated_coords[i]),
             'Remove'] <- 'Yes'
      withDups[((withDups$coordID == duplicated_coords[i]) &
                (withDups$upstream_geneID != '') &
                (withDups$downstream_geneID != '')),
             'Remove'] <- 'No'
    }
  }
  
  noDups <- withDups[withDups$Remove == 'No',]
  
  colsToKeep <- c("coordinate","direction","Non.template.strand",
                  "coord_local_region","coordID","coord_location",
                  "upstream_geneID","upstream_geneName",
                  "upstream_geneDesc",
                  "downstream_geneID","downstream_geneName",
                  "downstream_geneDesc")
  clean_DF <- noDups[,colsToKeep]
  
  print('Number of unique coordIDs:')
  print(length(unique(noDups$coordID)))
  print('')
  print('Number of DF rows:')
  print(nrow(noDups))
  return(noDups)
}

CFG_TTS_DF <- add_gene_annotations_TTS(CFG_TTS_sequences, ordered_geneDF)
inCell_TTS_DF <- add_gene_annotations_TTS(inCell_TTS_sequences, ordered_geneDF)

## Find CFG vs. in cellulo TTS overlap -----------------------------------------

find_coord_overlap <- function(str1, str2, ignore.case = TRUE, verbose = FALSE) {
  
  if(ignore.case) {
    str1 <- tolower(str1)
    str2 <- tolower(str2)
  }
  
  if(nchar(str1) < nchar(str2)) {
    x <- str2
    str2 <- str1
    str1 <- x
  }
  
  x <- strsplit(str2, "")[[1L]]
  n <- length(x)
  s <- sequence(seq_len(n))
  s <- split(s, cumsum(s == 1L))
  s <- rep(list(s), n)
  
  for(i in seq_along(s)) {
    s[[i]] <- lapply(s[[i]], function(x) {
      x <- x + (i-1L)
      x[x <= n]
    })
    s[[i]] <- unique(s[[i]])
  }
  
  s <- unlist(s, recursive = FALSE)
  s <- unique(s[order(-lengths(s))])
  
  i <- 1L
  len_s <- length(s)
  while(i < len_s) {
    lcs <- paste(x[s[[i]]], collapse = "")
    if(verbose) cat("now checking:", lcs, "\n")
    check <- grepl(lcs, str1, fixed = TRUE)
    if(check) {
      #      cat("the (first) longest common substring is:", lcs, "of length", nchar(lcs), "\n")
      return (lcs)
      break
    } else {
      i <- i + 1L 
    }
  }
}


# Identical TTSs
same_TTS <- CFG_TTS_DF[(CFG_TTS_DF$coord_local_region %in% inCell_TTS_DF$coord_local_region),]

# Extract regions with TTSs only in CFG
CFG_TTS_genes <- unique(CFG_TTS_DF$upstream_geneID)
inCell_TTS_genes <- unique(inCell_TTS_DF$upstream_geneID)
both_TTS_genes <- CFG_TTS_genes[(CFG_TTS_genes %in% inCell_TTS_genes)]
CFG_only_TTS_genes <- CFG_TTS_genes[!(CFG_TTS_genes %in% inCell_TTS_genes)]
inCell_only_TTS_genes <- inCell_TTS_genes[!(inCell_TTS_genes %in% CFG_TTS_genes)]

CFG_only_TTSs <- CFG_TTS_DF[CFG_TTS_DF$upstream_geneID %in% CFG_only_TTS_genes,]
inCell_only_TTSs <- inCell_TTS_DF[inCell_TTS_DF$upstream_geneID %in% inCell_only_TTS_genes,]

CFG_TTS_in_both <- CFG_TTS_DF[(CFG_TTS_DF$upstream_geneID %in% both_TTS_genes),]
inCell_TTS_in_both <- inCell_TTS_DF[(inCell_TTS_DF$upstream_geneID %in% both_TTS_genes),]

CFG_TTS_in_both$inCell_TTS <- 'None'
CFG_TTS_in_both$CFG_coord_region <- 'None'
CFG_TTS_in_both$inCell_distance <- 'No overlap'

inCell_TTS_in_both$CFG_TTS <- 'None'
inCell_TTS_in_both$CFG_coord_region <- 'None'
inCell_TTS_in_both$CFG_distance <- 'No overlap'

max_length_coord_region <- length(strsplit(inCell_TTS_in_both$coord_local_region[1],
                                           split='')[[1]])

for (i in 1:length(both_TTS_genes)) {
  
  CFG_gene_DF <- CFG_TTS_in_both[CFG_TTS_in_both$upstream_geneID == both_TTS_genes[i],]
  inCell_gene_DF <- inCell_TTS_in_both[inCell_TTS_in_both$upstream_geneID == both_TTS_genes[i],]
  
  for (j in 1:nrow(CFG_gene_DF)) {
    
    CFG_coord_region <- CFG_gene_DF$coord_local_region[j]
    CFG_coord <- CFG_gene_DF$coordID[j]
    
    for (k in 1:nrow(inCell_gene_DF)) {
      
      inCell_coord_region <- inCell_gene_DF$coord_local_region[k]
      
      largest_substring <- find_coord_overlap(CFG_coord_region, inCell_coord_region)
      substring_length <- length(strsplit(largest_substring, split='')[[1]])
      distance <- max_length_coord_region - substring_length
      
      # If the coordinates are in the same direction
      # and distance is calculated as being short:
      if ((inCell_gene_DF$direction[k] == CFG_gene_DF$direction[j]) & (distance <= 5)) {
        
        # Case where CFG TTS is slightly upstream of inCell TTS:
        end_of_inCell <- paste0(strsplit(inCell_coord_region, split='')[[1]][(max_length_coord_region - (substring_length-1)):max_length_coord_region],
                                collapse='')
        beginning_of_CFG <- paste0(strsplit(CFG_coord_region, split='')[[1]][1:(1+(substring_length-1))],
                                   collapse='')
        
        if ((end_of_inCell == largest_substring) & (beginning_of_CFG == largest_substring)) {
          
          inCell_TTS_in_both[inCell_TTS_in_both$coord_local_region == inCell_coord_region,
                             'CFG_TTS'] <- CFG_coord
          inCell_TTS_in_both[inCell_TTS_in_both$coord_local_region == inCell_coord_region,
                             'CFG_coord_region'] <- CFG_coord_region
          inCell_TTS_in_both[inCell_TTS_in_both$coord_local_region == inCell_coord_region,
                             'CFG_distance'] <- distance
          
        }
        
        # Case where CFG TTS is slightly downstream of inCell TTS:
        end_of_CFG <- tolower(paste0(strsplit(CFG_coord_region, split='')[[1]][(max_length_coord_region - (substring_length-1)):max_length_coord_region],
                                     collapse = ''))
        beginning_of_inCell <- tolower(paste0(strsplit(inCell_coord_region, split='')[[1]][1:(1+(substring_length-1))],
                                              collapse=''))
        
        if ((end_of_CFG == largest_substring) & (beginning_of_inCell == largest_substring)) {
          
          inCell_TTS_in_both[inCell_TTS_in_both$coord_local_region == inCell_coord_region,
                             'CFG_TTS'] <- CFG_coord
          inCell_TTS_in_both[inCell_TTS_in_both$coord_local_region == inCell_coord_region,
                             'CFG_coord_region'] <- CFG_coord_region
          inCell_TTS_in_both[inCell_TTS_in_both$coord_local_region == inCell_coord_region,
                             'CFG_distance'] <- distance
          
        }
        
      }
    }
  }
}

## Final numbers for Venn diagrams ---------------------------------------------

TTSs_in_common <- unique(inCell_TTS_in_both[inCell_TTS_in_both$CFG_TTS != 'None','coordID'])
TTSs_CFG_only <- unique(CFG_TTS_DF$coordID[!(CFG_TTS_DF$coordID %in% TTSs_in_common)])
TTSs_inCell_only <- unique(inCell_TTS_DF$coordID[!(inCell_TTS_DF$coordID %in% TTSs_in_common)])

num_common_TTSs <- length(TTSs_in_common)
num_CFG_only_TTSs <- length(TTSs_CFG_only)
num_inCell_only_TTSs <- length(TTSs_inCell_only)
neither_TTS <- 4411709 - (num_common_TTSs + num_CFG_only_TTSs + num_inCell_only_TTSs)

fisher_matrix_TTS <- matrix(c(num_common_TTSs, 
                              num_CFG_only_TTSs,
                              num_inCell_only_TTSs,
                              neither_TTS),
                            nrow = 2, ncol= 2)

fishers_TTS_DF <- fisher.test(fisher_matrix_TTS)

## TTSs per million calculations
# TSSs per million calculations
average_CFG_depths <- mean(c(6.674843,
                             6.673739,
                             6.675003,
                             6.667759,
                             6.668986,
                             6.672933,
                             6.670218,
                             6.676320,
                             6.673104,
                             6.671555,
                             6.671324,
                             6.677252))

average_inCell_nonrRNA_depths <- mean(c(11.407316,
                                        9.723713,
                                        9.797698))


TTS_num_DF <- data.frame(CFG_only_raw = num_CFG_only_TTSs,
                         both_raw = num_common_TTSs,
                         inCell_only_raw = num_inCell_only_TTSs,
                         CFG_only_TTSperMill = num_CFG_only_TTSs / average_CFG_depths,
                         both_TTSperMill = num_common_TTSs / average_CFG_depths,
                         inCell_only_TTSperMill = num_inCell_only_TTSs / average_inCell_nonrRNA_depths)

write.csv(TTS_num_DF, 'fig2_inCellulo_vs_CFG/fig2b_TTS_numForVenn.csv')


## Add CFG overlap info to inCell DFs ------------------------------------------

common_TTS_DF <- inCell_TTS_in_both[inCell_TTS_in_both$CFG_TTS != 'None',
                                    c("coordID","CFG_TTS","CFG_coord_region","CFG_distance")]
inCell_with_CFGoverlap_TTS <- left_join(inCell_TTS_DF,
                                        common_TTS_DF,
                                        by = c("coordID" = "coordID"))
inCell_with_CFGoverlap_TTS[is.na(inCell_with_CFGoverlap_TTS$CFG_TTS),'CFG_TTS'] <- 'No overlap'
inCell_with_CFGoverlap_TTS[is.na(inCell_with_CFGoverlap_TTS$CFG_coord_region),'CFG_coord_region'] <- 'No overlap'
inCell_with_CFGoverlap_TTS[is.na(inCell_with_CFGoverlap_TTS$CFG_distance),'CFG_distance'] <- 'No overlap'

## Histogram -------------------------------------------------------------------

TTS_overlapping_only <- inCell_with_CFGoverlap_TTS[inCell_with_CFGoverlap_TTS$CFG_TTS != 'No overlap',]
TTS_overlapping_only$CFG_distance <- as.double(TTS_overlapping_only$CFG_distance)

TTS_hist <- ggplot(TTS_overlapping_only, aes(x = CFG_distance)) +
  geom_histogram(bins = 6,
                 color = 'black',
                 fill = 'grey',
                 binwidth = 0.5) +
  theme_minimal()
ggsave(filename = 'fig2_inCellulo_vs_CFG/extData_fig4b_CFG_inCell_TTSdistances.png',
       plot = TTS_hist,
       width = 100,
       height = 100,
       units = 'mm',
       dpi = 300)

#### Generate TTS supplementary table ------------------------------------------

# Add in TTS ID 
# (with direction swapped to match bigWigs, rather than actual transcription direction)
CFG_TTS_DF$CFG_coord_plus <- CFG_TTS_DF$coordinate + 4641652

CFG_TTS_DF$direction2 <- CFG_TTS_DF$direction
CFG_TTS_DF[CFG_TTS_DF$direction == '+','direction2'] <- '-'
CFG_TTS_DF[CFG_TTS_DF$direction == '-','direction2'] <- '+'
CFG_TTS_DF$CFG_coordID <- paste(CFG_TTS_DF$CFG_coord_plus,
                                    CFG_TTS_DF$direction2,
                                    sep='')

combine_all_TTSs <- full_join(inCell_with_CFGoverlap_TTS,
                              CFG_TTS_DF,
                              by = c("CFG_TTS" = 'coordID'))

combine_all_TTSs[is.na(combine_all_TTSs)] <- ''
combine_all_TTSs[combine_all_TTSs$CFG_distance == '','CFG_distance'] <- 'No overlap'

## Format TTS regions ----------------------------------------------------------

combine_all_TTSs$CFG_TTS_local_region <- combine_all_TTSs$coord_local_region.y
combine_all_TTSs$inCell_TTS_local_region <- combine_all_TTSs$coord_local_region.x

for (i in 1:nrow(combine_all_TTSs)) {
  
  if (combine_all_TTSs$CFG_TTS_local_region[i] != '') {
    
    CFG_TTS_region <- strsplit(combine_all_TTSs$CFG_TTS_local_region[i],split='')[[1]]
    combine_all_TTSs$CFG_TTS_local_region[i] <- paste0(paste0(tolower(CFG_TTS_region[1:10]),
                                                              collapse=''),
                                                       paste0(toupper(CFG_TTS_region[11]),
                                                              collapse=''),
                                                       paste0(tolower(CFG_TTS_region[12:20]),
                                                              collapse=''),
                                                       collapse='')
  }
  
  if (combine_all_TTSs$inCell_TTS_local_region[i] != '') {
    
    inCell_TTS_region <- strsplit(combine_all_TTSs$inCell_TTS_local_region[i],split='')[[1]]
    
    combine_all_TTSs$inCell_TTS_local_region[i] <- paste0(paste0(tolower(inCell_TTS_region[1:10]),
                                                                 collapse=''),
                                                          paste0(toupper(inCell_TTS_region[11]),
                                                                 collapse=''),
                                                          paste0(tolower(inCell_TTS_region[12:20]),
                                                                 collapse=''),
                                                          collapse='')
  }
  
}

# Re-order so that overlapping TTSs come first
ordered_TTSs <- combine_all_TTSs[order(combine_all_TTSs$CFG_distance),]

# Create a single column with gene information  --------------------------------
# for both CFG and in cellulo TTSs ---------------------------------------------
ordered_TTSs$upstream_geneID <- ordered_TTSs$upstream_geneID.x
ordered_TTSs$upstream_geneName <- ordered_TTSs$upstream_geneName.x

ordered_TTSs$downstream_geneID <- ordered_TTSs$downstream_geneID.x
ordered_TTSs$downstream_geneName <- ordered_TTSs$downstream_geneName.x

ordered_TTSs[ordered_TTSs$upstream_geneID == '',
                 'upstream_geneID'] <- ordered_TTSs[ordered_TTSs$upstream_geneID == '',
                                               'upstream_geneID.y']
ordered_TTSs[ordered_TTSs$upstream_geneName == '',
                 'upstream_geneName'] <- ordered_TTSs[ordered_TTSs$upstream_geneName == '',
                                                        'upstream_geneName.y']

ordered_TTSs[ordered_TTSs$downstream_geneID == '',
                 'downstream_geneID'] <- ordered_TTSs[ordered_TTSs$downstream_geneID == '',
                                                        'downstream_geneID.y']
ordered_TTSs[ordered_TTSs$downstream_geneName == '',
                 'downstream_geneName'] <- ordered_TTSs[ordered_TTSs$downstream_geneName == '',
                                                          'downstream_geneName.y']

ordered_TTSs[is.na(ordered_TTSs)] <- ''


keep_columns <- c("upstream_geneID","upstream_geneName",
                  "downstream_geneID","downstream_geneName",
                  "CFG_TTS_local_region","inCell_TTS_local_region","CFG_distance",
                  "coordinate.y", "direction.y",
                  "coord_location.y","CFG_coordID","Non.template.strand.y",
                  "coordinate.x", "direction.x",
                  "coord_location.x","coordID","Non.template.strand.x")

all_TTSs_curated <- ordered_TTSs[,keep_columns]

all_TTSs_curated[((all_TTSs_curated$upstream_geneName == '-') |
                   (all_TTSs_curated$upstream_geneName == '- (antisense)')),
                'upstream_geneName'] <- ''
all_TTSs_curated[((all_TTSs_curated$downstream_geneName == '-') |
                   (all_TTSs_curated$downstream_geneName == '- (antisense)')),
                'downstream_geneName'] <- ''

colnames_for_TTS_combined <- c("Upstream geneID","Upstream gene name",
                               "Downstream geneID","Downstream gene name",
                               "CFG TTS region (-9 to +10)","In cellulo TTS region (-9 to +10)",
                               "CFG vs. in cellulo TTS distance",
                               "CFG TTS coordinate","CFG TTS strand","CFG TTS location",
                               "CFG TTS ID","CFG non-template strand (-50 to +20)",
                               "In cellulo TTS coordinate","In cellulo TTS strand","In cellulo TTS location",
                               "In cellulo TTS ID","In cellulo non-template strand (-50 to +20)")
colnames(all_TTSs_curated) <- colnames_for_TTS_combined


# Export 
write.csv(all_TTSs_curated, 'fig2_inCellulo_vs_CFG/table_S3_all_TTSs_CFG_inCellulo.csv',
          row.names = FALSE)

### Prep settings for relative entropy logos -----------------------------------

cs = make_col_scheme(chars=c('A', 'T', 'C', 'G'),
                     cols=c('#C70000','#008000',
                            '#0000C7','#FFAE00'))

# Calculate background model with fasta-get-markov
# Load in (code from Mike Wolfe)
prior <- read_delim('genome_files_misc/Mtb_background.tsv',
                    col_names = c("value", "prior"),
                    comment = "#",
                    delim = " ")

prior[prior$value == 'U','value'] <- 'T'

relative_entropy_calc <- function(input_seqs, csv_name) {
  write.csv(toupper(input_seqs),
            paste(csv_name, '.csv',sep=''),
            row.names = FALSE)
  d <- read.csv(paste(csv_name, '.csv',sep=''),
                col.names = "seq")
  
  seq_length <- length(strsplit(d$seq[1],split = '')[[1]])
  
  # Manually calculate the PWM
  bg_cor <- d %>% mutate(seq = str_split(seq, "", n = seq_length, simplify = TRUE)) %>% 
    as.matrix() %>% as_tibble() %>% pivot_longer(everything()) %>% filter(value != "N") %>%
    mutate(name = as.numeric(str_remove(name, "seq."))) %>%
    group_by(name, value) %>% summarize(n=n()) %>% 
    group_by(name) %>% mutate(freq = n/sum(n)) %>% 
    left_join(prior, by = "value") %>% 
    group_by(name) %>% 
    mutate(total_inf = sum(freq*log2(freq/prior))) %>% ungroup() %>% 
    mutate(rel_ent = pmax(freq* log2(freq/prior), 0)) %>%
    mutate(height = rel_ent) %>% select(name, value, height) %>% 
    pivot_wider(names_from = name, values_from = height) %>%
    column_to_rownames("value") %>%
    as.matrix()
  return (bg_cor)
}

### Generate sequence logos for TTS regions ------------------------------------

# CFG TTSs
CFG_TTS_sequences <- combine_all_TTSs$Non.template.strand.y
CFG_TTS_sequences <- CFG_TTS_sequences[CFG_TTS_sequences != '']

TTS_relEntropy <- relative_entropy_calc(unique(toupper(CFG_TTS_sequences)),
                                        paste('fig2_inCellulo_vs_CFG/CFG_TTS',sep=''))

p <- ggseqlogo(TTS_relEntropy, method = "custom", seq_type = 'dna',
               col_scheme = cs)
CFG_TTS_logo <- p + theme_classic() + labs(y = "Relative Entropy") + ylim(0, 0.7)

ggsave(filename = paste('fig2_inCellulo_vs_CFG/relEntropy_CFG_TTSlogo.png',sep=''),
       plot = CFG_TTS_logo,
       width = 300,
       height = 100,
       units = 'mm',
       dpi = 300)

# In cellulo TTSs
inCellulo_TTS_sequences <- combine_all_TTSs$Non.template.strand.x
inCellulo_TTS_sequences <- inCellulo_TTS_sequences[inCellulo_TTS_sequences != '']

TTS_relEntropy <- relative_entropy_calc(unique(toupper(inCellulo_TTS_sequences)),
                                        paste('fig2_inCellulo_vs_CFG/inCellulo_TTS',sep=''))

p <- ggseqlogo(TTS_relEntropy, method = "custom", seq_type = 'dna',
               col_scheme = cs)
inCellulo_TTS_logo <- p + theme_classic() + labs(y = "Relative Entropy") + ylim(0, 0.7)

ggsave(filename = paste('fig2_inCellulo_vs_CFG/relEntropy_inCellulo_TTSlogo.png',sep=''),
       plot = inCellulo_TTS_logo,
       width = 300,
       height = 100,
       units = 'mm',
       dpi = 300)

