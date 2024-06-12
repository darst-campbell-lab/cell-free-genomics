# Integrating cell-free and in cellulo profiling at single-NT resolution for TSSs
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

#### Read in TSS files ---------------------------------------------------------

CFG_TSS_all <- read.csv('5enrich_CRP/selectThreshold/DESeq2/noCRP_CRP_20counts_ends_results_Wald_local.csv')
colnames(CFG_TSS_all) <- c("coordID",colnames(CFG_TSS_all[2:ncol(CFG_TSS_all)]))
# For CFG data, split coordinate and strand
CFG_TSS_all$coordinate <- as.integer(stri_sub(CFG_TSS_all$coordID,
                                              from = 1, to = -2))
CFG_TSS_all$direction <- sapply(strsplit(CFG_TSS_all$coordID,""), tail, 1)
CFG_TSS <- CFG_TSS_all[,c("coordinate","direction")]
CFG_TSS <- CFG_TSS[order(CFG_TSS$coordinate),]

inCell_TSS <- read.csv('genome_files_misc/TSS_cortes_S1F.csv')

#### Read in reference genomes -------------------------------------------------

CFG_genome <- read.fasta('genome_files_misc/Eco_Mtb_genome.fasta', seqtype="DNA")[[1]]
inCell_TSS_genome <- read.fasta('genome_files_misc/TSS_cortes_refGenome.fasta', seqtype="DNA")[[1]]

#### Extract flanking sequences -------------------------------------------------

upstream_end <- 100
downstream_end <- 20
TSS <- upstream_end + 1

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

# Extracting sequences surrounding the TSS coordinate
CFG_TSS_sequences <- data.frame(extract_coordinate_environment(CFG_genome,
                                                               CFG_TSS,
                                                               upstream_end,
                                                               downstream_end))
# Need to get coordinate in correct range for H37Rv gene model in next step
CFG_TSS_sequences$coordinate <- CFG_TSS_sequences$coordinate - 4641652
CFG_TSS_sequences$coordID <- paste(CFG_TSS_sequences$coordinate,
                                   CFG_TSS_sequences$direction,
                                   sep='')

inCell_TSS_sequences <- data.frame(extract_coordinate_environment(inCell_TSS_genome,
                                                                  inCell_TSS,
                                                                  upstream_end,
                                                                  downstream_end))
inCell_TSS_sequences$coordID <- paste(inCell_TSS_sequences$coordinate,
                                      inCell_TSS_sequences$direction,
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

### Add gene annotations (TSS-specific) -------------------------------

add_gene_annotations_TSS = function(coordDF, geneDF) {
  
  # Join the promoter and gene DFs
  addGenes <- sqldf("select *
        from coordDF left join geneDF
          on (coordDF.coordinate between geneDF.start and geneDF.end)")

  addGenes$geneID_curated <- addGenes$geneID
  addGenes$strand_curated <- addGenes$geneStrand

  # Replace intergenic with downstream gene in _curated columns
  for (i in 1:nrow(addGenes)) {
  
    # For rows that correspond to intergenic regions:
    x <- stri_match(regex="intergenic", str=addGenes$geneID[i])
  
    if (!is.na(x)) {
    
      genes <- strsplit(addGenes$geneID[i], split="-")[[1]]
    
      # If the coordinate is going in the + direction, annotate the downstream gene
      if (addGenes$direction[i] == "+") {
      
        genes2 <- strsplit(genes[2], split=" ")[[1]]
        addGenes$geneID_curated[i] <- genes2[1]
        addGenes$strand_curated[i] <- paste0(gene_DF_curated[gene_DF_curated$geneID == genes2[1],
                                                             "geneStrand"],collapse='')
      
        # If the coordinate is going in the - direction, annotate the upstream gene  
      }else {
      
        addGenes$geneID_curated[i] <- genes[1]
        addGenes$strand_curated[i] <- paste0(gene_DF_curated[gene_DF_curated$geneID == genes[1],
                                                             "geneStrand"], collapse='')
      
      }
    }
  }

  addGenes <- addGenes[order(addGenes$coordinate),]

  # Annotate whether the TSS corresponds to a sense or antisense transcript, based on the transcript vs. gene orientation
  addGenes$transcriptOrientation <- "sense"

  addGenes[is.na(addGenes)] <- 0
  addGenes[(addGenes$direction != addGenes$strand_curated),
         "transcriptOrientation"] <- "antisense"
  
  # New column to add to the geneID whether the TSS is antisense or not
  addGenes$geneID_withAS <- addGenes$geneID_curated
  
  for (i in 1:nrow(addGenes)) {
    
    if (addGenes$transcriptOrientation[i] == 'antisense') {
      
      addGenes$geneID_withAS[i] <- paste(addGenes$geneID_curated[i],
                                         "_as",
                                         sep='')
      
    }
  }

  # Add TSS location
  addGenes$coord_location <- "intragenic"
  addGenes[addGenes$strand == 'inter','coord_location'] <- 'intergenic'

  return(unique(addGenes))
}

CFG_TSS_DF <- add_gene_annotations_TSS(CFG_TSS_sequences, ordered_geneDF)
inCell_TSS_DF <- add_gene_annotations_TSS(inCell_TSS_sequences, ordered_geneDF)

## Find CFG vs. in cellulo TSS overlap -----------------------------------------

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

# Identical TSSs
same_TSS <- CFG_TSS_DF[(CFG_TSS_DF$coord_local_region %in% inCell_TSS_DF$coord_local_region),]

# Extract regions with TSSs only in CFG
CFG_TSS_genes <- unique(CFG_TSS_DF$geneID)
inCell_TSS_genes <- unique(inCell_TSS_DF$geneID)
both_TSS_genes <- CFG_TSS_genes[CFG_TSS_genes %in% inCell_TSS_genes]
CFG_only_TSS_genes <- CFG_TSS_genes[!(CFG_TSS_genes %in% inCell_TSS_genes)]
inCell_only_TSS_genes <- inCell_TSS_genes[!(inCell_TSS_genes %in% CFG_TSS_genes)]

CFG_only_TSSs <- CFG_TSS_DF[CFG_TSS_DF$geneID %in% CFG_only_TSS_genes,]
inCell_only_TSSs <- inCell_TSS_DF[inCell_TSS_DF$geneID %in% inCell_only_TSS_genes,]

CFG_TSS_in_both <- CFG_TSS_DF[(CFG_TSS_DF$geneID %in% both_TSS_genes),]
inCell_TSS_in_both <- inCell_TSS_DF[(inCell_TSS_DF$geneID %in% both_TSS_genes),]

CFG_TSS_in_both$inCell_TSS <- 'None'
CFG_TSS_in_both$CFG_coord_region <- 'None'
CFG_TSS_in_both$inCell_distance <- 'No overlap'

inCell_TSS_in_both$CFG_TSS <- 'None'
inCell_TSS_in_both$CFG_coord_region <- 'None'
inCell_TSS_in_both$CFG_distance <- 'No overlap'

max_length_coord_region <- length(strsplit(inCell_TSS_in_both$coord_local_region[1],
                                    split='')[[1]])

for (i in 1:length(both_TSS_genes)) {

  CFG_gene_DF <- CFG_TSS_in_both[CFG_TSS_in_both$geneID == both_TSS_genes[i],]
  inCell_gene_DF <- inCell_TSS_in_both[inCell_TSS_in_both$geneID == both_TSS_genes[i],]
  
  for (j in 1:nrow(CFG_gene_DF)) {

    CFG_coord_region <- CFG_gene_DF$coord_local_region[j]
    CFG_coord <- CFG_gene_DF$coordID[j]
    
    for (k in 1:nrow(inCell_gene_DF)) {

      inCell_coord_region <- inCell_gene_DF$coord_local_region[k]
      
      largest_substring <- find_coord_overlap(CFG_coord_region, inCell_coord_region)
      substring_length <- length(strsplit(largest_substring, split='')[[1]])
      distance <- max_length_coord_region - substring_length
      
      if ((inCell_gene_DF$direction[k] == CFG_gene_DF$direction[j]) & (distance <= 5)) {
        # If the coordinates are in the same direction
        # and distance is calculated as being short:
        if ((inCell_gene_DF$direction[k] == CFG_gene_DF$direction[j]) & (distance <= 5)) {
          
          # Case where CFG TTS is slightly upstream of inCell TTS:
          end_of_inCell <- paste0(strsplit(inCell_coord_region, split='')[[1]][(max_length_coord_region - (substring_length-1)):max_length_coord_region],
                                  collapse='')
          beginning_of_CFG <- paste0(strsplit(CFG_coord_region, split='')[[1]][1:(1+(substring_length-1))],
                                     collapse='')
          
          if ((end_of_inCell == largest_substring) & (beginning_of_CFG == largest_substring)) {
            
            inCell_TSS_in_both[inCell_TSS_in_both$coord_local_region == inCell_coord_region,
                               'CFG_TSS'] <- CFG_coord
            inCell_TSS_in_both[inCell_TSS_in_both$coord_local_region == inCell_coord_region,
                               'CFG_coord_region'] <- CFG_coord_region
            inCell_TSS_in_both[inCell_TSS_in_both$coord_local_region == inCell_coord_region,
                               'CFG_distance'] <- distance
            
          }
          
          # Case where CFG TSS is slightly downstream of inCell TSS:
          end_of_CFG <- tolower(paste0(strsplit(CFG_coord_region, split='')[[1]][(max_length_coord_region - (substring_length-1)):max_length_coord_region],
                                       collapse = ''))
          beginning_of_inCell <- tolower(paste0(strsplit(inCell_coord_region, split='')[[1]][1:(1+(substring_length-1))],
                                                collapse=''))
          
          if ((end_of_CFG == largest_substring) & (beginning_of_inCell == largest_substring)) {
            
            inCell_TSS_in_both[inCell_TSS_in_both$coord_local_region == inCell_coord_region,
                               'CFG_TSS'] <- CFG_coord
            inCell_TSS_in_both[inCell_TSS_in_both$coord_local_region == inCell_coord_region,
                               'CFG_coord_region'] <- CFG_coord_region
            inCell_TSS_in_both[inCell_TSS_in_both$coord_local_region == inCell_coord_region,
                               'CFG_distance'] <- distance
            
          }
          
        }
      }
    }
  }
}

## Final numbers for Venn diagrams ---------------------------------------------

TSSs_in_common <- unique(inCell_TSS_in_both[inCell_TSS_in_both$CFG_TSS != 'None','coordID'])
TSSs_CFG_only <- unique(CFG_TSS_DF$coordID[!(CFG_TSS_DF$coordID %in% TSSs_in_common)])
TSSs_inCell_only <- unique(inCell_TSS_DF$coordID[!(inCell_TSS_DF$coordID %in% TSSs_in_common)])

num_common_TSSs <- length(TSSs_in_common)
num_CFG_only_TSSs <- length(TSSs_CFG_only)
num_inCell_only_TSSs <- length(TSSs_inCell_only)
neither_TSS <- 4411709 - (num_common_TSSs + num_CFG_only_TSSs + num_inCell_only_TSSs)

fisher_matrix_TSS <- matrix(c(num_common_TSSs, 
                          num_CFG_only_TSSs,
                          num_inCell_only_TSSs,
                          neither_TSS),
                        nrow = 2, ncol= 2)

fishers_TSS_DF <- fisher.test(fisher_matrix_TSS)

# TSSs per million calculations
average_CFG_depths <- mean(c(7.181966,
                        7.183429,
                        7.182328,
                        7.178550,
                        7.182842,
                        7.178846))

average_inCell_nonrRNA_depths <- mean(c(4.098506,4.295579,5.864551))

TSS_num_DF <- data.frame(CFG_only_raw = num_CFG_only_TSSs,
                         both_raw = num_common_TSSs,
                         inCell_only_raw = num_inCell_only_TSSs,
                         CFG_only_TSSperMill = num_CFG_only_TSSs / average_CFG_depths,
                         both_TSSperMill = num_common_TSSs / average_inCell_nonrRNA_depths,
                         inCell_only_TSSperMill = num_inCell_only_TSSs / average_inCell_nonrRNA_depths)

write.csv(TSS_num_DF, 'fig2_inCellulo_vs_CFG/fig2b_TSS_numForVenn.csv')
                               
## Function to remove duplicate values from a TSS dataframe --------------------
## Add in leaderless annotations; choose sense orientation over antisense; -----
## choose more downstream ORF if there are two overlapping ---------------------

# Preconditions:
# column names: coord_location, transcriptOrientation, start, coordinate, coordID
remove_coordinate_duplicates_TSS <- function(input_DF) {

  input_DF[((input_DF$coordinate == input_DF$start) &
              (input_DF$transcriptOrientation == 'sense')),
           'coord_location'] <- 'leaderless'

  input_DF$Remove <- 'No'

  unique_coordinates <- unique(input_DF$coordID)

  # Remove duplicated rows that arise due to leaderless TSSs
  for (i in 1:length(unique_coordinates)) {
  
    coordID <- unique_coordinates[i]
    coord_DF <- input_DF[input_DF$coordID == coordID,]
    leaderless <- coord_DF[coord_DF$coord_location == 'leaderless',]
  
    if ((nrow(coord_DF) > 1) & (nrow(leaderless) > 0)) {
    
      input_DF[((input_DF$coordID == coordID) & 
                  input_DF$coord_location != 'leaderless'),'Remove'] <- 'Yes'
    
    }
  
  }

  withLeaderless <- unique(input_DF[input_DF$Remove == 'No',])

  # Remove duplicated rows that arise due to both sense and antisense orientations
  # (choose sense orientation)
  for (i in 1:length(unique_coordinates)) {
  
    coordID <- unique_coordinates[i]
    coord_DF <- withLeaderless[withLeaderless$coordID == coordID,]
    sense <- coord_DF[coord_DF$transcriptOrientation == 'sense',]
    antisense <- coord_DF[coord_DF$transcriptOrientation == 'sense',]
  
    if ((nrow(coord_DF) > 1) & (nrow(sense) > 0) & (nrow(antisense) > 0)) {
    
      withLeaderless[((withLeaderless$coordID == coordID) &
                                          withLeaderless$transcriptOrientation == 'antisense'),
                                       'Remove'] <- 'Yes'
    
    }
  
  }

  # Remove duplicated rows that arise due to overlapping ORFs in same direction
  # (choose more downstream)
  chooseSense <- unique(withLeaderless[withLeaderless$Remove == 'No',])
  remaining_duplicates <- chooseSense[duplicated(chooseSense$coordID), 'coordID']

  for (i in 1:length(remaining_duplicates)) {
  
    coordID <- remaining_duplicates[i]
    coord_DF <- chooseSense[chooseSense$coordID == coordID,]
    direction <- unique(coord_DF$direction)
    transcriptOrientation <- unique(coord_DF$transcriptOrientation)
  
    # For TSSs in the positive direction,
  # choose more downstream gene (gene with higher start value)
    if (direction == '+') {
    
      chooseSense[((chooseSense$coordID == coordID) &
                                     chooseSense$start != max(coord_DF$start)),
                                  'Remove'] <- 'Yes'
    
    }
  
    # For TSSs in the negative direction,
    # choose more downstream gene (gene with lower start value)
    if (direction == '-') {
    
      chooseSense[((chooseSense$coordID == coordID) &
                                     chooseSense$start != min(coord_DF$start)),
                                  'Remove'] <- 'Yes'
    
    }
  
  }

  no_overlapping_ORFs <- unique(chooseSense[chooseSense$Remove == 'No',])
  
  print('Number of unique coordIDs:')
  print(length(unique(no_overlapping_ORFs$coordID)))
  print('')
  print('Number of DF rows:')
  print(nrow(no_overlapping_ORFs))
  return(no_overlapping_ORFs)
}


## Add CFG overlap info to inCell DFs ------------------------------------------

common_TSS_DF <- inCell_TSS_in_both[inCell_TSS_in_both$CFG_TSS != 'None',
                                    c("coordID","CFG_TSS","CFG_coord_region","CFG_distance")]
inCell_with_CFGoverlap_TSS <- left_join(inCell_TSS_DF,
                                    common_TSS_DF,
                                    by = c("coordID" = "coordID"))
inCell_with_CFGoverlap_TSS[is.na(inCell_with_CFGoverlap_TSS$CFG_TSS),'CFG_TSS'] <- 'No overlap'
inCell_with_CFGoverlap_TSS[is.na(inCell_with_CFGoverlap_TSS$CFG_coord_region),'CFG_coord_region'] <- 'No overlap'
inCell_with_CFGoverlap_TSS[is.na(inCell_with_CFGoverlap_TSS$CFG_distance),'CFG_distance'] <- 'No overlap'

## Histogram -------------------------------------------------------------------

inCell_CFGoverlap_TSS_noDup <- remove_coordinate_duplicates_TSS(inCell_with_CFGoverlap_TSS)
TSS_overlapping_only <- inCell_CFGoverlap_TSS_noDup[inCell_CFGoverlap_TSS_noDup$CFG_TSS != 'No overlap',]
TSS_overlapping_only$CFG_distance <- as.double(TSS_overlapping_only$CFG_distance)

TSS_hist <- ggplot(TSS_overlapping_only, aes(x = CFG_distance)) +
  geom_histogram(bins = 6,
                 color = 'black',
                 fill = 'grey',
                 binwidth = 0.5) +
  theme_minimal() +
  scale_y_continuous(limits = c(0, 600),
                     breaks = seq(0, 600, 100))
ggsave(filename = 'fig2_inCellulo_vs_CFG/extData_fig4a_CFG_inCell_TSSdistances.png',
       plot = TSS_hist,
       width = 100,
       height = 100,
       units = 'mm',
       dpi = 300)

#### Generate TSS supplementary table ------------------------------------------

CFG_TSSs_noDup <- remove_coordinate_duplicates_TSS(CFG_TSS_DF)
CFG_TSSs_noDup$CFG_coord_plus <- CFG_TSSs_noDup$coordinate + 4641652
CFG_TSSs_noDup$CFG_coordID <- paste(CFG_TSSs_noDup$CFG_coord_plus,
                              CFG_TSSs_noDup$direction,
                              sep='')

combine_all_TSSs <- full_join(inCell_CFGoverlap_TSS_noDup,
                              CFG_TSSs_noDup,
                              by = c("CFG_TSS" = 'coordID'))

# Create a single column with gene information 
# for both CFG and in cellulo TSSs
combine_all_TSSs$geneID <- combine_all_TSSs$geneID_withAS.x
combine_all_TSSs$geneName <- combine_all_TSSs$geneName.x
combine_all_TSSs$geneDesc <- combine_all_TSSs$geneDesc.x


combine_all_TSSs[is.na(combine_all_TSSs$geneID),
                 'geneID'] <- combine_all_TSSs[is.na(combine_all_TSSs$geneID),
                                               'geneID_withAS.y']
combine_all_TSSs[is.na(combine_all_TSSs$gene_start),
                 'geneName'] <- combine_all_TSSs[is.na(combine_all_TSSs$gene_start),
                                               'geneName.y']
combine_all_TSSs[is.na(combine_all_TSSs$gene_end),
                 'geneDesc'] <- combine_all_TSSs[is.na(combine_all_TSSs$gene_end),
                                                   'geneDesc.y']

combine_all_TSSs[is.na(combine_all_TSSs)] <- ''

# Add in gene names
combine_all_TSSs[is.na(combine_all_TSSs)] <- ''
combine_all_TSSs[combine_all_TSSs$CFG_distance == '','CFG_distance'] <- 'No overlap'

## Format TSS regions ----------------------------------------------------------

combine_all_TSSs$CFG_TSS_local_region <- combine_all_TSSs$coord_local_region.y
combine_all_TSSs$inCell_TSS_local_region <- combine_all_TSSs$coord_local_region.x

for (i in 1:nrow(combine_all_TSSs)) {
  
  if (combine_all_TSSs$CFG_TSS_local_region[i] != '') {

    CFG_TSS_region <- strsplit(combine_all_TSSs$CFG_TSS_local_region[i],split='')[[1]]
    combine_all_TSSs$CFG_TSS_local_region[i] <- paste0(paste0(tolower(CFG_TSS_region[1:10]),
                                                            collapse=''),
                                                     paste0(toupper(CFG_TSS_region[11]),
                                                            collapse=''),
                                                     paste0(tolower(CFG_TSS_region[12:20]),
                                                            collapse=''),
                                                     collapse='')
  }
  
  if (combine_all_TSSs$inCell_TSS_local_region[i] != '') {
    
    inCell_TSS_region <- strsplit(combine_all_TSSs$inCell_TSS_local_region[i],split='')[[1]]
    
    combine_all_TSSs$inCell_TSS_local_region[i] <- paste0(paste0(tolower(inCell_TSS_region[1:10]),
                                                               collapse=''),
                                                        paste0(toupper(inCell_TSS_region[11]),
                                                               collapse=''),
                                                        paste0(tolower(inCell_TSS_region[12:20]),
                                                               collapse=''),
                                                        collapse='')
  }
}

### Final formatting for export ------------------------------------------------

# Re-order so overlapping ones come first
ordered_TSSs <- combine_all_TSSs[order(combine_all_TSSs$CFG_distance),]

keep_columns <- c("geneID","geneName","geneDesc",
                  "CFG_TSS_local_region","inCell_TSS_local_region","CFG_distance",
                  "coordinate.y", "direction.y", "transcriptOrientation.y",
                  "coord_location.y","CFG_TSS","CFG_coordID","Non.template.strand.y",
                  "coordinate.x", "direction.x", "transcriptOrientation.x",
                  "coord_location.x","coordID","Non.template.strand.x")

all_TSSs_curated <- ordered_TSSs[,keep_columns]

colnames_for_TSS_combined <- c("geneID","geneName","geneDesc",
                               "CFG_TSS_local_region","inCell_TSS_local_region","TSS_distance",
                               "CFG_coord","CFG_direction","CFG_TSS_orientation","CFG_TSS_location",
                               "CFG_TSS_ID_lower","CFG_TSS_ID_higher","CFG_NT_strand",
                               "inCell_coord","inCell_direction","inCell_TSS_orientation","inCell_TSS_location",
                               "inCell_TSS_ID","inCell_NT_strand")
colnames(all_TSSs_curated) <- colnames_for_TSS_combined

# Export to combine with output of promoter_motifs R script (add promoter annotations)
write.csv(all_TSSs_curated, 'fig2_inCellulo_vs_CFG/all_TSSs_needMotifs.csv',
          row.names = FALSE)

## Export non-template strand for CRP motif scanning ---------------------------

# CFG TSSs
CFG_TSS_sequences <- all_TSSs_curated[all_TSSs_curated$CFG_NT_strand != '',]
CFG_ordered <- CFG_TSS_sequences[order(CFG_TSS_sequences$CFG_TSS_ID_higher),]

write.fasta(sequences = strsplit(CFG_ordered$CFG_NT_strand, split=''),
            name = CFG_ordered$CFG_TSS_ID_higher, open = "w",
            as.string = FALSE,
            file.out = 'fig2_inCellulo_vs_CFG/all_NT_CFG.fasta')

# In cellulo TSSs
inCell_TSS_sequences <- all_TSSs_curated[all_TSSs_curated$inCell_NT_strand != '',]
inCell_ordered <- inCell_TSS_sequences[order(as.integer(inCell_TSS_sequences$inCell_TSS_ID)),]

write.fasta(sequences = strsplit(inCell_ordered$inCell_NT_strand, split=''),
            name = inCell_ordered$inCell_TSS_ID, open = "w",
            as.string = FALSE,
            file.out = 'fig2_inCellulo_vs_CFG/all_NT_inCellulo.fasta')

## End -------------------------------------------------------------------------