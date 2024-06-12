# Promoter annotation script and motif/histogram generation
# Authors: Ruby Froom
# Date: May 24, 2024

#### Libraries -----------------------------------------------------------------

# Required to read in genome FASTA files
library(seqinr)
# Required to read DNAStrings to reverse-complement and complement in 
# extract_promoter and annotate_promoter
library(Biostrings)
# Required for full_join feature in annotate_GO_terms
library(dplyr)
# Required for str_extract and other string manipulation functions
library(stringr)
# For joining columns based on a range
library(sqldf)
# For regex pattern matching with overlapping motifs
library(stringi)
# For generating motif logos
library(universalmotif)
library(ggplot2)
library(ggseqlogo)
library(readr)
library(textshape)
library(tidyr)

## Input variable names for reference genomes and TSS coordinates --------------

csv_name <- 'genome_files_misc/TSS_cortes_S1F.csv'
sample_name <- "inCell_Cortes"
genome_name <- 'genome_files_misc/TSS_cortes_refGenome.fasta'

sample_name <- 'CFG'
csv_name <- '5enrich_CRP/selectThreshold/DESeq2/noCRP_CRP_20counts_ends_results_Wald_local.csv'
genome_name <- 'genome_files_misc/Eco_Mtb_genome.fasta'

# Loading in files
genome <- read.fasta(genome_name, seqtype="DNA")[[1]]
coordinateDF <- read.csv(csv_name)

# For in-cell data, turn the coordinate into an integer value
coordinateDF$coordinate <- as.double(coordinateDF$coordinate)

# For CFG data, split coordinate and strand
coordinateDF$coordinate <- as.integer(stri_sub(coordinateDF$coordID, from = 1, to = -2))
coordinateDF$direction <- sapply(strsplit(coordinateDF$coordID,""), tail, 1)
coordinateDF <- coordinateDF[,c("coordinate","direction")]

## Initializes acceptable length ranges ----------------------------------------
## How many basepairs of the promoter are included from the extract_promoter dataframe?
upstream_end <- 70
downstream_end <- 10
TSS <- upstream_end + 1
# The number of nucleotides acceptable between the -10 and -35 element
max_spacer_length <- 20
min_spacer_length <- 15
# The distance range acceptable between the -7T and +1
min_minus7_distance <- 4
max_minus7_distance <- 10
# How long are the -10, extended -10 and -35 motifs?
min10_consensus_bp <- 6
extmin10_consensus_bp <- 5
min35_consensus_bp <- 6
# In the absence of a minus10 element, how long should the discriminator be?
discriminator_length <- 7
flanking_min35_motif <- '..........ttg...'

# Calculates the upper bound of the -10 region to add to the local indices to make them global when looking for min35
minus10_upper_bound <- TSS - min10_consensus_bp - max_minus7_distance + 1

# Initializes the regex promoter motifs
cre_motif <- "g"
minus5_G_motif <- ".g"
minus10_motif <- ".a...t"
minus7T_motif <- ".....t"
minus11A_motif <- ".a...."
extended_minus10_motif <- paste0(c("..[^a]","[tg]","."),collapse="")
W <- "[at]"
D <- "[^c]"
H <- "[^g]"
K <- "[gt]"
minus35_motif <- paste0(c(W,D,D,H,".."),collapse="")
minus35_motif <- 'ttg...'
proximal_UP_element_motif <- "[ca]aaa[ta]a[ag].[ag]"
distal_UP_element_motif <- "..[^c][^c][^c][at][at][at][at][at][gt]t[^g]."
full_UP_element_motif <- "..a[at][^g][at][at][at][at]tt[gt][gt]...[^g][at]a..."
permissive_UP_element_motif <- "[at][at][at]"

proximal_UP_motif_length <- 9
distal_UP_motif_length <- 14


# Gets GC content from a sequence
get_gc_content = function(seq){
  return (mean(seq == "C" | seq == "G" | seq == "c" | seq == "g"))
}

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
                    'Template.strand',"IndexingCoordinate","TSS_local_region")
  
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
  
  sense_TSS_local_start <- sense_DF$IndexingCoordinate - 6
  sense_TSS_local_end <- sense_DF$IndexingCoordinate + 6
  sense_DF$TSS_local_region <- paste(str_sub(double_genome, sense_TSS_local_start, sense_DF$IndexCoordMin1),
                                     sense_TSS,
                                     str_sub(double_genome, sense_DF$IndexCoordPlus1, sense_TSS_local_end),sep='')
  
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
  
  antisense_extract <- str_sub(double_genome, antisense_DF$IndexingCoordinate - 6, antisense_DF$IndexingCoordinate + 6)
  antisense_DNAseq <- DNAStringSet(antisense_extract)
  antisense_DF$TSS_local_region <- as.character(reverseComplement(antisense_DNAseq))
  
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
sequences <- data.frame(extract_coordinate_environment(genome,
                                                       coordinateDF,
                                                       upstream_end,
                                                       downstream_end))

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

### Generate sequence logos for TSS region -------------------------------------

TSS_relEntropy <- relative_entropy_calc(unique(toupper(sequences$TSS_local_region)),
                      paste('fig2_inCellulo_vs_CFG/',sample_name,'_TSS',sep=''))

p <- ggseqlogo(TSS_relEntropy, method = "custom", seq_type = 'dna',
               col_scheme = cs)
TSS_logo <- p + theme_classic() + labs(y = "Relative Entropy")

ggsave(filename = paste('fig2_inCellulo_vs_CFG/extData_fig4d_relEntropy_',sample_name,
                        '_TSSlogo.png',sep=''),
       plot = TSS_logo,
       width = 300,
       height = 100,
       units = 'mm',
       dpi = 300)

## Annotate promoter features and assign weights to choose ---------------------
## most likely set of regulatory elements --------------------------------------
annotate_promoter = function(promoter_DF) {
  
  promoter_DF <- sequences
  # Initializing the output dataframe
  features <- c("TSS_nt","CRE","Min10","Disc","DiscLength","DiscGC","Min5G",
                "Min11A","Min12T","ExtMin10","Min35","flanking_min35","Spacer",
                "SpacerLength","SpacerGC","PermissiveUP","UPregion")
  promoter_feature_DF <- data.frame(matrix(ncol=(ncol(promoter_DF)+length(features)),
                                           nrow=nrow(promoter_DF)))
  col_names <- c(colnames(promoter_DF),features)
  colnames(promoter_feature_DF) <- col_names
  promoter_feature_DF[,1:ncol(promoter_DF)] <- promoter_DF
  promoter_feature_DF[,((ncol(promoter_DF)+1):(ncol(promoter_DF)+length(features)))] <- 0
  
  coordinates <- promoter_feature_DF$coordinate
  
  for (i in 1:length(coordinates)){
    # Turns the promoter string into a subset-selectable vector (non-template strand)
    promoter_gene_string <- promoter_feature_DF[promoter_feature_DF$coordinate == coordinates[i],
                                                "Non.template.strand"]
    promoter_gene_vector <- strsplit(promoter_gene_string,split='')
    promoter_gene <- promoter_gene_vector[[1]]
    
    # Generate different strings corresponding to regions to search for promoter motifs
    # to minimize detection of the motifs in other promoter regions when pulling out indices
    promoter <- promoter_gene[0:TSS+1]
    promoter_string <- paste0(promoter,collapse="")
    
    # Fills in starting nucleotide on template strand
    promoter_NT <- promoter_gene[1:TSS]
    NT_string <- paste0(promoter,collapse="")
    TSS_bp <- promoter_NT[TSS]
    promoter_feature_DF[promoter_feature_DF$coordinate == coordinates[i],"TSS_nt"] <- paste0(TSS_bp,collapse='')
    
    # Sets region for the CRE motif based on the TSS
    cre_region <- paste0(promoter[TSS+1],collapse="")
    
    # Searches for CRE motif and indicates if it is there
    cre <- str_extract(cre_region, cre_motif)
    boolean_cre <- str_detect(cre_region, cre_motif)
    
    if(boolean_cre == TRUE) {
      promoter_feature_DF[promoter_feature_DF$coordinate == coordinates[i],"CRE"] = 1
    }
    
    # Searches for a -10 element (default: NANNNT in the non-template strand)
    # Constrains the searchable region: default is -7T being between 4-10 bases away from +1
    minus10_region <- paste0(promoter[(TSS-(max_minus7_distance+min10_consensus_bp)):(TSS-min_minus7_distance-1)],collapse="")
    minus10 <- stri_match_all_regex(minus10_region, paste('(?=(',minus10_motif,'))',sep=''))[[1]][,2]
    
    # If there is a -10 element:
    # Fill in the sequence, 
    # Indicate how far away the -7T is from +1
    # Search for the discriminator, extended -10, -35, and UP elements based on -10 element location
    if(!is.na(minus10[1])) {
      
      min10_DF <- data.frame(matrix(data=NA, nrow=length(minus10), ncol=ncol(promoter_feature_DF)))
      colnames(min10_DF) <- colnames(promoter_feature_DF)
      
      # For each possible minus10 element found:
      for (j in 1:length(minus10)) {
        
        # Get the indices of rows that are empty to indicate where to start filling in
        IndexFirstNA <- which(is.na(min10_DF[,1]))[1]
        
        min10_DF[IndexFirstNA,] <- promoter_feature_DF[promoter_feature_DF$coordinate == coordinates[i],]
        min10_DF$"Min10"[IndexFirstNA] <- minus10[j]
        
        # Find the local indices of the minus10 element within the -10 region
        minus10_local_index <- str_locate_all(pattern=minus10[j], minus10_region)
        minus10_index <- minus10_local_index[[1]]
        
        # Set the minus10 upstream and downstream boundaries
        upstream_end_minus10 <- minus10_index[1] + minus10_upper_bound - 2
        downstream_end_minus10 <- minus10_index[2] + minus10_upper_bound - 2
        
        # Sets the discriminator region and where to look for -5G
        minus5_region <- paste0(promoter[(downstream_end_minus10+1):(downstream_end_minus10+2)],collapse="")
        discriminator_region <- promoter[(downstream_end_minus10+1):(TSS-1)]
        
        min10_DF$"Disc"[IndexFirstNA] <- paste0(discriminator_region,collapse='')
        min10_DF$"DiscLength"[IndexFirstNA] <- length(discriminator_region)
        
        # Searches for -5G in discriminator and adds if it is there
        # Also adds GC content of whole discriminator region (end of -10 motif to +1)
        minus5_G <- str_extract(minus5_region, minus5_G_motif)
        boolean_minus5_G <- str_detect(minus5_region, minus5_G_motif)
        discriminator_GC <- get_gc_content(discriminator_region)
        min10_DF$"DiscGC"[IndexFirstNA] <- discriminator_GC * 100
        
        if(boolean_minus5_G == TRUE) {
          min10_DF$"Min5G"[IndexFirstNA] = 1
        }
        
        # Searches for -12T relative to -10 position
        minus12_region <- paste0(promoter[(upstream_end_minus10-1):(upstream_end_minus10)],collapse="")
        minus12T <- str_extract(minus12_region, ".t")
        boolean_minus12T <- str_detect(minus12_region, ".t")
        
        if(boolean_minus12T == TRUE) {
          min10_DF$"Min12T"[IndexFirstNA] = 1
        }
        
        # Sets region for the extended -10 based on the upstream end of the -10
        extended_minus10_region <- paste0(promoter[(upstream_end_minus10-extmin10_consensus_bp):(upstream_end_minus10-1)],
                                          collapse="")
        
        # Adds extended -10 region
        min10_DF$"ExtMin10"[IndexFirstNA] = extended_minus10_region
        
        # Sets region for -35 based on the upstream end of the -10
        minus35_upstream_bound <- upstream_end_minus10 - (max_spacer_length+min35_consensus_bp)
        minus35_downstream_bound <- upstream_end_minus10 - min_spacer_length - 1
        minus35_region <- paste0(promoter[minus35_upstream_bound:minus35_downstream_bound],collapse='')
        flanking_min35_region <- paste0(promoter[(minus35_upstream_bound-10):minus35_downstream_bound],collapse='')
        
        # Searches for -35 element
        minus35 <- stri_match_all_regex(minus35_region, paste('(?=(',minus35_motif,'))',sep=''))[[1]][,2]
        flanking_min35 <- stri_match_all_regex(flanking_min35_region, paste('(?=(',flanking_min35_motif,'))',sep=''))[[1]][,2]
        
        # If there is a -35 element:
        # Calculates SpacerLength
        # Searches for UP element
        if(!is.na(minus10[1]) & !is.na(minus35[1])) {
          
          min35_DF <- data.frame(matrix(data=NA, nrow=length(minus35), ncol=ncol(promoter_feature_DF)))
          colnames(min35_DF) <- colnames(promoter_feature_DF)
          
          # For each possible minus10 element found:
          for (k in 1:length(minus35)) {
            
            min35_DF[k,] <- min10_DF[IndexFirstNA,]
            min35_DF$"Min35"[k] <- minus35[k]
            min35_DF$"flanking_min35" <- flanking_min35[k]
            
            # Locates indices of -35 element from within the whole promoter string
            minus35_index <- str_locate_all(pattern=minus35[k], paste0(promoter[(minus35_upstream_bound):(minus35_downstream_bound)],collapse=""))
            minus35_index <- minus35_index[[1]]
            
            # Initialize -35 boundaries so they can be updated from within "if" statements
            upstream_end_minus35 <- minus35_index[1] + minus35_upstream_bound
            downstream_end_minus35 <- minus35_index[2] + minus35_upstream_bound
            min35_DF$"SpacerLength"[k] <- upstream_end_minus10 - downstream_end_minus35
            
            spacer_region <- promoter[(downstream_end_minus35):(upstream_end_minus10-1)]
            min35_DF$"Spacer"[k] <- paste0(spacer_region,collapse='')
            spacer_GC <- get_gc_content(spacer_region)
            min35_DF$"SpacerGC"[k] <- spacer_GC * 100
            
            proximal_up_boundary <- upstream_end_minus35 - (proximal_UP_motif_length + 1)
            proximal_down_boundary <- upstream_end_minus35 + 1
            distal_up_boundary <- upstream_end_minus35 - (proximal_UP_motif_length + distal_UP_motif_length + 1)
            distal_down_boundary <- upstream_end_minus35 - (proximal_UP_motif_length + 1)
            
            permissive_UP_element_region <- paste0(promoter[distal_up_boundary:(proximal_down_boundary-3)],collapse="")
            min35_DF$"UPregion"[k] = paste0(permissive_UP_element_region,collapse="")
            permissive_UP_element <- str_extract_all(permissive_UP_element_region, permissive_UP_element_motif)
            boolean_permissiveUP <- str_detect(permissive_UP_element_region, permissive_UP_element_motif)  
            
            if(boolean_permissiveUP == TRUE) {
              min35_DF$"PermissiveUP"[k] = paste0(permissive_UP_element[[1]],collapse=",")
            }
          }
          
          min10_DF <- rbind(min35_DF, min10_DF)
          min10_DF <- min10_DF[!(min10_DF$Min35 == 0),]
          
        }
        
        # If there is no -35 element but there is a -10 element:
        # Searches for an UP element
        # Upstream bound assumptions: 20bp spacer, 6bp -35 motif stand-in, + UP element motif length
        # Downstream bound assumptions: 15bp spacer, 6bp -35 motif stand-in (for distal, + proximal UP element 9bp motif)
        if(!is.na(minus10[1]) & is.na(minus35[1])) {
          
          proximal_up_boundary <- upstream_end_minus10 - (max_spacer_length + min35_consensus_bp + proximal_UP_motif_length + 1)
          proximal_down_boundary <- upstream_end_minus10 - (min_spacer_length + min35_consensus_bp)
          distal_up_boundary <- upstream_end_minus10 - (max_spacer_length + min35_consensus_bp + proximal_UP_motif_length + distal_UP_motif_length + 1 - 6)
          distal_down_boundary <- upstream_end_minus10 - (min_spacer_length + min35_consensus_bp + proximal_UP_motif_length + 1)
          permissive_UP_element_region <- paste0(promoter[distal_up_boundary:(proximal_down_boundary-3)],collapse="")
          
          min10_DF$"UPregion"[IndexFirstNA] = paste0(permissive_UP_element_region,collapse="")
          permissive_UP_element <- str_extract_all(permissive_UP_element_region, permissive_UP_element_motif)
          boolean_permissiveUP <- str_detect(permissive_UP_element_region, permissive_UP_element_motif)   
          
          if(boolean_permissiveUP == TRUE) {
            min10_DF$"PermissiveUP"[IndexFirstNA] = paste0(permissive_UP_element[[1]],collapse=",")
          }
        }
      }
      
      promoter_feature_DF[promoter_feature_DF$coordinate == coordinates[i],] <- min10_DF[1,]
      promoter_feature_DF <- rbind(min10_DF[2:(nrow(min10_DF)),], promoter_feature_DF)
      promoter_feature_DF <- promoter_feature_DF[order(promoter_feature_DF$coordinate),]
      promoter_feature_DF <- promoter_feature_DF[!is.na(promoter_feature_DF$coordinate),]
    }
    
    minus11A <- stri_match_all_regex(minus10_region, paste('(?=(',minus11A_motif,'))',sep=''))[[1]][,2]
    
    # If there's no -10 element but there is a potential -11A:
    # Fill in the sequence, 
    # Indicate how far away the -7 position is from +1
    # Search for the discriminator, extended -10, -35, and UP elements based on -11A element location
    if(is.na(minus10[1]) & !is.na(minus11A[1])) {
      
      min11A_DF <- data.frame(matrix(data=NA, nrow=length(minus11A), ncol=ncol(promoter_feature_DF)))
      colnames(min11A_DF) <- colnames(promoter_feature_DF)
      
      # For each possible minus10 element found:
      for (j in 1:length(minus10)) {
        
        # Get the indices of rows that are empty to indicate where to start filling in
        IndexFirstNA <- which(is.na(min11A_DF[,1]))[1]
        
        min11A_DF[IndexFirstNA,] <- promoter_feature_DF[promoter_feature_DF$coordinate == coordinates[i],]
        min11A_DF$"Min11A"[IndexFirstNA] <- minus11A[j]
        
        # Find the local indices of the minus10 element within the -10 region
        minus11A_local_index <- str_locate_all(pattern=minus11A[j], minus10_region)
        minus11A_index <- minus11A_local_index[[1]]
        
        # Set the minus11A upstream and downstream boundaries
        upstream_end_minus11A <- minus11A_index[1] + minus10_upper_bound - 2
        downstream_end_minus11A <- minus11A_index[2] + minus10_upper_bound - 2
        
        # Sets the discriminator region and where to look for -5G
        minus5_region <- paste0(promoter[(downstream_end_minus11A+1):(downstream_end_minus11A+2)],collapse="")
        discriminator_region <- promoter[(downstream_end_minus11A+1):(TSS-1)]
        
        min11A_DF$"Disc"[IndexFirstNA] <- paste0(discriminator_region,collapse='')
        min11A_DF$"DiscLength"[IndexFirstNA] <- length(discriminator_region)
        
        # Searches for -5G in discriminator and adds if it is there
        # Also adds GC content of whole discriminator region (end of -10 motif to +1)
        minus5_G <- str_extract(minus5_region, minus5_G_motif)
        boolean_minus5_G <- str_detect(minus5_region, minus5_G_motif)
        discriminator_GC <- get_gc_content(discriminator_region)
        min11A_DF$"DiscGC"[IndexFirstNA] <- discriminator_GC * 100
        
        if(boolean_minus5_G == TRUE) {
          min11A_DF$"Min5G"[IndexFirstNA] = 1
        }
        
        # Searches for -12T relative to -10 position
        minus12_region <- paste0(promoter[(upstream_end_minus11A-1):(upstream_end_minus11A)],collapse="")
        minus12T <- str_extract(minus12_region, ".t")
        boolean_minus12T <- str_detect(minus12_region, ".t")
        
        if(boolean_minus12T == TRUE) {
          min11A_DF$"Min12T"[IndexFirstNA] = 1
        }
        
        # Sets region for the extended -10 based on the upstream end of the -10
        extended_minus10_region <- paste0(promoter[(upstream_end_minus11A-extmin10_consensus_bp):(upstream_end_minus11A-1)],
                                          collapse="")
        
        # Searches for extended -10 element and adds if it is there
        ext_minus10 <- str_extract(extended_minus10_region, extended_minus10_motif)
        boolean_ext_minus10 <- str_detect(extended_minus10_region, extended_minus10_motif)
        min11A_DF$"ExtMin10"[IndexFirstNA] = extended_minus10_region
        
        # Sets region for -35 based on the upstream end of the -10
        minus35_upstream_bound <- upstream_end_minus11A - (max_spacer_length+min35_consensus_bp)
        minus35_downstream_bound <- upstream_end_minus11A - min_spacer_length - 1
        minus35_region <- paste0(promoter[minus35_upstream_bound:minus35_downstream_bound],collapse='')
        flanking_min35_region <- paste0(promoter[(minus35_upstream_bound-10):minus35_downstream_bound],collapse='')
        
        # Searches for -35 element
        minus35 <- stri_match_all_regex(minus35_region, paste('(?=(',minus35_motif,'))',sep=''))[[1]][,2]
        flanking_min35 <- stri_match_all_regex(flanking_min35_region, paste('(?=(',flanking_min35_motif,'))',sep=''))[[1]][,2]
        
        # If there is a -35 element:
        # Calculates SpacerLength
        # Searches for UP element
        if(!is.na(minus11A[1]) & !is.na(minus35[1]) & is.na(minus10[1])) {
          
          min35_DF <- data.frame(matrix(data=NA, nrow=length(minus35), ncol=ncol(promoter_feature_DF)))
          colnames(min35_DF) <- colnames(promoter_feature_DF)
          
          # For each possible minus10 element found:
          for (k in 1:length(minus35)) {
            
            min35_DF[k,] <- min11A_DF[IndexFirstNA,]
            min35_DF$"Min35"[k] <- minus35[k]
            min35_DF$"flanking_min35" <- flanking_min35[k]
            
            # Locates indices of -35 element from within the whole promoter string
            minus35_index <- str_locate_all(pattern=minus35[k], paste0(promoter[(minus35_upstream_bound):(minus35_downstream_bound)],collapse=""))
            minus35_index <- minus35_index[[1]]
            
            # Initialize -35 boundaries so they can be updated from within "if" statements
            upstream_end_minus35 <- minus35_index[1] + minus35_upstream_bound
            downstream_end_minus35 <- minus35_index[2] + minus35_upstream_bound
            min35_DF$"SpacerLength"[k] <- upstream_end_minus11A - downstream_end_minus35
            
            spacer_region <- promoter[(downstream_end_minus35):(upstream_end_minus11A-1)]
            min35_DF$"Spacer"[k] <- paste0(spacer_region,collapse='')
            spacer_GC <- get_gc_content(spacer_region)
            min35_DF$"SpacerGC"[k] <- spacer_GC * 100
            
            proximal_up_boundary <- upstream_end_minus35 - (proximal_UP_motif_length + 1)
            proximal_down_boundary <- upstream_end_minus35 - 1
            distal_up_boundary <- upstream_end_minus35 - (proximal_UP_motif_length + distal_UP_motif_length + 1)
            distal_down_boundary <- upstream_end_minus35 - (proximal_UP_motif_length + 1)
            
            permissive_UP_element_region <- paste0(promoter[distal_up_boundary:(proximal_down_boundary-3)],collapse="")
            min35_DF$"UPregion"[k] = paste0(permissive_UP_element_region,collapse="")
            permissive_UP_element <- str_extract_all(permissive_UP_element_region, permissive_UP_element_motif)
            boolean_permissiveUP <- str_detect(permissive_UP_element_region, permissive_UP_element_motif)  
            
            if(boolean_permissiveUP == TRUE) {
              min35_DF$"PermissiveUP"[k] = paste0(permissive_UP_element[[1]],collapse=",")
            }
          }
          
          min11A_DF <- rbind(min35_DF, min11A_DF)
          min11A_DF <- min11A_DF[!(min11A_DF$Min35 == 0),]
          
        }
        
        # If there is no -35 element but there is a -11A element:
        # Searches for an UP element
        # Upstream bound assumptions: 20bp spacer, 6bp -35 motif stand-in, + UP element motif length
        # Downstream bound assumptions: 15bp spacer, 6bp -35 motif stand-in (for distal, + proximal UP element 9bp motif)
        if(!is.na(minus11A[1]) & is.na(minus35[1])) {
          
          proximal_up_boundary <- upstream_end_minus11A - (max_spacer_length + min35_consensus_bp + proximal_UP_motif_length + 1)
          proximal_down_boundary <- upstream_end_minus11A - (min_spacer_length + min35_consensus_bp + 1)
          distal_up_boundary <- upstream_end_minus11A - (max_spacer_length + min35_consensus_bp + proximal_UP_motif_length + distal_UP_motif_length + 1)
          distal_down_boundary <- upstream_end_minus11A - (min_spacer_length + min35_consensus_bp + proximal_UP_motif_length + 1)
          permissive_UP_element_region <- paste0(promoter[distal_up_boundary:(proximal_down_boundary-3)],collapse="")
          
          min11A_DF$"UPregion"[IndexFirstNA] = paste0(permissive_UP_element_region,collapse="")
          permissive_UP_element <- str_extract_all(permissive_UP_element_region, permissive_UP_element_motif)
          boolean_permissiveUP <- str_detect(permissive_UP_element_region, permissive_UP_element_motif)   
          
          if(boolean_permissiveUP == TRUE) {
            min11A_DF$"PermissiveUP"[IndexFirstNA] = paste0(permissive_UP_element[[1]],collapse=",")
          }
        }
      }
      
      promoter_feature_DF[promoter_feature_DF$coordinate == coordinates[i],] <- min11A_DF[1,]
      promoter_feature_DF <- rbind(min11A_DF[2:(nrow(min11A_DF)),], promoter_feature_DF)
      promoter_feature_DF <- promoter_feature_DF[order(promoter_feature_DF$coordinate),]
      promoter_feature_DF <- promoter_feature_DF[!is.na(promoter_feature_DF$coordinate),]
    }
  }
  return (promoter_feature_DF)
}

promoters <- annotate_promoter(sequences)

# Initializing weights for crude linear model
weights <- read.delim("genome_files_misc/promoter_weights.txt",sep=" ",skip=6,header=F)
rownames(weights) <- seq(from=-41, to=-1,by=1)
colnames(weights) <- c("a","c","g","t")

promoters$DiscLength <- promoters$DiscLength + 1
promoters$CombinedMin10_Min11A <- 0
promoters[promoters$Min10 == 0,"CombinedMin10_Min11A"] <- promoters[promoters$Min10 == 0,"Min11A"]
promoters[promoters$Min11A == 0,"CombinedMin10_Min11A"] <- promoters[promoters$Min11A == 0,"Min10"]

promoters$Min10weight <- max(weights)*min10_consensus_bp
promoters$ExtMin10weight <- max(weights)*extmin10_consensus_bp
promoters$Min35weight <- max(weights)*min35_consensus_bp

for (i in 1:nrow(promoters)) {
  if (promoters$CombinedMin10_Min11A[i] != 0){
    
    Min10vector <- strsplit(promoters$CombinedMin10_Min11A[i],split="")[[1]]
    Min12 <- Min10vector[1] 
    Min11 <- Min10vector[2]
    Min10 <- Min10vector[3]
    Min9 <- Min10vector[4]
    Min8 <- Min10vector[5]
    Min7 <- Min10vector[6]
    
    Min12value <- weights[rownames(weights) == -12,colnames(weights) == Min12]
    Min11value <- weights[rownames(weights) == -11,colnames(weights) == Min11]
    Min10value <- weights[rownames(weights) == -10,colnames(weights) == Min10]
    Min9value <- weights[rownames(weights) == -9,colnames(weights) == Min9]
    Min8value <- weights[rownames(weights) == -8,colnames(weights) == Min8]
    Min7value <- weights[rownames(weights) == -7,colnames(weights) == Min7]
    
    promoters$Min10weight[i] <- sum(Min12value,Min11value,Min10value,Min9value,Min8value,Min7value)
    
  }
  
  if (promoters$ExtMin10[i] != 0){
    
    ExtMin10vector <- strsplit(promoters$ExtMin10[i],split="")[[1]]
    Min17 <- ExtMin10vector[1] 
    Min16 <- ExtMin10vector[2]
    Min15 <- ExtMin10vector[3]
    Min14 <- ExtMin10vector[4]
    Min13 <- ExtMin10vector[5]
    
    Min17value <- weights[rownames(weights) == -17,colnames(weights) == Min15]
    Min16value <- weights[rownames(weights) == -16,colnames(weights) == Min15]
    Min15value <- weights[rownames(weights) == -15,colnames(weights) == Min15]
    Min14value <- weights[rownames(weights) == -14,colnames(weights) == Min14]
    Min13value <- weights[rownames(weights) == -13,colnames(weights) == Min13]
    
    promoters$ExtMin10weight[i] <- sum(Min17value,Min16value,Min15value,Min14value,Min13value)
    
  }
  
  if (promoters$Min35[i] != 0){
    
    Min35vector <- strsplit(promoters$Min35[i],split="")[[1]]
    Min36 <- Min35vector[1] 
    Min35 <- Min35vector[2]
    Min34 <- Min35vector[3]
    Min33 <- Min35vector[4]
    Min32 <- Min35vector[5]
    Min31 <- Min35vector[6]
    
    Min36value <- weights[rownames(weights) == -36,colnames(weights) == Min36]
    Min35value <- weights[rownames(weights) == -35,colnames(weights) == Min35]
    Min34value <- weights[rownames(weights) == -34,colnames(weights) == Min34]
    Min33value <- weights[rownames(weights) == -33,colnames(weights) == Min33]
    Min32value <- weights[rownames(weights) == -32,colnames(weights) == Min32]
    Min31value <- weights[rownames(weights) == -31,colnames(weights) == Min31]
    
    promoters$Min35weight[i] <- sum(Min36value,Min35value,Min34value,Min33value,Min32value,Min31value)
    
  }
  
}

# Weight minus10 based on discriminator length
promoters$DiscLength <- promoters$DiscLength - 1
promoters$DiscWeight <- 1
rawDiscpromoters <- na.omit(read.csv("genome_files_misc/%TSS_MASTER.csv"))
DiscLengthPercentages <- colSums(rawDiscpromoters)
promoters[promoters$DiscLength == 6,"DiscWeight"] <- DiscLengthPercentages["X6"] / DiscLengthPercentages["X7"]
promoters[promoters$DiscLength == 5,"DiscWeight"] <- DiscLengthPercentages["X5"] / DiscLengthPercentages["X7"]
promoters[promoters$DiscLength == 4,"DiscWeight"] <- DiscLengthPercentages["X4"] / DiscLengthPercentages["X7"]
promoters[promoters$DiscLength == 8,"DiscWeight"] <- DiscLengthPercentages["X8"] / DiscLengthPercentages["X7"]
promoters[promoters$DiscLength == 9,"DiscWeight"] <- DiscLengthPercentages["X9"] / DiscLengthPercentages["X7"]

promoters$Min10weightDiscAdjusted <- promoters$Min10weight
promoters[promoters$Min10weight < 0,"Min10weightDiscAdjusted"] <- promoters[promoters$Min10weight < 0,"Min10weight"] * promoters[promoters$Min10weight < 0,"DiscWeight"]
promoters[promoters$Min10weight > 0,"Min10weightDiscAdjusted"] <- promoters[promoters$Min10weight > 0,"Min10weight"] / promoters[promoters$Min10weight > 0,"DiscWeight"]

promoters$ExtMin10weightDiscAdjusted <- promoters$ExtMin10weight
promoters[promoters$ExtMin10weight < 0,"ExtMin10weightDiscAdjusted"] <- promoters[promoters$ExtMin10weight < 0,"ExtMin10weight"] * promoters[promoters$ExtMin10weight < 0,"DiscWeight"]
promoters[promoters$ExtMin10weight > 0,"ExtMin10weightDiscAdjusted"] <- promoters[promoters$ExtMin10weight > 0,"ExtMin10weight"] / promoters[promoters$ExtMin10weight > 0,"DiscWeight"]

promoters$SpacerWeight <- 0.0001
# KBK2 rates from Mulligan, Brosius & McClure 1985
# Extrapolations for 15, 19 and 20bp came from intro text of same paper
rate17bp <- 8.4
rate16bp <- 1.5
rate18bp <- 3.5
promoters[promoters$SpacerLength == 17,"SpacerWeight"] <- 1
promoters[promoters$SpacerLength == 16,"SpacerWeight"] <- rate16bp / rate17bp
promoters[promoters$SpacerLength == 18,"SpacerWeight"] <- rate18bp / rate17bp
promoters[promoters$SpacerLength == 15,"SpacerWeight"] <- (rate16bp * 0.02) / rate17bp
promoters[promoters$SpacerLength == 20,"SpacerWeight"] <- (rate18bp * 0.25) / rate17bp
promoters[promoters$SpacerLength == 19,"SpacerWeight"] <- (rate18bp * 0.25) / rate17bp

promoters$Min35weightSpacerAdjusted <- promoters$Min35weight
promoters[promoters$Min35weight < 0,"Min35weightSpacerAdjusted"] <- promoters[promoters$Min35weight < 0,"Min35weight"] * promoters[promoters$Min35weight < 0,"SpacerWeight"]
promoters[promoters$Min35weight > 0,"Min35weightSpacerAdjusted"] <- promoters[promoters$Min35weight > 0,"Min35weight"] / promoters[promoters$Min35weight > 0,"SpacerWeight"]

promoters$OverallWeightNoSpacing <- rowSums(promoters[,c("Min10weight","ExtMin10weight","Min35weight")])
promoters$OverallWeightSpacing <- rowSums(promoters[,c("Min10weightDiscAdjusted","ExtMin10weightDiscAdjusted","Min35weightSpacerAdjusted")])

# Pick promoter option with the lowest weights
promoters$LowestWeightNoSpacing <- "No"
promoters$LowestWeightSpacing <- "No"

for (i in 1:nrow(promoters)) {
  
  # Isolate rows where coordinates that are the same
  promoters_subset <- promoters[promoters$coordinate == promoters$coordinate[i],]
  promoters_subset[promoters_subset$OverallWeightSpacing == min(promoters_subset$OverallWeightSpacing),"LowestWeightSpacing"] <-   promoters_subset[promoters_subset$OverallWeightSpacing == min(promoters_subset$OverallWeightSpacing),"OverallWeightSpacing"]
  promoters_subset[promoters_subset$OverallWeightNoSpacing == min(promoters_subset$OverallWeightNoSpacing),"LowestWeightNoSpacing"] <-   promoters_subset[promoters_subset$OverallWeightNoSpacing == min(promoters_subset$OverallWeightNoSpacing),"OverallWeightNoSpacing"]
  promoters[promoters$coordinate == promoters$coordinate[i],] <- promoters_subset
  
}


#### Choose one set of regulatory elements per TSS as the most likely set ------
promoters$WeightConflict <- "No"
promoters[(promoters$LowestWeightNoSpacing == "No" & promoters$LowestWeightSpacing != "No") | (promoters$LowestWeightNoSpacing != "No" & promoters$LowestWeightSpacing == "No"),"WeightConflict"] <- "Yes"
promoters$TopPromoter <- 'Yes'
promoters[promoters$LowestWeightSpacing == "No","TopPromoter"] <- "No"

top_promoters <- unique(promoters[promoters$TopPromoter == 'Yes',])
top_promoters$coordID <- paste(top_promoters$coordinate,
                               top_promoters$direction,
                               sep='')

## Export supplementary table --------------------------------------------------

# Only keep columns with coordIDs for joining
# and annotations to add to other supplementary DF from inCell_TSS script
colsToKeep <- c("coordID","UPregion","Min35","Spacer","SpacerLength",
                "ExtMin10","CombinedMin10_Min11A","Disc","DiscLength","TSS_nt")

promoters_cleaned <- top_promoters[,colsToKeep]

write.csv(promoters_cleaned, paste('fig2_inCellulo_vs_CFG/',sample_name,"_annotated_promoters.csv",sep=''))

### Generate -10 consensus motif -----------------------------------------------
top_promoters$min10_region <- 0
top_promoters$min10_region_length <- 0

for (i in 1:nrow(top_promoters)) {
  
  NT_vector <- strsplit(top_promoters$Non.template.strand[i],
                        split='')[[1]]
  
  if (top_promoters$CombinedMin10_Min11A[i] != 0) {
    
    min10_location <- unlist(gregexpr(top_promoters$CombinedMin10_Min11A[i],
                                      top_promoters$Non.template.strand[i]))
    if (length(min10_location) > 1) {
      min10_location <- min10_location[2]
    }
    top_promoters$min10_region[i] <- paste0(NT_vector[(min10_location - 6):(min10_location+11)],
                                        collapse='')
    top_promoters$min10_region_length[i] <- length(strsplit(top_promoters$min10_region[i],split='')[[1]])
    
    if (top_promoters$min10_region_length[i] > 19) {

          top_promoters$min10_region[i] <- paste0(strsplit(top_promoters$min10_region[i],
                                                   split='')[[1]][1:18],collapse='')
    }
    
  }
  else {
    top_promoters$min10_region[i] <- paste0(NT_vector[53:70],
                                        collapse='')
    
  }
  
}

# There are a few top_promoters where the -10 motif occurs twice, close to each other
# For these, choose the first motif that occurs
for (i in 1:nrow(top_promoters)) {
  
  search_NAs <- unlist(gregexpr("NA",
                                top_promoters$min10_region[i]))
  
  if (length(search_NAs) > 1) {
    min10_location <- unlist(gregexpr(top_promoters$CombinedMin10_Min11A[i],
                                      top_promoters$Non.template.strand[i]))
    top_promoters$min10_region[i] <- paste0(NT_vector[(min10_location[1] - 6):(min10_location[1]+11)],
                                        collapse='')
  }
}

min10_relEntropy <- relative_entropy_calc(unique(top_promoters$min10_region),
                                        paste('fig2_inCellulo_vs_CFG/extData_fig4d_',sample_name,'_min10',sep=''))

p <- ggseqlogo(min10_relEntropy, method = "custom", seq_type = 'dna',
               col_scheme = cs)
min10_logo <- p + theme_classic() + labs(y = "Relative Entropy")

ggsave(filename = paste('fig2_inCellulo_vs_CFG/extData_fig4d_relEntropy_',sample_name,
                        '_min10logo.png',sep=''),
       plot = min10_logo,
       width = 300,
       height = 100,
       units = 'mm',
       dpi = 300)

## Plotting discriminator lengths ----------------------------------------------

no_zero_disc <- top_promoters %>%
  filter(DiscLength != 0) 

percent_with_disc <- (nrow(no_zero_disc) / nrow(top_promoters)) * 100

disc_hist <- ggplot(no_zero_disc, aes(x = DiscLength)) +
  geom_histogram(bins = 7, fill = "grey",
                 color = "black") +
  theme_minimal() +
  scale_x_continuous(n.breaks = 7)

ggsave(filename = paste('fig2_inCellulo_vs_CFG/extData_fig4d_',sample_name,
                        '_discHist.png',sep=''),
       plot = disc_hist,
       width = 300,
       height = 100,
       units = 'mm',
       dpi = 300)

### Generate -35 consensus motif -----------------------------------------------
top_promoters$min35_region <- 0
top_promoters$min35_region_length <- 0

for (i in 1:nrow(top_promoters)) {
  
  NT_vector <- strsplit(top_promoters$Non.template.strand[i],
                        split='')[[1]]
  
  # If a -35 was identified:
  if (top_promoters$Min35[i] != 0) {
    
    # Find the location of the -35
    min35_location <- unlist(gregexpr(top_promoters$Min35[i],
                                      top_promoters$Non.template.strand[i]))
    
    # If the motif occurs twice, choose the second option 
    # (empirically, this leads to the correct spacing for -35)
    if (length(min35_location) > 1) {
      min35_location <- min35_location[2]
    }
    top_promoters$min35_region[i] <- paste0(NT_vector[(min35_location - 6):(min35_location+11)],
                                        collapse='')
    top_promoters$min35_region_length[i] <- length(strsplit(top_promoters$min35_region[i],
                                                            split='')[[1]])
    
    if (top_promoters$min35_region_length[i] > 19) {
      top_promoters$min35_region[i] <- paste0(strsplit(top_promoters$min35_region[i],
                                                   split='')[[1]][1:18],collapse='')
    }
    
  }
  else {
    top_promoters$min35_region[i] <- paste0(NT_vector[36:53],
                                        collapse='')
    
  }
  
}

# There are a few top_promoters where the -35 motif occurs twice, close to each other
# For these, choose the first motif that occurs
for (i in 1:nrow(top_promoters)) {
  
  search_NAs <- unlist(gregexpr("NA",
                                top_promoters$min35_region[i]))
  
  if (length(search_NAs) > 1) {
    min35_location <- unlist(gregexpr(top_promoters$Min35[i],
                                      top_promoters$Non.template.strand[i]))
    top_promoters$min35_region[i] <- paste0(NT_vector[(min35_location[1] - 6):(min35_location[1]+11)],
                                        collapse='')
  }
}

min35_relEntropy <- relative_entropy_calc(unique(top_promoters$min35_region),
                                        paste('fig2_inCellulo_vs_CFG/extData_fig4d_',sample_name,'_min35',sep=''))

p <- ggseqlogo(min35_relEntropy, method = "custom", seq_type = 'dna',
               col_scheme = cs)
min35_logo_all <- p + theme_classic() + labs(y = "Relative Entropy") + ylim(0, 1.0)

ggsave(filename = paste('fig2_inCellulo_vs_CFG/extData_fig4d_relEntropy_',sample_name,
                        '_min35logoAll.png',sep=''),
       plot = min35_logo_all,
       width = 300,
       height = 100,
       units = 'mm',
       dpi = 300)

## Plotting spacer lengths -----------------------------------------------------

no_zero_spacer <- top_promoters %>%
  filter(SpacerLength != 0) 

percent_with_spacer <- (nrow(no_zero_spacer) / nrow(top_promoters)) * 100

spacer_hist <- ggplot(no_zero_spacer, aes(x = SpacerLength)) +
  geom_histogram(bins = 6, fill = "grey",
                 color = "black") +
  theme_minimal() +
  scale_x_continuous(n.breaks = 6)

ggsave(filename = paste('fig2_inCellulo_vs_CFG/extData_fig4d_',
                        sample_name,
                        '_spacerHist.png',sep=''),
       plot = spacer_hist,
       width = 300,
       height = 100,
       units = 'mm',
       dpi = 300)

## End -------------------------------------------------------------------------

