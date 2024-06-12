#### Comparing Mtb, Bsu, and Eco intrinsic terminators
#### Authors: Ruby Froom, Mike Wolfe (relative entropy calculations)
#### Date: June 11, 2024

library(seqinr)
library(stringr)
library(Biostrings)
library(readr)
library(textshape)
library(ggseqlogo)
library(tidyr)

# Y-axis maximum for sequence logos
sequence_ymax <- 1.7

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
  column_names <- c(colnames(coordinate_DF),'NT_strand',
                    'Template strand',"IndexingCoordinate",'TSS_local_region')
  
  # Initializes the output dataframe with the correct dimensions
  new_DF <- data.frame(matrix(ncol=length(column_names), nrow=nrow(coordinate_DF)))
  
  # Fills in the new dataframe with the coordinates and direction of transcription
  new_DF[,1:ncol(coordinate_DF)] <- coordinate_DF
  
  # Names the output dataframe columns using the vector column_names
  colnames(new_DF) <- column_names
  
  new_DF$IndexingCoordinate <- new_DF$coordinate
  new_DF[(new_DF$coordinate < upstream_end) & (new_DF$direction == "+"),"IndexingCoordinate"] <- new_DF[(new_DF$coordinate < upstream_end) & (new_DF$direction == "+"),"coordinate"] + length(genome_seq)
  new_DF[(new_DF$coordinate < downstream_end) & (new_DF$direction == "-"),"IndexingCoordinate"] <- new_DF[(new_DF$coordinate < downstream_end) & (new_DF$direction == "-"),"coordinate"] + length(genome_seq)
  
  double_genome <- paste0(c(genome_seq, genome_seq),collapse='')
  
  sense_DF <- new_DF[new_DF$direction == "+",]
  antisense_DF <- new_DF[new_DF$direction == "-",]
  
  # Initializes the range of the genome to extract
  sense_DF$start <- sense_DF$IndexingCoordinate - upstream_end
  sense_DF$end <- sense_DF$IndexingCoordinate + downstream_end
  sense_DF$IndexCoordMin1 <- sense_DF$IndexingCoordinate - 1
  sense_DF$IndexCoordPlus1 <- sense_DF$IndexingCoordinate + 1
  
  # Breaks up the extracted region to capitalize the TSS
  # Promoter is the area directly upstream of the coordinate
  sense_promoters <- str_sub(double_genome, sense_DF$start, sense_DF$IndexCoordMin1)
  sense_TSS <- toupper(str_sub(double_genome, sense_DF$IndexingCoordinate, sense_DF$IndexingCoordinate))
  sense_downstream <- str_sub(double_genome, sense_DF$IndexCoordPlus1, sense_DF$end)
  sense_DF$NT_strand <- paste(sense_promoters, sense_TSS, sense_downstream,sep='')
  
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
  sense_DF$`Template strand` <- sense_extracted_regions_T
  
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
  antisense_downstream <- tolower(str_sub(antisense_DNAseq_NT, upstream_end + 2, upstream_end + downstream_end+1))
  antisense_DF$NT_strand <- paste(antisense_promoter,antisense_TSS,antisense_downstream,sep='')
  
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
  antisense_DF$`Template strand` <- antisense_extracted_regions_T
  
  complete_DF <- rbind(sense_DF, antisense_DF)
  complete_DF_ordered <- complete_DF[order(complete_DF$coordinate),]
  final_DF <- complete_DF_ordered[,1:(ncol(coordinate_DF)+4)]
  
  return (final_DF)
}

## Bsu terminators (comprehensive atlas paper) ---------------------------------

# Loading in genome
Bsu_genome <- read.fasta('genome_files_misc/Bsu_genome.fasta', seqtype="DNA")[[1]]

# Read in DF output of motif_annotated
Bsu_annotated_all <- read.csv('3enrich_NusAG/termination_figures/terminator_annotations/Bsu_terms_annotated.csv')
Bsu_select <- Bsu_annotated_all[sample(nrow(Bsu_annotated_all), 86,
                                      replace = TRUE),]
Bsu_coordinateDF_ordered <- Bsu_select[order(Bsu_select$coordinate),]

# Extracting sequences surrounding the TSS coordinate
Bsu_seqs <- extract_coordinate_environment(Bsu_genome,
                                            Bsu_coordinateDF_ordered,
                                            51,
                                            38)

# Load output from fasta-get-markov (background model)
prior <- read_delim('genome_files_misc/Bsu_background.tsv',
                    col_names = c("value", "prior"),
                    comment = "#",
                    delim = " ")

relative_entropy_calc <- function(input_seqs, csv_name) {
  write.csv(toupper(na.omit(input_seqs)),
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

cs = make_col_scheme(chars=c('A', 'T', 'C', 'G'),
                     cols=c('#C70000','#008000',
                            '#0000C7','#FFAE00'))

bsu_relEntropy <- relative_entropy_calc(Bsu_seqs$NT_strand,
                                              '3enrich_NusAG/termination_figures/extData_fig9_bsu_terminators_output')

p <- ggseqlogo(bsu_relEntropy, method = "custom", seq_type = 'dna',
               col_scheme = cs)
bsu_terms <- p + theme_classic() + labs(y = "Relative Entropy") + ylim(0, sequence_ymax)
bsu_terms$scales$scales[[1]] <- scale_x_continuous(breaks = seq(1, 90, by = 1), 
                                                      labels= c(seq(-52, -52+89,
                                                                    by=1)))


ggsave(filename = '3enrich_NusAG/termination_figures/extData_fig9a_bsu_term_motif.png',
       plot = bsu_terms,
       width = 400,
       height = 100,
       units = 'mm',
       dpi = 300)

## Eco terminators (Peters et al) ----------------------------------------------

Eco_genome <- read.fasta('genome_files_misc/Eco_genome.fasta', seqtype="DNA")[[1]]

# Read in DF output of motif_annotated
Eco_annotated_all <- read.csv('3enrich_NusAG/termination_figures/terminator_annotations/Eco_terms_annotated.csv')
Eco_select <- Eco_annotated_all[sample(nrow(Eco_annotated_all), 86,
                                       replace = TRUE),]
Eco_coordinateDF_ordered <- Eco_select[order(Eco_select$coordinate),]

Eco_seqs <- extract_coordinate_environment(Eco_genome,
                                           Eco_coordinateDF_ordered,
                                           51,
                                           38)

# Load output from fasta-get-markov (background model)
prior <- read_delim('3enrich_NusAG/termination_figures/genomic_info/Eco_background.tsv',
                    col_names = c("value", "prior"),
                    comment = "#",
                    delim = " ")

relative_entropy_calc <- function(input_seqs, csv_name) {
  write.csv(toupper(na.omit(input_seqs)),
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

cs = make_col_scheme(chars=c('A', 'T', 'C', 'G'),
                     cols=c('#C70000','#008000',
                            '#0000C7','#FFAE00'))

eco_relEntropy <- relative_entropy_calc(Eco_seqs$NT_strand,
                                        '3enrich_NusAG/termination_figures/extData_fig9a_eco_terminators_output')

p <- ggseqlogo(eco_relEntropy, method = "custom", seq_type = 'dna',
               col_scheme = cs)
eco_terms <- p + theme_classic() + labs(y = "Relative Entropy") + ylim(0, sequence_ymax)
eco_terms$scales$scales[[1]] <- scale_x_continuous(breaks = seq(1, 90, by = 1), 
                                                   labels= c(seq(-52, -52+89,
                                                                 by=1)))


ggsave(filename = '3enrich_NusAG/termination_figures/extData_fig9a_eco_term_motif.png',
       plot = eco_terms,
       width = 400,
       height = 100,
       units = 'mm',
       dpi = 300)

## Mtb intergenic terminators --------------------------------------------------

intergenic_allCols <- read.csv('3enrich_NusAG/termination_figures/terminator_annotations/Mtb_terms_annotated.csv')

intergenic <- intergenic_allCols %>%
  filter(cluster != 'random') %>%
  filter(GeneLocation == 'inter')

Mtb_terms <- intergenic[,c("coordinate","strand")]
Mtb_terms$direction <- Mtb_terms$strand
Mtb_intergenic <- Mtb_terms[,c("coordinate","direction")]

Mtb_genome <- read.fasta('genome_files_misc/Eco_Mtb_genome.fasta', seqtype="DNA")[[1]]

Mtb_seqs <- extract_coordinate_environment(Mtb_genome,
                                           Mtb_intergenic,
                                           51,
                                           38)

# Load output from fasta-get-markov (background model)
prior <- read_delim('3enrich_NusAG/termination_figures/genomic_info/Mtb_background.tsv',
                    col_names = c("value", "prior"),
                    comment = "#",
                    delim = " ")
prior[prior$value == 'U','value'] <- 'T'


relative_entropy_calc <- function(input_seqs, csv_name) {
  write.csv(toupper(na.omit(input_seqs)),
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


Mtb_relEntropy <- relative_entropy_calc(toupper(intergenic$Nontemplate_strand),
                                        '3enrich_NusAG/termination_figures/extData_fig9a_mtb_terminators_output')

p <- ggseqlogo(Mtb_relEntropy, method = "custom", seq_type = 'dna',
               col_scheme = cs)
Mtb_terms <- p + theme_classic() + labs(y = "Relative Entropy") + ylim(0, sequence_ymax)
Mtb_terms$scales$scales[[1]] <- scale_x_continuous(breaks = seq(1, 90, by = 1), 
                                                   labels= c(seq(-52, -52+89,
                                                                 by=1)))

ggsave(filename = '3enrich_NusAG/termination_figures/extData_fig9a_Mtb_term_motif.png',
       plot = Mtb_terms,
       width = 400,
       height = 100,
       units = 'mm',
       dpi = 300)

## Compare hairpin features ----------------------------------------------------

Mtb_seqs <- intergenic

Bsu_seqs$species <- 'Bsu'
Eco_seqs$species <- 'Eco'
Mtb_seqs$species <- 'Mtb'

Bsu_seqs$AT_count_downstream <- Bsu_seqs$U_count_downstreamRegion + Bsu_seqs$A_count_downstreamRegion
Eco_seqs$AT_count_downstream <- Eco_seqs$U_count_downstreamRegion + Eco_seqs$A_count_downstreamRegion

colsToKeep <- c("species","AT_count_downstream","U_count_tractRegion","U_count",
                "RNAfold_hairpin_region_44_deltaG","RNAfold_hairpin_region_44_spacerLength",
                "RNAfold_hairpin_region_44_stemLength","RNAfold_hairpin_region_44_loopLength")

bsu <- Bsu_seqs[,colsToKeep]
eco <- Eco_seqs[,colsToKeep]
mtb <- Mtb_seqs[,colsToKeep]

species_terms <- rbind(bsu, eco, mtb)

## Distribution plots ----------------------------------------------------------

library(dunn.test)
library(ggplot2)

Mtb_col <- '#CBADC1'
Eco_col <- '#9eb5a3'
Bsu_col <- '#c47d72'

species_downstreamATs <- ggplot(species_terms, aes(x = species,
                                                y = AT_count_downstream,
                                                fill = species)) + 
  geom_violin(linewidth = 2) +
  geom_boxplot(width = 0.2,
               linewidth = 2) +
  scale_y_continuous(limits = c(0,10),
                     breaks = seq(0,10,2)) +
  scale_fill_manual(values = c(Bsu_col, Eco_col, Mtb_col)) +
  theme_minimal() +
  theme(legend.position = 'none')
ggsave(filename = '3enrich_NusAG/termination_figures/extData_fig9b_downstreamATs_species_boxplot.png',
       plot = species_downstreamATs,
       width = 200,
       height = 200,
       units = 'mm',
       dpi = 300)
padj_method = 'BH'
kruskal.test(AT_count_downstream ~ species,
             species_terms)
dsDNA_padj <- dunn.test(species_terms$AT_count_downstream,
                        species_terms$species,
                        method = padj_method)

# Hairpin deltaG differences
species_hairpinDeltaG <- ggplot(species_terms, aes(x = species,
                          y = RNAfold_hairpin_region_44_deltaG,
                          fill = species)) + 
  geom_violin(linewidth = 2) +
  geom_boxplot(width = 0.2,
               linewidth = 2) +
  scale_y_continuous(limits = c(-40,0),
                     breaks = seq(-40,0,10)) +
  scale_fill_manual(values = c(Bsu_col, Eco_col, Mtb_col)) +
  theme_minimal() +
  theme(legend.position = 'none')
ggsave(filename = '3enrich_NusAG/termination_figures/extData_fig9b_deltaG_species_boxplot.png',
       plot = species_hairpinDeltaG,
       width = 200,
       height = 200,
       units = 'mm',
       dpi = 300)
padj_method = 'BH'
kruskal.test(RNAfold_hairpin_region_44_deltaG ~ species,
             species_terms)
deltaG_padj <- dunn.test(species_terms$RNAfold_hairpin_region_44_deltaG,
                        species_terms$species,
                        method = padj_method)


padj_method = 'BH'
kruskal.test(RNAfold_hairpin_region_44_spacerLength ~ species,
             species_terms)
spacer_padj <- dunn.test(species_terms$RNAfold_hairpin_region_44_spacerLength,
                          species_terms$species,
                          method = padj_method)


## Hairpin spacer lengths normalized
spacerLength_percentage_DF <- data.frame(spacer_length = seq(0,10,1),
                              Mtb_spacer_percentage = c(nrow(species_terms[(species_terms$RNAfold_hairpin_region_44_spacerLength == 0) &
                                                                                         (species_terms$species == 'Mtb'),]) /
                                                                nrow(species_terms[(species_terms$species == 'Mtb'),]),
                                                              nrow(species_terms[(species_terms$RNAfold_hairpin_region_44_spacerLength == 1) &
                                                                                         (species_terms$species == 'Mtb'),]) /
                                                                nrow(species_terms[(species_terms$species == 'Mtb'),]),
                                                              nrow(species_terms[(species_terms$RNAfold_hairpin_region_44_spacerLength == 2) &
                                                                                         (species_terms$species == 'Mtb'),]) /
                                                                nrow(species_terms[(species_terms$species == 'Mtb'),]),
                                                              nrow(species_terms[(species_terms$RNAfold_hairpin_region_44_spacerLength == 3) &
                                                                                         (species_terms$species == 'Mtb'),]) /
                                                                nrow(species_terms[(species_terms$species == 'Mtb'),]),
                                                              nrow(species_terms[(species_terms$RNAfold_hairpin_region_44_spacerLength == 4) &
                                                                                         (species_terms$species == 'Mtb'),]) /
                                                                nrow(species_terms[(species_terms$species == 'Mtb'),]),
                                                              nrow(species_terms[(species_terms$RNAfold_hairpin_region_44_spacerLength == 5) &
                                                                                         (species_terms$species == 'Mtb'),]) /
                                                                nrow(species_terms[(species_terms$species == 'Mtb'),]),
                                                              nrow(species_terms[(species_terms$RNAfold_hairpin_region_44_spacerLength == 6) &
                                                                                         (species_terms$species == 'Mtb'),]) /
                                                                nrow(species_terms[(species_terms$species == 'Mtb'),]),
                                                              nrow(species_terms[(species_terms$RNAfold_hairpin_region_44_spacerLength == 7) &
                                                                                         (species_terms$species == 'Mtb'),]) /
                                                                nrow(species_terms[(species_terms$species == 'Mtb'),]),
                                                              nrow(species_terms[(species_terms$RNAfold_hairpin_region_44_spacerLength == 8) &
                                                                                         (species_terms$species == 'Mtb'),]) /
                                                                nrow(species_terms[(species_terms$species == 'Mtb'),]),
                                                              nrow(species_terms[(species_terms$RNAfold_hairpin_region_44_spacerLength == 9) &
                                                                                         (species_terms$species == 'Mtb'),]) /
                                                                nrow(species_terms[(species_terms$species == 'Mtb'),]),
                                                              nrow(species_terms[(species_terms$RNAfold_hairpin_region_44_spacerLength == 10) &
                                                                                         (species_terms$species == 'Mtb'),]) /
                                                                nrow(species_terms[(species_terms$species == 'Mtb'),])),
                              Bsu_spacer_percentage = c(nrow(species_terms[(species_terms$RNAfold_hairpin_region_44_spacerLength == 0) &
                                                                             (species_terms$species == 'Bsu'),]) /
                                                          nrow(species_terms[(species_terms$species == 'Bsu'),]),
                                                        nrow(species_terms[(species_terms$RNAfold_hairpin_region_44_spacerLength == 1) &
                                                                             (species_terms$species == 'Bsu'),]) /
                                                          nrow(species_terms[(species_terms$species == 'Bsu'),]),
                                                        nrow(species_terms[(species_terms$RNAfold_hairpin_region_44_spacerLength == 2) &
                                                                             (species_terms$species == 'Bsu'),]) /
                                                          nrow(species_terms[(species_terms$species == 'Bsu'),]),
                                                        nrow(species_terms[(species_terms$RNAfold_hairpin_region_44_spacerLength == 3) &
                                                                             (species_terms$species == 'Bsu'),]) /
                                                          nrow(species_terms[(species_terms$species == 'Bsu'),]),
                                                        nrow(species_terms[(species_terms$RNAfold_hairpin_region_44_spacerLength == 4) &
                                                                             (species_terms$species == 'Bsu'),]) /
                                                          nrow(species_terms[(species_terms$species == 'Bsu'),]),
                                                        nrow(species_terms[(species_terms$RNAfold_hairpin_region_44_spacerLength == 5) &
                                                                             (species_terms$species == 'Bsu'),]) /
                                                          nrow(species_terms[(species_terms$species == 'Bsu'),]),
                                                        nrow(species_terms[(species_terms$RNAfold_hairpin_region_44_spacerLength == 6) &
                                                                             (species_terms$species == 'Bsu'),]) /
                                                          nrow(species_terms[(species_terms$species == 'Bsu'),]),
                                                        nrow(species_terms[(species_terms$RNAfold_hairpin_region_44_spacerLength == 7) &
                                                                             (species_terms$species == 'Bsu'),]) /
                                                          nrow(species_terms[(species_terms$species == 'Bsu'),]),
                                                        nrow(species_terms[(species_terms$RNAfold_hairpin_region_44_spacerLength == 8) &
                                                                             (species_terms$species == 'Bsu'),]) /
                                                          nrow(species_terms[(species_terms$species == 'Bsu'),]),
                                                        nrow(species_terms[(species_terms$RNAfold_hairpin_region_44_spacerLength == 9) &
                                                                             (species_terms$species == 'Bsu'),]) /
                                                          nrow(species_terms[(species_terms$species == 'Bsu'),]),
                                                        nrow(species_terms[(species_terms$RNAfold_hairpin_region_44_spacerLength == 10) &
                                                                             (species_terms$species == 'Bsu'),]) /
                                                          nrow(species_terms[(species_terms$species == 'Bsu'),])),
                              Eco_spacer_percentage = c(nrow(species_terms[(species_terms$RNAfold_hairpin_region_44_spacerLength == 0) &
                                                                             (species_terms$species == 'Eco'),]) /
                                                          nrow(species_terms[(species_terms$species == 'Eco'),]),
                                                        nrow(species_terms[(species_terms$RNAfold_hairpin_region_44_spacerLength == 1) &
                                                                             (species_terms$species == 'Eco'),]) /
                                                          nrow(species_terms[(species_terms$species == 'Eco'),]),
                                                        nrow(species_terms[(species_terms$RNAfold_hairpin_region_44_spacerLength == 2) &
                                                                             (species_terms$species == 'Eco'),]) /
                                                          nrow(species_terms[(species_terms$species == 'Eco'),]),
                                                        nrow(species_terms[(species_terms$RNAfold_hairpin_region_44_spacerLength == 3) &
                                                                             (species_terms$species == 'Eco'),]) /
                                                          nrow(species_terms[(species_terms$species == 'Eco'),]),
                                                        nrow(species_terms[(species_terms$RNAfold_hairpin_region_44_spacerLength == 4) &
                                                                             (species_terms$species == 'Eco'),]) /
                                                          nrow(species_terms[(species_terms$species == 'Eco'),]),
                                                        nrow(species_terms[(species_terms$RNAfold_hairpin_region_44_spacerLength == 5) &
                                                                             (species_terms$species == 'Eco'),]) /
                                                          nrow(species_terms[(species_terms$species == 'Eco'),]),
                                                        nrow(species_terms[(species_terms$RNAfold_hairpin_region_44_spacerLength == 6) &
                                                                             (species_terms$species == 'Eco'),]) /
                                                          nrow(species_terms[(species_terms$species == 'Eco'),]),
                                                        nrow(species_terms[(species_terms$RNAfold_hairpin_region_44_spacerLength == 7) &
                                                                             (species_terms$species == 'Eco'),]) /
                                                          nrow(species_terms[(species_terms$species == 'Eco'),]),
                                                        nrow(species_terms[(species_terms$RNAfold_hairpin_region_44_spacerLength == 8) &
                                                                             (species_terms$species == 'Eco'),]) /
                                                          nrow(species_terms[(species_terms$species == 'Eco'),]),
                                                        nrow(species_terms[(species_terms$RNAfold_hairpin_region_44_spacerLength == 9) &
                                                                             (species_terms$species == 'Eco'),]) /
                                                          nrow(species_terms[(species_terms$species == 'Eco'),]),
                                                        nrow(species_terms[(species_terms$RNAfold_hairpin_region_44_spacerLength == 10) &
                                                                             (species_terms$species == 'Eco'),]) /
                                                          nrow(species_terms[(species_terms$species == 'Eco'),])))



spacerLength_species_percentage <- data.frame(spacerLength = rep(spacerLength_percentage_DF$spacer_length, 3),
                                          condition = c(rep('Mtb',11),
                                                        rep('Bsu',11),
                                                        rep('Eco',11)),
                                          spacer_percentage = c(spacerLength_percentage_DF$Mtb_spacer_percentage,
                                                                spacerLength_percentage_DF$Bsu_spacer_percentage,
                                                                spacerLength_percentage_DF$Eco_spacer_percentage))
spacerLength_species_percentage$spacerLength_min <- -(spacerLength_species_percentage$spacerLength)

hairpin_spacer_hist_species <- ggplot(spacerLength_species_percentage,
                                   aes(x = spacerLength_min,
                                       y = spacer_percentage,
                                       fill = condition)) +
  geom_col(position = 'dodge',
           alpha = 0.8, width = 0.8,
           color = 'black') +
  scale_x_continuous(limits = c(-11,1),
                     breaks = seq(-10, 0, 1)) +
  scale_fill_manual(values = c(Bsu_col, Eco_col, Mtb_col))+
  theme(legend.position = "none")+
  theme_minimal()
ggsave(filename = '3enrich_NusAG/termination_figures/extData_fig9b_spacerLength_hist_species.png',
       plot = hairpin_spacer_hist_species,
       width = 400,
       height = 300,
       units = 'mm',
       dpi = 300)


## Hairpin stem length differences
species_stemLength <- ggplot(species_terms, aes(x = species,
                                                   y = RNAfold_hairpin_region_44_stemLength,
                                                   fill = species)) + 
  geom_violin(linewidth = 2) +
  geom_boxplot(width = 0.2,
               linewidth = 2) +
  scale_y_continuous(limits = c(0,25),
                     breaks = seq(0,25,5)) +
  scale_fill_manual(values = c(Bsu_col, Eco_col, Mtb_col)) +
  theme_minimal() +
  theme(legend.position = 'none')
ggsave(filename = '3enrich_NusAG/termination_figures/extData_fig9b_stemLength_species_boxplot.png',
       plot = species_stemLength,
       width = 200,
       height = 200,
       units = 'mm',
       dpi = 300)
padj_method = 'BH'
kruskal.test(RNAfold_hairpin_region_44_stemLength ~ species,
             species_terms)
stem_padj <- dunn.test(species_terms$RNAfold_hairpin_region_44_stemLength,
                       species_terms$species,
                       method = padj_method)



utractTs <- ggplot(species_terms,
                   aes(x = U_count_tractRegion)) +
  geom_density(data = subset(species_terms, species == 'Bsu'),
               alpha = 0.5, fill = Bsu_col) +
  geom_vdensity(data = bsu$U_count_tractRegion,
                at = median(bsu$U_count_tractRegion),
                color = Bsu_col,
                linewidth = vline_width) +
  geom_density(data = subset(species_terms, species == 'Eco'),
               alpha = 0.5, fill = Eco_col) +
  geom_vdensity(data = eco$U_count_tractRegion,
                at = median(eco$U_count_tractRegion),
                color = Eco_col,
                linewidth = vline_width) +
  geom_density(data = subset(species_terms, species == 'Mtb'),
               alpha = 0.5, fill = Mtb_col) +
  geom_vdensity(data = mtb$U_count_tractRegion,
                at = median(mtb$U_count_tractRegion),
                color = Mtb_col,
                linewidth = vline_width) +
  theme_minimal() +
  theme(panel.grid.minor = element_blank())
ggsave(filename = '3enrich_NusAG/termination_figures/extData_fig9b_species_utractTs.png',
       plot = utractTs,
       width = 144,
       height = 100,
       units = 'mm',
       dpi = 300)
padj_method = 'BH'
kruskal.test(U_count_tractRegion ~ species,
             species_terms)
Utract_padj <- dunn.test(species_terms$U_count_tractRegion,
                         species_terms$species,
                         method = padj_method)


ggplot(species_terms,
       aes(x = U_count)) +
  geom_density(data = subset(species_terms, species == 'Bsu'),
               alpha = 0.5, fill = Bsu_col) +
  geom_vdensity(data = bsu$U_count,
                at = median(bsu$U_count),
                color = Bsu_col,
                linewidth = vline_width) +
  geom_density(data = subset(species_terms, species == 'Eco'),
               alpha = 0.5, fill = Eco_col) +
  geom_vdensity(data = eco$U_count,
                at = median(eco$U_count),
                color = Eco_col,
                linewidth = vline_width) +
  geom_density(data = subset(species_terms, species == 'Mtb'),
               alpha = 0.5, fill = Mtb_col) +
  geom_vdensity(data = mtb$U_count,
                at = median(mtb$U_count),
                color = Mtb_col,
                linewidth = vline_width) +
  theme_minimal() +
  theme(panel.grid.minor = element_blank())

# Hairpin loop length differences
hairpin_loopLength <- ggplot(species_terms,
                             aes(x = RNAfold_hairpin_region_44_loopLength)) +
  geom_density(data = subset(species_terms, species == 'Bsu'),
               alpha = 0.5, fill = Bsu_col) +
  geom_vdensity(data = bsu$RNAfold_hairpin_region_44_loopLength,
                at = median(bsu$RNAfold_hairpin_region_44_loopLength),
                color = Bsu_col,
                linewidth = vline_width) +
  geom_density(data = subset(species_terms, species == 'Eco'),
               alpha = 0.5, fill = Eco_col) +
  geom_vdensity(data = eco$RNAfold_hairpin_region_44_loopLength,
                at = median(eco$RNAfold_hairpin_region_44_loopLength),
                color = Eco_col,
                linewidth = vline_width) +
  geom_density(data = subset(species_terms, species == 'Mtb'),
               alpha = 0.5, fill = Mtb_col) +
  geom_vdensity(data = mtb$RNAfold_hairpin_region_44_loopLength,
                at = median(mtb$RNAfold_hairpin_region_44_loopLength),
                color = Mtb_col,
                linewidth = vline_width) +
  theme_minimal() +
  theme(panel.grid.minor = element_blank())
ggsave(filename = '3enrich_NusAG/termination_figures/extData_fig9b_species_loopLength.png',
       plot = hairpin_loopLength,
       width = 144,
       height = 100,
       units = 'mm',
       dpi = 300)
padj_method = 'BH'
kruskal.test(RNAfold_hairpin_region_44_loopLength ~ species,
             species_terms)
loop_padj <- dunn.test(species_terms$RNAfold_hairpin_region_44_loopLength,
                         species_terms$species,
                         method = padj_method)

posthoc_DF <- data.frame(comparisons = dsDNA_padj$comparisons,
                         dsDNA_padj = dsDNA_padj$P.adjusted,
                         dsDNA_sig = dsDNA_padj$P.adjusted <= 0.05,
                         deltaG_padj = deltaG_padj$P.adjusted,
                         deltaG_sig = deltaG_padj$P.adjusted <= 0.05,
                         spacer_padj = spacer_padj$P.adjusted,
                         spacer_sig = spacer_padj$P.adjusted <= 0.05,
                         loop_padj = loop_padj$P.adjusted,
                         loop_sig = loop_padj$P.adjusted <= 0.05,
                         stem_padj = stem_padj$P.adjusted,
                         stem_sig = stem_padj$P.adjusted <= 0.05)

write.csv(posthoc_DF, '3enrich_NusAG/termination_figures/extData_fig9b_species_posthoc_DF.csv')
