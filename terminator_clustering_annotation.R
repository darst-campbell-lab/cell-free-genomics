#### Identifying Nus-regulated terminators from cell-free genomics
#### Authors: Ruby Froom, Tom Carroll (k-means clustering and silhouette tests),
####          Mike Wolfe (relative entropy calculations)
#### Date: February 23, 2024

#### Libraries ----------------------------------------------------------------

library(dplyr)
library(DESeq2)
library(pheatmap)
library(NbClust)
library(viridis)
library(stringr)
library(stringi)
library(Biostrings)
library(seqinr)
library(universalmotif)
library(EnhancedVolcano)
library(ggseqlogo)
library(GenomicRanges)
library(seqinr)
library(Biostrings)
library(sqldf)

#### Prepare data and run LRT with reduced model using DESeq2 ------------------

# Adjust input csv files and formula as needed
counts <- '70'
sig_value <- 0.05
DESeq2_path <- 'selectThreshold/DESeq2'
design_matrix_string = paste(DESeq2_path,
                             'multifactor_LRT_int_design_matrix.csv',
                             sep='')
count_data_string = paste(DESeq2_path,
                          'multifactor_countdata_ends_',
                          counts,'counts.csv',
                          sep ='')
formula_string= '~condition'

# Read in design matrix
design_matrix = as.data.frame(read.table(design_matrix_string, sep=",",
                                         header=TRUE,colClasses = "factor"))
conditions <- levels(design_matrix[['condition']])

# Read in count data
count_data_all <- read.csv(count_data_string)

# Remove first extraneous column from count data
counts_only <- count_data_all[,c(2:ncol(count_data_all))]

# Add in the coordinate ID (location and strand) as row names
rownames(counts_only) <- count_data_all[,1]

# Generate DESeqDataSet
deseq_dataset <- DESeqDataSetFromMatrix(countData = counts_only,
                                        colData = design_matrix,
                                        design = formula(formula_string))

# The coordIDs with 'Eco_' correspond to genes from spike RNA (E. coli RNA)
# Isolate those rows
num_spike_coords <- nrow(counts_only %>%
                           filter(grepl("Eco_",rownames(counts_only))))

# Estimate size factors based on spike counts
# Note: the total amounts of RNA didn't happen to differ by very much,
# so the size factor differences aren't dramatic
spike_row_indices <- ((nrow(counts_only)-num_spike_coords)+1):(nrow(counts_only))
dds_sizeFactors <- estimateSizeFactors(deseq_dataset,
                                       controlGenes = spike_row_indices)

# Run LRT
dds <- DESeq(dds_sizeFactors, parallel=TRUE, test = 'LRT',
             reduced = ~1, fitType = 'local')
vsd <- varianceStabilizingTransformation(dds, fitType = 'local')

# Visually inspect model fits 
png("3enrich_NusAG/termination_figures/extdata_fig7a_3end_disp_ests.png",
    res = 300,
    width = 5,
    height = 5,
    units = 'in')
plotDispEsts(dds)
dev.off()

#### Prepare Nus-regulated hits for clustering ---------------------------------

# Isolate results of LRT test
across_groups_all <- results(dds)
across_groups_ordered <- across_groups_all[order(across_groups_all$pvalue),]

# Filter out any Eco coordinates (since they correspond to spike RNA sample)
across_groups_eco <- as.data.frame(across_groups_ordered) %>%
  filter(grepl("Eco_",rownames(as.data.frame(across_groups_ordered))))
across_groups_noEco <- across_groups_all[!(rownames(across_groups_all) %in%
                                             rownames(across_groups_eco)),]

# Isolate coordinates with at least 1 significant difference between means
# (Putative Nus-regulated hits)
sig_changes <- rownames(across_groups_noEco)[across_groups_noEco$padj <= sig_value &
                                               !is.na(across_groups_noEco$padj)]
normTF <- normTransform(dds)
rlogTF_all <- rlog(dds, fitType = 'local')
rlogTF <- rlogTF_all[!str_detect(rownames(rlogTF_all), 'Eco_'),]
rlog_matrix <- assay(rlogTF)

# Isolate Nus-regulated hits from transformed data
sigMat <- rlog_matrix[rownames(rlog_matrix) %in% sig_changes,]

#### Sequential rounds of k-means clustering -----------------------------------
# Note 1: here, I only proceed with sub-clustering if the data cleanly separate

# Note 2: In comments/variable names, I describe the cluster "phenotypes" of 
# Nus regulation based on visual impression of z-scores, 
# but see Fig. 3D for actual log2foldchanges associated with each cluster

## Cluster 1
set.seed(153)
cluster_number <- 3
k1 <- pheatmap(sigMat, scale = "row", kmeans_k = cluster_number)
clusterDF1 <- as.data.frame(factor(k1$kmeans$cluster))
colnames(clusterDF1) <- "Cluster"
clusterDF1[1:10, , drop = FALSE]
OrderByCluster1 <- sigMat[order(clusterDF1$Cluster),]

# Perform silhouette test to optimize cluster number
rowScaledMat1 <- t(scale(t(sigMat)))
clusterNum1 <- NbClust(rowScaledMat1, distance = "euclidean", min.nc = 2,
                       max.nc = 12,
                       method = "kmeans", index = "silhouette")
clusterNum1$Best.nc
orderedCluster1 <- sort(clusterNum1$Best.partition)

sigMat1 <- sigMat[match(names(orderedCluster1), rownames(sigMat)), ]

ColorCode  = list(cluster_nums = c("#5ec962",
                                   "#000000",
                                   "#FFA904"))

ColorCode  = list(cluster_nums = c("#5ec962",
                                   "#000000"))

first_heatmap <- pheatmap(sigMat1, scale = "row",
                          annotation_row = clusterDF1, show_rownames = FALSE,
                          cluster_rows = FALSE, border_color = "NA",
                          annotation_color = ColorCode)

ggsave(filename = '3enrich_NusAG/termination_figures/extData_fig7b_heatmap1_clusterAB_CD.png',
       plot = first_heatmap,
       width = 500,
       height = 500,
       units = 'mm',
       dpi = 300)

clusterB <- names(orderedCluster1[orderedCluster1 == 1])
clusterA <- names(orderedCluster1[orderedCluster1 == 3])

# Cluster C+D split
remove_cluster <- orderedCluster1[orderedCluster1 == 1]
remove_cluster2 <- orderedCluster1[orderedCluster1 == 3]

sigMat2 <- sigMat[!(rownames(sigMat) %in% names(remove_cluster)),]
sigMat2 <- sigMat2[!(rownames(sigMat2) %in% names(remove_cluster2)),]

pheatmap(sigMat2, scale = "row", show_rownames = FALSE)
set.seed(153)
cluster_number <- 2
k2 <- pheatmap(sigMat2, scale = "row", kmeans_k = cluster_number)
clusterDF2 <- as.data.frame(factor(k2$kmeans$cluster))
colnames(clusterDF2) <- "Cluster"
clusterDF2[1:10, , drop = FALSE]
OrderByCluster2 <- sigMat2[order(clusterDF2$Cluster),]
pheatmap(OrderByCluster2, scale = "row",
         annotation_row = clusterDF2, show_rownames = FALSE,
         cluster_rows = FALSE, border_color = "NA")

rowScaledMat2 <- t(scale(t(sigMat2)))
clusterNum2 <- NbClust(rowScaledMat2, distance = "euclidean", min.nc = 2, max.nc = 12,
                       method = "kmeans", index = "silhouette")
clusterNum2$Best.nc
orderedCluster2 <- sort(clusterNum2$Best.partition)

sigMat2 <- sigMat2[match(names(orderedCluster2), rownames(sigMat2)), ]

second_heatmap <- pheatmap(sigMat2, scale = "row", annotation_row = clusterDF2, show_rownames = FALSE,
         cluster_rows = FALSE, border_color = "NA")

ggsave(filename = '3enrich_NusAG/termination_figures/extData_fig7b_heatmap2_clusterCD.png',
       plot = second_heatmap,
       width = 500,
       height = 500,
       units = 'mm',
       dpi = 300)

clusterC <- names(orderedCluster2[orderedCluster2 == 1])
clusterD <- names(orderedCluster2[orderedCluster2 == 2])

#### Final heatmap -------------------------------------------------------------

row_order <- c(clusterA,
               clusterB,
               clusterC,
               clusterD)

order_Matrix <- as.data.frame(sigMat[match(row_order, rownames(sigMat)),])

cluster_nums <- c(rep(1,length(clusterA)),
                  rep(2,length(clusterB)),
                  rep(3,length(clusterC)),
                  rep(4,length(clusterD))
)


sub_anno <- as.data.frame(x = cluster_nums,
                          row.names = row_order)
as.factor(sub_anno$cluster_nums)

ColorCode  = list(cluster_nums = c("#FFA904",
                                   "#5ec962",
                                   "#21918c",
                                   "#3b528b"))

final_heatmap <- pheatmap(order_Matrix, scale = "row",
                          cluster_rows = FALSE, annotation_row = sub_anno,
                          annotation_color = ColorCode, border_color = "NA",
                          show_rownames = FALSE)

ggsave(filename = '3enrich_NusAG/termination_figures/fig4c_heatmap_final.png',
       plot = final_heatmap,
       width = 500,
       height = 500,
       units = 'mm',
       dpi = 300)

## NusA, NusG, and NusAG volcano plots -----------------------------------------

dds_Wald <- DESeq(dds_sizeFactors, parallel=TRUE, test = 'Wald',
             fitType = 'local')
vsd_Wald <- varianceStabilizingTransformation(dds_Wald, fitType = 'local')
res_NusG <- results(dds_Wald,contrast=c("condition","NusG","noTF"))
res_NusA <- results(dds_Wald,contrast=c("condition","NusA","noTF"))
res_NusAG <- results(dds_Wald,contrast=c("condition","NusA_NusG","noTF"))

NusG_eco <- as.data.frame(res_NusG) %>%
  filter(grepl("Eco_",rownames(as.data.frame(res_NusG))))
NusG_noEco <- res_NusG[!(rownames(res_NusG) %in%
                           rownames(NusG_eco)),]

NusA_eco <- as.data.frame(res_NusA) %>%
  filter(grepl("Eco_",rownames(as.data.frame(res_NusA))))
NusA_noEco <- res_NusA[!(rownames(res_NusA) %in%
                           rownames(NusA_eco)),]


NusAG_eco <- as.data.frame(res_NusAG) %>%
  filter(grepl("Eco_",rownames(as.data.frame(res_NusAG))))
NusAG_noEco <- res_NusAG[!(rownames(res_NusAG) %in%
                             rownames(NusAG_eco)),]

coords_to_highlight <- c('7196579-',
                         '8366386+',
                         '8366385+',
                         '6118706-',
                         '6118708-')

NusA_noEco$category <- 1
NusA_noEco[rownames(NusA_noEco) %in% coords_to_highlight,
           'category'] <- 2
NusA_noEco <- data.frame(NusA_noEco) %>%
  arrange(category)

nrow(NusA_noEco[(NusA_noEco$padj <= 0.05) & (NusA_noEco$log2FoldChange <= -1),])
nrow(NusA_noEco[(NusA_noEco$padj <= 0.05) & (NusA_noEco$log2FoldChange >= 1),])
nrow(NusA_noEco[(NusA_noEco$padj > 0.05) | (NusA_noEco$log2FoldChange < 1),])

keyvals_NusA <- ifelse(
  (rownames(NusA_noEco) == coords_to_highlight[1]), 'green',
  ifelse(
    (rownames(NusA_noEco) == coords_to_highlight[2]), 'blue',
    ifelse(
      (rownames(NusA_noEco) == coords_to_highlight[3]), 'black',
      ifelse(
        (rownames(NusA_noEco) == coords_to_highlight[4]) | (rownames(NusA_noEco) == coords_to_highlight[5]), 'red',
        ifelse(
          (!(rownames(NusA_noEco) %in% coords_to_highlight) & (NusA_noEco$padj < 0.05) & (abs(NusA_noEco$log2FoldChange) > 1)), '#6C98E6',
          'grey')))))
keyvals_NusA[is.na(keyvals_NusA)] <- 'grey'
names(keyvals_NusA)[keyvals_NusA == '#6C98E6'] <- 'NusA_affected'
names(keyvals_NusA)[keyvals_NusA == 'black'] <- 'icd1_coord2'
names(keyvals_NusA)[keyvals_NusA == 'green'] <- 'pitB'
names(keyvals_NusA)[keyvals_NusA == 'blue'] <- 'icd1_coord1'
names(keyvals_NusA)[keyvals_NusA == 'red'] <- 'rrf'

NusA_volcano <- EnhancedVolcano(NusA_noEco,
                                lab = rownames(NusA_noEco),
                                labSize = 0,
                                pointSize = 3,
                                x = 'log2FoldChange',
                                y = 'padj',
                                pCutoff = 0.05,
                                colCustom = keyvals_NusA,
                                legendPosition = 0,
                                colAlpha = 1,
                                xlim = c(-2,3.2),
                                ylim = c(0,32))

ggsave(filename = '3enrich_NusAG/termination_figures/extData_fig7d_NusA_volcano_Wald.png',
       plot = NusA_volcano,
       width = 300,
       height = 200,
       units = 'mm',
       dpi = 300)

NusG_noEco$category <- 1
NusG_noEco[rownames(NusG_noEco) %in% coords_to_highlight,
           'category'] <- 2
NusG_noEco <- data.frame(NusG_noEco) %>%
  arrange(category)

keyvals_NusG <- ifelse(
  (rownames(NusG_noEco) == coords_to_highlight[1]), 'green',
  ifelse(
    (rownames(NusG_noEco) == coords_to_highlight[2]), 'blue',
    ifelse(
      (rownames(NusG_noEco) == coords_to_highlight[3]), 'black',
      ifelse(
        (rownames(NusA_noEco) == coords_to_highlight[4]) | (rownames(NusA_noEco) == coords_to_highlight[5]), 'red',
        ifelse(
          (!(rownames(NusG_noEco) %in% coords_to_highlight) & (NusG_noEco$padj < 0.05) & (abs(NusG_noEco$log2FoldChange) > 1)), '#E26B7A',
          'grey')))))
keyvals_NusG[is.na(keyvals_NusG)] <- 'grey'
names(keyvals_NusG)[keyvals_NusG == '#E26B7A'] <- 'NusG_affected'
names(keyvals_NusG)[keyvals_NusG == 'black'] <- 'icd1_coord2'
names(keyvals_NusG)[keyvals_NusG == 'green'] <- 'pitB'
names(keyvals_NusG)[keyvals_NusG == 'blue'] <- 'icd1_coord1'
names(keyvals_NusG)[keyvals_NusG == 'red'] <- 'rrf'

nrow(NusG_noEco[(NusG_noEco$padj <= 0.05) & (NusG_noEco$log2FoldChange <= -1),])
nrow(NusG_noEco[(NusG_noEco$padj <= 0.05) & (NusG_noEco$log2FoldChange >= 1),])
nrow(NusG_noEco[(NusG_noEco$padj > 0.05) | (NusG_noEco$log2FoldChange < 1),])

NusG_volcano <- EnhancedVolcano(NusG_noEco,
                                lab = rownames(NusG_noEco),
                                labSize = 0,
                                pointSize = 3,
                                x = 'log2FoldChange',
                                y = 'padj',
                                pCutoff = 0.05,
                                colCustom = keyvals_NusG,
                                legendPosition = 0,
                                colAlpha = 1,
                                xlim = c(-2,3.2),
                                ylim = c(0,32))

ggsave(filename = '3enrich_NusAG/termination_figures/extData_fig7d_NusG_volcano_Wald.png',
       plot = NusG_volcano,
       width = 300,
       height = 200,
       units = 'mm',
       dpi = 300)

NusAG_noEco$category <- 1
NusAG_noEco[rownames(NusAG_noEco) %in% coords_to_highlight,
            'category'] <- 2
NusAG_noEco <- data.frame(NusAG_noEco) %>%
  arrange(category)

keyvals_NusAG <- ifelse(
  (rownames(NusAG_noEco) == coords_to_highlight[1]), 'green',
  ifelse(
    (rownames(NusAG_noEco) == coords_to_highlight[2]), 'blue',
    ifelse(
      (rownames(NusAG_noEco) == coords_to_highlight[3]), 'black',
      ifelse(
        (rownames(NusA_noEco) == coords_to_highlight[4]) | (rownames(NusA_noEco) == coords_to_highlight[5]), 'red',
        ifelse(
          (!(rownames(NusAG_noEco) %in% coords_to_highlight) & (NusAG_noEco$padj < 0.05) & (abs(NusAG_noEco$log2FoldChange) > 1)), 'purple',
          'grey')))))
keyvals_NusAG[is.na(keyvals_NusAG)] <- 'grey'
names(keyvals_NusAG)[keyvals_NusAG == 'purple'] <- 'NusAG_affected'
names(keyvals_NusAG)[keyvals_NusAG == 'black'] <- 'icd1_coord2'
names(keyvals_NusAG)[keyvals_NusAG == 'green'] <- 'pitB'
names(keyvals_NusAG)[keyvals_NusAG == 'blue'] <- 'icd1_coord1'
names(keyvals_NusAG)[keyvals_NusAG == 'red'] <- 'rrf'

nrow(NusAG_noEco[(NusAG_noEco$padj <= 0.05) & (NusAG_noEco$log2FoldChange <= -1),])
nrow(NusAG_noEco[(NusAG_noEco$padj <= 0.05) & (NusAG_noEco$log2FoldChange >= 1),])
nrow(NusAG_noEco[(NusAG_noEco$padj > 0.05) | (NusAG_noEco$log2FoldChange < 1),])

NusAG_volcano <- EnhancedVolcano(NusAG_noEco,
                                 lab = rownames(NusAG_noEco),
                                 labSize = 0,
                                 pointSize = 3,
                                 x = 'log2FoldChange',
                                 y = 'padj',
                                 pCutoff = 0.05,
                                 colCustom = keyvals_NusAG,
                                 legendPosition = 0,
                                 colAlpha = 1,
                                 xlim = c(-2,3.2),
                                 ylim = c(0,32))

ggsave(filename = '3enrich_NusAG/termination_figures/extData_fig7d_NusAG_volcano_Wald.png',
       plot = NusAG_volcano,
       width = 300,
       height = 200,
       units = 'mm',
       dpi = 300)

## Read in annotation and label with cluster  ----------------------------------
annotated_coords <- read.csv(paste('3enrich_NusAG/termination_figures/terminator_annotations/multifactor_',
                                   counts,
                                   'counts_annotated.csv',sep=''))

annotated_coords$cluster <- 'NE'
annotated_coords[annotated_coords$coordID %in% clusterA,'cluster'] <- 'A'
annotated_coords[annotated_coords$coordID %in% clusterB,'cluster'] <- 'B'
annotated_coords[annotated_coords$coordID %in% clusterC,'cluster'] <- 'C'
annotated_coords[annotated_coords$coordID %in% clusterD,'cluster'] <- 'D'

## Add in NusA and NusG stimulation info ---------------------------------------

NusG_noEco$NusG_padj <- NusG_noEco$padj
NusG_noEco$NusG_log2FC <- NusG_noEco$log2FoldChange
NusG_noEco$coordID <- rownames(NusG_noEco)

NusG_noEco$NusA_padj <- NusA_noEco$padj
NusA_noEco$NusA_log2FC <- NusA_noEco$log2FoldChange
NusA_noEco$coordID <- rownames(NusA_noEco)

NusAG_noEco$NusAG_padj <- NusAG_noEco$padj
NusAG_noEco$NusAG_log2FC <- NusAG_noEco$log2FoldChange
NusAG_noEco$coordID <- rownames(NusAG_noEco)

add_NusG <- left_join(annotated_coords, data.frame(NusG_noEco),
                      by = c("coordID" = "coordID"))
add_NusA <- left_join(add_NusG, data.frame(NusA_noEco),
                      by = c("coordID" = "coordID"))
add_NusAG <- left_join(add_NusA, data.frame(NusAG_noEco),
                       by = c("coordID" = "coordID"))

colsToKeep <- c("Nontemplate_strand",
                "coordID", "cluster", 'coordinate', 'strand', 'downstream_region',
                'Elemental_pause_match.','RNAfold_hairpin_region_44_structure',
                "RNAfold_hairpin_region_44_deltaG","U_count",
                "RNAfold_hairpin_region_44_spacerLength",
                "RNAfold_hairpin_region_44_stemLength",
                "RNAfold_hairpin_region_44_loopLength",
                "RNAfold_hairpin_region_44_deltaG",'U_count_tractRegion',
                "U_count_downstreamRegion","A_count_downstreamRegion",
                'NusA_log2FC','NusG_log2FC','NusAG_log2FC',
                'NusA_padj','NusG_padj','NusAG_padj')

observed_DF <- add_NusAG[,colsToKeep]
observed_DF[is.na(observed_DF$NusA_padj),'NusA_padj'] <- 10
observed_DF[is.na(observed_DF$NusG_padj),'NusG_padj'] <- 10
observed_DF[is.na(observed_DF$NusAG_padj),'NusAG_padj'] <- 10

observed_DF$NusA_log2FC <- as.double(observed_DF$NusA_log2FC)
observed_DF$NusG_log2FC <- as.double(observed_DF$NusG_log2FC)
observed_DF$NusAG_log2FC <- as.double(observed_DF$NusAG_log2FC)

observed_DF$NusG_stim <- 'No'
observed_DF[(observed_DF$NusG_padj <= 0.05) & (observed_DF$NusG_log2FC > 0),
            'NusG_stim'] <- 'Yes'
observed_DF$NusA_stim <- 'No'
observed_DF[(observed_DF$NusA_padj <= 0.05) & (observed_DF$NusA_log2FC > 0),
            'NusA_stim'] <- 'Yes'

## Plot l2fc for each cluster --------------------------------------------------

# For each pairwise comparison of interest:
# (each cluster in cluster_list compared to NE)
cluster_list <- c("NE","A","B","C","D")

l2fc_values <- c(observed_DF$NusA_log2FC,
                 observed_DF$NusG_log2FC,
                 observed_DF$NusAG_log2FC)

l2fc_labels <- c(rep("NusA",nrow(observed_DF)),
                 rep("NusG",nrow(observed_DF)),
                 rep("NusA + NusG",nrow(observed_DF)))

clusters <- c(rep(observed_DF$cluster,3))

l2fc_plot <- data.frame(l2fc = l2fc_values,
                        label = l2fc_labels,
                        cluster = clusters,
                        x_val = rep("x",length(l2fc_values)))

toPlot <- l2fc_plot %>%
  mutate(cluster = factor(cluster,
                          levels = c('NE','A','B','C','D'))) %>%
  mutate(label = factor(label,
                        levels = c('NusA','NusG','NusA + NusG'))) %>%
  filter(!(cluster == 'NE' & l2fc > 4))

for (i in 1:length(cluster_list)) {
  
  # Isolate the cluster comparison of interest
  clusterPlot <- toPlot[(toPlot$cluster == cluster_list[i]),]
  
  # Plot and save the NusG l2fc values
  l2fc_plot <- ggplot(clusterPlot,
                      aes(x = cluster, y = l2fc, fill = label)) +
    geom_violin(position = "dodge",
                linewidth = 2) +
    geom_boxplot(width = 0.2,
                 position = position_dodge(width = 0.9),
                 outlier.size = 3,
                 linewidth = 2) +
    theme_minimal() +
    geom_hline(yintercept=0, linetype = 'longdash',
               linewidth = 2) +
    theme(text=element_text(size=16, 
                            family="Arial"),
          legend.position = "none") +
    scale_fill_manual(values = c("#5F8ADC","#E26B7A","#BD67F3")) +
    ylim(c(-1.5,3.5))
  
  ggsave(filename = paste('3enrich_NusAG/termination_figures/extData_fig7c_cluster',cluster_list[i],'_l2fc_',
                          counts,'counts_padj',sig_value,'.png',sep=''),
         plot = l2fc_plot,
         width = 300,
         height = 300,
         units = 'mm',
         dpi = 300)
}

### Calculate L2FC statistical significance of clusters ------------------------

# Generate vectors for independent variable (cluster) 
# and dependent variable (L2FCs)
observed_DF$cluster_nusG <- paste(observed_DF$cluster, "_NusG", sep='')
observed_DF$cluster_nusA <- paste(observed_DF$cluster, "_NusA", sep='')
observed_DF$cluster_nusAG <- paste(observed_DF$cluster, "_NusAG", sep='')

all_l2fc <- c(observed_DF$NusG_log2FC, observed_DF$NusA_log2FC,
              observed_DF$NusAG_log2FC)
all_cluster <- c(observed_DF$cluster_nusG, observed_DF$cluster_nusA,
                 observed_DF$cluster_nusAG)

# Generate new DF for one-way ANOVA and Tukey posthoc testing
# Note: this will generate many more comparisons than required, 
# but will allow for more conservative p-value adjustment
anova_DF <- data.frame(cbind(all_l2fc, all_cluster))
anova_DF$all_l2fc <- as.double(anova_DF$all_l2fc)

# Run one-way ANOVA
model <- aov(all_l2fc ~ all_cluster, data = anova_DF)
# Run Tukey posthoc test
posthoc_info <- data.frame(TukeyHSD(model)$all_cluster)
posthoc_info$sig <- posthoc_info$p.adj <= 0.05

comparisons_of_interest <- c('NE_NusA-A_NusA',
                             'NE_NusG-A_NusG',
                             'NE_NusAG-A_NusAG',
                             'NE_NusA-B_NusA',
                             'NE_NusG-B_NusG',
                             'NE_NusAG-B_NusAG',
                             'NE_NusA-C_NusA',
                             'NE_NusG-C_NusG',
                             'NE_NusAG-C_NusAG',
                             'NE_NusA-D_NusA',
                             'NE_NusG-D_NusG',
                             'NE_NusAG-D_NusAG')

posthoc_DF_filt <- posthoc_info[rownames(posthoc_info) %in% comparisons_of_interest,]

write.csv(posthoc_DF_filt, '3enrich_NusAG/termination_figures/anova_cluster_posthoc_DF.csv')

## Add random coordinates ------------------------------------------------------

random_coords1 <- read.csv(paste('3enrich_NusAG/termination_figures/terminator_annotations/multifactor_random1_',
                                 counts,'counts_annotated.csv',sep=''))
random_coords2 <- read.csv(paste('3enrich_NusAG/termination_figures/terminator_annotations/multifactor_random2_',
                                 counts,'counts_annotated.csv',sep=''))
random_coords3 <- read.csv(paste('3enrich_NusAG/termination_figures/terminator_annotations/multifactor_random3_',
                                 counts,'counts_annotated.csv',sep=''))

random_coords1$coordID <- paste(random_coords1$coordinate,
                                random_coords1$strand,
                                sep='')
random_coords2$coordID <- paste(random_coords2$coordinate,
                                random_coords2$strand,
                                sep='')
random_coords3$coordID <- paste(random_coords3$coordinate,
                                random_coords3$strand,
                                sep='')

random_coords1$coordinate_category <- 'random'
random_coords2$coordinate_category <- 'random'
random_coords3$coordinate_category <- 'random'

random_coords1$random_replicate <- 1
random_coords2$random_replicate <- 2
random_coords3$random_replicate <- 3

random_coords1$NusG_log2FC <- "random"
random_coords2$NusG_log2FC <- "random"
random_coords3$NusG_log2FC <- "random"

random_coords1$NusA_log2FC <- "random"
random_coords2$NusA_log2FC <- "random"
random_coords3$NusA_log2FC <- "random"

random_coords1$NusAG_log2FC <- "random"
random_coords2$NusAG_log2FC <- "random"
random_coords3$NusAG_log2FC <- "random"

random_coords1$NusA_padj <- "random"
random_coords2$NusA_padj <- "random"
random_coords3$NusA_padj <- "random"

random_coords1$NusG_padj <- "random"
random_coords2$NusG_padj <- "random"
random_coords3$NusG_padj <- "random"

random_coords1$NusAG_padj <- "random"
random_coords2$NusAG_padj <- "random"
random_coords3$NusAG_padj <- "random"

random_coords1$NusG_stim <- "random"
random_coords2$NusG_stim <- "random"
random_coords3$NusG_stim <- "random"

random_coords1$NusA_stim <- "random"
random_coords2$NusA_stim <- "random"
random_coords3$NusA_stim <- "random"

random_coords1$NusAG_stim <- "random"
random_coords2$NusAG_stim <- "random"
random_coords3$NusAG_stim <- "random"

random_coords1$cluster<- "random"
random_coords2$cluster <- "random"
random_coords3$cluster <- "random"

random1 <- random_coords1[,c(colsToKeep, 'NusA_stim','NusG_stim')]
random2 <- random_coords2[,c(colsToKeep, 'NusA_stim','NusG_stim')]
random3 <- random_coords3[,c(colsToKeep, 'NusA_stim','NusG_stim')]

# Make sure that none of the random coords include observed
filtered_random1 <- data.frame(random1) %>%
  filter(!(coordID %in% observed_DF$coordID))
filtered_random2 <- data.frame(random2) %>%
  filter(!(coordID %in% observed_DF$coordID))
filtered_random3 <- data.frame(random3) %>%
  filter(!(coordID %in% observed_DF$coordID))

colsToRemove <- c("cluster_nusG","cluster_nusA","cluster_nusAG")

observed_DF$coordinate_category <- 'observed'
filtered_random1$coordinate_category <- 'random'
filtered_random2$coordinate_category <- 'random'
filtered_random3$coordinate_category <- 'random'

observed_random_DF <- rbind(observed_DF[,!(colnames(observed_DF) %in% colsToRemove)],
                            filtered_random1,
                            filtered_random2,
                            filtered_random3)

### Add inter- vs. intra-genic locations ---------------------------------------

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
                       "dot","gene_strand","dot2","metadata")

IDs = get_list_of_features(gene_DF$metadata, col=1)
Names = get_list_of_features(gene_DF$metadata, col=2)
Description = get_list_of_features(gene_DF$metadata, col=4)

gene_DF$geneID <- IDs
gene_DF$geneName <- Names
gene_DF$desc <- Description

# keep start, end, strand, gene_id, gene_name, and description columns
keep_cols <- c(4,5,7,10:12)
gene_DF_curated <- gene_DF[,keep_cols]
gene_DF_curated$start <- gene_DF_curated$start + 4641652
gene_DF_curated$end <- gene_DF_curated$end + 4641652
# Generate columns needed for the "inverse" dataframe
gene_DF_curated$startMin1 <- lead(gene_DF_curated$start, n=1)
gene_DF_curated$namesMin1 <- lead(gene_DF_curated$geneID, n=1)
# Generate the "inverse" dataframe- all of the intergenic spaces
inverseDF <- data.frame(geneID = paste(gene_DF_curated$geneID, "-",
                                       gene_DF_curated$namesMin1," intergenic region",sep=''),
                        start = gene_DF_curated$end, end = gene_DF_curated$startMin1,
                        gene_strand = "inter", geneName = gene_DF_curated$geneName,
                        desc = gene_DF_curated$desc)
fullDF <- na.omit(rbind(gene_DF_curated[,c(1:6)],inverseDF))
fullDF$start <- as.integer(fullDF$start)
fullDF$end <- as.integer(fullDF$end)

# Only include rows where the end is greater than the start 
#(because of overlapping genes, sometimes some weird rows get generated)
geneDF1 <- fullDF[(fullDF$end > (fullDF$start - 1)),]
geneDF <- geneDF1[(geneDF1$end != geneDF1$start),]

ordered_geneDF <- geneDF[order(geneDF$start),]

# For intergenic coordinates, add one ahead
ordered_geneDF$geneID_next <- lead(ordered_geneDF$geneID, n = 1)
ordered_geneDF$geneID_previous <- lag(ordered_geneDF$geneID, n = 1)

ordered_geneDF$geneName_next <- lead(ordered_geneDF$geneName, n = 1)
ordered_geneDF$geneName_previous <- lag(ordered_geneDF$geneName, n = 1)

ordered_geneDF$desc_next <- lead(ordered_geneDF$desc, n = 1)
ordered_geneDF$desc_previous <- lag(ordered_geneDF$desc, n = 1)

ordered_geneDF$gene_strand_next <- lead(ordered_geneDF$gene_strand, n = 1)
ordered_geneDF$gene_strand_previous <- lag(ordered_geneDF$gene_strand, n = 1)

ordered_geneDF$start_next <- lead(ordered_geneDF$start, n = 1)
ordered_geneDF$start_previous <- lag(ordered_geneDF$start, n = 1)

ordered_geneDF$end_next <- lead(ordered_geneDF$end, n = 1)
ordered_geneDF$end_previous <- lag(ordered_geneDF$end, n = 1)

# For intragenic coordinates, need 2 ahead
ordered_geneDF$geneID_next2 <- lead(ordered_geneDF$geneID, n = 2)
ordered_geneDF$geneID_previous2 <- lag(ordered_geneDF$geneID, n = 2)

ordered_geneDF$geneName_next2 <- lead(ordered_geneDF$geneName, n = 2)
ordered_geneDF$geneName_previous2 <- lag(ordered_geneDF$geneName, n = 2)

ordered_geneDF$desc_next2 <- lead(ordered_geneDF$desc, n = 2)
ordered_geneDF$desc_previous2 <- lag(ordered_geneDF$desc, n = 2)

ordered_geneDF$gene_strand_next2 <- lead(ordered_geneDF$gene_strand, n = 2)
ordered_geneDF$gene_strand_previous2 <- lag(ordered_geneDF$gene_strand, n = 2)

ordered_geneDF$start_next2 <- lead(ordered_geneDF$start, n = 2)
ordered_geneDF$start_previous2 <- lag(ordered_geneDF$start, n = 2)

ordered_geneDF$end_next2 <- lead(ordered_geneDF$end, n = 2)
ordered_geneDF$end_previous2 <- lag(ordered_geneDF$end, n = 2)

## Add gene annotations to coordinates -----------------------------------------

addGenes <- sqldf("select *
      from observed_random_DF inner join ordered_geneDF
        on (observed_random_DF.coordinate between ordered_geneDF.start and ordered_geneDF.end)")

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
addGenes[(addGenes$gene_strand == 'inter') & (addGenes$strand == '+') & (addGenes$gene_strand_previous == '+'),
         'upstream_geneID'] <- addGenes[(addGenes$gene_strand == 'inter') & (addGenes$strand == '+') & (addGenes$gene_strand_previous == '+'),
                                        'geneID_previous']
addGenes[(addGenes$gene_strand == 'inter') & (addGenes$strand == '+') & (addGenes$gene_strand_previous == '+'),
         'upstream_geneDistance'] <- addGenes[(addGenes$gene_strand == 'inter') & (addGenes$strand == '+') & (addGenes$gene_strand_previous == '+'),
                                              'coordinate'] - addGenes[(addGenes$gene_strand == 'inter') & (addGenes$strand == '+') & (addGenes$gene_strand_previous == '+'),
                                        'end_previous']
# If coordinate and upstream gene are in opposite orientations:
addGenes[(addGenes$gene_strand == 'inter') & (addGenes$strand == '+') & (addGenes$gene_strand_previous == '-'),
         'upstream_geneID'] <- paste(addGenes[(addGenes$gene_strand == 'inter') & (addGenes$strand == '+') & (addGenes$gene_strand_previous == '-'),
                                        'geneID_previous'], "_as", sep='')
addGenes[(addGenes$gene_strand == 'inter') & (addGenes$strand == '+') & (addGenes$gene_strand_previous == '-'),
         'upstream_geneDistance'] <- addGenes[(addGenes$gene_strand == 'inter') & (addGenes$strand == '+') & (addGenes$gene_strand_previous == '-'),
                                              'coordinate'] - addGenes[(addGenes$gene_strand == 'inter') & (addGenes$strand == '+') & (addGenes$gene_strand_previous == '-'),
                                                                       'end_previous']


addGenes[(addGenes$gene_strand == 'inter') & (addGenes$strand == '+') & (addGenes$gene_strand_previous == '+'),
         'upstream_geneName'] <- addGenes[(addGenes$gene_strand == 'inter') & (addGenes$strand == '+') & (addGenes$gene_strand_previous == '+'),
                                        'geneName_previous']
addGenes[(addGenes$gene_strand == 'inter') & (addGenes$strand == '+') & (addGenes$gene_strand_previous == '-'),
         'upstream_geneName'] <- paste(addGenes[(addGenes$gene_strand == 'inter') & (addGenes$strand == '+') & (addGenes$gene_strand_previous == '-'),
                                              'geneName_previous'], " (antisense)", sep='')

addGenes[(addGenes$gene_strand == 'inter') & (addGenes$strand == '+') & (addGenes$gene_strand_previous == '+'),
         'upstream_geneDesc'] <- addGenes[(addGenes$gene_strand == 'inter') & (addGenes$strand == '+') & (addGenes$gene_strand_previous == '+'),
                                        'desc_previous']
addGenes[(addGenes$gene_strand == 'inter') & (addGenes$strand == '+') & (addGenes$gene_strand_previous == '-'),
         'upstream_geneDesc'] <- paste(addGenes[(addGenes$gene_strand == 'inter') & (addGenes$strand == '+') & (addGenes$gene_strand_previous == '-'),
                                              'desc_previous'], " (antisense)", sep='')

addGenes[(addGenes$gene_strand == 'inter') & (addGenes$strand == '+') & (addGenes$gene_strand_next == '+'),
         'downstream_geneID'] <- addGenes[(addGenes$gene_strand == 'inter') & (addGenes$strand == '+') & (addGenes$gene_strand_next == '+'),
                                        'geneID_next']
addGenes[(addGenes$gene_strand == 'inter') & (addGenes$strand == '+') & (addGenes$gene_strand_next == '+'),
         'downstream_geneDistance'] <- addGenes[(addGenes$gene_strand == 'inter') & (addGenes$strand == '+') & (addGenes$gene_strand_next == '+'),
                                                'start_next'] - addGenes[(addGenes$gene_strand == 'inter') & (addGenes$strand == '+') & (addGenes$gene_strand_next == '+'),
                                              'coordinate']
addGenes[(addGenes$gene_strand == 'inter') & (addGenes$strand == '+') & (addGenes$gene_strand_next == '-'),
         'downstream_geneID'] <- paste(addGenes[(addGenes$gene_strand == 'inter') & (addGenes$strand == '+') & (addGenes$gene_strand_next == '-'),
                                              'geneID_next'], "_as", sep='')
addGenes[(addGenes$gene_strand == 'inter') & (addGenes$strand == '+') & (addGenes$gene_strand_next == '-'),
         'downstream_geneDistance'] <- addGenes[(addGenes$gene_strand == 'inter') & (addGenes$strand == '+') & (addGenes$gene_strand_next == '-'),
                                                                         'start_next'] - addGenes[(addGenes$gene_strand == 'inter') & (addGenes$strand == '+') & (addGenes$gene_strand_next == '-'),
                                                                                                  'coordinate']

addGenes[(addGenes$gene_strand == 'inter') & (addGenes$strand == '+') & (addGenes$gene_strand_next == '+'),
         'downstream_geneName'] <- addGenes[(addGenes$gene_strand == 'inter') & (addGenes$strand == '+') & (addGenes$gene_strand_next == '+'),
                                          'geneName_next']
addGenes[(addGenes$gene_strand == 'inter') & (addGenes$strand == '+') & (addGenes$gene_strand_next == '-'),
         'downstream_geneName'] <- paste(addGenes[(addGenes$gene_strand == 'inter') & (addGenes$strand == '+') & (addGenes$gene_strand_next == '-'),
                                                'geneName_next'], " (antisense)", sep='')

addGenes[(addGenes$gene_strand == 'inter') & (addGenes$strand == '+') & (addGenes$gene_strand_next == '+'),
         'downstream_geneDesc'] <- addGenes[(addGenes$gene_strand == 'inter') & (addGenes$strand == '+') & (addGenes$gene_strand_next == '+'),
                                          'desc_next']
addGenes[(addGenes$gene_strand == 'inter') & (addGenes$strand == '+') & (addGenes$gene_strand_next == '-'),
         'downstream_geneDesc'] <- paste(addGenes[(addGenes$gene_strand == 'inter') & (addGenes$strand == '+') & (addGenes$gene_strand_next == '-'),
                                                'desc_next'], " (antisense)", sep='')

# Negative strand intergenic terminators:
addGenes[(addGenes$gene_strand == 'inter') & (addGenes$strand == '-') & (addGenes$gene_strand_next == '-'),
         'upstream_geneID'] <- addGenes[(addGenes$gene_strand == 'inter') & (addGenes$strand == '-') & (addGenes$gene_strand_next == '-'),
                                        'geneID_next']
addGenes[(addGenes$gene_strand == 'inter') & (addGenes$strand == '-') & (addGenes$gene_strand_next == '-'),
         'upstream_geneDistance'] <- addGenes[(addGenes$gene_strand == 'inter') & (addGenes$strand == '-') & (addGenes$gene_strand_next == '-'),
                                              'start_next'] - addGenes[(addGenes$gene_strand == 'inter') & (addGenes$strand == '-') & (addGenes$gene_strand_next == '-'),
                                              'coordinate']
addGenes[(addGenes$gene_strand == 'inter') & (addGenes$strand == '-') & (addGenes$gene_strand_next == '+'),
         'upstream_geneID'] <- paste(addGenes[(addGenes$gene_strand == 'inter') & (addGenes$strand == '-') & (addGenes$gene_strand_next == '+'),
                                              'geneID_next'], "_as", sep='')
addGenes[(addGenes$gene_strand == 'inter') & (addGenes$strand == '-') & (addGenes$gene_strand_next == '+'),
         'upstream_geneDistance'] <- addGenes[(addGenes$gene_strand == 'inter') & (addGenes$strand == '-') & (addGenes$gene_strand_next == '+'),
                                              'start_next'] - addGenes[(addGenes$gene_strand == 'inter') & (addGenes$strand == '-') & (addGenes$gene_strand_next == '+'),
                                              'coordinate']

addGenes[(addGenes$gene_strand == 'inter') & (addGenes$strand == '-') & (addGenes$gene_strand_next == '-'),
         'upstream_geneName'] <- addGenes[(addGenes$gene_strand == 'inter') & (addGenes$strand == '-') & (addGenes$gene_strand_next == '-'),
                                          'geneName_next']
addGenes[(addGenes$gene_strand == 'inter') & (addGenes$strand == '-') & (addGenes$gene_strand_next == '+'),
         'upstream_geneName'] <- paste(addGenes[(addGenes$gene_strand == 'inter') & (addGenes$strand == '-') & (addGenes$gene_strand_next == '+'),
                                                'geneName_next'], " (antisense)", sep='')

addGenes[(addGenes$gene_strand == 'inter') & (addGenes$strand == '-') & (addGenes$gene_strand_next == '-'),
         'upstream_geneDesc'] <- addGenes[(addGenes$gene_strand == 'inter') & (addGenes$strand == '-') & (addGenes$gene_strand_next == '-'),
                                          'desc_next']
addGenes[(addGenes$gene_strand == 'inter') & (addGenes$strand == '-') & (addGenes$gene_strand_next == '+'),
         'upstream_geneDesc'] <- paste(addGenes[(addGenes$gene_strand == 'inter') & (addGenes$strand == '-') & (addGenes$gene_strand_next == '+'),
                                                'desc_next'], " (antisense)", sep='')

addGenes[(addGenes$gene_strand == 'inter') & (addGenes$strand == '-') & (addGenes$gene_strand_previous == '-'),
         'downstream_geneID'] <- addGenes[(addGenes$gene_strand == 'inter') & (addGenes$strand == '-') & (addGenes$gene_strand_previous == '-'),
                                          'geneID_previous']
addGenes[(addGenes$gene_strand == 'inter') & (addGenes$strand == '-') & (addGenes$gene_strand_previous == '-'),
         'downstream_geneDistance'] <- addGenes[(addGenes$gene_strand == 'inter') & (addGenes$strand == '-') & (addGenes$gene_strand_previous == '-'),
                                                'coordinate'] - addGenes[(addGenes$gene_strand == 'inter') & (addGenes$strand == '-') & (addGenes$gene_strand_previous == '-'),
                                                'end_previous']
addGenes[(addGenes$gene_strand == 'inter') & (addGenes$strand == '-') & (addGenes$gene_strand_previous == '+'),
         'downstream_geneID'] <- paste(addGenes[(addGenes$gene_strand == 'inter') & (addGenes$strand == '-') & (addGenes$gene_strand_previous == '+'),
                                                'geneID_previous'], "_as", sep='')
addGenes[(addGenes$gene_strand == 'inter') & (addGenes$strand == '-') & (addGenes$gene_strand_previous == '+'),
         'downstream_geneDistance'] <- addGenes[(addGenes$gene_strand == 'inter') & (addGenes$strand == '-') & (addGenes$gene_strand_previous == '+'),
                                                'coordinate'] - addGenes[(addGenes$gene_strand == 'inter') & (addGenes$strand == '-') & (addGenes$gene_strand_previous == '+'),
                                                'end_previous']
addGenes[(addGenes$gene_strand == 'inter') & (addGenes$strand == '-') & (addGenes$gene_strand_previous == '-'),
         'downstream_geneName'] <- addGenes[(addGenes$gene_strand == 'inter') & (addGenes$strand == '-') & (addGenes$gene_strand_previous == '-'),
                                            'geneName_previous']
addGenes[(addGenes$gene_strand == 'inter') & (addGenes$strand == '-') & (addGenes$gene_strand_previous == '+'),
         'downstream_geneName'] <- paste(addGenes[(addGenes$gene_strand == 'inter') & (addGenes$strand == '-') & (addGenes$gene_strand_previous == '+'),
                                                  'geneName_previous'], " (antisense)", sep='')

addGenes[(addGenes$gene_strand == 'inter') & (addGenes$strand == '-') & (addGenes$gene_strand_previous == '-'),
         'downstream_geneDesc'] <- addGenes[(addGenes$gene_strand == 'inter') & (addGenes$strand == '-') & (addGenes$gene_strand_previous == '-'),
                                            'desc_previous']
addGenes[(addGenes$gene_strand == 'inter') & (addGenes$strand == '-') & (addGenes$gene_strand_previous == '+'),
         'downstream_geneDesc'] <- paste(addGenes[(addGenes$gene_strand == 'inter') & (addGenes$strand == '-') & (addGenes$gene_strand_previous == '+'),
                                                  'desc_previous'], " (antisense)", sep='')

# Positive strand intragenic terminators:
addGenes[(addGenes$gene_strand == '+') & (addGenes$strand == '+'),
         'upstream_geneID'] <- addGenes[(addGenes$gene_strand == '+') & (addGenes$strand == '+'),
                                        'geneID']
addGenes[(addGenes$gene_strand == '+') & (addGenes$strand == '+'),
         'upstream_geneDistance'] <- 0
addGenes[(addGenes$gene_strand == '-') & (addGenes$strand == '+'),
         'upstream_geneID'] <- paste(addGenes[(addGenes$gene_strand == '-') & (addGenes$strand == '+'),
                                              'geneID'], "_as", sep='')
addGenes[(addGenes$gene_strand == '-') & (addGenes$strand == '+'),
         'upstream_geneDistance'] <- 0

addGenes[(addGenes$gene_strand == '+') & (addGenes$strand == '+'),
         'upstream_geneName'] <- addGenes[(addGenes$gene_strand == '+') & (addGenes$strand == '+'),
                                        'geneName']
addGenes[(addGenes$gene_strand == '-') & (addGenes$strand == '+'),
         'upstream_geneName'] <- paste(addGenes[(addGenes$gene_strand == '-') & (addGenes$strand == '+'),
                                              'geneName'], " (antisense)", sep='')

addGenes[(addGenes$gene_strand == '+') & (addGenes$strand == '+'),
         'upstream_geneDesc'] <- addGenes[(addGenes$gene_strand == '+') & (addGenes$strand == '+'),
                                        'desc']
addGenes[(addGenes$gene_strand == '-') & (addGenes$strand == '+'),
         'upstream_geneDesc'] <- paste(addGenes[(addGenes$gene_strand == '-') & (addGenes$strand == '+'),
                                              'desc'], " (antisense)", sep='')

addGenes[(addGenes$gene_strand != 'inter') & (addGenes$strand == '+') & (addGenes$gene_strand_next2 == '+'),
         'downstream_geneID'] <- addGenes[(addGenes$gene_strand != 'inter') & (addGenes$strand == '+') & (addGenes$gene_strand_next2 == '+'),
                                          'geneID_next2']
addGenes[(addGenes$gene_strand != 'inter') & (addGenes$strand == '+') & (addGenes$gene_strand_next2 == '+'),
         'downstream_geneDistance'] <- addGenes[(addGenes$gene_strand != 'inter') & (addGenes$strand == '+') & (addGenes$gene_strand_next2 == '+'),
                                          'start_next2'] - addGenes[(addGenes$gene_strand != 'inter') & (addGenes$strand == '+') & (addGenes$gene_strand_next2 == '+'),
                                                                    'coordinate']
addGenes[(addGenes$gene_strand != 'inter') & (addGenes$strand == '+') & (addGenes$gene_strand_next2 == '-'),
         'downstream_geneID'] <- paste(addGenes[(addGenes$gene_strand != 'inter') & (addGenes$strand == '+') & (addGenes$gene_strand_next2 == '-'),
                                                'geneID_next2'], " (antisense)", sep='')
addGenes[(addGenes$gene_strand != 'inter') & (addGenes$strand == '+') & (addGenes$gene_strand_next2 == '-'),
         'downstream_geneDistance'] <- addGenes[(addGenes$gene_strand != 'inter') & (addGenes$strand == '+') & (addGenes$gene_strand_next2 == '-'),
                                                'start_next2'] - addGenes[(addGenes$gene_strand != 'inter') & (addGenes$strand == '+') & (addGenes$gene_strand_next2 == '-'),
                                                                          'coordinate']

addGenes[(addGenes$gene_strand != 'inter') & (addGenes$strand == '+') & (addGenes$gene_strand_next2 == '+'),
         'downstream_geneName'] <- addGenes[(addGenes$gene_strand != 'inter') & (addGenes$strand == '+') & (addGenes$gene_strand_next2 == '+'),
                                          'geneName_next2']
addGenes[(addGenes$gene_strand != 'inter') & (addGenes$strand == '+') & (addGenes$gene_strand_next2 == '-'),
         'downstream_geneName'] <- paste(addGenes[(addGenes$gene_strand != 'inter') & (addGenes$strand == '+') & (addGenes$gene_strand_next2 == '-'),
                                                'geneName_next2'], " (antisense)", sep='')

addGenes[(addGenes$gene_strand != 'inter') & (addGenes$strand == '+') & (addGenes$gene_strand_next2 == '+'),
         'downstream_geneDesc'] <- addGenes[(addGenes$gene_strand != 'inter') & (addGenes$strand == '+') & (addGenes$gene_strand_next2 == '+'),
                                          'desc_next2']
addGenes[(addGenes$gene_strand != 'inter') & (addGenes$strand == '+') & (addGenes$gene_strand_next2 == '-'),
         'downstream_geneDesc'] <- paste(addGenes[(addGenes$gene_strand != 'inter') & (addGenes$strand == '+') & (addGenes$gene_strand_next2 == '-'),
                                                'desc_next2'], " (antisense)", sep='')

# Negative strand intragenic terminators:
addGenes[(addGenes$gene_strand == '-') & (addGenes$strand == '-'),
         'upstream_geneID'] <- addGenes[(addGenes$gene_strand == '-') & (addGenes$strand == '-'),
                                        'geneID']
addGenes[(addGenes$gene_strand == '-') & (addGenes$strand == '-'),
         'upstream_geneDistance'] <- 0
addGenes[(addGenes$gene_strand == '+') & (addGenes$strand == '-'),
         'upstream_geneID'] <- paste(addGenes[(addGenes$gene_strand == '+') & (addGenes$strand == '-'),
                                              'geneID'], "_as", sep='')
addGenes[(addGenes$gene_strand == '+') & (addGenes$strand == '-'),
         'upstream_geneDistance'] <- 0

addGenes[(addGenes$gene_strand == '-') & (addGenes$strand == '-'),
         'upstream_geneName'] <- addGenes[(addGenes$gene_strand == '-') & (addGenes$strand == '-'),
                                        'geneName']
addGenes[(addGenes$gene_strand == '+') & (addGenes$strand == '-'),
         'upstream_geneName'] <- paste(addGenes[(addGenes$gene_strand == '+') & (addGenes$strand == '-'),
                                              'geneName'], " (antisense)", sep='')

addGenes[(addGenes$gene_strand == '-') & (addGenes$strand == '-'),
         'upstream_geneDesc'] <- addGenes[(addGenes$gene_strand == '-') & (addGenes$strand == '-'),
                                        'desc']
addGenes[(addGenes$gene_strand == '+') & (addGenes$strand == '-'),
         'upstream_geneDesc'] <- paste(addGenes[(addGenes$gene_strand == '+') & (addGenes$strand == '-'),
                                              'desc'], " (antisense)", sep='')

addGenes[(addGenes$gene_strand != 'inter') & (addGenes$strand == '-') & (addGenes$gene_strand_previous2 == '-'),
         'downstream_geneID'] <- addGenes[(addGenes$gene_strand != 'inter') & (addGenes$strand == '-') & (addGenes$gene_strand_previous2 == '-'),
                                          'geneID_previous2']
addGenes[(addGenes$gene_strand != 'inter') & (addGenes$strand == '-') & (addGenes$gene_strand_next2 == '-'),
         'downstream_geneDistance'] <- addGenes[(addGenes$gene_strand != 'inter') & (addGenes$strand == '-') & (addGenes$gene_strand_next2 == '-'),
                                                'coordinate'] - addGenes[(addGenes$gene_strand != 'inter') & (addGenes$strand == '-') & (addGenes$gene_strand_next2 == '-'),
                                                'end_previous2']
addGenes[(addGenes$gene_strand != 'inter') & (addGenes$strand == '-') & (addGenes$gene_strand_previous2 == '+'),
         'downstream_geneID'] <- paste(addGenes[(addGenes$gene_strand != 'inter') & (addGenes$strand == '-') & (addGenes$gene_strand_previous2 == '+'),
                                                'geneID_previous2'], "_as", sep='')
addGenes[(addGenes$gene_strand != 'inter') & (addGenes$strand == '-') & (addGenes$gene_strand_next2 == '+'),
         'downstream_geneDistance'] <- addGenes[(addGenes$gene_strand != 'inter') & (addGenes$strand == '-') & (addGenes$gene_strand_next2 == '+'),
                                                'coordinate'] - addGenes[(addGenes$gene_strand != 'inter') & (addGenes$strand == '-') & (addGenes$gene_strand_next2 == '+'),
                                                                         'end_previous2']

addGenes[(addGenes$gene_strand != 'inter') & (addGenes$strand == '-') & (addGenes$gene_strand_previous2 == '-'),
         'downstream_geneName'] <- addGenes[(addGenes$gene_strand != 'inter') & (addGenes$strand == '-') & (addGenes$gene_strand_previous2 == '-'),
                                          'geneName_previous2']
addGenes[(addGenes$gene_strand != 'inter') & (addGenes$strand == '-') & (addGenes$gene_strand_previous2 == '+'),
         'downstream_geneName'] <- paste(addGenes[(addGenes$gene_strand != 'inter') & (addGenes$strand == '-') & (addGenes$gene_strand_previous2 == '+'),
                                                'geneName_previous2'], " (antisense)", sep='')

addGenes[(addGenes$gene_strand != 'inter') & (addGenes$strand == '-') & (addGenes$gene_strand_previous2 == '-'),
         'downstream_geneDesc'] <- addGenes[(addGenes$gene_strand != 'inter') & (addGenes$strand == '-') & (addGenes$gene_strand_previous2 == '-'),
                                          'desc_previous2']
addGenes[(addGenes$gene_strand != 'inter') & (addGenes$strand == '-') & (addGenes$gene_strand_previous2 == '+'),
         'downstream_geneDesc'] <- paste(addGenes[(addGenes$gene_strand != 'inter') & (addGenes$strand == '-') & (addGenes$gene_strand_previous2 == '+'),
                                                'desc_previous2'], " (antisense)", sep='')

# Specify whether a coordinate is inter- or intragenic
addGenes$GeneLocation <- 'intra'
addGenes[addGenes$gene_strand == 'inter',
             'GeneLocation'] <- 'inter'

### Enrichment of coordinates in intra vs. intergenic regions ------------------

# Filter subset dataframes for enrichment analysis
observed_coords_only <- addGenes %>%
  filter(coordinate_category == 'observed')

random_only <- addGenes %>%
  filter(coordinate_category == 'random')

equal_size_random <- random_only[c(1:nrow(observed_coords_only)),]

equal_size_random_observed <- rbind(observed_coords_only,
                                    equal_size_random)

# Functions to carry out Fisher's exact testing on location annotation
fishers_coordCategory_GeneLocation <- function(input_DF, location, category, num_tests) {
  
  category_location <- input_DF[(input_DF$coordinate_category == category) & (input_DF$GeneLocation == location),]
  category_nolocation <- input_DF[(!(input_DF$coordID %in% category_location$coordID)) & (input_DF$coordinate_category == category),]
  nocategory_location <- input_DF[(input_DF$coordinate_category != category) & (input_DF$GeneLocation == location),]
  nocategory_nolocation <- input_DF[(input_DF$coordinate_category != category) & (!(input_DF$coordID %in% nocategory_location$coordID)),]
  
  category_location_num <- nrow(category_location)
  category_nolocation_num <- nrow(category_nolocation)
  nocategory_location_num <- nrow(nocategory_location)
  nocategory_nolocation_num <- nrow(nocategory_nolocation)
  
  fisher_matrix <- matrix(c(category_location_num, 
                            category_nolocation_num,
                            nocategory_location_num,
                            nocategory_nolocation_num),
                          nrow = 2, ncol= 2)
  
  fishers_DF <- fisher.test(fisher_matrix, conf.level=(100-(5/num_tests))/100)
  return(fishers_DF)
}

fishers_NusA_GeneLocation <- function(input_DF, location, category, num_tests) {
  
  category_location <- input_DF[(input_DF$NusA_stim == category) & (input_DF$GeneLocation == location),]
  category_nolocation <- input_DF[(!(input_DF$coordID %in% category_location$coordID)) & (input_DF$NusA_stim == category),]
  nocategory_location <- input_DF[(input_DF$NusA_stim != category) & (input_DF$GeneLocation == location),]
  nocategory_nolocation <- input_DF[(input_DF$NusA_stim != category) & (!(input_DF$coordID %in% nocategory_location$coordID)),]
  
  category_location_num <- nrow(category_location)
  category_nolocation_num <- nrow(category_nolocation)
  nocategory_location_num <- nrow(nocategory_location)
  nocategory_nolocation_num <- nrow(nocategory_nolocation)
  
  fisher_matrix <- matrix(c(category_location_num, 
                            category_nolocation_num,
                            nocategory_location_num,
                            nocategory_nolocation_num),
                          nrow = 2, ncol= 2)
  
  fishers_DF <- fisher.test(fisher_matrix, conf.level=(100-(5/num_tests))/100)
  return(fishers_DF)
}

fishers_NusG_GeneLocation <- function(input_DF, location, category, num_tests) {
  
  category_location <- input_DF[(input_DF$NusG_stim == category) & (input_DF$GeneLocation == location),]
  category_nolocation <- input_DF[(!(input_DF$coordID %in% category_location$coordID)) & (input_DF$NusG_stim == category),]
  nocategory_location <- input_DF[(input_DF$NusG_stim != category) & (input_DF$GeneLocation == location),]
  nocategory_nolocation <- input_DF[(input_DF$NusG_stim != category) & (!(input_DF$coordID %in% nocategory_location$coordID)),]
  
  category_location_num <- nrow(category_location)
  category_nolocation_num <- nrow(category_nolocation)
  nocategory_location_num <- nrow(nocategory_location)
  nocategory_nolocation_num <- nrow(nocategory_nolocation)
  
  fisher_matrix <- matrix(c(category_location_num, 
                            category_nolocation_num,
                            nocategory_location_num,
                            nocategory_nolocation_num),
                          nrow = 2, ncol= 2)
  
  fishers_DF <- fisher.test(fisher_matrix, conf.level=(100-(5/num_tests))/100)
  return(fishers_DF)
}

location_list <- c('intra', 'inter')

# Initialize empty list
p_values <- c()
odds_ratios <- c()
CI_lower <- c()
CI_upper <- c()

for (i in 1:length(location_list)) {
  fishers_DF <- fishers_coordCategory_GeneLocation(equal_size_random_observed,
                                                   location_list[i],
                                                   "observed",
                                                   length(location_list))
  p_values <- c(p_values, fishers_DF$p.value)
  odds_ratios <- c(odds_ratios, fishers_DF$estimate)
  CI_lower <- c(CI_lower, fishers_DF$conf.int[1])
  CI_upper <- c(CI_upper, fishers_DF$conf.int[2])
}

OR_DF <- data.frame(location = location_list,
                    pvalues = p_values,
                    odds_ratios = odds_ratios,
                    CI_lower = CI_lower,
                    CI_upper = CI_upper)

OR_DF$Bonferroni_sig_cutoff <- 0.05/length(location_list)

OR_DF$Sig <- OR_DF$pvalues < OR_DF$Bonferroni_sig_cutoff
OR_DF$padj <- p.adjust(OR_DF$pvalues, method = 'bonferroni')

OR_DF$log10_OR <- log10(OR_DF$odds_ratios)
OR_DF$log10_CI_upper <- log10(OR_DF$CI_upper)
OR_DF$log10_CI_lower <- log10(OR_DF$CI_lower)

observed_coord_enrichment <- ggplot(OR_DF, aes(x = log10_OR, y = location)) + 
  geom_vline(aes(xintercept = 0), linewidth = .25, linetype = "dashed") + 
  geom_errorbarh(aes(xmax = log10_CI_upper, xmin = log10_CI_lower),
                 size = 2, height = .2, color = "gray") +
  geom_point(size = 10, color = "grey") +
  theme_bw() +
  theme(panel.grid.minor = element_blank())

ggsave(filename = paste('3enrich_NusAG/termination_figures/extData_fig6e_',
                        'location_random',counts,
                        'counts_padj',sig_value,'.png',sep=''),
       plot = observed_coord_enrichment,
       width = 200,
       height = 100,
       units = 'mm',
       dpi = 300)
write.csv(OR_DF, paste('3enrich_NusAG/termination_figures/extData_fig8d_',
                       'location_observed_random_',counts,
                       'counts_padj',sig_value,'.csv',sep=''))

equal_size_random_observed$coordinate_category <- factor(equal_size_random_observed$coordinate_category,
                                                         levels = c("random", "observed"))
equal_size_random_observed$GeneLocation <- factor(equal_size_random_observed$GeneLocation,
                                                  levels = c("intra", "inter"))

observed_random_bar_graph <- ggplot(equal_size_random_observed) +
  geom_bar(data = subset(equal_size_random_observed, coordinate_category == 'random'),
           aes(x = GeneLocation, y = after_stat(count / sum(count))),
           stat = 'count', just = 1, fill = 'black',
           width = 0.3) + 
  geom_bar(data = subset(equal_size_random_observed, coordinate_category == 'observed'),
           aes(x = GeneLocation, y = after_stat(count / sum(count))),
           stat = 'count', just = 0, fill = 'grey',
           width = 0.3) + 
  ylim(0,1) +
  theme_minimal()

ggsave(filename = paste('3enrich_NusAG/termination_figures/extData_fig6e_',
                        'barGraph_observed_vs_random.png',sep=''),
       plot = observed_random_bar_graph,
       width = 200,
       height = 200,
       units = 'mm',
       dpi = 300)

# Initialize empty list
p_values <- c()
odds_ratios <- c()
CI_lower <- c()
CI_upper <- c()

for (i in 1:length(location_list)) {
  fishers_DF <- fishers_NusG_GeneLocation(observed_coords_only, location_list[i],
                                          "Yes", length(location_list))
  p_values <- c(p_values, fishers_DF$p.value)
  odds_ratios <- c(odds_ratios, fishers_DF$estimate)
  CI_lower <- c(CI_lower, fishers_DF$conf.int[1])
  CI_upper <- c(CI_upper, fishers_DF$conf.int[2])
}

OR_DF <- data.frame(location = location_list,
                    pvalues = p_values,
                    odds_ratios = odds_ratios,
                    CI_lower = CI_lower,
                    CI_upper = CI_upper)

OR_DF$Bonferroni_sig_cutoff <- 0.05/length(location_list)

OR_DF$Sig <- OR_DF$pvalues < OR_DF$Bonferroni_sig_cutoff
OR_DF$padj <- p.adjust(OR_DF$pvalues, method = 'bonferroni')

OR_DF$log10_OR <- log10(OR_DF$odds_ratios)
OR_DF$log10_CI_upper <- log10(OR_DF$CI_upper)
OR_DF$log10_CI_lower <- log10(OR_DF$CI_lower)

NusG_coord_enrichment <- ggplot(OR_DF, aes(x = log10_OR, y = location)) + 
  geom_vline(aes(xintercept = 0), linewidth = .25, linetype = "dashed") + 
  geom_errorbarh(aes(xmax = log10_CI_upper, xmin = log10_CI_lower),
                 size = 2, height = .2, color = "#F26B7A") +
  geom_point(size = 10, color = "#F26B7A") +
  theme_bw() +
  theme(panel.grid.minor = element_blank())

ggsave(filename = paste('3enrich_NusAG/termination_figures/extData_fig8d_',
                        'location_NusG_stim_',counts,
                        'counts_padj',sig_value,'.png',sep=''),
       plot = NusG_coord_enrichment,
       width = 200,
       height = 100,
       units = 'mm',
       dpi = 300)
write.csv(OR_DF, paste('3enrich_NusAG/termination_figures/extData_fig8d_',
                       'location_NusG_stim_',counts,
                       'counts_padj',sig_value,'.csv',sep=''))

observed_coords_only$coordinate_category <- factor(observed_coords_only$NusG_stim,
                                                   levels = c("No", "Yes"))
observed_coords_only$GeneLocation <- factor(observed_coords_only$GeneLocation,
                                            levels = c("intra", "inter"))

NusG_stim_bar_graph <- ggplot(observed_coords_only) +
  geom_bar(data = subset(observed_coords_only, NusG_stim == 'No'),
           aes(x = GeneLocation, y = after_stat(count / sum(count))),
           stat = 'count', just = 1, fill = 'grey',
           width = 0.3) + 
  geom_bar(data = subset(observed_coords_only, NusG_stim == 'Yes'),
           aes(x = GeneLocation, y = after_stat(count / sum(count))),
           stat = 'count', just = 0, fill = '#F26B7A',
           width = 0.3) + 
  ylim(0,1) +
  theme_minimal()

ggsave(filename = paste('3enrich_NusAG/termination_figures/extData_fig8d_',
                        'barGraph_NusGstim.png',sep=''),
       plot = NusG_stim_bar_graph,
       width = 200,
       height = 100,
       units = 'mm',
       dpi = 300)

# Initialize empty list
p_values <- c()
odds_ratios <- c()
CI_lower <- c()
CI_upper <- c()

for (i in 1:length(location_list)) {
  fishers_DF <- fishers_NusA_GeneLocation(observed_coords_only, location_list[i],
                                          "Yes", length(location_list))
  p_values <- c(p_values, fishers_DF$p.value)
  odds_ratios <- c(odds_ratios, fishers_DF$estimate)
  CI_lower <- c(CI_lower, fishers_DF$conf.int[1])
  CI_upper <- c(CI_upper, fishers_DF$conf.int[2])
}

OR_DF <- data.frame(location = location_list,
                    pvalues = p_values,
                    odds_ratios = odds_ratios,
                    CI_lower = CI_lower,
                    CI_upper = CI_upper)

OR_DF$Bonferroni_sig_cutoff <- 0.05/length(location_list)

OR_DF$Sig <- OR_DF$pvalues < OR_DF$Bonferroni_sig_cutoff
OR_DF$padj <- p.adjust(OR_DF$pvalues, method = 'bonferroni')

OR_DF$log10_OR <- log10(OR_DF$odds_ratios)
OR_DF$log10_CI_upper <- log10(OR_DF$CI_upper)
OR_DF$log10_CI_lower <- log10(OR_DF$CI_lower)

NusA_coord_enrichment <- ggplot(OR_DF, aes(x = log10_OR, y = location)) + 
  geom_vline(aes(xintercept = 0), linewidth = .25, linetype = "dashed") + 
  geom_errorbarh(aes(xmax = log10_CI_upper, xmin = log10_CI_lower),
                 size = 2, height = .2, color = "#5F9AEC") +
  geom_point(size = 10, color = "#5F9AEC") +
  theme_bw() +
  theme(panel.grid.minor = element_blank())

ggsave(filename = paste('3enrich_NusAG/termination_figures/extData_fig8d_',
                        'location_NusA_stim_',counts,
                        'counts_padj',sig_value,'.png',sep=''),
       plot = NusA_coord_enrichment,
       width = 200,
       height = 100,
       units = 'mm',
       dpi = 300)
write.csv(OR_DF, paste('3enrich_NusAG/termination_figures/extData_fig8d_',
                       'location_NusA_stim_',counts,
                       'counts_padj',sig_value,'.csv',sep=''))

observed_coords_only$coordinate_category <- factor(observed_coords_only$NusA_stim,
                                                   levels = c("No", "Yes"))
observed_coords_only$GeneLocation <- factor(observed_coords_only$GeneLocation,
                                            levels = c("intra", "inter"))

NusA_stim_bar_graph <- ggplot(observed_coords_only) +
  geom_bar(data = subset(observed_coords_only, NusA_stim == 'No'),
           aes(x = GeneLocation, y = after_stat(count / sum(count))),
           stat = 'count', just = 1, fill = 'grey',
           width = 0.3) + 
  geom_bar(data = subset(observed_coords_only, NusA_stim == 'Yes'),
           aes(x = GeneLocation, y = after_stat(count / sum(count))),
           stat = 'count', just = 0, fill = '#5F9AEC',
           width = 0.3) + 
  ylim(0,1) +
  theme_minimal()

ggsave(filename = paste('3enrich_NusAG/termination_figures/extData_fig8d_',
                        'barGraph_NusAstim.png',sep=''),
       plot = NusA_stim_bar_graph,
       width = 200,
       height = 200,
       units = 'mm',
       dpi = 300)

### Add in location annotations ------------------------------------------------

addGenes$LocationAnnotation <- NA

annotateGenesLocations <- function(coordinateDF) {
  
  for (i in 1:nrow(coordinateDF)) {
    # For rows that correspond to intergenic regions:
    x <- stri_match(regex="intergenic", str=coordinateDF$geneID[i])
    coordinate_direction <- coordinateDF$strand[i]
    
    # If the 3' end is occurring in an intergenic region:
    if (!is.na(x)) {
      
      genes <- strsplit(coordinateDF$geneID[i], split="-")[[1]]
      gene1 <- genes[1]
      gene2 <- strsplit(genes[2], split = " intergenic")[[1]][1]
      
      gene1_direction <- na.omit(gene_DF_curated[gene_DF_curated$geneID == gene1,"gene_strand"])[1]
      gene2_direction <- na.omit(gene_DF_curated[gene_DF_curated$geneID == gene2,"gene_strand"])[1]
      
      gene1_end <- na.omit(gene_DF_curated[gene_DF_curated$geneID == gene1,"end"])[1]
      gene2_start <- na.omit(gene_DF_curated[gene_DF_curated$geneID == gene2,"start"])[1]
      
      gene1_distance <- coordinateDF$coordinate[i] - gene1_end
      gene2_distance <- gene2_start - coordinateDF$coordinate[i]
      
      # If the coordinate is closer to gene1 than gene2...
      if (gene1_distance < gene2_distance) {
        # If the type has not already been annotated (i.e. is not a bidirectional terminator)
        if (is.na(coordinateDF$LocationAnnotation[i])) {
          # If gene1 & coordinate are in the same direction...
          if (gene1_direction == coordinateDF$strand[i]) {
            
            # If gene1 & coordinate are in the + direction (coordinate occurs after gene end)...
            if (coordinateDF$strand[i] == "+") {
              # ... annotate as an Terminator
              coordinateDF$LocationAnnotation[i] <- "Terminator"
              
              # If the gene and coordinate are in the - direction (coordinate occurs before gene start)...
            }else{
              # .... annotate as an attenuator
              coordinateDF$LocationAnnotation[i] <- "Attenuator"
              
            }
            # If gene1 & coordinate are in opposite directions:
          }else{
            
            # If gene1 is + and coordinate is "-"...
            if (gene1_direction == "+") {
              
              # If gene2 is going in the "-" direction...
              if(gene2_direction == "-"){
                # ... then the coordinate is an Terminator for gene2
                coordinateDF$LocationAnnotation[i] <- "Terminator"
                # If gene2 is going in the + direction...
              }else{
                # ... then the coordinate is terminating antisense to the 5' end of gene2
                coordinateDF$LocationAnnotation[i] <- "Antisense"
              }
              # If gene1 is - and coordinate is +...
            }else{
              # ... then the coordinate is antisense to the 3' end of gene1
              coordinateDF$LocationAnnotation[i] <- "Antisense"
            }
          }
        }
      }else{
        # If the type has not already been annotated (i.e. is a bidirectional terminator)
        if (is.na(coordinateDF$LocationAnnotation[i])) {
          # If gene2 & coordinate are in the same direction...
          if (gene2_direction == coordinateDF$strand[i]) {
            
            # If the gene & coordinate are in the + direction (coordinate occurs before gene start)...
            if (coordinateDF$strand[i] == "+") {
              # ... annotate as an attenuator
              coordinateDF$LocationAnnotation[i] <- "Attenuator"
              # If the gene and coordinate are in the - direction (coordinate occurs after gene end)...
            }else{
              # .... annotate as an attenuator
              coordinateDF$LocationAnnotation[i] <- "Terminator"
            }
            # If gene2 & coordinate are in opposite directions:
          }else{
            
            # If gene2 is - and coordinate is +...
            if (gene2_direction == "-") {
              
              # If gene1 is going in the "+" direction...
              if(gene1_direction == "+"){
                # ... then the coordinate is an Terminator for gene
                coordinateDF$LocationAnnotation[i] <- "Terminator"
                # If gene1 is going in the - direction...
              }else{
                # ... then the coordinate is terminating antisense to the 5' end of gene2
                coordinateDF$LocationAnnotation[i] <- "Antisense"
              }
              # If gene2 is + and coordinate is -...
            }else{
              # ... then the coordinate is antisense to the 5' end of gene2
              coordinateDF$LocationAnnotation[i] <- "Antisense"
            }
          }
        }
      }
      # If the coordinate occurs within a gene body...
    }else{
      
      # If the coordinate and gene are in the same direction...
      if (coordinateDF$gene_strand[i] == coordinateDF$strand[i]) {
        # ... annotate as an Pause
        coordinateDF$LocationAnnotation[i] <- "Pause"
        # If the coordinate and gene are in opposite directions..
      }else{
        # ... annotate as antisense
        coordinateDF$LocationAnnotation[i] <- "Antisense"
      }
    }
  }
  
  return(coordinateDF)
}

fullyAnnotated <- annotateGenesLocations(addGenes)

## Identify coordinates where >2 coords per 10bp genomic region are present ------

geneIDs <- unique(fullyAnnotated$geneID)
num_IDs <- c()
multiples <- c()

for (i in 1:length(geneIDs)) {
  
  geneDF <- fullyAnnotated[fullyAnnotated$geneID == geneIDs[i],]
  num_IDs <- c(num_IDs, nrow(geneDF))
  
  if (nrow(geneDF) > 2) {
    
    multiples <- c(multiples,geneIDs[i])
  }
}

# Based on manual inspection of the geneIDs in the multiples variable
# If more than 2 coordinates called in a 10bp genomic window,
# only kept the 2 with the highest baseMeans as estimated in DESeq2
coordIDs_toRemove <- c("6115288-",
                       "6115292-",
                       "6115294-",
                       "8292065-",
                       "8292094-",
                       "8292101-",
                       "8292107-",
                       "8292109-",
                       "8292133-",
                       "8292135-",
                       "8979884+",
                       "8478754+",
                       "4937048+",
                       "6004544+",
                       "6113503-",
                       "6554036+",
                       "6554041+",
                       "7461339+",
                       "7461340+",
                       "7623125-",
                       "8742807-",
                       "8742808-",
                       "8742809-",
                       "8742811-",
                       "8742813-",
                       "8742814-",
                       "8742815-",
                       "8742816-",
                       "8980102+",
                       "8980106+",
                       "8980110+",
                       "9029476-",
                       "9029488-",
                       "9029490-",
                       "9029494-",
                       "9029496-",
                       "9029497-",
                       "9029498-",
                       "9029499-",
                       "9029503-",
                       "9029507-",
                       "9029513-")

fullyAnnotated$Removed <- "False"
fullyAnnotated[fullyAnnotated$coordID %in% coordIDs_toRemove,
                      'Removed'] <- "True"

## Prepare gene IDs to add vulnerability and PATRIC calls ----------------------

fullyAnnotated$upstream_geneID_noAS <- fullyAnnotated$upstream_geneID
fullyAnnotated$downstream_geneID_noAS <- fullyAnnotated$downstream_geneID

for (i in 1:nrow(fullyAnnotated)) {
  
  if (grepl("_as", fullyAnnotated$upstream_geneID[i])) {
    
    split_geneID <- strsplit(fullyAnnotated$upstream_geneID[i],
                             split = "_as")[[1]][1]
    fullyAnnotated$upstream_geneID_noAS[i] <- split_geneID
  }
  
  if (grepl("_as", fullyAnnotated$downstream_geneID[i])) {
    
    split_geneID <- strsplit(fullyAnnotated$downstream_geneID[i],
                             split = "_as")[[1]][1]
    fullyAnnotated$downstream_geneID_noAS[i] <- split_geneID
  }
  
}

## Add vulnerability calls -----------------------------------------------------

vulnerability <- read.csv("genome_files_misc/vulnerability_calls_CRISPRi.csv")
vulnerability_curated <- vulnerability[,c("locus_tag","tnseq_ess",
                                          "crispr_ess")]
colnames(vulnerability_curated) <- c("geneID_noAS","tnseq_ess",
                                     "crispr_ess")
upstream_vulnerability <- vulnerability_curated
colnames(upstream_vulnerability) <- paste("upstream_", 
                                          colnames(upstream_vulnerability),
                                          sep='')
downstream_vulnerability <- vulnerability_curated
colnames(downstream_vulnerability) <- paste("downstream_", 
                                          colnames(downstream_vulnerability),
                                          sep='')

add_upstream_vuln <- left_join(fullyAnnotated, upstream_vulnerability)
add_downstream_vuln <- left_join(add_upstream_vuln, downstream_vulnerability)

## Add PATRIC and KEGG annotations ---------------------------------------------

kegg_all <- read.csv('genome_files_misc/KEGG_Mtb_annotations.csv')
kegg_toKeep <- kegg_all[,c("gene","Category")]
colnames(kegg_toKeep) <- c("geneID_noAS", 'kegg_category')

upstream_kegg <- kegg_toKeep
colnames(upstream_kegg) <- paste("upstream_",
                                 colnames(upstream_kegg),
                                 sep='')

downstream_kegg <- kegg_toKeep
colnames(downstream_kegg) <- paste("downstream_",
                                 colnames(downstream_kegg),
                                 sep='')

add_upstream_kegg <- left_join(add_downstream_vuln, upstream_kegg)
add_downstream_kegg <- left_join(add_upstream_kegg, downstream_kegg)


patric_all <- read.csv('genome_files_misc/PATRIC_Mtb_annotations.csv')
patric_toKeep <- patric_all[,c("gene","Superclass_x","Class_x","Subclass_x")]
colnames(patric_toKeep) <- c("geneID_noAS","patric_superclass",
                             'patric_class','patric_subclass')
upstream_patric <- patric_toKeep
colnames(upstream_patric) <- paste("upstream_",
                                 colnames(upstream_patric),
                                 sep='')

downstream_patric <- patric_toKeep
colnames(downstream_patric) <- paste("downstream_",
                                   colnames(downstream_patric),
                                   sep='')
add_upstream_patric <- left_join(add_downstream_kegg, upstream_patric)
add_downstream_patric <- left_join(add_upstream_patric, downstream_patric)

add_downstream_patric[is.na(add_downstream_patric)] <- 'None found'

## Collapse KEGG and PATRIC terms into one row ---------------------------------

add_downstream_patric$upstream_kegg_collapse <- 0
add_downstream_patric$upstream_patric_superCollapse <- 0
add_downstream_patric$upstream_patric_classCollapse <- 0
add_downstream_patric$upstream_patric_subCollapse <- 0

add_downstream_patric$downstream_kegg_collapse <- 0
add_downstream_patric$downstream_patric_superCollapse <- 0
add_downstream_patric$downstream_patric_classCollapse <- 0
add_downstream_patric$downstream_patric_subCollapse <- 0

for (i in 1:nrow(add_downstream_patric)) {
  coord_DF <- add_downstream_patric[add_downstream_patric$coordID == add_downstream_patric$coordID[i],]
  
  upstream_kegg_categories <- unique(coord_DF$upstream_kegg_category)
  upstream_kegg_categories <- upstream_kegg_categories[!upstream_kegg_categories == 'None found']
  
  upstream_patric_superclasses <- unique(coord_DF$upstream_patric_superclass)
  upstream_patric_superclasses <- upstream_patric_superclasses[!upstream_patric_superclasses == 'None found']
  
  upstream_patric_classes <- unique(coord_DF$upstream_patric_class)
  upstream_patric_classes <- upstream_patric_classes[!upstream_patric_classes == 'None found']
  
  upstream_patric_subclasses <- unique(coord_DF$upstream_patric_subclass)
  upstream_patric_subclasses <- upstream_patric_subclasses[!upstream_patric_subclasses == 'None found']
  
  add_downstream_patric$upstream_kegg_collapse[i] <- paste0(upstream_kegg_categories, collapse = ', ')
  add_downstream_patric$upstream_patric_superCollapse[i] <- paste0(upstream_patric_superclasses, collapse = ', ')
  add_downstream_patric$upstream_patric_classCollapse[i] <- paste0(upstream_patric_classes, collapse = ', ')
  add_downstream_patric$upstream_patric_subCollapse[i] <- paste0(upstream_patric_subclasses, collapse = ', ')
  
  downstream_kegg_categories <- unique(coord_DF$downstream_kegg_category)
  downstream_kegg_categories <- downstream_kegg_categories[!downstream_kegg_categories == 'None found']
  
  downstream_patric_superclasses <- unique(coord_DF$downstream_patric_superclass)
  downstream_patric_superclasses <- downstream_patric_superclasses[!downstream_patric_superclasses == 'None found']
  
  downstream_patric_classes <- unique(coord_DF$downstream_patric_class)
  downstream_patric_classes <- downstream_patric_classes[!downstream_patric_classes == 'None found']
  
  downstream_patric_subclasses <- unique(coord_DF$downstream_patric_subclass)
  downstream_patric_subclasses <- downstream_patric_subclasses[!downstream_patric_subclasses == 'None found']
  
  add_downstream_patric$downstream_kegg_collapse[i] <- paste0(downstream_kegg_categories, collapse = ', ')
  add_downstream_patric$downstream_patric_superCollapse[i] <- paste0(downstream_patric_superclasses, collapse = ', ')
  add_downstream_patric$downstream_patric_classCollapse[i] <- paste0(downstream_patric_classes, collapse = ', ')
  add_downstream_patric$downstream_patric_subCollapse[i] <- paste0(downstream_patric_subclasses, collapse = ', ')
  
}

old_KEGG_PATRIC_cols <- c("upstream_kegg_category", "upstream_patric_superclass",
                          "upstream_patric_class", "upstream_patric_subclass",
                          "downstream_kegg_category", "downstream_patric_superclass",
                          "downstream_patric_class", "downstream_patric_subclass")

no_old_cols <- add_downstream_patric[,!(colnames(add_downstream_patric) %in% old_KEGG_PATRIC_cols)]

collapsed_DF <- unique(no_old_cols)

## Adjust coordinate to match Mtb Broad ref genome -----------------------------

len_eco_spike_genome <- 4641652
collapsed_DF$norm_coordinate <- collapsed_DF$coordinate - len_eco_spike_genome

## Add in cluster descriptions -------------------------------------------------

collapsed_DF$cluster_description <- 0
collapsed_DF[collapsed_DF$cluster == 'NE','cluster_description'] <- 'No Nus effect'
collapsed_DF[collapsed_DF$cluster == 'A','cluster_description'] <- 'Nus suppressed'
collapsed_DF[collapsed_DF$cluster == 'B','cluster_description'] <- 'Solely NusG stimulated'
collapsed_DF[collapsed_DF$cluster == 'C','cluster_description'] <- 'Solely NusA stimulated'
collapsed_DF[collapsed_DF$cluster == 'D','cluster_description'] <- 'NusA + NusG stimulated'

## Calculate downstream DNA AT % -----------------------------------------------

collapsed_DF$AT_count_downstream <- collapsed_DF$U_count_downstreamRegion + collapsed_DF$A_count_downstreamRegion
length_downstream_region <- length(strsplit(collapsed_DF$downstream_region, split = '')[[1]])
collapsed_DF$AT_percent_downstream <- (collapsed_DF$AT_count_downstream / length_downstream_region) * 100

## Look for elemental pause motif ----------------------------------------------

pause_regex <- 'gg........[ct]g'
collapsed_DF$Elemental_pause_region <- 'Not transcribed in CFG'
collapsed_DF$Elemental_pause <- 'Not transcribed in CFG'

for (i in 1:nrow(collapsed_DF)) {
  
  if (collapsed_DF$Nontemplate_strand[i] != 'Not transcribed in CFG') {
    
    pause_region <- strsplit(collapsed_DF$Nontemplate_strand[i], split='')[[1]][42:53]
    collapsed_DF$Elemental_pause_region[i] <- paste0(pause_region,collapse='')
    collapsed_DF$Elemental_pause[i] <- grepl(pause_regex,
                                             paste0(tolower(pause_region),
                                                    collapse=''))
  }
  
  
}

## Final supplementary table organization --------------------------------------

coordinate_info <- c("norm_coordinate", "strand",
                     "upstream_geneID","upstream_geneName",
                     "downstream_geneID","downstream_geneName",
                     "GeneLocation","LocationAnnotation","Elemental_pause",
                     "Removed","coordID")
Nus_effect_info <- c("cluster","cluster_description",
                     "NusA_log2FC","NusG_log2FC","NusAG_log2FC",
                     "NusA_padj","NusG_padj","NusAG_padj")
gene_class_info <- c("upstream_geneDesc","upstream_crispr_ess","upstream_tnseq_ess",
                     "downstream_geneDesc","downstream_crispr_ess","downstream_tnseq_ess",
                     "upstream_kegg_collapse", "upstream_patric_superCollapse",
                     "upstream_patric_classCollapse","upstream_patric_subCollapse",
                     "downstream_kegg_collapse", "downstream_patric_superCollapse",
                     "downstream_patric_classCollapse","downstream_patric_subCollapse")
sequence_info <- c("RNAfold_hairpin_region_44_structure","RNAfold_hairpin_region_44_deltaG",
                   'RNAfold_hairpin_region_44_spacerLength',
                   'AT_percent_downstream','RNAfold_hairpin_region_44_stemLength',
                   "RNAfold_hairpin_region_44_loopLength",
                   'U_count_tractRegion','Nontemplate_strand')

col_order <- c(coordinate_info, Nus_effect_info,
               gene_class_info, sequence_info)
cleaned_DF <- collapsed_DF[,col_order]
sorted_DF <- cleaned_DF %>%
  arrange(GeneLocation, cluster)

sorted_DF[sorted_DF == 'None found'] <- ''
sorted_DF[sorted_DF == 'None'] <- ''
sorted_DF[sorted_DF$upstream_geneName == '-', 'upstream_geneName'] <- ''
sorted_DF[sorted_DF$upstream_geneName == '- (antisense)', 'upstream_geneName'] <- ''
sorted_DF[sorted_DF$downstream_geneName == '-', 'downstream_geneName'] <- ''
sorted_DF[sorted_DF$downstream_geneName == '- (antisense)', 'downstream_geneName'] <- ''

new_colnames <- c("TTS coordinate","TTS strand",
                  "Upstream gene ID", "Upstream gene name","Downstream gene ID", "Downstream gene name",
                  "TTS location","Predicted regulatory class",'Elemental pause match',
                  "Removed","CFG TTS ID",
                  "Cluster","Cluster description",
                  "NusA log2FC", "NusG log2FC","NusA + NusG log2FC",
                  "NusA p-adj.","NusG p-adj.","NusA + NusG p-adj.",
                  "Upstream gene desc.","Upstream gene CRISPRi call","Upstream gene Tn-seq call",
                  "Downstream gene desc.","Downstream gene CRISPRi call", "Downstream gene Tn-seq call",
                  "Upstream gene KEGG class", "Upstream gene PATRIC superclass",
                  "Upstream gene PATRIC class","Upstream gene PATRIC subclass",
                  "Downstream gene KEGG class","Downstream gene PATRIC superclass",
                  "Downstream gene PATRIC class","Downstream gene PATRIC subclass",
                  "Predicted hairpin structure (-44 to -1)","Hairpin folding free energy change",
                  'Hairpin spacer length','dsDNA AT %','Hairpin stem length','Hairpin loop length',
                  'Number of Us in U-tract (-8 to -1)',
                  "Non-template strand (-50 to +38)")

colnames(sorted_DF) <- new_colnames

no_random <- sorted_DF[sorted_DF$Cluster != 'random',]

write.csv(no_random, '3enrich_NusAG/termination_figures/supp_table6_terminators_genes.csv')

## Remove coordinates where >2 coords per 10bp genomic region are present ------

analysis_DF <- collapsed_DF[!(collapsed_DF$coordID %in% coordIDs_toRemove),]

## Sequence logos for NusA vs. NusG stim ---------------------------------------

library(tidyverse)

cs = make_col_scheme(chars=c('A', 'T', 'C', 'G'),
                     cols=c('#C70000','#008000',
                            '#0000C7','#FFAE00'))

location_list <- c("inter", "intra")

stim_list <- c("Yes", "No", 'random')

# Calculate background model with fasta-get-markov
# Load in (code from Mike Wolfe)
prior <- read_delim('genome_files_misc/Mtb_background.tsv',
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

## Rest of code
for (i in 1:length(location_list)) {
  
  location_DF <- analysis_DF %>%
    filter(GeneLocation == location_list[i])
  
  for (j in 1:length(stim_list)) {
    
    location_NusGstim_DF <- location_DF %>%
      filter(NusG_stim == stim_list[j])
    
    print(paste("Location: ",location_list[i],sep=''))
    print(paste("NusG_stim: ",stim_list[j],sep=''))
    print(paste("Num rows: ",nrow(location_NusGstim_DF),sep=''))
    print(" ")
    
    NusG_relEntropy <- relative_entropy_calc(location_NusGstim_DF$Nontemplate_strand,
                                             paste('3enrich_NusAG/termination_figures/extDat_fig8b_',
                                                   location_list[i],'_NusGstim_',
                                                   stim_list[i],sep=''))
    
    p <- ggseqlogo(NusG_relEntropy, method = "custom", seq_type = 'dna',
                   col_scheme = cs)
    NusG_logo <- p + theme_classic() + labs(y = "Relative Entropy") + ylim(0,2.2)
    NusG_logo$scales$scales[[1]] <- scale_x_continuous(breaks = seq(1, 90, by = 1), 
                                                       labels= c(seq(-52, 37, by=1)))
    
    ggsave(filename = paste('3enrich_NusAG/termination_figures/extDat_fig8b_',location_list[i],
                            '_NusGstim_',stim_list[j],'_logo_',
                            counts,'counts_padj',sig_value,'.png',sep=''),
           plot = NusG_logo,
           width = 500,
           height = 100,
           units = 'mm',
           dpi = 300)
    
    location_NusAstim_DF <- location_DF %>%
      filter(NusA_stim == stim_list[j])
    
    print(paste("Location: ",location_list[i],sep=''))
    print(paste("NusA_stim: ",stim_list[j],sep=''))
    print(paste("Num rows: ",nrow(location_NusAstim_DF),sep=''))
    print(" ")
    
    NusA_relEntropy <- relative_entropy_calc(location_NusAstim_DF$Nontemplate_strand,
                                             paste('3enrich_NusAG/termination_figures/extDat_fig8b_',
                                                   location_list[i],'_NusAstim_',
                                                   stim_list[i],sep=''))
    
    p <- ggseqlogo(NusA_relEntropy, method = "custom", seq_type = 'dna',
                   col_scheme = cs)
    NusA_logo <- p + theme_classic() + labs(y = "Relative Entropy") + ylim(0,2.2)
    NusA_logo$scales$scales[[1]] <- scale_x_continuous(breaks = seq(1, 90, by = 1), 
                                                       labels= c(seq(-52, 37, by=1)))
    
    ggsave(filename = paste('3enrich_NusAG/termination_figures/extDat_fig8b_',location_list[i],
                            '_NusAstim_',stim_list[j],'_logo_',
                            counts,'counts_padj',sig_value,'.png',sep=''),
           plot = NusA_logo,
           width = 500,
           height = 100,
           units = 'mm',
           dpi = 300)
  }
}

## Export intergenic terminators for further sequence analysis -----------------

intergenic_allCols <- analysis_DF %>%
  filter(GeneLocation == 'inter')

write.csv(intergenic_allCols, '3enrich_NusAG/termination_figures/intergenic_terminators.csv')

## End -------------------------------------------------------------------------

