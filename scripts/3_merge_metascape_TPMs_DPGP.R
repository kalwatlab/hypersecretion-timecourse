# script for merging the RNAseq TPMs and read counts with DPGP clustering results and metascape analysis just for genes signficant by edgeR

suppressPackageStartupMessages(library(readxl))
suppressPackageStartupMessages(library(tidyverse))

#load in multiple tables: metascape; TPMs from total RNAseq; DPGP cluster identities

# Need to rename the gene symbol header to be the same across all tables in order to merge them. 
# Relabel the clusters to indicate its from SW or Tg.

RNASeq_TPMs <- read.delim2("data/ILMN_989_Kalwat_IBRI_TotalRNASeq_TPMs.txt", header = T)
RNASeq_TPMs <- dplyr::rename(RNASeq_TPMs, "Gene.Symbol" = "GeneSymbol")


Metascape <- read_excel("data/metascape_result.xlsx", sheet = "Annotation") # goes with result_ID kalwat.rerun.all.twz_eh1i6
Metascape <- Metascape %>% 
  dplyr::rename(
  SW016789_sig = SW016789,
  Tg_sig = Thapsigargin
  )

Metascape$Gene <- NULL
colnames(Metascape)
Metascape[,c(13:35)] <- NULL # can drop most of the messy columns

# Metascape already has Gene.Symbol as the header for identifier
DPGPclusters_SW <- read.table("data/output_SW_optimal_clustering.txt", header = T)
DPGPclusters_SW <- dplyr::rename(DPGPclusters_SW, "Gene.Symbol" = "gene")
DPGPclusters_SW <- dplyr::rename(DPGPclusters_SW, "SW_DPGP_cluster" = "cluster")

DPGPclusters_Tg <- read.table("data/output_Tg_optimal_clustering.txt", header = T)
DPGPclusters_Tg <- dplyr::rename(DPGPclusters_Tg, "Gene.Symbol" = "gene")
DPGPclusters_Tg <- dplyr::rename(DPGPclusters_Tg, "Tg_DPGP_cluster" = "cluster")

DPGP_clusters <- full_join(DPGPclusters_SW,DPGPclusters_Tg, by = "Gene.Symbol")
DPGP_clusters <- DPGP_clusters[,c("Gene.Symbol","SW_DPGP_cluster","Tg_DPGP_cluster")] #reorder the columns to put gene symbol first

DPGP_clusters <- DPGP_clusters %>%
  mutate(SW_DPGP_cluster = ifelse(is.na(SW_DPGP_cluster), 0, SW_DPGP_cluster),
         Tg_DPGP_cluster = ifelse(is.na(Tg_DPGP_cluster), 0, Tg_DPGP_cluster))

# write this merged reformatted DPGP result file to the data import folder for downstream use
write.table(DPGP_clusters, "data/DPGP_clusters_SW_Tg.tsv", sep = '\t', row.names = F)

# merge tables
merged_list <- full_join(DPGP_clusters,Metascape, by = c("Gene.Symbol" = "original_id"))

# Some gene names were mismatched due to renaming by Metascape. Manually repair these.
which(merged_list$Gene.Symbol == 'BC003331')
merged_list[202,6] <- "Ord4"
which(merged_list$Gene.Symbol == 'BC030870')
merged_list[576,6] <- "Smim31"

merged_list[864,1] <- "Phf10"
merged_list[864,2] <- 3
merged_list[864,3] <- 1

which(merged_list$Gene.Symbol == 'Phf10')
merged_list <- merged_list[-99,] # remove the row with gene Phf10 because it was fixed above
which(merged_list$Gene.Symbol == 'LOC106740')
merged_list <- merged_list[-87,] # remove the LOC106740 gene, which is just a lncRNA at the PHF10 locus

View(merged_list)

merged_list_2 <- full_join(merged_list, RNASeq_TPMs, by = "Gene.Symbol")

write.table(merged_list_2, "output/2024_MIN6_SW_Tg_timecourse_DPGP_metascape_TPMs.tsv", sep = '\t', row.names = F)
