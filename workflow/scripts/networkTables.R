################################################################################
## This script generates interaction and vertex tables for gene networks based on
## the top MI (Mutual Information) cutoff. It annotates interactions using a 
## BioMart-derived annotation file. The script expects as input:
##   - An MI matrix (gene x gene, with MI values)
##   - A gene annotation file (from BioMart)
##   - Parameters for MI cutoff and condition (from Snakemake)
##
## Output:
##   - Interactions table: Top MI gene pairs with annotations and distances
##   - Vertices table: Unique genes with annotation details
################################################################################

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library(vroom)
library(tidyr)
library(dplyr)

CUTOFF <- as.numeric(snakemake@params[["cutoff"]])
COND <- snakemake@params[["cond"]]

cat("Loading files\n")
load(snakemake@params[["annot"]])

MImatrix <- vroom::vroom(snakemake@input[["mi_matrix"]])
genes <- colnames(MImatrix)

cat("MI matrix with ", nrow(MImatrix), " rows and ", ncol(MImatrix), " columns loaded \n")
MImatrix <- MImatrix %>% 
  as.matrix()

MImatrix[is.na(MImatrix)] <- 0
MImatrix[lower.tri(MImatrix, diag = T)] <- NA

MImatrix <- as_tibble(MImatrix)
MImatrix$source <- genes

cat("Annotating interactions\n")
MIvals <- MImatrix %>% pivot_longer(cols = starts_with("ENSG"), 
                                     names_to = "target",
                                     values_to = "mi",
                                     values_drop_na = TRUE) %>% 
  filter(source != target) %>%
  arrange(desc(mi)) %>% 
  mutate(row_num = row_number()) %>% 
  filter(row_num <= CUTOFF)

cat("Merging annotations\n")
MIvals$cond <- COND

annot <- annot %>% select(ensembl_id, chr, start, end, gene_name)
colnames(annot) <-  c("source", "source_chr", "source_start", "source_end","source_name")
MIvals <- merge(annot, MIvals)
colnames(annot) <-  c("target", "target_chr", "target_start", "target_end",  "target_name")
MIvals <- merge(annot, MIvals)
MIvals <- MIvals %>% mutate(inter = if_else(source_chr == target_chr,  F, T), 
                      interaction_type = if_else(inter == T, "Inter", "Intra"),
                      distance = if_else(inter == F, as.integer(pmax(source_start, target_start) - 
                                           pmin(source_start, target_start)), as.integer(-1)))

targets <- MIvals %>% select(target, target_chr, target_start, target_end, target_name)
sources <- MIvals %>% select(source, source_chr, source_start, source_end, source_name)
colnames(targets) <- c("ensembl_id", "chr", "start", "end",  "symbol")
colnames(sources)  <-  c("ensembl_id", "chr", "start", "end", "symbol")
vertices <- bind_rows(targets, sources)
vertices <- vertices[!duplicated(vertices$ensembl_id), ] 

MIvals <- MIvals %>% 
  select(source, target, mi, distance, row_num, interaction_type, cond) %>%
  arrange(row_num) 

cat("Saving files\n")
vroom_write(vertices, file = snakemake@output[["vertices"]])
vroom_write(MIvals, file = snakemake@output[["interactions"]])
