################################################################################
## This script calculates the fraction of intra-chromosomal and inter-chromosomal
## interactions across different Mutual Information (MI) thresholds, using both
## accumulated and chunked binning strategies. 
##
## Inputs:
##   - MI matrix (gene x gene, with MI values)
##   - Gene annotation file (from BioMart)
##
## Outputs:
##   - Fraction tables for intra/inter interactions in 1k-sized bins (TSV)
##   - Fraction tables for intra/inter interactions in log-sized bins (TSV)
################################################################################

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library(vroom)
library(magrittr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggthemes)
library(furrr)

Sys.setenv(VROOM_CONNECTION_SIZE=500072)

MCCORES <- snakemake@threads[[1]]
future::plan(future::multisession, workers = MCCORES)

cat("Loading annotation\n")
load(snakemake@params[["annot"]])

cat("Reading matrix\n")
MImatrix <- vroom::vroom(snakemake@input[["mi_matrix"]])
COND <- snakemake@params[["cond"]]

MImatrix <- MImatrix %>% 
  as.matrix()

MImatrix[is.na(MImatrix)] <- 0
MImatrix[lower.tri(MImatrix, diag = T)] <- NA

genes <- colnames(MImatrix)

cat("Annotating interactions\n")

#Extratct upper triangle indices and MI values without pivot_longer
upper_indices <- which(
       upper.tri(as.matrix(MImatrix[, -ncol(MImatrix)])),
       arr.ind = TRUE
)

# Create a tibble with source, target, and MI values directly
mi_vals <- tibble(
       source = genes[upper_indices[, 1]],
       target = genes[upper_indices[, 2]],
       mi = MImatrix[upper_indices]
       ) 
cat("Done building tibble\n")
rm(MImatrix)
cat("Done deleting matrix\n")
mi_vals <- mi_vals |>
  arrange(desc(mi))    
cat("Done arrenging mi vals\n")

mi_vals <- mi_vals |>
  left_join(annot %>% dplyr::select(ensembl_id, chr), 
                                 by = c("source" = "ensembl_id"), multiple = "first") %>% 
  rename("source_chr" = "chr") %>%
  left_join(annot %>% dplyr::select(ensembl_id, chr), 
                                         by = c("target" = "ensembl_id"),  multiple = "first") %>% 
  rename("target_chr" = "chr") %>%
  mutate(interaction = ifelse(source_chr == target_chr, "intra", "inter"),
         nrow = row_number())

cat("Creating onek and log chunks\n")
### bins by a thousand interactions

onek_chunks <- mi_vals %>%  
  mutate(bin = floor((nrow-1)/1000))

cat("Getting chunks\n")
onek_chunks <- onek_chunks |>
  group_by(bin) %>%
  count(interaction) %>% 
  pivot_wider(id_cols = bin, names_from = interaction, values_from = n, 
              values_fill = 0)

cat("Saving chunks\n")
onek_chunks <- onek_chunks |>
  mutate(inter_fraction = round(inter/(intra+inter), 4),
         intra_fraction = round(intra/(intra+inter), 4),
         cond = COND) %>%
  ungroup() %T>% 
  vroom_write(file = snakemake@output[["onek_chunks"]])

cat("Saving bins\n")
bins <- unique(onek_chunks$bin)
names(bins) <- unique(onek_chunks$bin)

furrr::future_map_dfr(bins, ~ onek_chunks %>% 
                        filter(bin <= .x) %>% 
                        select(intra, inter) %>% colSums(),
                      .id = "bin") %>%
  mutate(bin = as.numeric(bin),
         inter_fraction = round(inter/(intra+inter), 4),
         intra_fraction = round(intra/(intra+inter), 4),
         cond = COND) %>%
  vroom_write(file = snakemake@output[["onek_bins"]])

rm(onek_chunks)

cat("Getting log chunks\n")
log_chunks <- mi_vals %>%
  mutate(bin_size = 10^floor(log10(nrow)))

cat("Deleting mi_vals\n")

rm(mi_vals)

log_chunks <- log_chunks |>
  	mutate(bin_size = ifelse(bin_size < 1000, 1000, bin_size))
log_chunks <- log_chunks |>
         mutate(bin = (floor((nrow-1)/bin_size)+1)*bin_size) 

cat("Saving log chunks\n")
log_chunks <- log_chunks |>
  group_by(bin) %>%
  count(interaction) %>% 
  pivot_wider(id_cols = bin, names_from = interaction, values_from = n, 
              values_fill = 0) %>%
  mutate(inter_fraction = round(inter/(intra+inter), 4),
         intra_fraction = round(intra/(intra+inter), 4),
         cond = COND) %>%
  ungroup() %T>%
  vroom_write(file = snakemake@output[["log_chunks"]])

bins <- unique(log_chunks$bin)
names(bins) <- unique(log_chunks$bin)

cat("Saving log bins\n")
furrr::future_map_dfr(bins, ~ log_chunks %>% 
                        filter(bin <= .x) %>% 
                        select(intra, inter) %>% colSums(),
                      .id = "bin") %>%
  mutate(bin = as.numeric(bin),
         inter_fraction = round(inter/(intra+inter), 4),
         intra_fraction = round(intra/(intra+inter), 4),
         cond = COND) %>%
  vroom_write(file = snakemake@output[["log_bins"]])
