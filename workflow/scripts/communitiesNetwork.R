## #############################################################
## This file builds a network of communities per condition
## by counting the links between communities
## It requires interactions and vertices files from networkTables.R 
## and membership file from communities.R as input
################################################################

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library(readr)
library(dplyr)
library(tidyr)

COND <- snakemake@params[["cond"]]
TISSUE <- snakemake@params[["tissue"]]
cat("Reading files\n")
interactions <- read_tsv(snakemake@input[["interactions"]])
vertices <-  read_tsv(snakemake@input[["vertices"]],  col_types = cols(chr=col_character())) %>%
  select(ensembl_id, symbol) 
membership <- read_tsv(snakemake@input[["membership"]])

vertices <- vertices %>% inner_join(membership, by ="ensembl_id")

cat("Getting interactions among communities\n")
colnames(vertices) <- c("source", "source_symbol", "source_comm")
interactions <- interactions %>% inner_join(vertices, by="source")
colnames(vertices) <- c("target", "target_symbol",  "target_comm")
interactions <- interactions %>% inner_join(vertices, by="target")

interactions <- interactions %>% select(source_comm, target_comm) %>%
  mutate(id = paste0(pmin(source_comm, target_comm), "_", pmax(source_comm, target_comm)))

interactions_by_communities <-interactions %>% group_by(id) %>% tally()
interactions_by_communities <- interactions_by_communities %>% 
  separate(id, into = c("source", "target")) %>%
  mutate(source = paste0(TISSUE, "_", source), target = paste0(TISSUE, "_", target))

cat("Saving files\n")
write_tsv(interactions_by_communities, file = snakemake@output[[1]])  
