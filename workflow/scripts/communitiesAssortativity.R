################################################################################
## This script calculates assortativity values for network communities based on
## chromosomal location and gene expression. For each community, it computes
## chromosomal assortativity, and for cancer networks, also computes expression
## and differential expression assortativity. 
##
## Inputs:
##   - Interactions table (from networkTables.R)
##   - Vertices table (from networkTables.R)
##   - Expression table (with log2FC and baseMean)
##   - Community membership table
##   - Condition (normal or cancer; from Snakemake params)
##
## Outputs:
##   - Chromosomal assortativity table (TSV)
##   - Expression assortativity table (TSV, cancer only)
##   - Differential expression summary table (TSV, cancer only)
################################################################################

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library(readr)
library(dplyr)
library(igraph)

getAssortativityByAttr <- function(gr, n_attr) {
  vs <- as_data_frame(gr, what = "vertices")
  comm_assort <- lapply(unique(vs$community), function(com){
    vc <- vs %>% filter(community == com) %>% 
      select(name) %>% unlist()
    com_graph <- induced_subgraph(gr, vc)
    v_attr <- as_data_frame(com_graph, what = "vertices")[, n_attr]
    newman_a <- assortativity_nominal(com_graph,
                                      types = as.numeric(as.factor(v_attr)), 
                                      directed = T)
    A <- as.matrix(get.adjacency(com_graph))
    l_times_A <- outer(v_attr, v_attr, `==`) * A
    l_dif_A <- outer(v_attr, v_attr, `!=`) * A
    bynode <- rowSums(l_times_A) / rowSums(A)
    bynode_frac <- mean(bynode)
    total_frac <- sum(l_times_A)/sum(A)
    dif_fraction <- (sum(l_times_A)-sum(l_dif_A))/sum(A)
    return(list(community_id = com, newman = newman_a, bynode = bynode_frac, 
                totalfrac = total_frac, diffraction=dif_fraction))
  })
  comm_assort <- bind_rows(comm_assort)
  return(comm_assort)
}
 
COND <- snakemake@params[["cond"]]
interactions <- read_tsv(snakemake@input[["interactions"]])
vertices <-  read_tsv(snakemake@input[["vertices"]])
expression <- read_tsv(snakemake@input[["expression"]])
membership <- read_tsv(snakemake@input[["membership"]])

vertices <- vertices %>% 
  left_join(membership, by = "ensembl_id") %>%
  left_join(expression, by = c("ensembl_id" = "ensembl_id"))

if(COND == "cancer") {
  vertices <- vertices %>% mutate(exp = case_when(log2FC > 0 ~ "up",
                                                  log2FC < 0 ~ "down",
                                                  log2FC == 0 ~ "NA"),
                                  exp = as.factor(exp),
                                  diff_exp = case_when(log2FC >= 1 ~ "up",
                                                  log2FC <= -1 ~ "down",
                                                  TRUE ~ "no"),
                                  diff_exp = as.factor(diff_exp),)
}

g <- graph_from_data_frame(interactions, vertices = vertices,
                           directed = FALSE)  

chr_assortativity <- getAssortativityByAttr(g, "chr")
write_tsv(chr_assortativity, file = snakemake@output[["chr_assortativity"]])

## No expression assortativity for normal tissue
if(COND == "cancer") {
  exp_summary <- vertices %>% group_by(community) %>%
    summarise(mean_log_fc = mean(log2FC),
              mean_avg_exp  =  mean(baseMean))
  
  diff_exp_summary <- vertices %>% group_by(community, diff_exp) %>% 
    summarise(n = n(),  mean_log_fc = mean(log2FC)) %>% 
    mutate(freq = n/sum(n))
  
  exp_assortativity <- getAssortativityByAttr(g, "exp")
  
  exp_summary <- exp_assortativity %>%
    inner_join(exp_summary,  by = c("community_id" = "community")) 
  
  write_tsv(exp_summary, file = snakemake@output[["expr_assortativity"]])  
  write_tsv(diff_exp_summary, file = snakemake@output[["diff_expr_summary"]])  
} else {
  write_lines(c("NA"), file = snakemake@output[["expr_assortativity"]])
  write_lines(c("NA"), file = snakemake@output[["diff_expr_summary"]])  
}
