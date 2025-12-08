log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library(vroom)
library(dplyr)
library(tidyr)
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)

cat("Reading files\n")

mem <- vroom(file=snakemake@input[["membership"]])
gene_universe <- vroom(file=snakemake@input[["universe"]], col_select = "name")

gene_universe <- gene_universe |>
  pull(name)

all_enrichments <- lapply(X = unique(mem$community),
                          FUN = function(com){
			
		cat("Working with community: ", com, "\n")

    gene_list <- mem |>
      filter(community == com) |>
      pull(name)
    
    if(length(gene_list) >= 5) {
      
      terms <- enrichGO(gene          = gene_list,
                      universe      = gene_universe,
                      OrgDb         = org.Hs.eg.db,
                      keyType       = "ENSEMBL",
                      ont           = "BP",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.005,
                      qvalueCutoff  = 0.01,
                      minGSSize     = 10,
                      readable      = FALSE)
      
      tryCatch({
        bp <- pairwise_termsim(terms)
        bp2 <- simplify(bp, cutoff=0.7, by="p.adjust", select_fun=min)
        simple_results <- bp2@result
        if(nrow(simple_results) > 0) {
          simple_results$commun <- com 
        }
        return(simple_results)
      }, error = function(cond) {
        return(NULL)
      })
     
    }
    return(NULL)
})

bind_rows(all_enrichments) |>
  janitor::clean_names() |>
  vroom_write(snakemake@output[[1]])