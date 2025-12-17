log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library(vroom)
library(dplyr)

log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

boots_files <- unlist(snakemake@input)

all_bins <- lapply(1:length(boots_files), function(i) {
  cat("Reading file ", boots_files[i], "\n")
  df <- vroom::vroom(paste0(boots_files[i]))
  df <- df |>
    select(bin, intra_fraction, cond) |>
    mutate(sample = i)
  return(df)
})

all_bins <- bind_rows(all_bins)
bin_summ <- all_bins |>
  group_by(bin, cond) |>
  summarise(mean_intra = mean(intra_fraction),
            sd_intra = sd(intra_fraction)) |>
  ungroup()

vroom::vroom_write(bin_summ, snakemake@output[[1]], delim = "\t")
