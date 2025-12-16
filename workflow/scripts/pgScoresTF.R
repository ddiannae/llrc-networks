library(igraph)
library(vroom)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(stringr)

TISSUE <- snakemake@params[["tissue"]]
TISSUE <- str_to_title(TISSUE)

tfs <- vroom(snakemake@params[["tfs"]], col_names = FALSE, delim = ",") %>%
  pull(X1)

all_pgscores <- lapply(c("normal", "cancer"), function(cond) {
  interactions <- vroom(snakemake@input[[paste0(cond, "_interactions")]]) %>%
    dplyr::rename("from" = "source", "to" = "target") |>
    dplyr::select(from, to, everything())

  vertices <- vroom(snakemake@input[[paste0(cond, "_vertices")]]) %>%
    dplyr::select(name, everything())

  net <- graph_from_data_frame(interactions, directed = F, vertices = vertices)
  pgscores <- page_rank(net)
  pgscores <- tibble(name = names(pgscores$vector), pg_score = pgscores$vector) %>%
    inner_join(
      vertices |>
        select(name, symbol),
      by = "name"
    ) |>
    mutate(is_tf = ifelse(name %in% tfs, TRUE, FALSE),
      cond = cond)
  return(pgscores)
})

all_pgscores <- bind_rows(all_pgscores) |>
  mutate(cond = str_to_title(cond))
cond_colors <- c("Normal" = "#e3a098", "Cancer" = "#a32e27")
ggplot(all_pgscores, aes(y = pg_score, x = is_tf)) +
  geom_boxplot(aes(color = cond)) +
  geom_jitter(aes(color = cond), width = 0.2, alpha = 0.3) +
  stat_compare_means(
    method = "wilcox.test",
    label.y.npc = 0.9,
    label = "p.format"
  ) +
  labs(x = "", y = "Page Rank score") +
  scale_x_discrete(
    labels = as_labeller(c(
      `FALSE` = "Non Transcription\nFactors",
      `TRUE` = "Transcription\nFactors"
    ))
  ) +
  scale_color_manual(name = "Condition", values = cond_colors) +
  theme_minimal(base_size = 8) +
  facet_wrap(
    ~cond
  ) +
  ggtitle(TISSUE)

ggsave(filename = snakemake@output[[1]],
       width = 5, height = 3, units = "in")
