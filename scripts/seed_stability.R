#Significant term count per seed - diagnostic for 04's consensus approach, reads 04's output only

source(here::here("scripts", "_paths.R"))
library(ggplot2)

counts <- read.csv(file.path(DIR_TABLES, "seed_stability_counts.csv"))

consensus_size <- function(database, dataset) {
  file <- paste0("GSEA_consensus_", ifelse(database == "KEGG", "KEGG", "GO"), "_", dataset, ".csv")
  tab <- read.csv(file.path(DIR_TABLES, file))
  sum(tab$seed_fraction >= 0.8)
}
consensus <- unique(counts[, c("database", "dataset")])
consensus$n_consensus <- mapply(consensus_size, consensus$database, consensus$dataset)

ranges <- aggregate(n_significant ~ database + dataset, counts, function(x) paste0(min(x), "-", max(x)))
names(ranges)[3] <- "range_label"
label_y <- aggregate(n_significant ~ database, counts, max)
names(label_y)[2] <- "y_pos"
ranges <- merge(ranges, label_y, by = "database")
ranges$label <- paste0(ranges$dataset, ": ", ranges$range_label)
ranges <- ranges[order(ranges$database, ranges$dataset), ]
ranges$y_pos <- ave(ranges$y_pos, ranges$database, FUN = function(y) y * seq(1, 0.85, length.out = length(y)))

p <- ggplot(counts, aes(x = factor(seed), y = n_significant, color = dataset)) +
  geom_point(size = 2.5) +
  geom_line(aes(group = dataset), alpha = 0.5) +
  geom_hline(data = consensus, aes(yintercept = n_consensus, color = dataset),
             linetype = "dashed", linewidth = 0.6) +
  geom_text(data = ranges, aes(x = 1, y = y_pos, label = label, color = dataset),
            hjust = 0, size = 3.3, fontface = "bold", show.legend = FALSE) +
  facet_wrap(~ database, scales = "free_y") +
  labs(title = "Significant gene sets per seed (p.adjust < 0.05)",
       subtitle = "Dashed line = consensus (>=80% of seeds significant) - single seeds swing, consensus doesn't",
       x = "Seed", y = "Significant count", color = "Dataset") +
  theme_bw(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(DIR_FIGURES, "seed_stability_counts.png"), plot = p,
       width = 11, height = 6, dpi = 300, bg = "white")
