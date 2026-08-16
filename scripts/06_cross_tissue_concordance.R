#Cross-tissue logFC concordance - GSE56500 (Spinal Cord) & GSE68605 (Motor Cortex)
#Tests: Spearman correlation with permutation p | threshold sensitivity | shared-gene expression

source(here::here("scripts", "_paths.R"))
source(here::here("scripts", "_helpers.R"))
library(ggplot2)

#Load DEG results
top_56500 <- read.csv(file.path(DIR_TABLES, "ALS_vs_Control_all.csv"))
top_68605 <- read.csv(file.path(DIR_TABLES, "GSE68605_ALS_vs_Control_all.csv"))
sig_56500 <- read.csv(file.path(DIR_TABLES, "ALS_vs_Control_significant.csv"))
sig_68605 <- read.csv(file.path(DIR_TABLES, "GSE68605_ALS_vs_Control_significant.csv"))
overlap_3 <- read.csv(file.path(DIR_TABLES, "DEG_overlap_GSE56500_GSE68605.csv"))

universe <- intersect(top_56500$GeneSymbol, top_68605$GeneSymbol)
merged <- merge(top_56500[top_56500$GeneSymbol %in% universe, c("GeneSymbol", "logFC")],
                 top_68605[top_68605$GeneSymbol %in% universe, c("GeneSymbol", "logFC")],
                 by = "GeneSymbol", suffixes = c("_56500", "_68605"))
cat("Universe size:", length(universe), "\n")

#Spearman correlation, all shared genes, with permutation test
sp <- cor.test(merged$logFC_56500, merged$logFC_68605, method = "spearman")
set.seed(42)
nperm <- 10000
perm_rhos <- replicate(nperm, cor(merged$logFC_56500, sample(merged$logFC_68605), method = "spearman"))
emp_p <- (sum(abs(perm_rhos) >= abs(sp$estimate)) + 1) / (nperm + 1)
cat("Spearman rho (all shared genes):", sp$estimate, "| permutation p:", emp_p, "\n")

#Spearman correlation, genes significant in either dataset
sig_either <- intersect(union(sig_56500$GeneSymbol, sig_68605$GeneSymbol), universe)
merged_sig <- merged[merged$GeneSymbol %in% sig_either, ]
sp_sig <- cor.test(merged_sig$logFC_56500, merged_sig$logFC_68605, method = "spearman")
cat("Spearman rho (significant in either, n =", nrow(merged_sig), "):", sp_sig$estimate,
    "| p:", sp_sig$p.value, "\n")

#Concordance stats - the numbers behind the plot subtitle, as a table
concordance_stats <- data.frame(
  metric = c("universe_n", "rho_all_genes", "perm_p_all_genes",
             "n_significant_either", "rho_significant_either", "p_significant_either"),
  value  = c(length(universe), sp$estimate, emp_p,
             nrow(merged_sig), sp_sig$estimate, sp_sig$p.value))
write.csv(concordance_stats, file.path(DIR_TABLES, "cross_tissue_concordance_stats.csv"), row.names = FALSE)

#Concordance scatter plot
merged$Overlap3 <- merged$GeneSymbol %in% overlap_3$GeneSymbol
p_concordance <- ggplot(merged, aes(x = logFC_56500, y = logFC_68605)) +
  geom_point(data = subset(merged, !Overlap3), color = "grey70", alpha = 0.4, size = 1) +
  geom_point(data = subset(merged, Overlap3), color = "#E63946", size = 3.5) +
  geom_text(data = subset(merged, Overlap3), aes(label = GeneSymbol),
            vjust = -0.9, color = "#E63946", fontface = "bold", size = 3.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_smooth(method = "lm", se = FALSE, color = "steelblue", linewidth = 0.6) +
  labs(title = "Cross-tissue logFC concordance (14,614 shared genes)",
       subtitle = sprintf("Spearman rho = %.3f, permutation p = %.4f", sp$estimate, emp_p),
       x = "logFC, GSE56500 (Spinal Cord)", y = "logFC, GSE68605 (Motor Cortex)") +
  theme_bw(base_size = 13)
ggsave(file.path(DIR_FIGURES, "cross_tissue_logFC_concordance.png"), plot = p_concordance,
       width = 8, height = 6.5, dpi = 300, bg = "white")

#Threshold sensitivity - hypergeometric overlap at increasing adj.P.Val cutoffs
sensitivity <- lapply(c(0.01, 0.05, 0.10, 0.20), function(cutoff) {
  sA <- top_56500$GeneSymbol[top_56500$adj.P.Val < cutoff & abs(top_56500$logFC) > 0.5 & top_56500$GeneSymbol %in% universe]
  sB <- top_68605$GeneSymbol[top_68605$adj.P.Val < cutoff & abs(top_68605$logFC) > 0.5 & top_68605$GeneSymbol %in% universe]
  shared <- intersect(sA, sB)
  N <- length(universe); K <- length(sA); n <- length(sB); k <- length(shared)
  expected <- K * n / N
  p <- phyper(k - 1, K, N - K, n, lower.tail = FALSE)
  data.frame(cutoff = cutoff, sig_56500 = K, sig_68605 = n, shared = k,
             expected = round(expected, 3), hypergeometric_p = signif(p, 4))
})
sensitivity <- do.call(rbind, sensitivity)
print(sensitivity, row.names = FALSE)
write.csv(sensitivity, file.path(DIR_TABLES, "cross_tissue_threshold_sensitivity.csv"), row.names = FALSE)

#Shared-gene mean expression by group, both datasets
rds_56500 <- readRDS(file.path(DIR_INTERIM, "ALS_DEG_final.rds"))
rds_68605 <- readRDS(file.path(DIR_INTERIM, "GSE68605_DEG_final.rds"))

mean_expr <- function(expr, pheno, gene) {
  g <- expr[gene, ]
  c(Control = mean(g[pheno$ALS_status == "Control"]), ALS = mean(g[pheno$ALS_status == "ALS"]))
}
genes3 <- overlap_3$GeneSymbol
shared_expr <- do.call(rbind, lapply(genes3, function(g) {
  m56500 <- mean_expr(rds_56500$expr_mat, rds_56500$pheno, g)
  m68605 <- mean_expr(rds_68605$expr_mat, rds_68605$pheno, g)
  data.frame(Gene = g,
             GSE56500_Control = round(m56500["Control"], 3), GSE56500_ALS = round(m56500["ALS"], 3),
             GSE68605_Control = round(m68605["Control"], 3), GSE68605_ALS = round(m68605["ALS"], 3))
}))
print(shared_expr, row.names = FALSE)
write.csv(shared_expr, file.path(DIR_TABLES, "cross_tissue_shared_genes_expression.csv"), row.names = FALSE)

cat("Shared genes in SCFA panel:", any(genes3 %in% all_scfa), "\n")
