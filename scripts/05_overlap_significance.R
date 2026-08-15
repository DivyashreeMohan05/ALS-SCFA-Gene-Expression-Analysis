#Cross-tissue KEGG/GO BP overlap significance - GSE56500 (Spinal Cord) & GSE68605 (Motor Cortex)
#Tests: hypergeometric on shared-term count | binomial on NES sign concordance
#Uses 04's consensus tables (>=80% of 20 seeds significant), not a single seed -
#single-seed GSEA counts aren't reproducible enough to test overlap on directly.

source(here::here("scripts", "_paths.R"))

overlap_test <- function(tabA, tabB, id_col = "ID") {
  universe <- intersect(tabA[[id_col]], tabB[[id_col]])
  N <- length(universe)

  sig_A  <- tabA[[id_col]][tabA$seed_fraction >= 0.8 & tabA[[id_col]] %in% universe]
  sig_B  <- tabB[[id_col]][tabB$seed_fraction >= 0.8 & tabB[[id_col]] %in% universe]
  shared <- intersect(sig_A, sig_B)

  k <- length(shared); K <- length(sig_A); n <- length(sig_B)
  expected  <- K * n / N
  p_overlap <- phyper(k - 1, K, N - K, n, lower.tail = FALSE)

  nes_A <- setNames(tabA$median_NES, tabA[[id_col]])
  nes_B <- setNames(tabB$median_NES, tabB[[id_col]])

  agree  <- if (k > 0) sum(sign(nes_A[shared]) == sign(nes_B[shared])) else 0
  pA     <- if (K > 0) mean(nes_A[sig_A] > 0) else NA
  pB     <- if (n > 0) mean(nes_B[sig_B] > 0) else NA
  p_null <- pA * pB + (1 - pA) * (1 - pB)
  bt     <- if (k > 0) binom.test(agree, k, p_null, alternative = "greater") else list(p.value = NA)

  concordance <- data.frame(
    Term          = tabA$Description[match(shared, tabA[[id_col]])],
    NES_GSE56500  = round(nes_A[shared], 3),
    NES_GSE68605  = round(nes_B[shared], 3),
    SameDirection = sign(nes_A[shared]) == sign(nes_B[shared]),
    row.names     = NULL)

  stats <- data.frame(
    metric = c("universe_N", "sig_GSE56500", "sig_GSE68605", "shared",
               "expected_by_chance", "hypergeometric_p",
               "sign_agreement", "null_agreement_prob", "binomial_p"),
    value  = c(N, K, n, k, round(expected, 3), signif(p_overlap, 4),
               paste0(agree, "/", k), round(p_null, 4), signif(bt$p.value, 4)))

  list(N = N, K = K, n = n, k = k, expected = expected, p_overlap = p_overlap,
       agree = agree, p_null = p_null, binom_p = bt$p.value,
       concordance = concordance, stats = stats)
}

report <- function(label, res) {
  cat("\n===", label, "===\n")
  cat("Universe N:", res$N, "\n")
  cat("Consensus-significant - GSE56500:", res$K, "| GSE68605:", res$n, "\n")
  cat("Shared:", res$k, "| expected:", round(res$expected, 3), "\n")
  cat("Hypergeometric p =", signif(res$p_overlap, 4), "\n")
  cat("Sign agreement:", res$agree, "/", res$k, "| null prob:", round(res$p_null, 3), "\n")
  cat("Binomial p =", signif(res$binom_p, 4), "\n")
  if (res$k > 0) print(res$concordance$Term)
}

kegg_A <- read.csv(file.path(DIR_TABLES, "GSEA_consensus_KEGG_GSE56500.csv"))
kegg_B <- read.csv(file.path(DIR_TABLES, "GSEA_consensus_KEGG_GSE68605.csv"))
go_A   <- read.csv(file.path(DIR_TABLES, "GSEA_consensus_GO_GSE56500.csv"))
go_B   <- read.csv(file.path(DIR_TABLES, "GSEA_consensus_GO_GSE68605.csv"))

kegg_result <- overlap_test(kegg_A, kegg_B)
go_result   <- overlap_test(go_A, go_B)

report("KEGG consensus overlap", kegg_result)
report("GO BP consensus overlap", go_result)

write.csv(kegg_result$concordance, file.path(DIR_TABLES, "KEGG_overlap_concordance.csv"), row.names = FALSE)
write.csv(go_result$concordance, file.path(DIR_TABLES, "GO_overlap_concordance.csv"), row.names = FALSE)

stats <- rbind(kegg_result$stats,
               transform(go_result$stats, metric = paste0("GO_", metric)))
write.csv(stats, file.path(DIR_TABLES, "overlap_significance_stats.csv"), row.names = FALSE)
