#Cross-tissue KEGG overlap significance - GSE56500 (Spinal Cord) & GSE68605 (Motor Cortex)
#Tests: hypergeometric on shared pathway count | binomial on NES sign concordance

source(here::here("scripts", "_paths.R"))
library(clusterProfiler)
library(org.Hs.eg.db)
set.seed(42)

if (!exists("make_ranked_list") || !exists("symbol_to_entrez")) {
  source(here::here("scripts", "04_Pathway_enrichment.R"))
}

#Load DEG results
top_56500 <- read.csv(file.path(DIR_TABLES, "ALS_vs_Control_all.csv"))
top_68605 <- read.csv(file.path(DIR_TABLES, "GSE68605_ALS_vs_Control_all.csv"))

ranked_A <- symbol_to_entrez(make_ranked_list(top_56500))
ranked_B <- symbol_to_entrez(make_ranked_list(top_68605))

#GSEA with pvalueCutoff = 1 to recover the full tested set
gsea_all <- function(ranked_entrez) {
  as.data.frame(
    gseKEGG(
      geneList     = ranked_entrez,
      organism     = "hsa",
      minGSSize    = 15,
      maxGSSize    = 500,
      pvalueCutoff = 1,
      nPermSimple  = 10000,
      seed         = TRUE,
      verbose      = FALSE))
}

all_A <- gsea_all(ranked_A)
all_B <- gsea_all(ranked_B)

write.csv(all_A, file.path(DIR_TABLES, "GSEA_KEGG_GSE56500_alltested.csv"), row.names = FALSE)
write.csv(all_B, file.path(DIR_TABLES, "GSEA_KEGG_GSE68605_alltested.csv"), row.names = FALSE)

#Overlap count - hypergeometric
universe <- intersect(all_A$ID, all_B$ID)
N <- length(universe)

sig_A  <- all_A$ID[all_A$p.adjust < 0.05 & all_A$ID %in% universe]
sig_B  <- all_B$ID[all_B$p.adjust < 0.05 & all_B$ID %in% universe]
shared <- intersect(sig_A, sig_B)

k <- length(shared); K <- length(sig_A); n <- length(sig_B)
expected  <- K * n / N
p_overlap <- phyper(k - 1, K, N - K, n, lower.tail = FALSE)

cat("\nPathways tested in both:", N, "\n")
cat("Significant - GSE56500:", K, "| GSE68605:", n, "\n")
cat("Shared:", k, "| expected:", round(expected, 1), "\n")
cat("Hypergeometric p =", signif(p_overlap, 3), "\n")

#Direction concordance - binomial against marginal sign proportions
nes_A <- setNames(all_A$NES, all_A$ID)
nes_B <- setNames(all_B$NES, all_B$ID)

agree  <- sum(sign(nes_A[shared]) == sign(nes_B[shared]))
pA     <- mean(nes_A[sig_A] > 0)
pB     <- if (n > 0) mean(nes_B[sig_B] > 0) else NA
p_null <- pA * pB + (1 - pA) * (1 - pB)
bt     <- if (k > 0) binom.test(agree, k, p_null, alternative = "greater") else list(p.value = NA)

cat("\nActivated fraction - GSE56500:", round(pA, 2), "| GSE68605:", round(pB, 2), "\n")
cat("Expected agreement:", round(p_null, 3), "\n")
cat("Observed agreement:", agree, "/", k, "\n")
cat("Binomial p =", signif(bt$p.value, 3), "\n")

#Save results
concordance <- data.frame(
  Pathway       = all_A$Description[match(shared, all_A$ID)],
  NES_GSE56500  = round(nes_A[shared], 3),
  NES_GSE68605  = round(nes_B[shared], 3),
  SameDirection = sign(nes_A[shared]) == sign(nes_B[shared]),
  row.names     = NULL)

stats <- data.frame(
  metric = c("universe_N", "sig_GSE56500", "sig_GSE68605", "shared",
             "expected_by_chance", "hypergeometric_p",
             "sign_agreement", "null_agreement_prob", "binomial_p"),
  value  = c(N, K, n, k, round(expected, 2), signif(p_overlap, 4),
             paste0(agree, "/", k), round(p_null, 4), signif(bt$p.value, 4)))

write.csv(concordance, file.path(DIR_TABLES, "KEGG_overlap_concordance.csv"), row.names = FALSE)

#Same two tests, GO Biological Process terms instead of KEGG pathways
gsea_all_go <- function(ranked_entrez) {
  as.data.frame(
    gseGO(
      geneList     = ranked_entrez,
      OrgDb        = org.Hs.eg.db,
      ont          = "BP",
      minGSSize    = 15,
      maxGSSize    = 500,
      pvalueCutoff = 1,
      nPermSimple  = 10000,
      seed         = TRUE,
      verbose      = FALSE))
}

all_A_go <- gsea_all_go(ranked_A)
all_B_go <- gsea_all_go(ranked_B)

universe_go <- intersect(all_A_go$ID, all_B_go$ID)
N_go <- length(universe_go)

sig_A_go <- all_A_go$ID[all_A_go$p.adjust < 0.05 & all_A_go$ID %in% universe_go]
sig_B_go <- all_B_go$ID[all_B_go$p.adjust < 0.05 & all_B_go$ID %in% universe_go]
shared_go <- intersect(sig_A_go, sig_B_go)

k_go <- length(shared_go); K_go <- length(sig_A_go); n_go <- length(sig_B_go)
expected_go  <- K_go * n_go / N_go
p_overlap_go <- phyper(k_go - 1, K_go, N_go - K_go, n_go, lower.tail = FALSE)

cat("\nGO BP terms tested in both:", N_go, "\n")
cat("Significant - GSE56500:", K_go, "| GSE68605:", n_go, "\n")
cat("Shared:", k_go, "| expected:", round(expected_go, 2), "\n")
cat("Hypergeometric p =", signif(p_overlap_go, 3), "\n")

nes_A_go <- setNames(all_A_go$NES, all_A_go$ID)
nes_B_go <- setNames(all_B_go$NES, all_B_go$ID)

agree_go  <- sum(sign(nes_A_go[shared_go]) == sign(nes_B_go[shared_go]))
pA_go     <- if (K_go > 0) mean(nes_A_go[sig_A_go] > 0) else NA
pB_go     <- if (n_go > 0) mean(nes_B_go[sig_B_go] > 0) else NA
p_null_go <- pA_go * pB_go + (1 - pA_go) * (1 - pB_go)
bt_go     <- if (k_go > 0) binom.test(agree_go, k_go, p_null_go, alternative = "greater") else list(p.value = NA)

cat("\nActivated fraction - GSE56500:", round(pA_go, 2), "| GSE68605:", round(pB_go, 2), "\n")
cat("Expected agreement:", round(p_null_go, 3), "\n")
cat("Observed agreement:", agree_go, "/", k_go, "\n")
cat("Binomial p =", signif(bt_go$p.value, 3), "\n")

concordance_go <- data.frame(
  Term          = all_A_go$Description[match(shared_go, all_A_go$ID)],
  NES_GSE56500  = round(nes_A_go[shared_go], 3),
  NES_GSE68605  = round(nes_B_go[shared_go], 3),
  SameDirection = sign(nes_A_go[shared_go]) == sign(nes_B_go[shared_go]),
  row.names     = NULL)

write.csv(concordance_go, file.path(DIR_TABLES, "GO_overlap_concordance.csv"), row.names = FALSE)

stats_go <- data.frame(
  metric = paste0("GO_", c("universe_N", "sig_GSE56500", "sig_GSE68605", "shared",
                            "expected_by_chance", "hypergeometric_p",
                            "sign_agreement", "null_agreement_prob", "binomial_p")),
  value  = c(N_go, K_go, n_go, k_go, round(expected_go, 2), signif(p_overlap_go, 4),
             paste0(agree_go, "/", k_go), round(p_null_go, 4), signif(bt_go$p.value, 4)))

stats <- rbind(stats, stats_go)
write.csv(stats, file.path(DIR_TABLES, "overlap_significance_stats.csv"), row.names = FALSE)
