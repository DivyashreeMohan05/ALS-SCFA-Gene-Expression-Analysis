
# Gene Set Enrichment Analysis
# ALS vs Control Transcriptomic Datasets- GSE56500 (Spinal Cord)& GSE68605 (Motor Cortex)
## Method: clusterProfiler GSEA, KEGG pathway enrichment, GO Biological Process enrichment, Ranking metric = logFC × -log10(P-value)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(enrichplot)
source(here::here("scripts", "_paths.R"))
source(here::here("scripts", "_helpers.R"))

#seed=TRUE/nPermSimple are silently dropped by gseKEGG()/gseGO() (never
#forwarded to enrichit::gsea_gson). A single set.seed() per script only
#anchors the first permutation call, not each one - results depend on call
#order/count. Fixed by seeding immediately before every individual call.
#Even seeded, sampleSize (the actual multilevel-method precision knob) is
#unreachable through the public API and stuck at its default (101), so a
#single seed's result still isn't representative - hence consensus across
#seeds below, kept as the primary result.
CONSENSUS_SEEDS <- c(42, 1, 7, 123, 2024, 5, 11, 17, 23, 29,
                      31, 37, 41, 43, 47, 53, 59, 61, 67, 71)

run_kegg_once <- function(ranked_entrez, seed) {
  set.seed(seed)
  gseKEGG(geneList = ranked_entrez, organism = "hsa", minGSSize = 15, maxGSSize = 500,
          pvalueCutoff = 1, nPerm = 10000, verbose = FALSE)
}
run_go_once <- function(ranked_entrez, seed) {
  set.seed(seed)
  gseGO(geneList = ranked_entrez, OrgDb = org.Hs.eg.db, ont = "BP", minGSSize = 15, maxGSSize = 500,
        pvalueCutoff = 1, nPerm = 10000, verbose = FALSE)
}

#Per gene set: fraction of seeds where p.adjust < 0.05, plus median NES/p.adjust
consensus_table <- function(runs) {
  runs_df <- lapply(runs, as.data.frame)
  all_ids <- unique(unlist(lapply(runs_df, function(d) d$ID)))
  n <- length(runs_df)
  rows <- lapply(all_ids, function(id) {
    padj <- sapply(runs_df, function(d) { i <- which(d$ID == id); if (length(i)) d$p.adjust[i[1]] else NA })
    nes  <- sapply(runs_df, function(d) { i <- which(d$ID == id); if (length(i)) d$NES[i[1]] else NA })
    owner <- runs_df[[which(sapply(runs_df, function(d) id %in% d$ID))[1]]]
    data.frame(ID = id, Description = owner$Description[owner$ID == id][1],
               seed_fraction = mean(!is.na(padj) & padj < 0.05),
               median_NES = median(nes, na.rm = TRUE),
               median_p.adjust = median(padj, na.rm = TRUE))
  })
  do.call(rbind, rows)
}

#Consensus (primary) + single-seed (supplementary) output for one dataset
run_gsea <- function(ranked_entrez, label) {
  kegg_runs <- lapply(CONSENSUS_SEEDS, function(s) run_kegg_once(ranked_entrez, s))
  go_runs   <- lapply(CONSENSUS_SEEDS, function(s) run_go_once(ranked_entrez, s))

  kegg_cons <- consensus_table(kegg_runs)
  go_cons   <- consensus_table(go_runs)
  kegg_counts <- sapply(kegg_runs, function(r) sum(as.data.frame(r)$p.adjust < 0.05, na.rm = TRUE))
  go_counts   <- sapply(go_runs,   function(r) sum(as.data.frame(r)$p.adjust < 0.05, na.rm = TRUE))

  write.csv(kegg_cons[order(-kegg_cons$seed_fraction), ],
            file.path(DIR_TABLES, paste0("GSEA_consensus_KEGG_", label, ".csv")), row.names = FALSE)
  write.csv(go_cons[order(-go_cons$seed_fraction), ],
            file.path(DIR_TABLES, paste0("GSEA_consensus_GO_", label, ".csv")), row.names = FALSE)

  #Supplementary: single-seed (first of CONSENSUS_SEEDS) output, same shape as before the consensus fix
  gsea_kegg <- kegg_runs[[1]]
  gsea_go   <- go_runs[[1]]
  if (nrow(as.data.frame(gsea_kegg)) > 0) {
    gsea_kegg@result <- gsea_kegg@result[!is.na(gsea_kegg@result$pvalue), ]}
  if (nrow(as.data.frame(gsea_go)) > 0) {
    gsea_go@result <- gsea_go@result[!is.na(gsea_go@result$pvalue), ]}
  write.csv(as.data.frame(gsea_kegg), file.path(DIR_TABLES, paste0("GSEA_KEGG_", label, ".csv")), row.names = FALSE)
  if (nrow(as.data.frame(gsea_go)) > 0) {
    write.csv(as.data.frame(gsea_go), file.path(DIR_TABLES, paste0("GSEA_GO_BP_", label, ".csv")), row.names = FALSE)}

  cat("\n===", label, "(consensus across", length(CONSENSUS_SEEDS), "seeds) ===\n")
  cat("KEGG: count range", min(kegg_counts), "-", max(kegg_counts),
      "| consensus (>=80% of seeds):", sum(kegg_cons$seed_fraction >= 0.8), "\n")
  cat("GO BP: count range", min(go_counts), "-", max(go_counts),
      "| consensus (>=80% of seeds):", sum(go_cons$seed_fraction >= 0.8), "\n")

  return(list(kegg = gsea_kegg, go = gsea_go, kegg_consensus = kegg_cons, go_consensus = go_cons))}
#Generate and save GSEA dotplots (from the supplementary single-seed result)
save_dotplot <- function(gsea_obj, title, filename, n = 15) {
  if (nrow(as.data.frame(gsea_obj)) == 0) {
    cat("No significant terms for:", title, "\n")
    return(NULL)}
  p <- dotplot(gsea_obj, showCategory = n, split = ".sign") +
    facet_grid(. ~ .sign) +
    labs(title = title) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
      axis.text.y = element_text(size = 8))
  ggsave(
    file.path(DIR_FIGURES, filename),
    plot = p,
    width = 12,
    height = 8,
    dpi = 300,
    bg = "white")
   return(p)}
#GSEA analysis: GSE56500
top_56500 <- read.csv(file.path(DIR_TABLES, "ALS_vs_Control_all.csv"))
ranked_56500 <- make_ranked_list(top_56500)
ranked_56500_entrez <- symbol_to_entrez(ranked_56500)
cat("Genes mapped to Entrez (GSE56500):", length(ranked_56500_entrez), "\n")
gsea_56500 <- run_gsea(ranked_56500_entrez, "GSE56500")
save_dotplot(
  gsea_56500$kegg,
  "GSEA KEGG — ALS vs Control (GSE56500)",
  "GSEA_KEGG_dotplot_GSE56500.png")
save_dotplot(
  gsea_56500$go,
  "GSEA GO BP — ALS vs Control (GSE56500)",
  "GSEA_GOBP_dotplot_GSE56500.png")
#GSEA analysis: GSE68605
top_68605 <- read.csv(file.path(DIR_TABLES, "GSE68605_ALS_vs_Control_all.csv"))
ranked_68605 <- make_ranked_list(top_68605)
ranked_68605_entrez <- symbol_to_entrez(ranked_68605)
cat("Genes mapped to Entrez (GSE68605):", length(ranked_68605_entrez), "\n")
gsea_68605 <- run_gsea(ranked_68605_entrez, "GSE68605")
save_dotplot(
  gsea_68605$kegg,
  "GSEA KEGG — ALS vs Control (GSE68605)",
  "GSEA_KEGG_dotplot_GSE68605.png")
save_dotplot(
  gsea_68605$go,
  "GSEA GO BP — ALS vs Control (GSE68605)",
  "GSEA_GOBP_dotplot_GSE68605.png")
#Cross dataset comparison - consensus-significant sets (seed_fraction >= 0.8)
kegg_56500 <- gsea_56500$kegg_consensus[gsea_56500$kegg_consensus$seed_fraction >= 0.8, ]
kegg_68605 <- gsea_68605$kegg_consensus[gsea_68605$kegg_consensus$seed_fraction >= 0.8, ]
shared_kegg <- intersect(kegg_56500$Description, kegg_68605$Description)
cat("\nShared KEGG pathways (both datasets, consensus):",
    length(shared_kegg), "\n")
if (length(shared_kegg) > 0) {
  cat("Shared KEGG pathways replicated across spinal cord and motor cortex:\n")
  print(shared_kegg)
} else {
  cat("No shared KEGG pathways — tissue-specific enrichment patterns\n")
}
write.csv(
  data.frame(Pathway = shared_kegg),
  file.path(DIR_TABLES, "GSEA_shared_KEGG_pathways.csv"),
  row.names = FALSE)
