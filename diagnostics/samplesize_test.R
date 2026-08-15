#Diagnostic only - not part of the pipeline, see diagnostics/README.md
#Shows sampleSize (not nPerm/maxPerm) is the precision lever for
#gseKEGG()/gseGO()'s default method, unreachable through their public API

source(here::here("scripts", "_paths.R"))
library(clusterProfiler)
library(org.Hs.eg.db)
source(here::here("scripts", "_helpers.R"))

top_56500 <- read.csv(file.path(DIR_TABLES, "ALS_vs_Control_all.csv"))
ranked_A <- symbol_to_entrez(make_ranked_list(top_56500))

set.seed(42)
r0 <- gseKEGG(geneList = ranked_A, organism = "hsa", minGSSize = 15, maxGSSize = 500,
              pvalueCutoff = 1, verbose = FALSE)
gene_sets <- r0@geneSets
geneList  <- r0@geneList

#enrichit:::gsea() called directly - diagnostic only, unexported internal
run_direct <- function(seed, sampleSize) {
  set.seed(seed)
  res <- enrichit:::gsea(geneList = geneList, gene_sets = gene_sets,
                          minGSSize = 15, maxGSSize = 500,
                          method = "multilevel", sampleSize = sampleSize, verbose = FALSE)
  res$p.adjust <- p.adjust(res$pvalue, method = "BH")
  sum(res$p.adjust < 0.05, na.rm = TRUE)
}

seeds <- c(42, 1, 7)
cat("=== default sampleSize (101), via gseKEGG's actual default ===\n")
for (s in seeds) cat("seed", s, ":", run_direct(s, 101), "\n")

cat("\n=== sampleSize = 10001 (unreachable through gseKEGG()/gseGO()) ===\n")
for (s in seeds) cat("seed", s, ":", run_direct(s, 10001), "\n")
