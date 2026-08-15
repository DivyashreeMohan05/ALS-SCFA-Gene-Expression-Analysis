#Full pipeline - run from repo root with: Rscript scripts/00_run_all.R

source(here::here("scripts", "_paths.R"))
source(here::here("scripts", "_fetch_geo.R"))
source(here::here("scripts", "_helpers.R"))

scripts <- c(
  "01_DEG_Analysis.R",
  "02_SCFA_Analysis.R",
  "03_DEG_SCFA.R",
  "04_Pathway_enrichment.R",
  "seed_stability.R",
  "05_overlap_significance.R",
  "06_cross_tissue_concordance.R")

for (s in scripts) {
  message("=== ", s, " ===")
  source(here::here("scripts", s), echo = FALSE)
}

writeLines(capture.output(sessionInfo()), here::here("results", "sessionInfo.txt"))
