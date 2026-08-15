#Snapshot KEGG pathway definitions - gseKEGG pulls these live and they drift
source(here::here("scripts", "_paths.R"))
library(clusterProfiler)

kegg <- download_KEGG("hsa")
saveRDS(kegg, file.path(DIR_INTERIM, "kegg_hsa_snapshot.rds"))

writeLines(
  paste("KEGG hsa retrieved:", format(Sys.Date()),
        "| clusterProfiler", as.character(packageVersion("clusterProfiler"))),
  here::here("results", "kegg_version.txt"))
