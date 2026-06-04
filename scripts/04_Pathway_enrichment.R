
# Gene Set Enrichment Analysis 
# ALS vs Control Transcriptomic Datasets- GSE56500 (Spinal Cord)& GSE68605 (Motor Cortex)
## Method: clusterProfiler GSEA, KEGG pathway enrichment, GO Biological Process enrichment, Ranking metric = logFC × -log10(P-value)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(enrichplot)
#Creating a ranked gene list
make_ranked_list <- function(df) {
  df <- df[!is.na(df$GeneSymbol) & !is.na(df$logFC) & !is.na(df$P.Value), ]
  df <- df[!duplicated(df$GeneSymbol), ]
  df$P.Value <- pmax(df$P.Value, 1e-300)
  df$P.Value <- pmin(df$P.Value, 1 - 1e-6)
  ranked        <- df$logFC * -log10(df$P.Value)
  names(ranked) <- df$GeneSymbol
  sort(ranked, decreasing = TRUE)}
#Convert gene symbols to Entrez IDs for further analysis
symbol_to_entrez <- function(ranked) {
  # Removing non-standard names before mapping
  ranked <- ranked[!is.na(names(ranked))]
  ranked <- ranked[names(ranked) != ""]
  ranked <- ranked[!grepl("^[0-9]", names(ranked))]   # removes probe IDs starting with digits
  ranked <- ranked[!grepl("\\.", names(ranked))]  # removes names with ambiguous probes   
  ids <- bitr(names(ranked),
              fromType = "SYMBOL",
              toType   = "ENTREZID",
              OrgDb    = org.Hs.eg.db)
  ranked <- ranked[names(ranked) %in% ids$SYMBOL]
  names(ranked) <- ids$ENTREZID[match(names(ranked), ids$SYMBOL)]
  sort(ranked, decreasing = TRUE)}
# Run GSEA
run_gsea <- function(ranked_entrez, label) {
#Perform KEGG pathway enrichment analysis
  gsea_kegg <- gseKEGG(
    geneList = ranked_entrez,
    organism = "hsa",
    minGSSize = 15,
    maxGSSize = 500,
    pvalueCutoff = 0.05,
    nPermSimple  = 10000,
    verbose = FALSE)
  if (nrow(as.data.frame(gsea_kegg)) > 0) {
    gsea_kegg@result <- gsea_kegg@result[!is.na(gsea_kegg@result$pvalue), ]}
#Perform GO Biological Process enrichment analysis
  gsea_go <- gseGO(
    geneList = ranked_entrez,
    OrgDb = org.Hs.eg.db,
    ont = "BP",
    minGSSize = 15,
    maxGSSize = 500,
    pvalueCutoff = 0.05,
    nPermSimple  = 10000,
    verbose = FALSE)
  if (nrow(as.data.frame(gsea_go)) > 0) {
    gsea_go@result <- gsea_go@result[!is.na(gsea_go@result$pvalue), ]}
  write.csv(
      as.data.frame(gsea_kegg),
      paste0("GSEA_KEGG_", label, ".csv"),
      row.names = FALSE)
  if (nrow(as.data.frame(gsea_go)) > 0) {
    write.csv(
      as.data.frame(gsea_go),
      paste0("GSEA_GO_BP_", label, ".csv"),
      row.names = FALSE) }
  cat("\n===", label, "===\n")
  cat("KEGG pathways:", nrow(as.data.frame(gsea_kegg)), "\n")
  cat("GO BP terms: ", nrow(as.data.frame(gsea_go)), "\n")
  return(list(
    kegg = gsea_kegg,
    go = gsea_go
  ))}
#Generate and save GSEA dotplots
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
    filename,
    plot = p,
    width = 12,
    height = 8,
    dpi = 300,
    bg = "white")
   return(p)}
#GSEA analysis: GSE56500
top_56500 <- read.csv("../DEG/ALS_vs_Control_all.csv")
clean_symbol <- function(x) {
  if (is.na(x)) return(NA)
  if (grepl(" /// ", x)) {
    parts <- strsplit(x, " /// ")[[1]]
  } else if (grepl(" // ", x)) {
    parts <- strsplit(x, " // ")[[1]]
  } else {
    return(x)
  }
  for (i in 2:length(parts)) {
    gene <- trimws(parts[i])
    if (!is.na(gene) && gene != "" && gene != "---" &&
        !grepl("^[0-9]", gene) && !grepl("^NM_|^NR_|^XM_|^ENST|^uc|^AF|^AK|^BC", gene)) {
      return(gene)
    }
  }
  return(trimws(parts[length(parts)]))
}
top_56500$GeneSymbol <- sapply(top_56500$GeneSymbol, clean_symbol)
top_56500 <- top_56500[!is.na(top_56500$GeneSymbol), ]
top_56500 <- top_56500[!duplicated(top_56500$GeneSymbol), ]
ranked_56500 <- make_ranked_list(top_56500)
ranked_56500_entrez <- symbol_to_entrez(ranked_56500)
cat("Genes mapped to Entrez (GSE56500):", length(ranked_56500_entrez), "\n")
gsea_56500 <- run_gsea(ranked_56500_entrez, "GSE56500")
if (nrow(as.data.frame(gsea_56500$kegg)) == 0) {
  cat("Warning: No significant KEGG pathways for GSE56500\n")}
if (nrow(as.data.frame(gsea_56500$go)) == 0) {
  cat("Warning: No significant GO BP terms for GSE56500\n")
  }
save_dotplot(
  gsea_56500$kegg,
  "GSEA KEGG — ALS vs Control (GSE56500)",
  "GSEA_KEGG_dotplot_GSE56500.png")
save_dotplot(
  gsea_56500$go,
  "GSEA GO BP — ALS vs Control (GSE56500)",
  "GSEA_GOBP_dotplot_GSE56500.png")
#GSEA analysis: GSE68605
top_68605 <- read.csv("../DEG+SCFA/GSE68605_ALS_vs_Control_all.csv")
ranked_68605 <- make_ranked_list(top_68605)
ranked_68605_entrez <- symbol_to_entrez(ranked_68605)
cat("Genes mapped to Entrez (GSE68605):", length(ranked_68605_entrez), "\n") 
gsea_68605 <- run_gsea(ranked_68605_entrez, "GSE68605")
if (nrow(as.data.frame(gsea_68605$kegg)) == 0) {
  cat("Warning: No significant KEGG pathways for GSE68605\n")}
if (nrow(as.data.frame(gsea_68605$go)) == 0) {
  cat("Warning: No significant GO BP terms for GSE68605\n")}
save_dotplot(
  gsea_68605$kegg,
  "GSEA KEGG — ALS vs Control (GSE68605)",
  "GSEA_KEGG_dotplot_GSE68605.png")
save_dotplot(
  gsea_68605$go,
  "GSEA GO BP — ALS vs Control (GSE68605)",
  "GSEA_GOBP_dotplot_GSE68605.png")
#Cross dataset comparison
kegg_56500 <- as.data.frame(gsea_56500$kegg)
kegg_68605 <- as.data.frame(gsea_68605$kegg)
if (nrow(kegg_56500) > 0 & nrow(kegg_68605) > 0) {
  shared_kegg <- intersect(
    kegg_56500$Description,
    kegg_68605$Description)
  cat("\nShared KEGG pathways (both datasets):",
      length(shared_kegg), "\n")
  if (length(shared_kegg) > 0) {
    cat("Shared KEGG pathways replicated across spinal cord and motor cortex:\n")
    print(shared_kegg)
    write.csv(
      data.frame(Pathway = shared_kegg),
      "GSEA_shared_KEGG_pathways.csv",
      row.names = FALSE)} else {
     cat("No shared KEGG pathways — tissue-specific enrichment patterns\n")
    }}
#Save complete GSEA workspace
save.image("GSEA_analysis.RData")