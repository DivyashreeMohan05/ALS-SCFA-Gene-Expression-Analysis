#Shared GSEA helper functions - sourced by 04 and 05, requires
#clusterProfiler and org.Hs.eg.db already loaded by the caller

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
