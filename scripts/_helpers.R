#Shared helpers - ranking/mapping functions need clusterProfiler and org.Hs.eg.db loaded by caller

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

#SCFA gene panel
scfa_genes <- list(
  "FFA Receptors"       = c("FFAR2", "FFAR3", "FFAR4", "GPR109A"),
  "Transporters"        = c("SLC5A8", "SLC16A1", "SLC16A3"),
  "Butyrate Metabolism" = c("ACSS2", "ACAT1", "HADHA", "HADHB"),
  "HDAC Targets"        = c("HDAC1", "HDAC2", "HDAC3", "HDAC4",
                            "HDAC5", "HDAC6", "HDAC7", "HDAC8",
                            "SIRT1", "SIRT3"),
  "NF-kB"               = c("NFKB1", "RELA", "IKBKB", "NFKBIA"),
  "NLRP3"               = c("NLRP3", "CASP1", "IL1B", "IL18"),
  "Gut-Brain"           = c("TLR4", "MYD88", "TREM2", "CX3CR1")
)
all_scfa      <- unlist(scfa_genes, use.names = FALSE)
gene_category <- rep(names(scfa_genes), lengths(scfa_genes))
names(gene_category) <- all_scfa

#Write or update one row of preprocessing_summary.csv, keyed by accession
update_preprocessing_summary <- function(accession, ...) {
  file <- file.path(DIR_TABLES, "preprocessing_summary.csv")
  new_vals <- list(...)
  tab <- if (file.exists(file)) {
    read.csv(file, stringsAsFactors = FALSE)
  } else {
    data.frame(accession = character(0), stringsAsFactors = FALSE)
  }
  i <- if (accession %in% tab$accession) which(tab$accession == accession) else nrow(tab) + 1
  tab[i, "accession"] <- accession
  for (col in names(new_vals)) {
    if (!col %in% names(tab)) tab[[col]] <- NA
    tab[i, col] <- new_vals[[col]]
  }
  write.csv(tab, file, row.names = FALSE)
}
