#SCFA Gene Filtering – GSE56500 
#curated SCFA-related genes from DEG 

# Load DEG results
library(ggplot2)
library(pheatmap)
load("../DEG/ALS_DEG_final.RData")
clean_rownames <- function(x) {
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
        !grepl("^[0-9]", gene) && !grepl("^NM_|^NR_|^XM_|^ENST|^uc", gene)) {
      return(gene)
    }}
  return(trimws(parts[length(parts)]))}
rownames(expr_mat) <- sapply(rownames(expr_mat), clean_rownames)
top_ALS$GeneSymbol <- sapply(top_ALS$GeneSymbol, clean_rownames)
iqr_vals        <- apply(expr_mat, 1, IQR)
expr_mat        <- expr_mat[order(iqr_vals, decreasing = TRUE), ]
expr_mat        <- expr_mat[!duplicated(rownames(expr_mat)), ]
top_ALS  <- top_ALS[order(top_ALS$P.Value), ]
top_ALS  <- top_ALS[!duplicated(top_ALS$GeneSymbol), ]
required_objects <- c("top_ALS", "expr_mat", "pheno")
missing_objects <- required_objects[!sapply(required_objects, exists)]
if (length(missing_objects) > 0) {
  stop(
    paste("Missing objects in RData:",
          paste(missing_objects, collapse = ", ")))
}
# Define curated SCFA gene list by functions
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
#Filtering DEG table for SCFA genes
scfa_df          <- top_ALS[top_ALS$GeneSymbol %in% all_scfa, ]
scfa_df$Category <- gene_category[scfa_df$GeneSymbol]
scfa_df$sig <- ifelse(scfa_df$adj.P.Val < 0.05, "Significant", "NS")
scfa_df          <- scfa_df[order(scfa_df$Category, scfa_df$P.Value), ]
cat("SCFA genes found in dataset:", nrow(scfa_df), "of", length(all_scfa), "\n")
cat(
  "Significant SCFA genes:",
  sum(scfa_df$adj.P.Val < 0.05),"\n")
print(scfa_df[, c("GeneSymbol", "Category", "logFC", "P.Value", "adj.P.Val", "sig")])
write.csv(scfa_df, "SCFA_DEG_results.csv", row.names = FALSE)
category_colors <- c(
  "FFA Receptors"       = "#E63946",
  "Transporters"        = "#F4A261",
  "Butyrate Metabolism" = "#2A9D8F",
  "HDAC Targets"        = "#457B9D",
  "NF-kB"               = "#6A0572",
  "NLRP3"               = "#E76F51",
  "Gut-Brain"           = "#264653"
)
top_ALS$Category <- ifelse(top_ALS$GeneSymbol %in% all_scfa,
                           gene_category[top_ALS$GeneSymbol], "Other")
# Volcano plot- SCFA genes highlighted by category
p_volcano <- ggplot() +
  geom_point(data = subset(top_ALS, Category == "Other"),
             aes(x = logFC, y = -log10(P.Value)),
             color = "grey80", alpha = 0.4, size = 1.2) +
  geom_point(data = subset(top_ALS, Category != "Other"),
             aes(x = logFC, y = -log10(P.Value), color = Category),
             size = 3.5, shape = 18) +
  geom_text(data = subset(top_ALS, Category != "Other"),
            aes(x = logFC, y = -log10(P.Value),
                label = GeneSymbol, color = Category),
            size = 2.8, vjust = -0.8, fontface = "bold") +
  scale_color_manual(values = category_colors, name = "Gene Category") +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed",
             color = "black", alpha = 0.4) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed",
             color = "black", alpha = 0.4) +
  labs(title    = "SCFA-Related Genes in ALS vs Control (GSE56500)",
       subtitle = paste("SCFA genes found:", nrow(scfa_df), "| Significant:", sum(scfa_df$sig == "Significant")),
       x = "Log2 Fold Change",
       y = "-log10(P-value)") +
  theme_bw(base_size = 13) +
  theme(plot.title      = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle   = element_text(hjust = 0.5, color = "grey40"),
        legend.position = "right",
        legend.text     = element_text(size = 9))
ggsave("volcano_SCFA_highlighted.png", plot = p_volcano,
       width = 10, height = 6, dpi = 300, bg = "white")
found_genes <- rownames(expr_mat)[rownames(expr_mat) %in% all_scfa]
if (length(found_genes) == 0) {
  stop("No SCFA genes found in expression matrix.")
}
scfa_expr   <- expr_mat[found_genes, ]
row_annot <- data.frame(
  Category  = gene_category[rownames(scfa_expr)],
  row.names = rownames(scfa_expr)
)
col_annot <- data.frame(
  Group     = pheno$ALS_status,
  row.names = colnames(scfa_expr)
)
annot_colors <- list(
  Group    = c("ALS" = "#E63946", "Control" = "#457B9D"),
  Category = category_colors
)
row_order <- order(row_annot$Category)
scfa_expr <- scfa_expr[row_order, ]
row_annot <- row_annot[row_order, , drop = FALSE]

# Heatmap- row-scaled expression across samples
pheatmap(
  scfa_expr,
  annotation_col    = col_annot,
  annotation_row    = row_annot,
  annotation_colors = annot_colors,
  scale             = "row",
  cluster_cols      = TRUE,
  cluster_rows      = FALSE,      # keep category grouping intact
  show_colnames     = FALSE,
  fontsize_row      = 9,
  color             = colorRampPalette(c("#457B9D", "white", "#E63946"))(100),
  main              = "SCFA Gene Expression: ALS vs Control (GSE56500)",
  filename          = "heatmap_SCFA_genes.png",
  width             = 10,
  height            = 8)
save.image("SCFA_GSE56500_analysis.RData")
