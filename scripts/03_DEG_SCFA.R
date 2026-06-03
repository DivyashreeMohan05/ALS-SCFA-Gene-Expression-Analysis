#DEG Analysis + SCFA Overlay – GSE68605 (Motor Cortex)
#8 ALS + 3 Control
#Method used- limma | Threshold: adj.P.Val < 0.05, |logFC| > 0.5

library(GEOquery)
library(limma)
library(ggplot2)
library(pheatmap)
#Loading GEO dataset
gse      <- getGEO("GSE68605", GSEMatrix = TRUE)[[1]]
expr_mat <- exprs(gse)
pheno    <- pData(gse)
fdata    <- fData(gse)
#Assign ALS/Control labels
colnames(pheno)
print(table(pheno$`patient group:ch1`))
pheno$ALS_status <- ifelse(
  grepl("control", pheno$`patient group:ch1`, ignore.case = TRUE),
  "Control",
  "ALS")
table(pheno$ALS_status)
#Probe annotation
fdata$GeneSymbol <- trimws(as.character(fdata$`Gene Symbol`))
fdata$GeneSymbol[fdata$GeneSymbol == "" | fdata$GeneSymbol == "---"] <- NA
valid    <- !is.na(fdata$GeneSymbol)
expr_mat <- expr_mat[valid, ]
fdata    <- fdata[valid, ]
iqr_order          <- order(apply(expr_mat, 1, IQR), decreasing = TRUE)
expr_mat           <- expr_mat[iqr_order, ]
fdata              <- fdata[iqr_order, ]
keep               <- !duplicated(fdata$GeneSymbol)
expr_mat           <- expr_mat[keep, ]
rownames(expr_mat) <- fdata$GeneSymbol[keep]
cat("Genes retained after probe collapse:", nrow(expr_mat), "\n")
#DEG using limma
group            <- factor(pheno$ALS_status, levels = c("Control", "ALS"))
design           <- model.matrix(~0 + group)
colnames(design) <- c("Control", "ALS")
fit  <- lmFit(expr_mat, design)
fit2 <- contrasts.fit(fit, makeContrasts(ALS_vs_Control = ALS - Control,
                                         levels = design))
fit2 <- eBayes(fit2, trend = TRUE)
top_68605            <- topTable(fit2, coef = "ALS_vs_Control",
                                 adjust = "BH", number = Inf)
top_68605$GeneSymbol <- rownames(top_68605)
sig_68605 <- subset(top_68605, adj.P.Val< 0.05 & abs(logFC) > 0.5)
print(paste("GSE68605 significant DEGs:", nrow(sig_68605)))
write.csv(top_68605, "GSE68605_ALS_vs_Control_all.csv",         row.names = FALSE)
write.csv(sig_68605, "GSE68605_ALS_vs_Control_significant.csv", row.names = FALSE)
#Volcano plot
top_68605$sig <- ifelse(top_68605$adj.P.Val < 0.05 & abs(top_68605$logFC) > 0.5,
                        "Significant", "NS")
p_volcano <- ggplot(top_68605, aes(x = logFC, y = -log10(P.Value), color = sig)) +
  geom_point(alpha = 0.6, size = 1.8) +
  scale_color_manual(values = c("NS" = "grey70", "Significant" = "red"), name = "") +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed",
             color = "blue", alpha = 0.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed",
             color = "blue", alpha = 0.5) +
  annotate("text", x = Inf, y = Inf,
           label = paste("Significant n =", sum(top_68605$sig == "Significant")),
           hjust = 1.1, vjust = 1.5, size = 4.5, color = "red") +
  labs(title    = "ALS vs Control (GSE68605)",
       subtitle = "Affymetrix HG-U133 Plus 2.0",
       x = "Log2 Fold Change", y = "-log10(P-value)") +
  theme_bw(base_size = 13) +
  theme(plot.title      = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle   = element_text(hjust = 0.5, color = "grey40"),
        legend.position = "top")
ggsave("volcano_GSE68605.png", plot = p_volcano,
       width = 8, height = 6, dpi = 300, bg = "white")

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
cat("SCFA genes in expression matrix:",sum(all_scfa %in% rownames(expr_mat)),"of", length(all_scfa),"\n")
#Extract SCFA genes 
scfa_df          <- top_68605[top_68605$GeneSymbol %in% all_scfa, ]
scfa_df$Category <- gene_category[scfa_df$GeneSymbol]
scfa_df$sig      <- ifelse(scfa_df$adj.P.Val < 0.05, "Significant", "NS")
scfa_df          <- scfa_df[order(scfa_df$Category, scfa_df$P.Value), ]
cat("SCFA genes found in GSE68605:", nrow(scfa_df), "of", length(all_scfa), "\n")
print(scfa_df[, c("GeneSymbol", "Category", "logFC", "P.Value", "adj.P.Val", "sig")])
write.csv(scfa_df, "GSE68605_SCFA_results.csv", row.names = FALSE)
#SCFA-focused volcano plot
category_colors <- c(
  "FFA Receptors"       = "#E63946",
  "Transporters"        = "#F4A261",
  "Butyrate Metabolism" = "#2A9D8F",
  "HDAC Targets"        = "#457B9D",
  "NF-kB"               = "#6A0572",
  "NLRP3"               = "#E76F51",
  "Gut-Brain"           = "#264653"
)
top_68605$Category <- ifelse(top_68605$GeneSymbol %in% all_scfa,
                             gene_category[top_68605$GeneSymbol], "Other")
p_scfa_volcano <- ggplot() +
  geom_point(data = subset(top_68605, Category == "Other"),
             aes(x = logFC, y = -log10(P.Value)),
             color = "grey80", alpha = 0.4, size = 1.2) +
  geom_point(data = subset(top_68605, Category != "Other"),
             aes(x = logFC, y = -log10(P.Value), color = Category),
             size = 3.5, shape = 18) +
  geom_text(data = subset(top_68605, Category != "Other"),
            aes(x = logFC, y = -log10(P.Value),
                label = GeneSymbol, color = Category),
            size = 2.8, vjust = -0.8, fontface = "bold") +
  scale_color_manual(values = category_colors, name = "Gene Category") +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed",
             color = "black", alpha = 0.4) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed",
             color = "black", alpha = 0.4) +
  labs(title    = "SCFA-Related Genes in ALS vs Control (GSE68605)",
       subtitle = paste(
         "SCFA genes found:", nrow(scfa_df),"| Significant:", sum(scfa_df$sig == "Significant")),
       x = "Log2 Fold Change", y = "-log10(P-value)") +
  theme_bw(base_size = 13) +
  theme(plot.title      = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle   = element_text(hjust = 0.5, color = "grey40"),
        legend.position = "right",
        legend.text     = element_text(size = 9))
ggsave("volcano_SCFA_GSE68605.png", plot = p_scfa_volcano,
       width = 10, height = 6, dpi = 300, bg = "white")
#heatmap
found_genes <- rownames(expr_mat)[rownames(expr_mat) %in% all_scfa]
if (length(found_genes) == 0) {
  stop("No SCFA genes found in expression matrix.")}
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
pheatmap(
  scfa_expr,
  annotation_col    = col_annot,
  annotation_row    = row_annot,
  annotation_colors = annot_colors,
  scale             = "row",
  cluster_cols      = TRUE,
  cluster_rows      = FALSE,
  show_colnames     = FALSE,
  fontsize_row      = 9,
  color             = colorRampPalette(c("#457B9D", "white", "#E63946"))(100),
  main              = "SCFA Gene Expression: ALS vs Control (GSE68605)",
  filename          = "heatmap_SCFA_GSE68605.png",
  width             = 10, height = 8
)
#Comparing GSE68605 with GSE56500 to look for overlap
sig_56500 <- read.csv("../DEG/ALS_vs_Control_significant.csv")
overlap_genes <- intersect(sig_68605$GeneSymbol, sig_56500$GeneSymbol)
cat("\nOverlapping DEGs (GSE56500 ∩ GSE68605):", length(overlap_genes), "\n")
overlap_df <- sig_68605[sig_68605$GeneSymbol %in% overlap_genes,
                        c("GeneSymbol", "logFC", "P.Value", "adj.P.Val")]
colnames(overlap_df)[2:4] <- c("logFC_68605", "P.Value_68605", "adj.P.Val_68605")
overlap_56500 <- sig_56500[sig_56500$GeneSymbol %in% overlap_genes,
                           c("GeneSymbol", "logFC", "P.Value")]
overlap_final <- merge(overlap_df, overlap_56500, by = "GeneSymbol")
overlap_final <- overlap_final[order(overlap_final$P.Value_68605), ]
write.csv(overlap_final, "DEG_overlap_GSE56500_GSE68605.csv", row.names = FALSE)
print(head(overlap_final, 20))
#Identifies SCFA-related genes within DEG overlap
scfa_overlap <- overlap_genes[overlap_genes %in% all_scfa]
cat("SCFA genes in overlap:", length(scfa_overlap), "\n")
if (length(scfa_overlap) > 0) print(scfa_overlap)
save.image("GSE68605_analysis.RData")
