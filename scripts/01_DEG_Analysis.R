# DEG Analysis – GSE56500 (Spinal Cord, RNA-seq)
#ALS subtypes: csALS + c9ALS combined vs Control
#Methods used: limma | Threshold: adj.P.Val < 0.05, |logFC| > 0.5
library(GEOquery)
library(limma)
library(ggplot2)
# Loading GEO dataset
gse       <- getGEO("GSE56500", GSEMatrix = TRUE)[[1]]
expr_mat  <- exprs(gse) #expression matrix of probes × samples 
pheno     <- pData(gse)
fdata     <- fData(gse) #feature annotation
#Creating a ALS/Control label
pheno$ALS_status <- ifelse(pheno$`patient group:ch1` == "control", "Control", "ALS")
#Probe to Gene symbol mapping
extract_gene <- function(x) {
  if (is.na(x) || x == "---") return(NA)
  # trying /// separator, then // 
  if (grepl(" /// ", x)) {
    parts <- strsplit(x, " /// ")[[1]]
  } else if (grepl(" // ", x)) {
    parts <- strsplit(x, " // ")[[1]]
  } else {
    return(trimws(x))
  }
  #trying out part 2 
  for (i in 2:length(parts)) {
    gene <- trimws(parts[i])
    if (!is.na(gene) && gene != "" && gene != "---" && 
        !grepl("^[0-9]", gene) && !grepl("^NM_|^NR_|^XM_|^ENST|^uc", gene)) {
      return(gene)
    }}
  gene <- trimws(parts[1])
  if (!is.na(gene) && gene != "" && gene != "---") return(gene)
  return(NA)
}
fdata$GeneSymbol  <- sapply(fdata$gene_assignment, extract_gene)
valid             <- !is.na(fdata$GeneSymbol)
expr_mat          <- expr_mat[valid,]
fdata             <- fdata[valid,]
#Probes combined to one row per gene (IQR-based selection)
iqr_order  <- order(apply(expr_mat, 1, IQR), decreasing = TRUE)
expr_mat   <- expr_mat[iqr_order, ]
fdata      <- fdata[iqr_order, ]
keep       <- !duplicated(fdata$GeneSymbol)
expr_mat   <- expr_mat[keep, ]
rownames(expr_mat) <- fdata$GeneSymbol[keep]
cat("Genes retained after probe collapse:", nrow(expr_mat), "\n")
# Inspecting for non-standard names 
suspicious <- grep("^[0-9]|\\.|^-", rownames(expr_mat), value = TRUE)
if (length(suspicious) > 0) {
  cat("WARNING: Non-standard gene names detected:", length(suspicious), "\n")
  print(head(suspicious, 10))
}
#Differential gene expression
group   <- factor(pheno$ALS_status, levels = c("Control", "ALS"))
design  <- model.matrix(~0 + group)
colnames(design) <- c("Control", "ALS")
fit  <- lmFit(expr_mat, design)
fit2 <- contrasts.fit(fit, makeContrasts(ALS_vs_Control = ALS - Control, levels = design))
fit2 <- eBayes(fit2, trend = TRUE)
#Results table
top_ALS <- topTable(fit2, coef = "ALS_vs_Control", adjust = "BH", number = Inf)
top_ALS$GeneSymbol <- rownames(top_ALS)
#significance threshold (Using BH-adjusted p-values to control fdr)
sig_ALS <- subset(top_ALS,adj.P.Val < 0.05 & abs(logFC) > 0.5)
print(paste("Total significant DEGs:", nrow(sig_ALS)))
#Saving results
write.csv(top_ALS, "ALS_vs_Control_all.csv",row.names = FALSE)
write.csv(sig_ALS, "ALS_vs_Control_significant.csv", row.names = FALSE)
#Volcano plot
top_ALS$sig <- ifelse(top_ALS$adj.P.Val < 0.05 & abs(top_ALS$logFC) > 0.5,
                      "Significant", "NS")
p <- ggplot(top_ALS, aes(x = logFC, y = -log10(P.Value), color = sig)) +
  geom_point(alpha = 0.6, size = 1.8) +
  scale_color_manual(values = c("NS" = "grey70", "Significant" = "red"), name = "") +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "blue", alpha = 0.5) +
  geom_hline(yintercept = -log10(0.05),  linetype = "dashed", color = "blue", alpha = 0.5) +
  annotate("text", x = Inf, y = Inf,
           label = paste("Significant n =", sum(top_ALS$sig == "Significant")),
           hjust = 1.1, vjust = 1.5, size = 4.5, color = "red") +
  labs(title = "ALS vs Control (GSE56500)",
       subtitle = "csALS + c9ALS combined",
       x = "Log2 Fold Change", y = "-log10(P-value)") +
  theme_bw(base_size = 13) +
  theme(plot.title    = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, color = "grey40"),
        legend.position = "top")
ggsave("volcano_ALS_combined.png", plot = p,
       width = 8, height = 6, dpi = 300, bg = "white")
#Saving the entire workspace 
save.image("ALS_DEG_final.RData")
