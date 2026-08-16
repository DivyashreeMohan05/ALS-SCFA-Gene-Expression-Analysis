# DEG Analysis – GSE56500 (Spinal Cord, GPL5188 Affymetrix Human Exon 1.0 ST)
#ALS subtypes: csALS + c9ALS combined vs Control
#Methods used: limma | Threshold: adj.P.Val < 0.05, |logFC| > 0.5
library(GEOquery)
library(limma)
library(ggplot2)
library(ggrepel)
source(here::here("scripts", "_paths.R"))
source(here::here("scripts", "_fetch_geo.R"))
source(here::here("scripts", "_helpers.R"))
# Loading GEO dataset
gse       <- getGEO("GSE56500", GSEMatrix = TRUE, destdir = DIR_RAW)[[1]]
gpl       <- annotation(gse)
expr_mat  <- exprs(gse) #expression matrix of probes × samples
pheno     <- pData(gse)
fdata     <- fData(gse) #feature annotation
probes_total <- nrow(fdata)
#Creating a ALS/Control label
pheno$ALS_status <- ifelse(pheno$`patient group:ch1` == "control", "Control", "ALS")
#Probe to gene symbol mapping - /// separated, keep symbol only if unambiguous
resolve_symbol <- function(symbols) {
  symbols <- unique(symbols[!is.na(symbols) & symbols != ""])
  if (length(symbols) == 1) return(symbols)
  NA
}
extract_gene <- function(x) {
  if (is.na(x)) return(NA)
  blocks <- strsplit(x, " /// ")[[1]]
  fields2 <- character(0)
  for (b in blocks) {
    f <- strsplit(b, " // ")[[1]]
    if (length(f) >= 2) fields2 <- c(fields2, trimws(f[2]))
  }
  resolve_symbol(fields2)
}
#Diagnostic classification only - mirrors extract_gene, doesn't feed the pipeline
classify_probe <- function(x) {
  if (is.na(x)) return("no_symbol")
  blocks <- strsplit(x, " /// ")[[1]]
  fields2 <- character(0)
  for (b in blocks) {
    f <- strsplit(b, " // ")[[1]]
    if (length(f) >= 2) fields2 <- c(fields2, trimws(f[2]))
  }
  symbols <- unique(fields2[!is.na(fields2) & fields2 != ""])
  if (length(symbols) == 0) "no_symbol" else if (length(symbols) == 1) "resolved" else "multimapped"
}
probe_class <- table(factor(sapply(fdata$gene_assignment, classify_probe),
                             levels = c("no_symbol", "resolved", "multimapped")))
fdata$GeneSymbol  <- sapply(fdata$gene_assignment, extract_gene)
valid             <- !is.na(fdata$GeneSymbol)
expr_mat          <- expr_mat[valid,]
fdata             <- fdata[valid,]
#Probes combined to one row per gene: highest IQR wins, ties broken by probe ID
iqr_vals   <- apply(expr_mat, 1, IQR)
probe_id   <- rownames(fdata)
keep       <- tapply(seq_along(iqr_vals), fdata$GeneSymbol,
                      function(i) i[order(-iqr_vals[i], probe_id[i])[1]])
expr_mat   <- expr_mat[keep, ]
fdata      <- fdata[keep, ]
rownames(expr_mat) <- fdata$GeneSymbol
cat("Probes:", length(iqr_vals), "-> Genes retained after collapse:", nrow(expr_mat), "\n")
update_preprocessing_summary("GSE56500",
  gpl = gpl, platform = "Affymetrix Human Exon 1.0 ST",
  n_als = sum(pheno$ALS_status == "ALS"), n_control = sum(pheno$ALS_status == "Control"),
  probes_total = probes_total,
  probes_dropped_no_symbol = unname(probe_class["no_symbol"]),
  probes_dropped_multimapped = unname(probe_class["multimapped"]),
  probes_retained = length(iqr_vals), genes_after_collapse = nrow(expr_mat))
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
#Significance threshold - BH-adjusted
sig_ALS <- subset(top_ALS,adj.P.Val < 0.05 & abs(logFC) > 0.5)
print(paste("Total significant DEGs:", nrow(sig_ALS)))
#Saving results
write.csv(top_ALS, file.path(DIR_TABLES, "ALS_vs_Control_all.csv"), row.names = FALSE)
write.csv(sig_ALS, file.path(DIR_TABLES, "ALS_vs_Control_significant.csv"), row.names = FALSE)
#Volcano plot
top_ALS$sig <- ifelse(top_ALS$adj.P.Val < 0.05 & abs(top_ALS$logFC) > 0.5,
                      "Significant", "NS")
top10 <- top_ALS[order(top_ALS$adj.P.Val), ][1:10, ]
p <- ggplot(top_ALS, aes(x = logFC, y = -log10(adj.P.Val), color = sig)) +
  geom_point(alpha = 0.6, size = 1.8) +
  geom_text_repel(data = top10, aes(label = GeneSymbol), color = "black",
                   size = 3, max.overlaps = Inf) +
  scale_color_manual(values = c("NS" = "grey70", "Significant" = "red"), name = "") +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "blue", alpha = 0.5) +
  geom_hline(yintercept = -log10(0.05),  linetype = "dashed", color = "blue", alpha = 0.5) +
  annotate("text", x = Inf, y = Inf,
           label = paste("Significant n =", sum(top_ALS$sig == "Significant")),
           hjust = 1.1, vjust = 1.5, size = 4.5, color = "red") +
  labs(title = "ALS vs Control (GSE56500)",
       subtitle = "csALS + c9ALS combined",
       x = "Log2 Fold Change", y = "-log10(adj.P.Val)") +
  theme_bw(base_size = 13) +
  theme(plot.title    = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, color = "grey40"),
        legend.position = "top")
set.seed(42)
ggsave(file.path(DIR_FIGURES, "volcano_ALS_combined.png"), plot = p,
       width = 8, height = 6, dpi = 300, bg = "white")
#Saving objects needed by 02_SCFA_Analysis.R
saveRDS(list(expr_mat = expr_mat, top_ALS = top_ALS, pheno = pheno),
        file.path(DIR_INTERIM, "ALS_DEG_final.rds"))
