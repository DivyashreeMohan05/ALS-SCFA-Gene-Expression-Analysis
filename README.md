# In Silico Gene Expression Analysis of SCFA-Related Genes and NF-κB/NLRP3 Signaling in ALS Neuroinflammation

**Divyashree Mohan**  
MSc Bioinformatics, Saarland University  
*(Project initiated during BSc Hons Biomedical Sciences, SRIHER, India)*

---

## Overview

This project investigates whether disrupted short-chain fatty acid (SCFA) signaling is molecularly linked to NF-κB/NLRP3-driven neuroinflammation in amyotrophic lateral sclerosis (ALS). Using two publicly available GEO datasets from distinct CNS tissues, this analysis performs differential gene expression, SCFA gene filtering, and GSEA pathway enrichment entirely in R.

---

## Research Question

> Can in silico gene expression analysis reveal a molecular link between disrupted SCFA signaling and NF-κB/NLRP3-driven neuroinflammation in ALS?

---

## Datasets
Datasets were downloaded from GEO
GSE56500 | Spinal cord |  6 ALS (csALS + c9ALS) + 6 Control 

GSE68605 | Motor cortex |  8 ALS + 3 Control 

---

## Methods

- **DEG analysis:** limma with eBayes (trend = TRUE); BH-adjusted p-values; threshold adj.P.Val < 0.05, |logFC| > 0.5
- **SCFA gene panel:** 33 curated genes across 7 functional categories — FFA Receptors, Transporters, Butyrate Metabolism, HDAC Targets, NF-κB, NLRP3, Gut-Brain signaling
- **Pathway enrichment:** clusterProfiler gseKEGG + gseGO (Biological Process); ranking metric = logFC × −log10(P-value); nPermSimple = 10,000
- **Visualisation:** ggplot2, pheatmap

---

## Key Results

### Differential Gene Expression

**GSE56500** — broad downregulation pattern across both csALS and c9ALS subtypes.  
**GSE68605** — highly upregulated outliers pointing to active neuroinflammatory gene programs.


### SCFA Gene Dysregulation

32 of 33 curated SCFA genes were detected in GSE56500 (spinal cord); 0 reached significance after BH correction, consistent with limited detection power in a small RNA-seq cohort. In GSE68605 (motor cortex), 30 of 33 SCFA genes were detected with 6 significantly dysregulated.

Key findings from GSE68605:
- **MYD88, CX3CR1, TREM2** upregulated → gut-brain neuroinflammation
- **HDAC4, SIRT1** upregulated → impaired butyrate-mediated HDAC inhibition
- **NFKBIA** ↑ + **SLC5A8** ↓ → NF-κB activation + reduced SCFA transport
- **FFAR2, FFAR4** dysregulated → reduced SCFA receptor signalling
- **HADHA, HADHB, ACAT1** dysregulated → impaired mitochondrial butyrate oxidation

---

### GSEA KEGG Pathway Enrichment

**GSE56500 — Spinal Cord**
Suppressed: Oxidative phosphorylation, Parkinson disease, ALS pathway, TCA cycle  
Activated: Cytokine-cytokine receptor interaction, Complement and coagulation cascades, Hematopoietic cell lineage

**GSE68605 — Motor Cortex**
Suppressed: Parkinson disease, Sphingolipid metabolism, Motor proteins  
Activated: Complement and coagulation cascades, Cytokine signalling, PI3K-Akt signalling

---

### Cross-tissue Replicated KEGG Pathways

15 KEGG pathways were significantly enriched in both spinal cord and motor cortex, providing cross-tissue validation:


 Complement and coagulation cascade - Neuroinflammation — directly linked to SCFA suppression 
 Cytokine-cytokine receptor interaction - NF-κB/NLRP3 driven inflammatory signalling 
 PI3K-Akt signaling pathway - SCFA receptor (FFAR2/4) downstream signalling 
 Parkinson disease - Shared neurodegeneration signature 
 Motor proteins - ALS motor neuron degeneration 
> Suppressed mitochondrial function + activated neuroinflammation = fingerprint of SCFA disruption in ALS

---

## Conclusions

- SCFA-mediated HDAC inhibition and NF-κB signaling are disrupted in ALS motor cortex
- Cross-tissue replication of complement cascades, cytokine signalling, and PI3K-Akt pathways provides converging evidence linking SCFA receptor signalling disruption to neuroinflammation
- Findings support gut-brain axis involvement in ALS disease neuroinflammation
- Restoring SCFA signalling via butyrate supplementation or microbiota modulation is a rational therapeutic avenue requiring experimental validation

---


## Limitations

- Small sample sizes (GSE56500: n=12; GSE68605: n=11) limit statistical power
- GSE56500 spinal cord RNA-seq has limited probe coverage for gut-enriched SCFA receptor genes
- Cross-tissue comparison is exploratory; batch effects between platforms not corrected
- In silico analysis only — experimental validation required

---

## References

1. The Links between ALS and NF-κB 
https://pmc.ncbi.nlm.nih.gov/articles/PMC8070122/
2. Elevated NLRP3 Inflammasome Activation Is Associated with Motor Neuron Degeneration in ALS
https://pmc.ncbi.nlm.nih.gov/articles/PMC11202041/
3. The Role of Short-Chain Fatty Acids in Microbiota–Gut–Brain Cross-Talk with a Focus on ALS https://pmc.ncbi.nlm.nih.gov/articles/PMC10606032/
4. The emerging role of microbiota derived SCFAs in neurodegenerative disorders https://pmc.ncbi.nlm.nih.gov/articles/PMC12152874/
5. Mechanisms of Blood–Brain Barrier Protection by Microbiota-Derived Short-Chain Fatty Acids
https://pmc.ncbi.nlm.nih.gov/articles/PMC9954192/

---

*Presented at HIPS Young Investigators Symposium, May 2026*
