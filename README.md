# SCFA-Related Genes and Neuroinflammation in ALS

Reproducible R pipeline testing whether short-chain fatty acid (SCFA) signalling is
linked to NF-κB/NLRP3-driven neuroinflammation in amyotrophic lateral sclerosis,
using two independent GEO datasets from different parts of the central nervous
system.

## Overview

The SCFA hypothesis is not supported at these sample sizes. A cross-tissue signal
I was not looking for turned out to be the main positive finding.

**Datasets:** GSE56500 (spinal cord, Affymetrix Human Exon 1.0 ST, 6 ALS +
6 control) and GSE68605 (motor cortex, Affymetrix HG-U133 Plus 2.0, 8 ALS +
3 control).

## Pipeline

1. **Data retrieval** — GEO files are downloaded once and stored locally with
   recorded checksums, so every run uses the same input data
2. **Gene annotation** — probe IDs are matched to gene symbols; probes that map to
   more than one gene are dropped as ambiguous
3. **Probe collapse** — where several probes measure the same gene, the most
   variable one is kept
4. **Differential expression** — `limma`, BH-corrected p-values, thresholds
   adj.P.Val < 0.05 and |logFC| > 0.5
5. **SCFA panel** — 33 curated genes across 7 categories (FFA receptors,
   transporters, butyrate metabolism, HDAC targets, NF-κB, NLRP3, gut-brain
   signalling)
6. **Pathway enrichment** — `clusterProfiler` GSEA on KEGG and GO, run 20 times
   with different random seeds and reported as a consensus
7. **Cross-tissue comparison** — testing whether the two datasets agree, at both
   gene and pathway level

## Key results

**No SCFA gene is significant in either tissue.** 32 of 33 panel genes were
detected in spinal cord and 30 of 33 in motor cortex, but none passed the
significance threshold. The closest were SIRT1 in spinal cord (adj.P.Val 0.105)
and HDAC4 in motor cortex (0.566). This held even under an earlier version of the
code that was biased toward calling genes significant
(`results/figures/SCFA_panel_null.png`).

**The two tissues change in the same direction.** Across 14,614 shared genes the
fold-changes correlate weakly but reliably (Spearman rho = 0.113, permutation
p = 1 × 10⁻⁴). Among the 63 genes significant in either dataset the correlation is
strong (rho = 0.740) — a shared disease signal on top of ordinary tissue
differences (`results/figures/cross_tissue_logFC_concordance.png`).

**Three genes are significant in both tissues, all moving the same way.**
Differential expression found 41 significant genes in spinal cord and 28 in motor
cortex. Three appear in both: SERPINA3 (+3.63 spinal cord / +4.70 motor cortex),
AQP1 (+3.23 / +3.43) and PRUNE2 (−2.82 / −2.01). Chance alone predicts 0.072
shared genes, so finding 3 is unlikely (hypergeometric p = 4.9 × 10⁻⁵), and the
overlap gets stronger as the threshold is relaxed to 0.10 and 0.20. At the
stricter 0.01 cutoff no overlap is possible, since spinal cord has only 3
significant genes in total.

**Pathway enrichment is asymmetric, and nothing is shared.** Spinal cord shows 27
KEGG pathways and 190 GO terms in the consensus set; motor cortex shows 0 and 5.
No KEGG pathway or GO term is significant in both tissues, so the agreement
between datasets is at gene level only.

All numbers are in `results/tables/`.

## Why this pipeline was rebuilt

I could not reproduce my own results from a fresh copy of the repository.
Rebuilding it to run end to end turned up five errors:

1. **Gene name parsing.** The old code decided whether a string was a gene symbol
   by checking it against a list of accession prefixes. That list was incomplete,
   so entries like `BC001082 // RP11-529I10.4` were kept as gene names. Valid
   Entrez mappings in GSE56500 went from 2,807 to 14,503.
2. **Probes matching several genes** were assigned to the first one that passed
   the filter, instead of being dropped as ambiguous.
3. **Pathway significance.** `gseKEGG()` filters its output on the raw p-value,
   not the corrected one, and the old code counted those rows as significant.
   Shared KEGG pathways went from 15 to 0.
4. **Pathway figures came from a single random seed.** Across 20 seeds the number
   of significant pathways ranges from 20 to 37 (KEGG, spinal cord) and 5 to 41
   (GO, motor cortex), which is why results are now reported as a consensus
   (`results/figures/seed_stability_counts.png`).
5. **The previous README** reported 6 significant SCFA genes in motor cortex.
   Those genes pass the uncorrected p-value but not the corrected one — the code
   had always used the corrected threshold, so this was a writing error.

One consequence: the cross-tissue gene overlap, previously reported as zero, is
actually 3. The old figure was an artefact of malformed gene names failing to
match between datasets.

The version presented at the HIPS Symposium is kept at tag
[`v1.0`](../../tree/v1.0).

## Limitations

The cohorts are small (n = 12 and n = 11), so the SCFA result means "not
detectable here", not "not there". The exon array measures several probes per gene
and keeping only the most variable one discards some information. The two datasets
come from different platforms and no batch correction was applied, though
comparing fold-changes rather than raw values reduces this problem. This is a
computational analysis only and would need experimental follow-up.

## Reproducibility

Dependencies are managed with `renv`. Run `renv::restore()` after cloning, then
`Rscript scripts/00_run_all.R` to regenerate every table and figure.

All 24 result tables come out byte-identical from a fresh clone, including a full
20-seed GSEA sweep run from scratch.

R 4.6.0, package versions pinned in `renv.lock` (limma 3.68.2, clusterProfiler
4.20.0, GEOquery 2.80.0, org.Hs.eg.db 3.23.1, ggplot2 4.0.3). KEGG data retrieved
2026-08-15. Input file checksums in `results/data_checksums.txt` — the pipeline
stops if an input file has changed.

**Note on clusterProfiler 4.20.0:** `gseKEGG()` and `gseGO()` accept a `seed`
argument but never pass it on, so setting it has no effect — `set.seed()` has to be
called immediately before each call instead. `pvalueCutoff` also filters on the raw
p-value rather than the corrected one.

## Repository structure

```
scripts/    00_run_all.R    runs the whole pipeline
            01–06           DEG, SCFA panel, GSEA, cross-tissue analysis
            _*.R            paths, helper functions, GEO download
            seed_stability.R
results/    tables/         all numeric results
            figures/        volcano, SCFA panel, GSEA, concordance plots
renv.lock   pinned package versions
```

## Tools

R, limma, GEOquery, clusterProfiler, org.Hs.eg.db, ggplot2, ggrepel, pheatmap

## References

1. The Links between ALS and NF-κB —
   https://pmc.ncbi.nlm.nih.gov/articles/PMC8070122/
2. Elevated NLRP3 Inflammasome Activation Is Associated with Motor Neuron
   Degeneration in ALS — https://pmc.ncbi.nlm.nih.gov/articles/PMC11202041/
3. The Role of Short-Chain Fatty Acids in Microbiota–Gut–Brain Cross-Talk with a
   Focus on ALS — https://pmc.ncbi.nlm.nih.gov/articles/PMC10606032/
4. The emerging role of microbiota-derived SCFAs in neurodegenerative disorders —
   https://pmc.ncbi.nlm.nih.gov/articles/PMC12152874/
5. Mechanisms of Blood–Brain Barrier Protection by Microbiota-Derived Short-Chain
   Fatty Acids — https://pmc.ncbi.nlm.nih.gov/articles/PMC9954192/

Presented at HIPS Symposium, Young Investigators Talk, May 2026.
