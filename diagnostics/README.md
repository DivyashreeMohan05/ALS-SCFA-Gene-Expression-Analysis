# Diagnostics

Scripts here are investigation, not pipeline. Nothing in `scripts/00_run_all.R`
sources this directory, and nothing here is required to reproduce `results/`.

## samplesize_test.R

`gseKEGG()`/`gseGO()` (via `enrichit::gsea_gson`, v0.1.4, as shipped with
clusterProfiler 4.20.0) silently drop `seed` and `nPermSimple` - neither is a
named formal of `gseKEGG`/`gseGO`, so both land in `...`, which is never
forwarded to `gsea_gson()`. `nPerm` *is* forwarded, but only matters for
`method %in% c("sample","permute")`; the default `method="multilevel"` never
reads it. The actual precision knob for the multilevel method is
`sampleSize`, passed straight into the internal C++ engine
(`enrichit:::gsea_multilevel_cpp`) - but `sampleSize` isn't a named formal of
`gsea_gson()` either, so it's unreachable from `gseKEGG()`/`gseGO()` the same
way `seed` is.

`samplesize_test.R` confirms this by calling the unexported `enrichit:::gsea()`
directly (bypassing `gseKEGG()`) with `sampleSize` raised from the default 101
to 10001: GSE56500 KEGG goes from a 29-37 spread across seeds 42/1/7 down to a
flat 30 at every seed.

This is not something to fix by calling `enrichit:::gsea()` in the pipeline -
`:::` reaches into an unexported internal of a dependency, which can change or
disappear on any package update with no warning, undermining the whole point
of a reproducible pipeline. It's also enrichit-version-specific: nothing
guarantees `sampleSize` stays unreachable-but-present in a future release, or
that the internal function keeps the same name.

**This is why `04` reports a consensus (significant in >=80% of 20 seeds)
instead of a single-seed count**: the public API gives no way to raise
precision enough to make a single-seed count stable, so the pipeline reports
across-seed agreement instead of pretending one seed's number is the answer.
