
<!-- README.md is generated from README.Rmd. Please edit that file -->

# mgx.tools

Helper functions for metagenomics and amplicon sequencing processing in
R. Goes from raw reads to phyloseq object, with parallelized rewrites of
common functions and utilities for informed filtering decisions.

## Dependencies

This package relies on several **Bioconductor** packages that must be
installed manually before installing `mgx.tools`:

``` r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c(
    "Biostrings",
    "dada2",
    "DESeq2",
    "phyloseq",
    "ShortRead",
    "SummarizedExperiment"
))
```

## Installation

Once Bioconductor dependencies are installed:

``` r
# install.packages("pak")
pak::pkg_install("jorondo1/mgx.tools")
# pak::pkg_install("jorondo1/mgx.tools?reinstall", upgrade = TRUE)
```

## Overview

| Function              | Description                                   |
|-----------------------|-----------------------------------------------|
| `primer_occurence()`  | Check primer orientations across reads        |
| `run_cutadapt()`      | Cutadapt wrapper for use with `mclapply`      |
| `dropped_samples()`   | Identify samples lost after filtering         |
| `track_dada()`        | Compile read counts across DADA2 steps        |
| `plot_track_change()` | Visualise read loss across pipeline steps     |
| `gen_qplots()`        | Generate quality plots in parallel            |
| `save_qplots()`       | Export quality plots to PDF                   |
| `chimera_report()`    | Report chimera rates                          |
| `minimum_ASV_count()` | Plot ASV counts by minimum sequence threshold |
| `drop_rare_asvs()`    | Remove rare ASVs and empty samples            |
| `asv_to_fasta()`      | Export ASV sequences as FASTA                 |
