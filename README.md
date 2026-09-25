
<!-- README.md is generated from README.Rmd. Please edit that file -->

# mgx.tools

Helper functions for amplicon and shotgun metagenomics in R, from raw
reads to a [phyloseq](https://joey711.github.io/phyloseq/) object and on
to diversity analyses. Most functions are convenience wrappers,
parallelised rewrites of existing functions, or small utilities that
help make informed filtering decisions.

The package is under active development and function names and arguments
may still change.

## What’s in it

### DADA2 pipeline helpers.

Support for each step of the [DADA2](https://benjjneb.github.io/dada2/)
workflow: checking primer orientation, running cutadapt in parallel,
quality plots, chimera and rare-ASV reports, and formatting DECIPHER
taxonomy for DADA2. Includes a read-tracking system that logs counts to
disk after each step and plots where reads were lost, and checkpoints to
save and resume the pipeline.

### Building and tidying phyloseq objects.

Assembling a phyloseq object from shotgun taxonomic profiles
(e.g. MetaPhlAn, Kraken), subsetting samples and taxa, dropping
ultra-rare taxa, inspecting sequencing depth and sparsity.

### Phyloseq function rewrites.

Alternatives to `rarefy_even_depth()` (parallelised and reproducible),
`psmelt()` (much faster, drops zero-abundance rows) and `tax_glom()`
(names agglomerated taxa after their rank instead of an arbitrary ASV).

### Diversity.

Iterative rarefaction to estimate alpha diversity (Hill numbers, Faith
PD) and beta diversity (Bray-Curtis, robust Aitchison, UniFrac) in a
single call, plus rarefaction curves, variance-stabilizing
transformation and PCoA.

### Taxonomy and phylogeny.

Building ASV phylogenetic trees, clustering ASVs by phylogenetic
distance to bridge the gap between genus and ASV, and computing
classification rates across taxonomic ranks.

### Community composition plots.

Preparing data for stacked barplots of relative abundances averaged over
selected metadata

## Installation

``` r
# install.packages("pak")
pak::pkg_install("jorondo1/mgx.tools")
```

`pak` installs the required dependencies automatically, including the
Bioconductor ones. Be aware that mgx.tools has many dependencies (around
30, including DADA2 and phyloseq), so the first installation can take a
while.

Some functions rely on additional packages that are not installed by
default: tree building (DECIPHER, phangorn), phylogenetic diversity
(phangorn, GUniFrac, picante), variance-stabilizing transformation
(DESeq2), and `psflashmelt()` (data.table). If one is missing when you
call such a function, you will be offered to install it. To install them
all upfront:

``` r
pak::pkg_install("jorondo1/mgx.tools", dependencies = TRUE)
```

Some DADA2 helpers also require
[cutadapt](https://cutadapt.readthedocs.io/) to be installed on your
system.

## Citation

If you use mgx.tools in your work, please cite:

Rondeau-Leclaire J. mgx.tools: Helper functions for metagenomics and
amplicon sequencing processing. <https://github.com/jorondo1/mgx.tools>
<!-- TODO: add Zenodo DOI once v0.1.1 is archived -->
