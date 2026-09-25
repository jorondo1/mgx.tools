# mgx.tools 0.1.0

* First tagged release.
* Helpers for going from raw amplicon/shotgun reads to a phyloseq object:
  primer checks, cutadapt wrapper, DADA2 read tracking, quality plots,
  chimera reports, rare-ASV filtering, ASV trees and clustering.
* Parallelised rewrites: `rarefy_even_depth2()`, `tax_glom2()`, `psflashmelt()`,
  `rcurve()`, plus iterative rarefaction of alpha/beta diversity with
  `rarefy_diversity()`.
* `btools` (Faith PD in `rarefy_diversity()`), `data.table` (`psflashmelt()`)
  and `kableExtra` are optional (Suggests).
