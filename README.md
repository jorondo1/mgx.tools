## Helper functions for microbiome bioinformatics. 

This package in development contains functions to:

- Modularize the DADA2 pipeline, essentially wrapping up some steps into functions like `primer_occurence()`, `.allOrients()`; a wrapper function to run cutadapt to facilitate use with `mclapply()`; and cleaner ways to track reads and ASVs through the pipeline
- Do postprocessing on ASV/sequence tables when finishing the DADA2 pipeline, before creating a phyloseq object
- `samdat_as_tibble()` extracts a phyloseq object's sample_data component as a tibble (should be simple, right?)
- Do phyloseq stuff faster (multithreaded `psmelt()` `rarefy_even_depth()`)
- `tax_glom2()` : a rewritten `phyloseq::tax_glom()` that updates taxa names to reflect the agglomeration rank and removes counts from taxa with missing classifications at that rank.
- `rarecurve()` : relatively fast rarefaction curves
- `assemblePhyloseq()`: Helper functions to assemble phyloseq objects from the output of metaphlan and kraken
- `community_viz()`: a plotting function for sample-wise visualisation of microbiome compositions in a bar chart at desired taxonomic levels
- a bunch of functions to compute diversity in different ways using a phyloseq object
- `compute_pcoa` takes a phyloseq object and leverages vegan to produce an object containing a distance matricx, an sample data table with the first two PCo, and a data.frame of eigenvalues.
  
Some code is shamelessly stolen from other repositories to be improved upon. I know, I should open an issue or branch or whatever, I just did not get around to it. 
