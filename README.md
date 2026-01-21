Helper functions for microbiome research. 

This package in development contains functions to:
- Do postprocessing on ASV/sequence tables when finishing the DADA2 pipeline, before creating a phyloseq object
- `samdat_as_tibble()` extracts a phyloseq object's sample_data component as a tibble (should be simple, right?)
- Do phyloseq stuff faster (multithreaded `psmelt()` `rarefy_even_depth()`)
- `tax_glom2()` : a rewritten `phyloseq::tax_glom()` that updates taxa names to reflect the agglomeration rank and removes counts from taxa with missing classifications at that rank.
- `rarecurve()` : relatively fast rarefaction curves
- `assemblePhyloseq()`: Helper functions to assemble phyloseq objects from the output of metaphlan and kraken
- `community_viz()`: a plotting function for sample-wise visualisation of microbiome compositions in a bar chart at desired taxonomic levels
- a bunch of functions to compute diversity in different ways using a phyloseq object

Some code is shamelessly stolen from other repositories to be improved upon. I know, I should open an issue or branch or whatever, I just did not get around to it. 
