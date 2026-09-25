#' @keywords internal
"_PACKAGE"

# Functions called without a namespace prefix throughout the package.
# Importing them here means they work without the user attaching
# phyloseq, dplyr or ggplot2 first (and `filter()` resolves to dplyr, not stats).
#' @importFrom phyloseq access ntaxa otu_table phy_tree phyloseq_to_deseq2
#'   rank_names sample_data sample_data<- sample_sums sample_variables
#'   tax_table tax_table<- taxa_are_rows taxa_names taxa_sums
#' @importFrom dplyr across all_of case_when desc filter mutate n select where
#' @importFrom ggplot2 aes theme_minimal xlim
#' @importFrom dada2 getSequences getUniques
#' @importFrom purrr map
#' @importFrom rlang sym
#' @importFrom stringr str_detect
#' @importFrom methods as is
#' @importFrom stats aggregate as.dist cmdscale cutree na.omit setNames
#' @importFrom utils head setTxtProgressBar txtProgressBar
#' @importFrom graphics axis hist
#' @importFrom grDevices colorRampPalette
NULL

# Column names used in tidy evaluation, plus objects that run_cutadapt() and
# track_dada() deliberately read from the calling environment.
utils::globalVariables(c(
  ".", "Abundance", "aggTaxo", "asv_count", "Class", "Classification_rate",
  "classified", "denoisedF", "denoisedR", "Family", "filtered", "Genus", "i",
  "input", "Kingdom", "label", "nonchim", "Order", "OTU", "Phylum",
  "raw_seqtab", "reads", "relAb", "removeNs", "Sample", "sequence_sum",
  "Species_cluster", "step", "tax1", "taxColour", "taxRank", "values",
  "variable", "x", "xend", "y", "yend",
  "fnFs.cut", "fnRs.cut", "fnFs.filtN", "fnRs.filtN", "seqtab"
))
