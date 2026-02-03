#' utility tools to work with phyloseq objects
#' 1. samdat_as_tibble()
#' 2.
#'

#' Keep samples with metadata
#' @param seqtab matrix sequence table with samples as rows
#' @param samples vector of sample names
#' @export
subset_samples <- function(seqtab, samples) {
  seqtab[rownames(seqtab) %in% samples, ] %>% # subset
    .[, colSums(.) > 0] # Remove taxa with no hits
}

#' ASVs classified at the kingdom level and present in seqtab
#' @param taxonomy taxonomic table with unique ID as rows (e.g. OTU, ASV, Species)
#' @param seqtab matrix sequence table with samples as rows
#' @param min_seq taxa with fewer than min_seq overall (sum of all samples) are removed
#' @export
subset_tax_table <- function(taxonomy, seqtab, min_seq) {
  if (!is.data.frame(taxonomy)) {
    taxonomy <- as.data.frame(taxonomy)
  }

  keep <- colSums(seqtab) >= min_seq
  taxa_phylum <- subset(taxonomy, Phylum != "Unclassified") %>% # subset needs the input to be a df
    rownames()
  taxa <- taxa_phylum
    intersect(
      colnames(seqtab)[keep] # only keep taxa still present in seqtab
    )

  message(paste(length(keep)-length(taxa), 'taxa were removed'))
  message(paste('Of these,',length(taxa_phylum),'were unclassified at the Phylum level'))
  message('Consider using remove_ultra_rare() next.')
  as.matrix(taxonomy[taxa, ])
}

#' Remove samples with fewer than n sequences once taxa removed
#' @param seqtab matrix sequence table with samples as rows
#' @param taxonomy taxonomic table with unique ID as rows (e.g. OTU, ASV, Species)
#' @param n samples with fewer than n sequences are removed
#' @export
remove_ultra_rare <- function(seqtab, taxonomy, n) {
  message('Make sure to use subset_taxa() first.')
  if (!is.data.frame(taxonomy)) {
    taxonomy <- as.data.frame(taxonomy)
  }

  result <- seqtab[, rownames(taxonomy), drop = FALSE]
  result <- result[rowSums(result) > n, , drop = FALSE]  # Filter rows (samples). n = sum across ASVs in a sample
  result <- result[,colSums(result) > 1, drop = FALSE] # Remove singleton ASVs

  removed <- setdiff(rownames(seqtab), rownames(result))
  num_removed <- length(removed)
  if(num_removed>0){
    message("Samples removed:")
    print(removed)
  } else {
    message("No samples removed.")
  }
  message('Consider using viz_seqdepth() next.')
  return(result)
}

#' Visualise sequence count distribution across samples
#' @param seqtab matrix sequence table with samples as rows
#' @export
viz_seqdepth <- function(seqtab, breaks = 100, log_count = FALSE) {
  sums <- rowSums(seqtab)
  if(log_count){
    sums <- log10(sums)
    lab <- "Sample size (log10)"}
  else {lab <- "Sample size"}
  hist(sums, breaks = breaks,
       xlab = lab, xaxt = "n", main = 'Distribution of sequence count per sample')
  axis(1, at = pretty(sums, n = 40))  # adding ticks
}

#' 1. Get a phyloseq object's sample data as tibble
#' Creates a Sample column with the sample names
#' @params ps a phyloseq object
#' @params strings_as_factors convert strings to factors
#' @export
samdat_as_tibble <- function(ps, strings_as_factors = TRUE){
  phyloseq::sample_data(ps) %>%
    data.frame() %>%
    tibble::rownames_to_column('Sample') %>%
    tibble::tibble() %>%
    {
      if (strings_as_factors) {
        dplyr::mutate(., dplyr::across(where(is.character), as.factor))
      } else {
        .
      }
    }
}
