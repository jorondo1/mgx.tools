
#' Variance-stabilizing transformation
#' @param ps phyloseq object
#' @export

vst_ps_to_mx <- function(ps) {
  phyloseq_to_deseq2(ps, ~ 1) %>% # DESeq2 object
    DESeq2::estimateSizeFactors(., geoMeans = apply(
      DESeq2::counts(.), 1, function(x) exp(sum(log( x[x>0] )) / length(x)))) %>%
    DESeq2::varianceStabilizingTransformation(blind=TRUE, fitType = 'mean') %>% # VST
    SummarizedExperiment::assay(.) %>% t() %>%
    { .[. < 0] <- 0; . } # replace negatives by zeros
}
