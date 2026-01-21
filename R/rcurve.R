#' stolen from https://github.com/Russel88/MicEco/blob/master/R/rcurve.R
#' minor adjustments to trnsform otu as matrix
#' @export
rcurve <- function (physeq, subsamp = 10^c(1:5), trim = TRUE, add_sample_data = TRUE)
{
  otu <- otu_table(physeq)
  if (!taxa_are_rows(physeq)) {
    otu <- t(otu)
  }
  otu <- as(otu, 'matrix') # This is required otherwise vegan says: Error in as(x, "matrix")[i, j, drop = FALSE] :
  otu <- round(otu) # Sourmash may have non-integer abundances
  colS <- colSums(otu)
  pb <- txtProgressBar(min = 0, max = length(subsamp), style = 3)
  rars <- list()
  for (i in seq_along(subsamp)) {
    setTxtProgressBar(pb, i)
    rars[[i]] <- vegan::rarefy(otu, sample = subsamp[i],
                               MARGIN = 2)
  }
  mat <- do.call(cbind, rars)
  if (trim) {
    mat_bool <- sapply(subsamp, function(i) sapply(colS,
                                                   function(j) i <= j))
    mat_new <- mat * mat_bool
    mat_new[mat_new == 0] <- NA
  }
  else {
    mat_new <- mat
  }
  colnames(mat_new) <- subsamp
  df <- as.data.frame.table(mat_new)
  colnames(df) <- c("Sample", "Reads", "Richness")
  df$Reads <- as.numeric(as.character(df$Reads))
  if (trim) {
    df <- na.omit(df)
  }
  if (add_sample_data) {
    samp <- sample_data(physeq)
    df2 <- merge(df, samp, by = "Sample", by.y = "row.names")
    return(df2)
  }
  else {
    return(df)
  }
}
