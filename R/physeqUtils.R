#' utility tools to work with phyloseq objects
#' 1. samdat_as_tibble()
#' 2.


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
