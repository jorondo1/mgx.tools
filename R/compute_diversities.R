#' Takes a phyloseq object as input and creates a tibble containing all the sample data and one column per index.
#' @param physeq a phyloseq object
#' @export

diversity_tibble <- function(ps) {
  result <- mgx.tools::samdat_as_tibble(ps)
  
  result <- result %>%
    dplyr::left_join(tibble::enframe(mgx.tools::estimate_Hill(ps, q = 0), name = "Sample", value = "Hill_0"), by = "Sample") %>%
    dplyr::left_join(tibble::enframe(mgx.tools::estimate_Hill(ps, q = 1), name = "Sample", value = "Hill_1"), by = "Sample") %>%
    dplyr::left_join(tibble::enframe(mgx.tools::estimate_Hill(ps, q = 2), name = "Sample", value = "Hill_2"), by = "Sample") %>%
    dplyr::left_join(tibble::enframe(mgx.tools::estimate_diversity(ps, index = "Shannon"), name = "Sample", value = "Shannon"), by = "Sample") %>%
    dplyr::left_join(tibble::enframe(mgx.tools::estimate_diversity(ps, index = "Simpson"), name = "Sample", value = "Simpson"), by = "Sample") %>%
    dplyr::left_join(tibble::enframe(mgx.tools::estimate_diversity(ps, index = "Tail"), name = "Sample", value = "Tail"), by = "Sample")
  
  return(result)
}