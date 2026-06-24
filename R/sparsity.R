#' Compute sparsity from a matrix where samples are column
#' @param input matrix, tibble, data.frame or phyloseq object
#' @export
#' 
sparsity <- function(input){
  
  if( is.matrix(input) ) {
    seqtab <- as.data.frame(input)
  }
  
  else if( tibble::is_tibble(input) | is.data.frame(input) ) {
    seqtab <- input %>% 
      dplyr::select(where(is.double))
  } 
  
  else if ( class(input) == "phyloseq" ) {
    seqtab <- phyloseq::otu_table(input) %>% 
      as("matrix")
    if( ! phyloseq::taxa_are_rows(input) ) {
      seqtab <- t(seqtab)
    }
  }
  
  else stop('unsupported input.')
  
  message('Proportion of 0 in matrix:')
  sum(seqtab==0)/(nrow(seqtab)*ncol(seqtab))
}
