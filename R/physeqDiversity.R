#' Functions to compute various diversity indices from phyloseq object
#'

#' Hill numbers
#' @param ps phyloseq object
#' @param q Hill number order q (i-e 0, 1 or 2)
#' @export
estimate_Hill <- function(ps, q) {
  x <- phyloseq::otu_table(ps) %>% as("matrix")
  if (phyloseq::taxa_are_rows(ps)) {
    x <- t(x)
  }
  total <- rowSums(x)
  x <- sweep(x, 1, total, "/")

  if (q == 0) {  # Species richness
    div <- rowSums(x > 0)
  } else if (q == 1) { # Shannon diversity (exponential of Shannon entropy)
    div <- exp(-rowSums(x * log(x, base = exp(1)), na.rm = TRUE))
  } else {  # Hill number formula for q ≠ 0 and q ≠ 1
    div <- rowSums(x^q)^(1 / (1 - q))
  }
  return(div)
}

#' Esitmate diversity (Shannon, Simpson, Tail)
#' @param ps phyloseq object
#' @param index Richness, Shannon, Simpson or Tail
#' @export
estimate_diversity <- function(ps, index = 'Shannon') {
  x <- phyloseq::otu_table(ps) %>% as("matrix")
  if (phyloseq::taxa_are_rows(ps)) {
    x <- t(x)
  }
  total <- apply(x, 1, sum)
  x <- sweep(x, 1, total, "/")

  if(index == 'Tail') {
    tail_stat <- function(row) {
      values <- sort(row, decreasing = TRUE)
      sqrt(sum(values * ((seq_along(values)-1)^2)))
    }
    div <- apply(x, 1, tail_stat)
  }
  if(index == 'Shannon') {
    x <- -x * log(x, exp(1))
    div <- apply(x, 1, sum, na.rm = TRUE)
  }
  if(index == 'Simpson') {
    div <- apply((x * x), 1, sum, na.rm = TRUE)
  }
  if(index == 'Richness') {
    div <- apply(x, 1, function(x) sum(x != 0))
  }
  return(div)
}

#' Compute Hill 0, 1 and 2, as well as Tail
#' Returns a list
#' @param ps phyloseq object
#' @param idx a vector of index names
#' @export
div.fun <- function(ps, idx) {
  message('This function is obsolete, use diversity_tibble() instead.')
  div_estimate <- list()
  for (i in seq_along(idx)) { # compute Hill numbers
    H_q=paste0("H_",i-1) # format H_0, H_1...
    div_estimate[[H_q]] <- estimate_Hill(ps, idx[i])
  }
  div_estimate[["Tail"]] <- estimate_diversity(ps, index = "Tail")
  return(div_estimate)
}


