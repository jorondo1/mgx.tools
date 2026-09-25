#' Compute classification rates 
#' 
#' @description
#' Computes classification rates of amplicon taxonomy table using a
#' phyloseq object (or a list of them) as input. Produces 
#' @param ps a phyloseq object or a named list of phyloseq objects
#' @param Dataset String; creates a Dataset column in the output with that string as value. Managed automatically by compute_classification_rates() when named list of phylosq objects are provided.
#' @param ranks vector of taxonomic ranks to compute rate by

.compute_classrate <- function(ps, Dataset = NA, ranks) {
  
  taxranks <- phyloseq::tax_table(ps) %>% colnames()
  miss_ranks <- setdiff(ranks, taxranks)
  
  if (length(miss_ranks)>0) {
    stop(paste('Specified rank(s) missing from input data:', '\n', miss_ranks, collapse = '\n'))
  }
  
  ps.melted <- psflashmelt(ps) %>%
    dplyr::filter(Abundance > 0) %>%
    dplyr::group_by(Sample) %>%
    dplyr::mutate(relAb = Abundance / sum(Abundance))
  
  map(ranks, function(rank) {
    ps.melted %>%
      dplyr::select(Sample, !!rlang::sym(rank), relAb) %>%
      dplyr::mutate(classified = case_when(
        !!rlang::sym(rank) == 'Unclassified' | is.na(!!rlang::sym(rank)) |
          grepl(
            "Incertae|Unclassified|uncultured|unknown",
            !!rlang::sym(rank),
            ignore.case = TRUE
          ) ~ 0,
        TRUE ~ 1)) %>%
      dplyr::summarise(
        asv_prop   = sum(classified) / n(),
        relAb_prop = sum(classified * relAb)
      ) %>%
      tidyr::pivot_longer(
        cols = c('relAb_prop','asv_prop'),
                          names_to = 'proportion_type',
        values_to = 'Classification_rate') %>%
      dplyr::mutate(
        taxRank = factor(rank, levels = ranks),
        proportion_type = dplyr::case_when(
          proportion_type == 'asv_prop' ~ 'Proportion of ASVs',
          proportion_type == 'relAb_prop' ~ 'Proportion of ASV reads'  
        )
      )
  }) %>%
    purrr::list_rbind() %>%
    {if (is.na(Dataset)) {.} else dplyr::mutate(., Dataset = Dataset)}
}

#' Call the .compute_classrate() function over lists of ps objects if pertinent
#' @param ps a phyloseq object, or a named list of them (names become the `Dataset` column)
#' @param ranks taxonomic ranks to compute classification rates for
#' @export
compute_classification_rates <- function(
    ps,
    ranks = c('Phylum','Class','Order','Family','Genus')) {
  
  if (is.list(ps) && !inherits(ps, "phyloseq")) {
    if (is.null(names(ps)) || any(names(ps) == "")) {
      stop("`ps` must be a *named* list (names used as Datasets). ",
           "Either name the list elements or pass a single phyloseq object.")
    }
    purrr::imap(ps, \(ps_obj, Dataset) .compute_classrate(ps_obj, Dataset, ranks)) %>%
      purrr::list_rbind()
  } else {
    .compute_classrate(ps, Dataset = NA, ranks)
  }
}

#' Plot classification rates
#' @param classrates output of `compute_classification_rates()`
#' @export 

plot_class_rates <- function(classrates){
  
  p <- classrates %>% 
    ggplot2::ggplot(aes(y = Classification_rate, x = taxRank, colour = taxRank)) +
    ggplot2::geom_boxplot() +
    ggplot2::ylim(0,NA)+
    ggplot2::scale_colour_brewer(palette = 'Set2') +
    ggplot2::labs(y = 'Proportion of taxonomically labelled ASVs',
         colour = 'Taxonomic rank') +
    ggplot2::theme(
      axis.ticks.x = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_blank(),
      legend.position = 'none') 
  
  if ("Dataset" %in% names(classrates)) {
    p <- p + ggplot2::facet_grid(Dataset ~ proportion_type)
  }else{
    p <- p + ggplot2::facet_grid(. ~ proportion_type)
  }
  return(p)

}
