#' find top N taxa and generate a tibble of taxa with Others category
#' @param psmelt melted phyloseq data with a `relAb` column
#' @param taxLvl taxonomic rank to aggregate at (e.g. "Genus")
#' @param topN number of taxa to keep before lumping into "Others"
#' @export
topTaxa <- function(psmelt, taxLvl, topN) {
  psmelt |>
    dplyr::group_by(!!rlang::sym(taxLvl)) |> # group by tax level
    dplyr::filter(relAb != 'NaN') |>
    dplyr::summarise(relAb = mean(relAb)) |> # find top abundant taxa
    dplyr::arrange(desc(relAb)) |>
    dplyr::mutate(aggTaxo = as.factor(
      dplyr::case_when( # aggTaxo will become the plot legend
        row_number() <= topN ~ !!sym(taxLvl), #+++ We'll need to manually order the species!
        row_number() > topN ~ 'Others'))) # +1 to include the Others section!
}

#' Community-plot data from a melted phyloseq object
#' 
#' @description
#' long df prepared for barchart visualisation of community members
#' generates an "others" category for a specified number of taxa
#' at desired taxonomic rank
#' @param melted_ps output of `phyloseq::psmelt()` or `psflashmelt()`
#' @param nTaxa number of taxa to show before lumping into "Others"
#' @param taxRank taxonomic rank to aggregate at
#' @param grouping_vars sample data columns to keep (e.g. for faceting)
#' @param order_by legend order: "alphabetical" or "abundance"
#' @param seed random seed for colour assignment
#' @export
gen_comm_plot_data <- function(
    melted_ps,
    nTaxa = 20,
    taxRank = "Genus",
    grouping_vars,
    order_by = "alphabetical",
    seed = 1346134564) {
  
  order_by_possible <- c("alphabetical", "abundance")
  if(!(order_by %in% order_by_possible)) {
    stop("Invalid ordering method, choose from:",
         paste(order_by_possible, collapse = ", "))
  }
  melted_ps %<>%
    dplyr::group_by(Sample) |>
    dplyr::mutate(relAb = Abundance / sum(Abundance))
  
  # Compute top taxa and create "Others" category
  (top_taxa <- topTaxa(melted_ps, taxRank, nTaxa))
  
  top_taxa_lvls <- top_taxa |>
    dplyr::group_by(aggTaxo) %>%
    aggregate(relAb ~ aggTaxo, data = ., FUN = sum) %>%
    {
      if (order_by == 'alphabetical') {
        dplyr::arrange(., aggTaxo)
      } else if (order_by == 'abundance') {
        dplyr::arrange(., desc(relAb))
      }
    } %$% aggTaxo |>
    as.character() %>% # Others first:
    setdiff(., c('Others', 'Unclassified')) %>% c('Others', 'Unclassified', .)
  
  # Add colours as a column
  set.seed(seed)
  colour_matching <- data.frame(
    aggTaxo = top_taxa_lvls,
    taxColour = sample(colorRampPalette(RColorBrewer::brewer.pal(12, 'Paired'))(length(top_taxa_lvls)))
  )
  
  # Merge data, aggregate taxonomy, and colours
  melted_ps %>%
    dplyr::left_join(top_taxa %>% select(-relAb), by = taxRank) %>%
    dplyr::left_join(colour_matching, by = 'aggTaxo') %>%
    dplyr::mutate(aggTaxo = factor(aggTaxo, levels = top_taxa_lvls)) %>%
    dplyr::group_by(Sample, aggTaxo, taxColour, across(all_of(grouping_vars))) %>%
    dplyr::summarise(
      relAb = sum(relAb),
      Abundance = sum(Abundance),
      .groups = 'drop')
}

