#' Assemble phyloseq object 
#' 
#' @description
#' from a table that contains both taxonomy and sample identifiers in the columns
#' Tailored to work with the tibble-format output of parse_*() functions (mgx.tools package under development)

#' @param abundTable wide tibble with taxonomic ranks and sample IDs as column; one species (or ASV) per row
#' @param sampleData data frame with sample data and matching sample IDs as row names
#' @param filtering default FALSE; deprecated
#' @param min_sample_size default 100; samples with fewer sequences are discarded
#' @param unique_taxa_ID default "Species"; Column name containing unique taxa IDs. For shotgun (Metaphlan, Kraken, etc.), usually "Species".
#' @param min_taxa_count default 1; taxa with fewer sequences overall (sum of all samples) are removed
#' @param justBacteria default FALSE; filter to keep only bacteria counts
#' @export
assemble_phyloseq <- function(
    abunTable, sampleData,
    min_sample_size = 100,
    min_taxa_count = 1,
    unique_taxa_ID = 'Species',
    filtering = FALSE,
    justBacteria = FALSE,
    onlySpecies = FALSE,
    seqdepth_plot = TRUE
) {

  # Check formats
  if(tibble::is_tibble(sampleData)) {
    stop('sampleData must be data frame with sample IDs as rownames.')
  }

  # Cleanup the taxonomy
  abunTable %<>%
    {if(justBacteria) (.) %>%
        dplyr::filter(Kingdom %in% c("Bacteria","Archaea") | is.na(Kingdom))
      else .} %>%
    mutate(across(where(is.character), \(x) {
      stringr::str_replace_all(x,'_', ' ') %>%
        stringr::str_replace('Candidatus ', '') %>%
        stringr::str_remove(" [A-Z]$")  #https://gtdb.ecogenomic.org/faq#why-do-some-family-and-higher-rank-names-end-with-an-alphabetic-suffix
    }))

  # Extract abundance table with Species as identifier
  abund <- abunTable %>%
    dplyr::select(where(is.double), all_of(unique_taxa_ID)) %>%
    dplyr::group_by(across(all_of(unique_taxa_ID))) %>%
    dplyr::summarise(across(where(is.numeric), sum)) %>%
    tibble::column_to_rownames(unique_taxa_ID) %>%
    dplyr::select(where(~ sum(.) >= min_sample_size)) %>%
    dplyr::filter(rowSums(dplyr::select(., where(is.numeric))) >= min_taxa_count)

  # Visualise sequence depth per sample
  if(seqdepth_plot){
    abund_t <- as.data.frame(abund) %>% t()
    viz_seqdepth(abund_t, log_count = FALSE)
  }

  # Extract taxonomy
  tax <- abunTable %>%
    { if (onlySpecies) {
      dplyr::select(., !!sym(unique_taxa_ID))
    } else {
      dplyr::select(., where(is.character))
    }} %>%
    unique() %>% # because of renaming above, some species will be duplicate
    dplyr::mutate(tmp = !!sym(unique_taxa_ID)) %>%
    tibble::column_to_rownames('tmp') %>%
    as.matrix()

  # Some datasets may end up with very low read counts and lose samples.
  # We subset the sample dataset, but we add a check if all samples are lost:
  keep_samples <- which(rownames(sampleData) %in% colnames(abund))

  if (length(keep_samples)==0) {
    return(NULL)
    message('no samples left')
  } else {
    sampleData_subset <- sampleData[keep_samples,, drop = FALSE]

    # Build phyloseq
    ps <- phyloseq::phyloseq(
      phyloseq::otu_table(abund, taxa_are_rows = TRUE),
      phyloseq::sample_data(sampleData_subset),
      phyloseq::tax_table(tax)
    )
    #remove any empty samples
    phyloseq::prune_samples(phyloseq::sample_sums(ps) > 0, ps) %>%
      # remove taxa absent from all (may happen if you end up using not all the samples you parse, e.g. metadata missing so sample dropped in the process)
      phyloseq::prune_taxa(phyloseq::taxa_sums(.) > 0,.) %>%
      return()
  }
}
