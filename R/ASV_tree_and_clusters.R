#' Build phylogenetic tree using NJ 
#' @param physeq a phyloseq object
#' @param ncores number of cores to use
#' @export

ASV_tree_for_physeq <- function(
    physeq, 
    ncores = parallel::detectCores() ) {
  
  # Get ASVs
  original_asv_names <- phyloseq::taxa_names(physeq)  # save original names
  seqs <- Biostrings::DNAStringSet(colnames(phyloseq::otu_table(physeq)))
  names(seqs) <- original_asv_names  # assign names to DNAStringSet
  
  # Align with DECIPHER
  aligned <- DECIPHER::AlignSeqs(seqs, processors = ncores)
  
  # Convert to phangorn format & build tree
  phang_align <- phangorn::as.phyDat(aligned)
  dm <- phangorn::dist.ml(phang_align)  # pairwise distances
  tree <- phangorn::NJ(dm)  # neighbor-joining 
  
  # Add tree to phyloseq object
  phyloseq::phy_tree(physeq) <- ape::as.phylo(tree)
  
  return(physeq)
}

#' Create ASV clusters based on phylogenetic distance
#' 
#' @description
#' Adds a column to taxonomy, bridging the Genus-ASV gap for amplicons where the taxonomy is usually unresolved at Species level
#' Allows coarser-grain taxonomic analyses without losing unidentified ASVs
#' based on https://github.com/benjjneb/dada2/issues/947
#' @param physeq a phyloseq object with tree data
#' @param threshold the phylogenetic distance threshold passed to 
#' @param export_dir directory where to export dendogram as well as taxonomy-cluster consistency plot
#' @export

cluster_ASVs_physeq <- function(
    physeq,
    hclust.method = 'complete',
    threshold = 0.03){
  
  if (is.null(phyloseq::phy_tree(physeq, errorIfNULL = FALSE))) {
    stop("Provided phyloseq object doesn't have a tree! Aborting.")
  }  

  similarity=round(100*(1-threshold),0)
  
  # Extract tree
  tree <- phyloseq::phy_tree(physeq)
  
  # compute cophenetic distance for hierarchical clustering
  coph_dist <- stats::cophenetic(tree)
  
  # Hierarchical clustering:
  hc <- stats::hclust(stats::as.dist(coph_dist), method = hclust.method)
  clusters <- cutree(hc, h = threshold)
  
  # HC plot 
  # This should be a dendrogram with taxonomy superposed
  
  # hc.plot <- .dendro_plot(hc, threshold)
  # ggplot2::ggsave(
  #   file.path(export_dir, paste0("hc_",similarity,".pdf")), 
  #   bg = 'white', width = 2400, height = 1600, units = 'px', dpi = 200
  # )
  
  # build dataframe and Map back to phyloseq taxonomy table
  cluster_df <- data.frame(
    ASV = names(clusters),
    Species_cluster = as.factor(clusters)
  ) 
  
  # Merge into tax table
  tax <- data.frame(phyloseq::tax_table(physeq), check.names = FALSE) %>%
    tibble::rownames_to_column("ASV") %>%
    dplyr::left_join(cluster_df, by = "ASV") %>%
    tibble::column_to_rownames("ASV")
  
  .check_genus_consistency(tax)
  
  phyloseq::tax_table(physeq) <- as.matrix(tax)
  
  return(physeq)
}

#' Plots the tree
#' @param hc the output of hclust
#' @param threshold the phylogenetic distance threshold passed to 
#' @keywords internal
.dendro_plot <- function(hc, threshold) {
  
  dendro_data <- ggdendro::dendro_data(hc, type = "rectangle")
  
  p <- ggplot2::ggplot(dendro_data$segments) +
    ggplot2::geom_segment(
      ggplot2::aes(x = x, y = y, xend = xend, yend = yend)) +
    ggplot2::geom_vline(
      xintercept = NA, 
      linetype = "dashed", 
      color = "red", 
      linewidth = 0.5) +  # add threshold line
    ggplot2::geom_text(
      data = dendro_data$labels, 
      ggplot2::aes(x = x, y = -0.01, label = label), 
      angle = 90, hjust = 1, vjust = 0.5, size = 2
    ) +
    ggplot2::geom_hline(
      yintercept = threshold, 
      linetype = "dashed", 
      color = "red", 
      linewidth = 0.5) +
    ggplot2::labs(x = "ASV", y = "Branch length (distance)") +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank())
  
  return(p)
}


#' Evaluates the consistency between ASV clusters and Genus-level classification.
#' 
#' @description
#' Ideally, each cluster would have the same Genus-level taxonomy.
#' A cluster spanning multple genera suggests the clustering should be done at shorter phylogenetic distances (threshold)
#' @param tax_table a taxonomy data.frame with ASVs as rownames and up to Genus taxonomy column, including a Species_cluster column
#' @keywords internal

.check_genus_consistency <- function(tax_table) {
  
  tax_table <- tax_table %>% 
    dplyr::mutate(Genus = dplyr::case_when(
      is.na(Genus) ~ "Unclassified", TRUE~Genus
    )) 
  
  # Unique Genera by cluster
  clusters <- tax_table %>% 
    dplyr::group_by(Species_cluster) %>% 
    dplyr::distinct(Genus) %>% 
    dplyr::arrange(Species_cluster)
  
  # Cluster with more than one defined genus
  inconsistent_clusters <- clusters %>% 
    dplyr::filter(Genus != "Unclassified") %>% 
    dplyr::group_by(Species_cluster) %>% 
    dplyr::summarise(n=dplyr::n()) %>% 
    dplyr::filter(n>1) %>% 
    dplyr::pull(Species_cluster) %>% 
    unique()
  
  
  message(paste(length(unique(clusters$Species_cluster)), 
                'clusters created out of',
                nrow(tax_table), 'ASVs.'))
  
  if (rlang::is_empty(inconsistent_clusters)) {
    
    message('All clusters are consistent with identified genera.')
    
  } else {
    
    # print full taxonomy of inconstintent clusters
    inconst_clust_taxonomy <- clusters %>% 
      dplyr::filter(Species_cluster %in% inconsistent_clusters) %>% 
      dplyr::left_join(tax_table, by = join_by(Species_cluster, Genus)) %>% 
      dplyr::group_by(Species_cluster, Class, Order, Family, Genus) %>% 
      dplyr::summarise(n_ASVs = dplyr::n(), .groups = 'drop') 
    
    message(paste(length(inconsistent_clusters),'clusters span multiple genera. Excerpt:'))
    
    print(kableExtra::kable(head(inconst_clust_taxonomy, n = 20)))
  }
}






