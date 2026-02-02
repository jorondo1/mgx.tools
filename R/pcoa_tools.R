#' Distance-based analyses

#' Variance-stabilizing transformation
#' @param ps phyloseq object
vst_ps_to_mx <- function(ps) {
  phyloseq_to_deseq2(
    ps, ~ 1) %>% # DESeq2 object
    DESeq2::estimateSizeFactors(., geoMeans = apply(
      DESeq2::counts(.), 1, function(x) exp(sum(log(x[x>0]))/length(x)))) %>%
    DESeq2::varianceStabilizingTransformation(blind=T) %>% # VST
    SummarizedExperiment::assay(.) %>% t %>%
    { .[. < 0] <- 0; . } # replace negatives by zeros
}

#' PCOA
#' Return a list of 3 elements :
#' 1. the phyloseq sample_data with 2 first PCo added
#' 2. the eigenvalues
#' 3. the distance/dissimilarity matrix
#' @param ps phyloseq object
#' @param dist distance name (choose from vegan)
#' @param vst default FALSE, whether to use variance-stabilizing transformation
#' @param all_coordinates default FALSE; whether to include a table with all coordinates instead of just adding the first 2 to the sample_data table
#' @export
compute_pcoa <- function(ps, dist,
                         all_coordinates = FALSE, # adds a table with all coordinates
                         vst = FALSE # add variance-stabilizing transformation
) {

  seqtab <- phyloseq::otu_table(ps)

  # Validate distance
  unifrac_names <- c("unifrac.u","unifrac.w")
  dist_list <- c(unifrac_names,"manhattan", "euclidean", "canberra", "clark", "bray", "kulczynski", "jaccard", "gower", "altGower", "morisita", "horn", "mountford", "raup", "binomial", "chao", "cao", "mahalanobis", "chisq", "chord", "hellinger", "aitchison", "robust.aitchison")
  if (!dist %in% dist_list) {
    stop(paste(c("Distance must be one of the following:", dist_list), collapse = ", "))
  }

  # Validate tree if distance is UniFrac
  if (dist %in% unifrac_names & is.null(phyloseq::phy_tree(ps, errorIfNULL = FALSE))) {
    stop(paste("The provided phyloseq object does not contain a tree.", dist, "requires a reference tree."))
  }

  dist.mx <- if (dist == 'unifrac.w') {
    phyloseq::UniFrac(ps, weighted = TRUE, parallel = TRUE)

  } else if (dist == 'unifrac.u') {
    phyloseq::UniFrac(ps, weighted = FALSE, parallel = TRUE)

  } else {
    ps %>%
      { counts <- if (vst) vst_ps_to_mx(.) else seqtab
      if (phyloseq::taxa_are_rows(ps) & !vst) t(counts) else counts
      } %>%
      vegan::vegdist(method = dist)
  }

  num_k <- if(taxa_are_rows(seqtab)){
    ncol(seqtab) - 1
  } else {nrow(seqtab) - 1}

  PCoA <- cmdscale(dist.mx, k = num_k, eig = TRUE)
  eig <- round(PCoA$eig[1:3]/sum(PCoA$eig),2)
  message(paste("Variance explained by first PCo's:",eig[1], ',', eig[2], ',', eig[3]))
  # create output list
  out <- data.frame(sample_data(ps))
  out$PCo1 <- PCoA$points[,1]
  out$PCo2 <- PCoA$points[,2]

  out <- list(metadata = out, eig = PCoA$eig, dist.mx = dist.mx)

  if(all_coordinates){
    out[['coordinates']] <- PCoA$points
  }
  return(out)
}

#' Pivot a dist object to a long dataframe
#' Possible to use only a subset of samples, provided their name
#' @param dist.mx distance matrix (dist object from vegan::)
#' @param sample_subset a vector of sample names to subset the distance matrix
#' @export
compile_dist_pairs <- function(dist.mx, sample_subset = NULL) {

  dist_matrix <- as.matrix(dist.mx)

  if (!is.null(sample_subset)) {
    dist_matrix <- dist_matrix[sample_subset,sample_subset]
  }

  upper_indices <- which(upper.tri(dist_matrix), arr.ind = TRUE)

  data.frame(
    Sample1 = rownames(dist_matrix)[upper_indices[, 1]],
    Sample2 = colnames(dist_matrix)[upper_indices[, 2]],
    Distance = dist_matrix[upper_indices]
  )
}
