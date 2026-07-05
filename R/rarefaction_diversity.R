#' A function to repeat rarefy and return an updated phyloseq object
#' Takes in a ps object and outputs sample_data as a tibble with one column per index
#' If phy_tree, computes Faith and both Unifrac metrics
#' 
#' Note: mc.cores > 1 relies on parallel::mcmapply (fork-based), which is
#' not available on Windows -- falls back to mc.cores = 1 there.
#' @param ps a phyloseq object
#' @param depth sequencing depth to rarefy to (default: min sample size)
#' @param n_iter number of rarefaction iterations
#' @param mc.cores cores to use
#' @param vst use variance-stabilizing transformation for bray-curtis and unifrac distances
#' @export
#' 
rarefy_diversity <- function(
    ps,
    depth = NULL, 
    n_iter = 100, 
    mc.cores = parallel::detectCores(),
    vst = FALSE) {
  
  # Root tree if necessary
  
  if (!ape::is.rooted(phyloseq::phy_tree(ps))) {
    phyloseq::phy_tree(ps) <- phangorn::midpoint(phyloseq::phy_tree(ps))
  }
  
  # extract counts, transpose if needed
  
  count_table <- as(phyloseq::otu_table(ps), 'matrix')
  if (phyloseq::taxa_are_rows(ps)) count_table <- t(count_table)
  
  depth <- depth %||% min(rowSums(count_table))
  
  ## ---- Drop samples below depth & now-empty taxa ------------------------
  # Filtering here (once) keeps ps, count_table, and n_samples consistent for
  # the rest of the function.
  
  seed_start   <- sample.int(1e6, 1)
  keep_samples <- rowSums(count_table) >= depth
  keep_taxa    <- colSums(count_table[keep_samples, , drop = FALSE]) > 0
  
  if (all(keep_samples)) {
    ps_depth <- ps
  } else {
    ps_depth <- phyloseq::prune_samples(rownames(count_table)[keep_samples], ps)
  }
  if (!all(keep_taxa)) {
    ps_depth <- phyloseq::prune_taxa(colnames(count_table)[keep_taxa], ps_depth)
  }
  
  count_table <- count_table[keep_samples, keep_taxa]
  n_samples   <- nrow(count_table)
  
  ## ---- 3. Per-iteration function: rarefy + alpha + beta ----------------------
  # Defined inside rarefy_diversity so it closes over ps, count_table, depth,
  # seed_start, n_samples, and vst.
  
  rarefy_iter <- function(i) {
    set.seed(seed_start + i - 1)
    
    # -- Rarefy:
    taxon_rare <- vegan::rrarefy(count_table, sample = depth)
    
    # -- re-filter
    # a taxon can rarefy down to zero reads across all samples in a given draw.
    keep_taxa    <- colSums(taxon_rare) > 0
    
    if (all(keep_taxa)) {
      ps_rare <- ps
    } else {
      ps_rare <- phyloseq::prune_taxa(colnames(taxon_rare)[keep_taxa], ps_depth)
      taxon_rare <- taxon_rare[, keep_taxa, drop = FALSE]
    }
    
    n_taxa          <- ncol(taxon_rare)
    n_samples_iter  <- nrow(taxon_rare)
    tail_multiplier <- (seq_len(n_taxa) - 1)^2
    
    ## -- Alpha diversity --
    row_sums <- rowSums(taxon_rare)
    richness <- rowSums(taxon_rare > 0)
    # Shannon
    p        <- taxon_rare / row_sums
    logp     <- ifelse(p > 0, log(p), 0)
    shannon  <- -rowSums(p * logp)
    simpson  <- rowSums(p^2)
    
    # Tail : 
    sorted_mat <- apply(taxon_rare, 1, sort, decreasing = TRUE)  # n_taxa x n_samples
    tail_vals  <- sqrt(colSums(sorted_mat * tail_multiplier))
    
    ## -- Rebuild rarefied ps object for Faith PD and UniFrac --
    
    phyloseq::otu_table(ps_rare) <- phyloseq::otu_table(
      taxon_rare, taxa_are_rows = FALSE
    )
    
    # Faith PD: 
    faith_res <- suppressWarnings(suppressMessages(btools::estimate_pd(ps_rare)))
    
    ## -- Beta diversity --
    
    counts_vst <-
      if (vst) {
        suppressMessages(mgx.tools::vst_ps_to_mx(ps_rare))
      } else { 
        phyloseq::otu_table(ps_rare) 
      }
    
    # Unifrac 
    gu <- suppressWarnings(
      GUniFrac::GUniFrac(
        taxon_rare, phyloseq::phy_tree(ps_rare), alpha = c(1), verbose = FALSE))

    # Return: 
    
    list(
      alpha = list(
        richness = richness, 
        shannon  = shannon, 
        simpson  = simpson, 
        tail     = tail_vals,
        faith    = faith_res$PD),
      
      beta = list(
        unifrac.w   = as.dist(gu$unifracs[, , "d_1"]),
        unifrac.u   = as.dist(gu$unifracs[, , "d_UW"]) ,
        bray        = vegan::vegdist(counts_vst, method = 'bray'),
        # note: r.aitchison uses the non-vst transformed, but rarefied table
        # motivated by https://pubmed.ncbi.nlm.nih.gov/38251877/
        r.aitchison = vegan::vegdist(taxon_rare, method = 'robust.aitchison')
      )
    )
  }
  
  ## ---- 4. Run all iterations (parallel) ------------------------------------
  
  # Exception if only one core
  plan_type <- if (mc.cores > 1) future::multisession else future::sequential
  future::plan(plan_type, workers = mc.cores)
  on.exit(future::plan(future::sequential), add = TRUE)
  
  # Loop 
  results <- furrr::future_map(
    1:n_iter, rarefy_iter,
    .options = furrr::furrr_options(seed = TRUE)
  )
  ## ---- 5. Aggregate alpha diversity across iterations -----------------------
  
  extract_alpha_metric <- function(metric) {
    matrix(
      unlist(lapply(results, function(x) x$alpha[[metric]])),
      nrow = n_samples, ncol = n_iter, byrow = FALSE
    )
  }
  
  richness_mat <- extract_alpha_metric("richness")
  shannon_mat  <- extract_alpha_metric("shannon")
  simpson_mat  <- extract_alpha_metric("simpson")
  tail_mat     <- extract_alpha_metric("tail")
  faith_mat    <- extract_alpha_metric("faith")
  
  row_means <- function(x) matrixStats::rowMeans2(x)
  
  alpha_average <- data.frame(
    Richness = row_means(richness_mat),
    Shannon  = row_means(shannon_mat),
    Hill_1   = exp(row_means(shannon_mat)),
    Simpson  = row_means(simpson_mat),
    Hill_2   = row_means(1 / simpson_mat),
    Tail     = row_means(tail_mat),
    Faith    = row_means(faith_mat),
    row.names = rownames(count_table)
  )
  
  sd <- phyloseq::sample_data(ps_depth)
  for (metric in names(alpha_average)) {
    sd[[metric]] <- alpha_average[[metric]]
  }
  
  alpha_tibble <- sd %>% 
    data.frame() %>% 
    tibble::rownames_to_column('Sample') %>% 
    tibble::tibble()
  
  ## ---- 6. Aggregate beta diversity across iterations (fuse) ----------------
  
  extract_beta_metric <- function(metric) {
    dist_list <- lapply(results, function(x) x$beta[[metric]])
    fused_dist <- do.call(analogue::fuse, dist_list)
    
    # Clean call and wheights because heavy
    attr(fused_dist, "call") <- NULL
    attr(fused_dist, "weights") <- NULL

    return(fused_dist)
  }  
  ## ---- 7. Return -------------------------------------------------------
  
  list(
    alpha = alpha_tibble,
    beta = list(
      bray        = extract_beta_metric("bray"),
      unifrac_u   = extract_beta_metric("unifrac.u"),
      unifrac_w   = extract_beta_metric("unifrac.w"),
      r_aitchison = extract_beta_metric("r.aitchison")
    )
  )
}