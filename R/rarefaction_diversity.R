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
  
  ## ---- 1. Extract count table, set rarefaction depth ---------------------
  
  count_table <- as(phyloseq::otu_table(ps), 'matrix')
  if (phyloseq::taxa_are_rows(ps)) count_table <- t(count_table)
  
  depth <- depth %||% min(rowSums(count_table))
  
  ## ---- 2. Drop samples below depth & now-empty taxa ------------------------
  # Filtering here (once) keeps ps, count_table, and n_samples consistent for
  # the rest of the function.
  
  seed_start   <- sample.int(1e6, 1)
  keep_samples <- rowSums(count_table) >= depth
  keep_taxa    <- colSums(count_table[keep_samples, , drop = FALSE]) > 0
  
  ps <- phyloseq::prune_samples(rownames(count_table)[keep_samples], ps)
  ps <- phyloseq::prune_taxa(colnames(count_table)[keep_taxa], ps)
  
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
    keep_samples <- rowSums(taxon_rare) >= depth
    keep_taxa    <- colSums(taxon_rare[keep_samples, , drop = FALSE]) > 0
    
    ps_rare <- phyloseq::prune_samples(rownames(taxon_rare)[keep_samples], ps)
    ps_rare <- phyloseq::prune_taxa(colnames(taxon_rare)[keep_taxa], ps_rare)
    
    taxon_rare <- taxon_rare[keep_samples, keep_taxa]
    
    n_taxa          <- ncol(taxon_rare)
    n_samples_iter  <- nrow(taxon_rare)
    tail_multiplier <- (seq_len(n_taxa) - 1)^2
    
    ## -- Alpha diversity (row-wise, pre-allocated) --
    
    richness  <- integer(n_samples_iter)
    shannon   <- numeric(n_samples_iter)
    simpson   <- numeric(n_samples_iter)
    tail_vals <- numeric(n_samples_iter)
    
    for (j in 1:n_samples_iter) {
      row <- taxon_rare[j, ]
      non_zero_idx  <- row > 0
      non_zero_vals <- row[non_zero_idx]
      
      if (length(non_zero_vals) > 0) {
        row_sum <- sum(non_zero_vals)
        richness[j] <- length(non_zero_vals)
        
        # Faster Shannon using precomputed log
        p <- non_zero_vals / row_sum
        shannon[j] <- -sum(p * log(p))
        simpson[j] <- sum(p * p)
      }
      
      # Tail statistic - sort entire row
      # Could be optimized with partial sort, but likely not worth it
      tail_vals[j] <- sqrt(sum(sort(row, decreasing = TRUE) * tail_multiplier))
    }
    
    ## -- Rebuild rarefied ps object for Faith PD and UniFrac --
    
    phyloseq::otu_table(ps_rare) <- phyloseq::otu_table(
      taxon_rare, taxa_are_rows = FALSE
    )
    
    # Faith PD: 
    faith_res <- suppressWarnings(suppressMessages(btools::estimate_pd(ps_rare)))
    
    ## -- Beta diversity --
    
    counts_vst <-
      if (vst) {
        suppressMessages(mgx.tools:::vst_ps_to_mx(ps_rare))
      } else { 
        phyloseq::otu_table(ps_rare) 
      }
    
    list(
      alpha = list(
        richness = richness, 
        shannon  = shannon, 
        simpson  = simpson, 
        tail     = tail_vals,
        faith    = faith_res$PD),
      
      beta = list(
        unifrac.w   = suppressWarnings(phyloseq::UniFrac(ps_rare, weighted = TRUE)),
        unifrac.u   = suppressWarnings(phyloseq::UniFrac(ps_rare, weighted = FALSE)),
        bray        = vegan::vegdist(counts_vst, method = 'bray'),
        # note: r.aitchison uses the non-vst transformed, but rarefied table
        # motivated by https://pubmed.ncbi.nlm.nih.gov/38251877/
        r.aitchison = vegan::vegdist(taxon_rare, method = 'robust.aitchison')
      )
    )
  }
  
  ## ---- 4. Run all iterations (parallel) ------------------------------------
  
  results <- parallel::mcmapply(
    rarefy_iter, 1:n_iter, SIMPLIFY = FALSE, mc.cores = mc.cores
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
  
  sd <- phyloseq::sample_data(ps)
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
    do.call(analogue::fuse, dist_list)
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