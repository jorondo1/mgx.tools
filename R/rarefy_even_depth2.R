#' Parallelize phyloseq's rarefaction function with reproducible results
#' @import foreach
#' @importFrom doParallel registerDoParallel
#' @importFrom doRNG %dorng%
#' @export
rarefy_even_depth2 <- function (
    physeq, sample.size = min(phyloseq::sample_sums(physeq)), rngseed = FALSE,
    replace = TRUE, trimOTUs = TRUE, verbose = TRUE, ncores = parallel::detectCores()-1
)
{
  # Check for required packages for parallel processing
  if (!requireNamespace("doParallel", quietly = TRUE)) {
    stop("The package 'doParallel' is required for this function.")
  }
  if (!requireNamespace("doRNG", quietly = TRUE)) {
    stop("The package 'doRNG' is required for reproducible parallel randomization.")
  }

  # --- Seed Management ---
  # This is the crucial part. We set the seed here, ONCE.
  # The `%dorng%` operator will handle distributing this correctly to the workers.
  if (is.numeric(rngseed)) {
    set.seed(rngseed)
    if (verbose) {
      message("`set.seed(", rngseed, ")` was used to initialize repeatable random subsampling.")
    }
  } else if (verbose) {
    message("You set `rngseed` to FALSE or a non-numeric value. Subsampling will not be reproducible.")
  }

  # --- Input Validation and Sample Pruning (same as original) ---
  if (length(sample.size) > 1) {
    warning("`sample.size` had more than one value. Using only the first.")
    sample.size <- sample.size[1]
  }
  if (sample.size <= 0) {
    stop("sample.size must be positive.")
  }

  if (min(sample_sums(physeq)) < sample.size) {
    rmsamples <- phyloseq::sample_names(physeq)[phyloseq::sample_sums(physeq) < sample.size]
    if (verbose) {
      message(length(rmsamples), " samples removed because they contained fewer reads than `sample.size`.")
    }
    physeq <- phyloseq::prune_samples(setdiff(phyloseq::sample_names(physeq), rmsamples), physeq)
  }

  # --- Parallel Rarefaction ---
  newsub <- physeq
  # Ensure taxa are rows for easy column-wise operation
  if (!taxa_are_rows(newsub)) {
    newsub <- phyloseq::t(newsub)
  }

  # Set up the parallel backend
  doParallel::registerDoParallel(cores = ncores)

  # Extract the OTU table for processing
  otu_tab <- as(phyloseq::otu_table(newsub), "matrix")

  # The parallel loop using `%dorng%` for reproducibility
  rarefied_list <- foreach::foreach(i = 1:ncol(otu_tab)) %dorng% {
    # Each iteration of this loop is a task sent to a worker.
    # `%dorng%` ensures the RNG state is handled correctly.
    phyloseq:::rarefaction_subsample(otu_tab[, i], sample.size, replace)
  }

  # Re-assemble the results from the list into a new OTU table matrix
  newotu <- do.call(cbind, rarefied_list)
  rownames(newotu) <- rownames(otu_tab)
  colnames(newotu) <- colnames(otu_tab)

  # Update the phyloseq object
  phyloseq::otu_table(newsub) <- phyloseq::otu_table(newotu, taxa_are_rows = TRUE)

  # --- OTU Trimming and Finalization (same as original) ---
  if (trimOTUs) {
    rmtaxa <- taxa_names(newsub)[taxa_sums(newsub) <= 0]
    if (length(rmtaxa) > 0) {
      if (verbose) {
        message(length(rmtaxa), " OTUs were removed because they are no longer ",
                "present in any sample after random subsampling.")
      }
      newsub <- phyloseq::prune_taxa(setdiff(taxa_names(newsub), rmtaxa), newsub)
    }
  }

  # Return to original orientation if necessary
  if (!taxa_are_rows(physeq)) {
    newsub <- phyloseq::t(newsub)
  }

  return(newsub)
}
