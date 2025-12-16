#' Paralellize phyloseq's rarefaction function
#' Shamelessly stolen and fed to Gemini
#' @import foreach
#' @importFrom doParallel registerDoParallel
#' @export
rarefy_even_depth2 <- function (
    physeq, sample.size = min(phyloseq::sample_sums(physeq)), rngseed = FALSE,
    replace = TRUE, trimOTUs = TRUE, verbose = FALSE, ncores = parallel::detectCores()-1
)
{
  if (!requireNamespace("doParallel", quietly = TRUE)) {
    stop("The package 'doParallel' is required but not installed.")
  }

  if (as(rngseed, "logical")) {
    set.seed(rngseed)
    if (verbose) {
      message("`set.seed(", rngseed, ")` was used to initialize repeatable random subsampling.")
      message("Please record this for your records so others can reproduce.")
      message("Try `set.seed(", rngseed, "); .Random.seed` for the full vector",
              sep = "")
      message("...")
    }
  } else if (verbose) {
    message("You set `rngseed` to FALSE. Make sure you've set & recorded\n",
            " the random seed of your session for reproducibility.\n",
            "See `?set.seed`\n")
    message("...")
  }

  if (length(sample.size) > 1) {
    warning("`sample.size` had more than one value. ", "Using only the first. \n ... \n")
    sample.size <- sample.size[1]
  }

  if (sample.size <= 0) {
    stop("sample.size less than or equal to zero. ", "Need positive sample size to work.")
  }

  if (min(sample_sums(physeq)) < sample.size) {
    rmsamples = phyloseq::sample_names(physeq)[
      phyloseq::sample_sums(physeq) < sample.size]
    if (verbose) {
      message(length(rmsamples), " samples removed", "because they contained fewer reads than `sample.size`.")
      message("Up to first five removed samples are: \n")
      message(paste(rmsamples[1:min(5, length(rmsamples))],
                    sep = "\t"))
      message("...")
    }
    physeq = phyloseq::prune_samples(
      setdiff(phyloseq::sample_names(physeq),
              rmsamples),
      physeq)
  }

  newsub <- physeq

  if (!taxa_are_rows(newsub)) {
    newsub <- phyloseq::t(newsub)
  }

  # Parallelization
  applyByCol <- function(df, FUN, ...) {
    res = foreach::`%dopar%`( # this is how we call package operators
      foreach::foreach(i=1:ncol(df)),
      FUN(as.vector(df[,i]), ...)
    )
    res = as.data.frame(res)
    colnames(res) = colnames(df)
    rownames(res) = rownames(df)
    return(res)
  }
  # apply through each sample, and replace
  threads = ncores
  doParallel::registerDoParallel(cores=threads)

  newotu <- applyByCol(
    phyloseq::otu_table(newsub),
    phyloseq:::rarefaction_subsample,
    sample.size=sample.size,
    replace=replace)
  # newotu <- apply(otu_table(newsub), 2, phyloseq:::rarefaction_subsample,
  #                 sample.size = sample.size, replace = replace)

  rownames(newotu) <- phyloseq::taxa_names(physeq)

  phyloseq::otu_table(newsub) <- phyloseq::otu_table(newotu, TRUE)

  if (trimOTUs) {
    rmtaxa = taxa_names(newsub)[taxa_sums(newsub) <= 0]
    if (length(rmtaxa) > 0) {
      if (verbose) {
        message(length(rmtaxa), "OTUs were removed because they are no longer \n",
                "present in any sample after random subsampling\n")
        message("...")
      }
      newsub = phyloseq::prune_taxa(
        setdiff(phyloseq::taxa_names(newsub),
                rmtaxa),
        newsub
      )
    }
  }
  if (!taxa_are_rows(physeq)) {
    newsub <- phyloseq::t(newsub)
  }
  return(newsub)
}
