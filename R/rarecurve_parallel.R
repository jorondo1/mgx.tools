#' Paralellized version of vegan::rarecurve()
#' Shamelessly adapted by Gemini
#' @params x count table with sample as rownames and taxa as colnames
#' @export
rarecurve_parallel <- function (x, step = 1, sample, xlab = "Sample Size", ylab = "Species",
                                label = TRUE, col, lty, tidy = FALSE, ncores = NULL, ...)
{
  x <- as.matrix(x)
  if (!identical(all.equal(x, round(x)), TRUE))
    stop("function accepts only integers (counts)")
  minobs <- min(x[x > 0])
  if (minobs > 1)
    warning(gettextf("most observed count data have counts 1, but smallest count is %d",
                     minobs))
  if (missing(col))
    col <- par("col")
  if (missing(lty))
    lty <- par("lty")
  tot <- rowSums(x)
  S <- specnumber(x)
  if (any(S <= 0)) {
    message("empty rows removed")
    x <- x[S > 0, , drop = FALSE]
    tot <- tot[S > 0]
    S <- S[S > 0]
  }
  nr <- nrow(x)
  col <- rep(col, length.out = nr)
  lty <- rep(lty, length.out = nr)

  # Parallel processing setup
  if (is.null(ncores)) {
    ncores <- parallel::detectCores() - 1  # Leave one core free
  }
  ncores <- min(ncores, nr)  # Don't use more cores than rows

  if (ncores > 1) {
    # Use parallel processing
    cl <- parallel::makeCluster(ncores)
    on.exit(parallel::stopCluster(cl))

    # Export necessary functions to clusters
    parallel::clusterExport(cl, varlist = c("rarefy"),
                            envir = environment())

    out <- parallel::parLapply(cl, seq_len(nr), function(i) {
      n <- seq(1, tot[i], by = step)
      if (n[length(n)] != tot[i]) {
        n <- c(n, tot[i], use.names = FALSE)
      }
      drop(suppressWarnings(rarefy(x[i, ], n)))
    })
  } else {
    # Use sequential processing (original code)
    out <- lapply(seq_len(nr), function(i) {
      n <- seq(1, tot[i], by = step)
      if (n[length(n)] != tot[i]) {
        n <- c(n, tot[i], use.names = FALSE)
      }
      drop(suppressWarnings(rarefy(x[i, ], n)))
    })
  }

  if (tidy) {
    len <- sapply(out, length)
    nm <- rownames(x)
    df <- data.frame(Site = factor(rep(nm, len), levels = nm),
                     Sample = unlist(lapply(out, attr, which = "Subsample")),
                     Species = unlist(out))
    return(df)
  }
  Nmax <- sapply(out, function(x) max(attr(x, "Subsample")))
  Smax <- sapply(out, max)
  plot(c(1, max(Nmax)), c(1, max(Smax)), xlab = xlab, ylab = ylab,
       type = "n", ...)
  if (!missing(sample)) {
    abline(v = sample)
    rare <- sapply(out, function(z) approx(x = attr(z, "Subsample"),
                                           y = z, xout = sample, rule = 1)$y)
    abline(h = rare, lwd = 0.5)
  }
  for (ln in seq_along(out)) {
    N <- attr(out[[ln]], "Subsample")
    lines(N, out[[ln]], col = col[ln], lty = lty[ln], ...)
  }
  if (label) {
    vegan::ordilabel(cbind(tot, S), labels = rownames(x), ...)
  }
  invisible(out)
}


#' convert rarecurve output to data frame
#' @params rare_output rarefaction output
#' @params ps_object a phyloseq object
#' @export
rarecurve_to_df <- function(rare_output, ps_object) {
  # Extract sample names
  if (phyloseq::taxa_are_rows(ps_object)) {
    # If taxa are rows, then samples are columns
    sample_names <- phyloseq::sample_names(ps_object)
  } else {
    # If samples are rows, use row names of OTU table
    sample_names <- rownames(phyloseq::otu_table(ps_object))
  }

  # Convert rarecurve output to data frame
  rare_df <- map_dfr(seq_along(rare_output), function(i) {
    data.frame(
      Sample = sample_names[i],
      Subsample = attr(rare_output[[i]], "Subsample"),
      Species = as.numeric(rare_output[[i]])
    )
  })

  return(rare_df)
}
