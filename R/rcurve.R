#' Rarefaction curves, faster
#' 
#' @description
#' stolen from https://github.com/Russel88/MicEco/blob/master/R/rcurve.R
#' minor adjustments to transform otu as matrix
#' 
#' @export
rcurve <- function (physeq, subsamp = 10^c(1:5), trim = TRUE, add_sample_data = TRUE) {
  message('Not quite efficient, try quickRareCurve() instead !') 
  
  otu <- phyloseq::otu_table(physeq)
  
  if (!phyloseq::taxa_are_rows(physeq)) {
    otu <- t(otu)
  }
  otu <- as(otu, 'matrix') # This is required otherwise vegan says: Error in as(x, "matrix")[i, j, drop = FALSE] :
  otu <- round(otu) # Sourmash may have non-integer abundances
  colS <- colSums(otu)
  pb <- txtProgressBar(min = 0, max = length(subsamp), style = 3)
  rars <- list()
  for (i in seq_along(subsamp)) {
    setTxtProgressBar(pb, i)
    rars[[i]] <- vegan::rarefy(otu, sample = subsamp[i],
                               MARGIN = 2)
  }
  mat <- do.call(cbind, rars)
  if (trim) {
    mat_bool <- sapply(subsamp, function(i) sapply(colS,
                                                   function(j) i <= j))
    mat_new <- mat * mat_bool
    mat_new[mat_new == 0] <- NA
  }
  else {
    mat_new <- mat
  }
  colnames(mat_new) <- subsamp
  df <- as.data.frame.table(mat_new)
  colnames(df) <- c("Sample", "Reads", "Richness")
  df$Reads <- as.numeric(as.character(df$Reads))
  if (trim) {
    df <- na.omit(df)
  }
  if (add_sample_data) {
    samp <- phyloseq::sample_data(physeq)
    df2 <- merge(df, samp, by = "Sample", by.y = "row.names")
    return(df2)
  }
  else {
    return(df)
  }
}

#' Quick Rarefaction Curves
#' 
#' @description
#' Accelerated version of Vegan::rarecurve, original code by Dave Clark https://dave-clark.github.io/post/speeding-up-rarefaction-curves-for-microbial-community-ecology/
#' 
#' Modified to handle matrix as well as phyloseq object, it has more dependencies (furrr/future). 
#' Optimized using Claude Sonnet 5 Medium
#' 
#' Tidy option has been added back.
#' @param x A samples-by-taxa count matrix, data frame, or a `phyloseq` object.
#' @param step Step size for the rarefaction subsampling sequence (number of
#'   individuals between subsample points). Smaller values give smoother
#'   curves at higher computational cost.
#' @param sample Optional. A rarefaction depth at which to draw a reference
#'   line and interpolate species richness for each sample (plot mode only).
#' @param xlab Label for the x-axis (plot mode only). Default `"Sample Size"`.
#' @param ylab Label for the y-axis (plot mode only). Default `"Species"`.
#' @param label Logical. If `TRUE`, label each curve with its sample name via
#'   `vegan::ordilabel()` (plot mode only).
#' @param col Line color(s) for the rarefaction curves (plot mode only).
#'   Recycled to the number of samples. Defaults to `par("col")`.
#' @param lty Line type(s) for the rarefaction curves (plot mode only).
#'   Recycled to the number of samples. Defaults to `par("lty")`.
#' @param tidy Logical. If `TRUE`, return a tidy `tibble` with columns
#'   `Site`, `Sample`, and `Species` instead of producing a plot, matching
#'   `vegan::rarefy(..., tidy = TRUE)` output.
#' @param max.cores Logical. If `TRUE`, use all available cores
#'   (`future::availableCores()`) for parallelization. Overrides `nCores`.
#' @param nCores Number of cores to use when `max.cores = FALSE`.
#' @param ... Additional arguments passed to the underlying plotting
#'   functions (`plot()`, `lines()`, `ordilabel()`); ignored when
#'   `tidy = TRUE`.
#' @export

quickRareCurve <- function(x, step = 1, sample, xlab = "Sample Size",
                           ylab = "Species", label = TRUE, col, lty,
                           tidy = FALSE, max.cores = TRUE, nCores = 1, ...) {
  
  # --- Accept phyloseq objects -----------------------------------------
  if (methods::is(x, "phyloseq")) {
    otu <- phyloseq::otu_table(x)
    # vegan wants samples as rows, taxa as columns
    x <- if (phyloseq::taxa_are_rows(otu)) {
      t(methods::as(otu, "matrix"))
    } else {
      methods::as(otu, "matrix")
    }
  }
  
  x <- as.matrix(x, rownames.force = TRUE)
  if (!isTRUE(all.equal(x, round(x))))
    stop("function accepts only integers (counts)")
  x <- round(x)
  
  minobs <- min(x[x > 0])
  if (minobs > 1)
    warning(gettextf("most observed count data have counts 1, but smallest count is %d", minobs))
  
  if (missing(col)) col <- graphics::par("col")
  if (missing(lty)) lty <- graphics::par("lty")
  
  tot <- rowSums(x)
  S   <- vegan::specnumber(x)
  
  if (any(S <= 0)) {
    message("empty rows removed")
    x   <- x[S > 0, , drop = FALSE]
    tot <- tot[S > 0]
    S   <- S[S > 0]
  }
  
  nr  <- nrow(x)
  col <- rep(col, length.out = nr)
  lty <- rep(lty, length.out = nr)
  
  # --- Parallel setup -----------------------------------------------
  mc <- getOption("mc.cores", ifelse(max.cores, future::availableCores(), nCores))
  
  # For tens-hundreds of rows, per-task process overhead in multisession
  # can rival the actual rarefy() cost, especially at low `step`. Chunking
  # rows into ~mc groups (vs. one future per row) amortizes that.
  # scheduling > 1 tells furrr to bundle roughly that many tasks per chunk.
  message(paste("Using", mc, "cores"))
  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)
  future::plan(future::multisession, workers = mc)
  
  row_list <- asplit(x, 1)
  
  out <- furrr::future_map(seq_len(nr), function(i) {
    n <- seq(1, tot[i], by = step)
    if (n[length(n)] != tot[i]) {
      n <- c(n, tot[i])
    }
    drop(suppressWarnings(vegan::rarefy(row_list[[i]], n)))
  }, .options = furrr::furrr_options(seed = TRUE, scheduling = 2))
  
  # --- Tidy output ----------------------------------------------------
  if (tidy) {
    len <- sapply(out, length)
    nm  <- rownames(x)
    df <- tibble::tibble(
      Site    = factor(rep(nm, len), levels = nm),
      Sample  = unlist(lapply(out, attr, which = "Subsample")),
      Species = unlist(out)
    )
    return(df)
  }
  
  # --- Plot output ------------------------------------------------
  Nmax <- sapply(out, function(z) max(attr(z, "Subsample")))
  Smax <- sapply(out, max)
  
  graphics::plot(c(1, max(Nmax)), c(1, max(Smax)), xlab = xlab, ylab = ylab,
                 type = "n", ...)
  
  if (!missing(sample)) {
    graphics::abline(v = sample)
    rare <- sapply(out, function(z) {
      stats::approx(x = attr(z, "Subsample"), y = z, xout = sample, rule = 1)$y
    })
    graphics::abline(h = rare, lwd = 0.5)
  }
  
  for (ln in seq_along(out)) {
    N <- attr(out[[ln]], "Subsample")
    graphics::lines(N, out[[ln]], col = col[ln], lty = lty[ln], ...)
  }
  
  if (label) {
    vegan::ordilabel(cbind(tot, S), labels = rownames(x), ...)
  }
  
  invisible(out)
}