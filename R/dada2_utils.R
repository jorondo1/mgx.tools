#' Find all orientation of primers in the data
#' stolen from DADA2 tutorial https://benjjneb.github.io/dada2/tutorial.html
#' @params primer primer sequence as string
#' @keywords internal
.allOrients <- function(primer) {
  # Create all orientations of the input sequence
  dna_string <- Biostrings::DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna_string,
               Complement = Biostrings::complement(dna_string),
               Reverse = Biostrings::reverse(dna_string),
               RevComp = Biostrings::reverseComplement(dna_string))
  return(sapply(orients, toString))  # Convert back to character vector
}

#' Primer counter
#' Count the number of times the primers appear in the forward and reverse read,
#' while considering all possible primer orientations
#' @keywords internal
.primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- Biostrings::vcountPattern(
    primer, ShortRead::sread(ShortRead::readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}


#' Compute primer occurence across all orientations
#' @param fnFs Forward read files
#' @param fnRs Reverse read files
#' @param FWD Forward primer sequence
#' @param REV Reverse primer sequence
#' @param ncores Number of cores to use (default: detectCores() - 1)
#' @export
primer_occurence <- function(fnFs, fnRs, FWD, REV, ncores = NULL){

  # Set up parallel processing
  if (is.null(ncores)) {
    ncores <- max(1, parallel::detectCores() - 1)
  }

  future::plan(future::multisession, workers = ncores)
  on.exit(future::plan(future::sequential), add = TRUE)

  # Get all primer orientations
  FWD.orients <- .allOrients(FWD)
  REV.orients <- .allOrients(REV)

  # Parallel computation
  FWD.ForwardReads <- future.apply::future_sapply(FWD.orients, .primerHits, fn = fnFs[[1]])
  FWD.ReverseReads <- future.apply::future_sapply(FWD.orients, .primerHits, fn = fnRs[[1]])
  REV.ForwardReads <- future.apply::future_sapply(REV.orients, .primerHits, fn = fnFs[[1]])
  REV.ReverseReads <- future.apply::future_sapply(REV.orients, .primerHits, fn = fnRs[[1]])

  # Return results
  rbind(FWD.ForwardReads = FWD.ForwardReads,
        FWD.ReverseReads = FWD.ReverseReads,
        REV.ForwardReads = REV.ForwardReads,
        REV.ReverseReads = REV.ReverseReads)
}

#' CUTADAPT wrapper function
#' Adapted cutadapt in function to allow using it within mclapply
#' @export
run_cutadapt <- function(i, cutadapt_path, R1.flags, R2.flags) {
  # Check cutadapt exists
  cutadapt_check <- suppressWarnings(
    system2(cutadapt_path, args = "--version", stdout = TRUE, stderr = TRUE)
  )

  if (length(cutadapt_check) == 0) {
    stop(
      "cutadapt not found. Please install it first.\n",
      "pass the \n",
      "Or make sure it's in your PATH",
      call. = FALSE
    )
  }

  system2(
    cutadapt_path, args = c(
      R1.flags, R2.flags,
      "-n", 2,
      "-m", 21, '-M', 300, # see https://github.com/benjjneb/dada2/issues/2045#issuecomment-2449416862
      "-o", fnFs.cut[i],
      "-p", fnRs.cut[i],
      fnFs.filtN[i], fnRs.filtN[i])
  )
}

#' Check if an expected file exists
#' Filtering with minlen may yield empty samples (e.g. neg. controls);
#' list files that did survive filtering:
#' @export

dropped_samples <- function(path_list) {
  dropped <- path_list[!file.exists(path_list)]
  num_dropped <- length(dropped)
  if(num_dropped>0) {
    message(paste(num_dropped,"sample(s) didn't pass filtering:"))
    message(paste(basename(dropped), collapse="\n"))
    message('Sample(s) removed from the list of expected output paths.')
  } else {
    message("All samples passed filtering!")
  }
  path_list[file.exists(path_list)]
}

#' Report chimera rate and reads
#' @export

chimera_report <- function(seqtab, seqtab.nochim) {
  ASV1 <- ncol(seqtab)
  ASV2 <- ncol(seqtab.nochim)
  SEQ1 <- sum(seqtab)
  SEQ2 <- sum(seqtab.nochim)
  chim <- ASV1-ASV2
  seqloss <- round(100*(SEQ1-SEQ2)/SEQ1,2)
  message(paste(chim, "chimeras were found out of", ASV1, "ASVs."))
  message(paste(seqloss,"% of sequences were lost."))
}


#' Plot ASV counts with minimum sequence count
#' @export
minimum_ASV_count <- function(seqtab.nochim){
  # Loop through n from 0 to k
  column_counts <- numeric()
  for (n in 1:100) {
    # Filter columns where colSums > n
    filtered_table <- seqtab.nochim[, colSums(seqtab.nochim) >= n, drop = FALSE]

    # Count columns and store result
    column_counts[n] <- ncol(filtered_table)  # n+1 because R indices start at 1
  }

  df <- data.frame(sequence_sum = seq(1:100),
                   asv_count = column_counts)

  ggplot(df, aes(x = sequence_sum, y = asv_count)) +
    geom_point() +
    labs(title = 'ASV with minimum total sequence count',
         x = 'Minimum number of sequences across all samples',
         y = 'ASV count',
         caption = 'Use this to find a conservative minimum below which to drop ASVs with fewer total sequence counts')
}

#' Remove ultra rare asvs
#' Remove samples with fewer than n sequences
#' as well as singletons, if any (shouldn't)
#' @export
drop_rare_asvs <- function(seqtab.nochim, at_least_n) {
  filtered_ASV <- seqtab.nochim[ , colSums(seqtab.nochim) >= at_least_n , drop = FALSE]
  ASV2 <- ncol(seqtab.nochim)
  ASV3 <- ncol(filtered_ASV)
  rare <- ASV2-ASV3
  message(paste(rare, "ASVs dropped because they had fewer than\n", at_least_n, "sequences across all datasets."))
  filtered_samples <- filtered_ASV[rowSums(filtered_ASV) >0 , , drop = FALSE]
  dropped_sample <- filtered_ASV[rowSums(filtered_ASV) == 0 , , drop = FALSE]
  samples1 <- nrow(filtered_ASV)
  samples2 <- nrow(filtered_samples)
  empty <- samples1-samples2
  if(empty>0) {
    message(paste(empty, "sample(s) dropped because they no longer had sequences:"))
    message(paste(rownames(dropped_sample), collapse="\n"))
  }
  message(paste("Sequence table has",nrow(filtered_samples), "samples and", ncol(filtered_samples), 'ASVs.'))
  return(filtered_samples)
}

#' Read counter
#' @keywords internal
.getN <- function(x) sum(getUniques(x))

#' Read tracker
#' compile reads across samples in long format
#' also computes changes between steps (last column flags this since format is long)
#' @export
track_dada <- function(out.N, out, sample.names,
                       dadaFs, dadaRs = NULL,
                       mergers_pooled = mergers_pooled,
                       seqtab.nochim = seqtab.nochim) {

  track <- cbind(out.N, out[,2])
  rownames(track) <- sample.names
  track <- track[which(track[,3]>0),]
  track <- cbind(track, sapply(dadaFs, .getN))

  # Conditionally add dadaRs column if it exists
  if(!is.null(dadaRs)) {
    track <- cbind(track, sapply(dadaRs, .getN))
    column_names <- c("input", "removeNs", "filtered", "denoisedF", "denoisedR", "raw_seqtab", "nonchim")
  } else {
    column_names <- c("input", "removeNs", "filtered", "denoisedF", "raw_seqtab", "nonchim")
  }

  # Add the remaining columns
  track <- cbind(track,
                 rowSums(seqtab),
                 rowSums(seqtab.nochim))

  colnames(track) <- column_names

  # Full tibble with counts and changes
  tibble_out <- data.frame(track) %>%
    tibble::rownames_to_column('Sample') %>%
    tibble::tibble() %>%
    dplyr::mutate(Nfilt_lost = (input-removeNs)/input,
                  filterAndTrim_lost = (removeNs-filtered)/removeNs,
                  denoising_lost = (filtered-denoisedF)/filtered,
                  merging_lost = (denoisedF-raw_seqtab)/denoisedF, # Proportion of reads lost to merging
                  bimera_lost = (raw_seqtab-nonchim)/raw_seqtab) %>%
    tidyr::pivot_longer(where(is.numeric), names_to = 'variable', values_to = 'values')

  list(
    counts_per_step = tibble_out %>% filter(!str_detect(variable, '_lost')),
    lost_per_step = tibble_out %>% filter(str_detect(variable, '_lost'))
  )
}


#' Plot track changes
#' @export
plot_track_change <- function(track_change) {

  p1 <- track_change[['counts_per_step']] %>%
    dplyr::mutate(variable = factor(variable, levels = c("input", "removeNs", "filtered", "denoisedF", "raw_seqtab", "nonchim"))) %>%
    ggplot2::ggplot(aes(x = variable, y = values)) +
    ggplot2::geom_line(aes(group = Sample), linewidth = 0.1) +
    ggrain::geom_rain() +
    ggplot2::labs(
      title = 'Read counts across main DADA2 pipeline steps.',
      x = "", y = 'Sequence count')

  p2 <- track_change[['lost_per_step']] %>%
    dplyr::mutate(variable = factor(variable, level = c('bimera_lost', 'merging_lost', 'denoising_lost', 'filterAndTrim_lost', 'Nfilt_lost'))) %>%
    ggplot2::ggplot(aes(y = variable, x = values)) +
    ggplot2::geom_jitter(height = 0.2, width =0) +
    ggplot2::labs(
      title = 'Proportion of reads lost at a specific pipeline step.',
      x = "Proportion reads lost",
      y = "") +
    xlim(0,1)

  p1/p2 &
    theme_minimal()

}


#' Export ASVs as fasta
#' @param seqtab matrix sequence table with samples as rows
#' @param path.out path where to save fasta
#' @export
asv_to_fasta <- function(seqtab, path.out) {
  seqs <- colnames(seqtab)
  fasta <- Biostrings::DNAStringSet(seqs)
  names(fasta) <- paste0("ASV_", seq_along(seqs))
  Biostrings::writeXStringSet(fasta, path.out)
}










