###
# SBC figures for the SBC study: the 9 skyline rates, plus a row of the
# per-taxon latents (b, d, tau1).
#
# Adapted from skyfbdr_SBC.Rmd (David Cerny) for SBC 0.5 and this repo's layout:
#  - default_diagnostics / messages / warnings must be one per replicate, not empty/NULL
#  - thin_ranks = 1 (there are ~750 retained draws per rep; adaptive thinning wants more)
#  - true values are aligned to the reps that actually produced a complete log
#
# The rates are one value per replicate, so they go through SBC as usual. b, d and
# tau1 are per-taxon, so one extinct taxon is drawn per replicate: taxa within a
# replicate share a history, and pooling them would judge the ranks against a band
# built for far more independent observations than there are. Every panel is therefore
# one rank per replicate. The latents still get their own SBC_results and are stacked
# underneath, because a replicate with no extinct taxon drops out and their n is lower.
#
# At rho = 1 the latent row covers extinct taxa only. Scoring b, d and tau1 together needs
# one rank per taxon for each, and d is pinned at the present for an extant taxon rather
# than estimated, so its rank would be meaningless. Extinction status is read off the
# data, so restricting to it does not disturb the ranks of b or tau1. Below rho = 1 it is
# not data, so every taxon is kept, a survivor's true d is the present, and SBC breaks
# the tie with the draws pinned there at random.
#
# tau1 is internal to the distribution, so master.Rev exposes it as a deterministic node
# (bd.getAugmentedFirstAges()) and mnModel logs it with the rates, on the same generation
# grid as the matrix monitor. tau_K is logged too (tauK[...]) but is not plotted yet.
#
# usage:  simfbdr summarize <batch>/runs/<name> ...
#
# Every path and hyperparameter comes from the run's manifest, so the scoring uses the
# same values the analysis used.

suppressMessages({
  library(SBC)
  library(posterior)
  library(gtools)   # mixedsort is essential for matching output to true values
  library(ggplot2)
  library(parallel)
})

# replicates are read in parallel, by SUMMARIZE_JOBS processes (simfbdr summarize --jobs)
JOBS <- max(1L, as.integer(Sys.getenv("SUMMARIZE_JOBS", detectCores())))

# a worker's error stops the run, and a worker that died returns NULL, which is one too
par_map <- function(x, f, ...) {
  out <- mclapply(x, f, ..., mc.cores = JOBS)
  for (o in out) if (is.null(o) || inherits(o, "try-error")) stop(if (is.null(o)) "a worker died" else o)
  out
}

# a log's column names exactly as RevBayes wrote them
log_header <- function(f) strsplit(readLines(f, n = 1), "\t", fixed = TRUE)[[1]]

# the columns at positions cols, parsing no others. colClasses recycles, so it is built from
# this file's own header
read_columns <- function(f, cols, hdr = log_header(f)) {
  read.table(f, header = TRUE, sep = "\t", check.names = FALSE,
             colClasses = ifelse(seq_along(hdr) %in% cols, "numeric", "NULL"))
}

# a finished log carries one header plus GENS/PRINTGEN + 1 samples
complete_lines_of <- function(mf) as.integer(mf$GENS) / as.integer(mf$PRINTGEN) + 2L

latent_order <- c("d", "tau1", "b")   # youngest to oldest

# Nothing here is inferred from a directory name. simfbdr writes the manifest when it
# runs the analysis, recording the batch, the reporting model and every hyperprior it fed
# to rb, so the scoring reads the same values the analysis used and cannot be pointed at the
# wrong truth.
# A manifest is read through manifest.py, which answers by config key. It is found beside this
# script, so a batch can be anywhere and this can be run from anywhere.
BIN <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
read_manifest <- function(dir) {
  out <- system2(Sys.getenv("PY", "python3"),
                 c(shQuote(file.path(BIN, "manifest.py")), "flat", shQuote(dir)), stdout = TRUE)
  kv <- strsplit(out, "\t", fixed = TRUE)
  setNames(lapply(kv, function(x) paste(x[-1], collapse = "\t")), sapply(kv, `[`, 1))
}

manifest_for <- function(run) {
  m <- read_manifest(run)
  if (!length(m)) stop("no manifest in ", run, "; simfbdr analyze it first")
  m
}

# a non-extended tree marginalizes the extinction times, so there is no d latent to score;
# tau1 and b are still estimated and still worth checking
is_nonextended <- function(mf) identical(mf$EXTENDED, "false")
latent_order_for <- function(mf) if (is_nonextended(mf)) c("tau1", "b") else latent_order

rates_of <- function(mf) {
  ni <- as.integer(mf$NINTERVALS)
  as.vector(t(outer(c("lambda", "mu", "psi"), seq_len(ni),
                    function(a, b) paste0(a, "[", b, "]"))))
}
pair_re <- "\\[ *[-0-9.eE+]+, *[-0-9.eE+]+ *\\]"

post_burnin <- function(n) (round(n * 0.25) + 1):n

# sim.R has reversed each rate triple before writing since 2026-08-03, so a batch whose
# manifest says youngest_first is already aligned to RevBayes. A batch with no such key
# predates the change and holds true_vals oldest-first.
truth_order_of <- function(mf) {
  if (!is.null(mf$true_vals_order)) return(mf$true_vals_order)
  v <- read_manifest(mf$batch_dir)$true_vals_order
  if (length(v)) v else "oldest_first"
}

align_skyline_index <- function(df, ni, order) {
  if (identical(order, "youngest_first")) return(df)
  for (r in c("lambda", "mu", "psi")) {
    cols <- paste0(r, "[", seq_len(ni), "]")
    if (all(cols %in% colnames(df))) df[, cols] <- df[, rev(cols)]
  }
  df
}

# A wrong mapping is indistinguishable from a depth-dependent defect, and it has been wrong
# in both directions, so check it against the data rather than trusting the key: a posterior
# mean correlates better with its own truth column than with the mirrored one.
check_alignment <- function(draws, true_vals, ni, ids) {
  # one interval has no ordering to mirror, and a single column collapses the matrices below
  if (length(draws) < 10 || ni < 2) return(invisible())
  for (r in c("lambda", "mu", "psi")) {
    cols <- paste0(r, "[", seq_len(ni), "]")
    if (!all(cols %in% colnames(true_vals))) next
    pm <- t(sapply(draws, function(d) colMeans(d[, cols, drop = FALSE])))
    kept <- apply(pm, 2, function(x) all(is.finite(x) & x > 0))
    cols <- cols[kept]; pm <- pm[, kept, drop = FALSE]
    tr <- log(true_vals[ids, cols, drop = FALSE])
    own <- mean(mapply(function(a, b) cor(log(pm[, a]), tr[[b]]), cols, cols))
    mir <- mean(mapply(function(a, b) cor(log(pm[, a]), tr[[b]]), cols, rev(cols)))
    if (mir > own)
      cat(sprintf("   WARNING: %s truth looks mirrored (corr %.2f against %.2f reversed);",
                  r, own, mir), "check true_vals_order\n")
  }
}

results_obj <- function(draws) {
  SBC_results(stats = data.frame(sim_id = seq_along(draws)), fits = draws,
              backend_diagnostics = NULL,
              default_diagnostics = data.frame(sim_id = seq_along(draws)),
              outputs = NULL,
              messages = vector("list", length(draws)),
              warnings = vector("list", length(draws)),
              errors = vector("list", length(draws)))
}

# thin_ranks = 1 keeps all ~750 retained draws (SBC's adaptive thinning wants more than
# we have). Autocorrelation within a chain inflates both rank tails, so SBC_THIN_RANKS
# is exposed to tell that artifact apart from a real miscalibration.
thin_ranks_opt <- as.integer(Sys.getenv("SBC_THIN_RANKS", "1"))

stats_for <- function(draws, true_vals) {
  datasets <- SBC_datasets(true_vals, generated = vector("list", length(draws)))
  recompute_SBC_statistics(results_obj(draws), datasets, backend = NULL,
                           thin_ranks = thin_ranks_opt)
}

# lambda_a epsilon_a / (lambda_a epsilon_a + lambda epsilon): a new species differs from its parent in
# each character with chance epsilon_a after an anagenesis and epsilon after a budding, so this is the
# share of all character change that is anagenetic. This is S, 0 for change only at budding. A is
# lambda_a / (lambda_a + lambda), the share of speciations that are anagenetic, and needs no characters
ana_share <- function(d, ni) {
  if (!"lambda_a" %in% names(d)) return(d)
  for (i in seq_len(ni)) {
    l <- d[[paste0("lambda[", i, "]")]]
    d[[paste0("A[", i, "]")]] <- d$lambda_a / (d$lambda_a + l)
    if (!all(c("epsilon", "epsilon_a") %in% names(d))) next
    a <- d$lambda_a * d$epsilon_a
    d[[paste0("S[", i, "]")]] <- a / (a + l * d$epsilon)
  }
  d
}

# the names the rates have always been read under: read.table's, with RevBayes's brackets put back
rate_names <- function(hdr, ni) {
  nm <- gsub("\\[$", "]", gsub("\\.", "[", make.names(hdr, unique = TRUE)))
  # one interval is not a skyline, so the Rev template declares scalar rates and the log
  # names them bare. Everything downstream is indexed, so index them here.
  if (ni == 1) for (a in c("lambda", "mu", "psi")) nm[nm == a] <- paste0(a, "[1]")
  nm
}

# what is scored once per replicate, and its truth, from one log's header
rate_plan <- function(mf, hdr) {
  sim_dir <- mf$batch_dir
  rates <- rates_of(mf)

  # the simulator records the origin as `age`; score it alongside the rates when it was logged
  probe <- make.names(hdr, unique = TRUE)
  if ("org" %in% probe) rates <- c(rates, "org")
  # a gamma rate parameter the batch left out is drawn per replicate and estimated, so it is a
  # parameter like any other and its ranks have to be uniform too
  for (h in c("gamma_rate_lambda", "gamma_rate_mu")) {
    if (h %in% probe) rates <- c(rates, h)
  }
  # so are the psi mixture's shape, the characters' chance to change and the anagenetic rate, when
  # the run estimated them
  for (h in c("alpha", "epsilon", "epsilon_a", "lambda_a", "beta", "retention", "char_rate", "turnover", "mix_r")) {
    if (h %in% probe) rates <- c(rates, h)
  }
  # the invisible mixture's weights, a simplex, ranked component by component
  rates <- c(rates, sub("^mix_w\\.(\\d+)\\.$", "mix_w[\\1]", grep("^mix_w\\.\\d+\\.$", probe, value = TRUE)))

  # with anagenesis, A, and with both chances to change too, S, in each interval
  ni <- as.integer(mf$NINTERVALS)
  share <- character(0)
  if ("lambda_a" %in% probe) share <- paste0("A[", seq_len(ni), "]")
  if (all(c("epsilon", "epsilon_a", "lambda_a") %in% probe)) share <- c(share, paste0("S[", seq_len(ni), "]"))

  true_vals <- read.table(file.path(sim_dir, "true_vals.tsv"), header = TRUE)
  colnames(true_vals) <- gsub("(\\D+)(\\d+)$", "\\1[\\2]", colnames(true_vals))
  names(true_vals)[names(true_vals) == "age"] <- "org"
  true_vals <- align_skyline_index(true_vals, ni, truth_order_of(mf))
  if (length(share)) true_vals <- ana_share(true_vals, ni)
  # a run under another character model than the batch's estimates parameters that have no truth
  rates <- intersect(rates, names(true_vals)); share <- intersect(share, names(true_vals))

  list(cols = c(rates, share), true_vals = true_vals, ni = ni)
}

# one log's posterior check and its post-burn-in draws of the scored columns, in one pass
read_rates <- function(f, plan, complete_lines) {
  if (length(readLines(f)) != complete_lines) return(list(status = "short"))
  hdr <- log_header(f)
  nm <- rate_names(hdr, plan$ni)
  keep <- which(nm %in% c("Posterior", plan$cols))
  d <- read_columns(f, keep, hdr)
  colnames(d) <- nm[keep]
  # a chain stranded at -inf accepts everything and its rates diffuse, but its log is still
  # full length, so length alone does not catch it
  if (!all(is.finite(d$Posterior))) return(list(status = "nonfinite"))
  d <- ana_share(d, plan$ni)
  list(status = "ok", draws = posterior::as_draws_matrix(d[post_burnin(nrow(d)), plan$cols]))
}

# the skyline rates, one value per replicate
rate_stats <- function(plan, ids, draws) {
  check_alignment(draws, plan$true_vals, plan$ni, ids)
  stats_for(draws, posterior::as_draws_matrix(plan$true_vals[ids, plan$cols]))
}

# b, d and tau1, one extinct taxon per replicate. Returns NULL if nothing usable was logged.
latent_stats <- function(mf, ids) {
  sim_dir <- mf$batch_dir
  present <- suppressWarnings(as.numeric(mf$PRESENT))   # "none" when the batch sets no present
  if (!length(present) || is.na(present)) present <- 0
  status_is_data <- is.null(mf$RHO) || as.numeric(mf$RHO) == 1
  nonext <- is_nonextended(mf)
  load(file.path(sim_dir, "sim_list.RData"))   # -> sims
  n_draw <- complete_lines_of(mf) - 1L          # every log in ids is complete

  # a replicate's scorable taxa, each with its log columns and its truth, read from the header
  # alone unless b's pinning has to be checked
  candidates <- function(rep) {
    spec <- file.path(sim_dir, mf$occurrences, sprintf("taxa_%d.tsv", rep))
    rng <- file.path(sim_dir, "true_ranges", sprintf("ranges_%d.tsv", rep))
    log <- file.path(mf$run_dir, "output", sprintf("rep_%d.log", rep))
    if (!all(file.exists(spec, rng, log))) return(list())

    # RevBayes orders taxa by byte order of name (t10 precedes t2); radix sorts in the C locale
    taxa <- sort(unique(as.character(read.table(spec, header = TRUE)$taxon)), method = "radix")
    ranges <- read.table(rng, header = TRUE)
    s <- sims[[rep]]
    n_taxa <- length(taxa)
    sim_idx <- as.integer(sub("t", "", taxa))   # simulator indices of the sampled taxa

    hdr <- log_header(log)
    tcols <- grep("^tau1\\[", hdr)
    if (length(tcols) != n_taxa || n_draw < 50) return(list())
    # b and d are deterministic nodes in the mnModel log, so they work for a tree model too
    bcols <- grep("^bb\\[", hdr)
    dcols <- grep("^dd\\[", hdr)
    if (length(bcols) != n_taxa || (length(dcols) != n_taxa && !nonext))
      return(list())   # the run did not monitor the per-taxon latents, so there are no latent ranks

    # In a non-extended tree an I[i] taxon's b is pinned to its sampled ancestor's range end, with
    # the speciation that separates them integrated out by the I-set term, so b is a structural
    # node age there rather than an estimated birth. Detect it as b coinciding with some taxon's
    # tau_K. This is a posterior summary, hence a function of the data, so the split is legitimate.
    kcols <- grep("^tauK\\[", hdr)
    b_mat <- k_mat <- NULL
    if (nonext && length(kcols) == n_taxa) {
      ll <- read_columns(log, c(bcols, kcols), hdr)
      keep <- post_burnin(nrow(ll))
      b_mat <- as.matrix(ll[keep, hdr[bcols], drop = FALSE])
      k_mat <- as.matrix(ll[keep, hdr[kcols], drop = FALSE])
    }

    # Which start the run logged decides which truth scores it, and the two differ by the chain
    # of unsampled ancestors. A divergence is where this lineage last branched off a SAMPLED
    # one: the stem above o_i uses the unsampled propagator q, which permits any number of
    # species changes between the branching and the first fossil, so it takes a walk up through
    # the unsampled ancestors. An origination is the taxon's own TS, with no walk. Scoring one
    # against the other would read the unsampled chain as bias.
    divergence <- identical(mf$RANGE_START, "divergence")
    true_b_all <- numeric(length(taxa))
    founder <- logical(length(taxa))
    off_root <- numeric(length(taxa))
    for (j in seq_along(taxa)) {
      k <- as.integer(sub("t", "", taxa[j]))
      if (!divergence) {
        true_b_all[j] <- s$TS[k]
        next
      }
      cur <- k
      par <- s$PAR[cur]
      while (!is.na(par) && !(par %in% sim_idx)) { cur <- par; par <- s$PAR[cur] }
      true_b_all[j] <- s$TS[cur]
      founder[j] <- is.na(par)        # the walk ran out of ancestors, so the chain reaches species 1
      # where that chain left species 1's lineage, for the rows that cannot hold the origin
      ch <- k
      pp <- s$PAR[ch]
      while (!is.na(pp) && pp != 1) { ch <- pp; pp <- s$PAR[ch] }
      off_root[j] <- s$TS[ch]
    }

    # The density allows one zero-gamma row and it must start at the origin, so only one chain
    # that reaches species 1 may claim it. Budding keeps species 1's lineage going, so the chain
    # that left it LAST is its continuation and holds the origin; the rest branched off that
    # lineage at off_root. Claiming it for all of them puts the truth outside the support.
    fi <- which(founder)
    if (length(fi) > 1) {
      rest <- setdiff(fi, fi[which.min(off_root[fi])])
      true_b_all[rest] <- off_root[rest]
    }

    cands <- list()
    for (j in seq_along(taxa)) {
      k <- as.integer(sub("t", "", taxa[j]))
      true_d <- max(if (is.na(s$TE[k])) 0 else s$TE[k], present)
      # at rho = 1 survival is data and pins d; below it, skipping survivors would select on the truth
      if (true_d <= present && status_is_data) next

      # any pinning at all makes b a mixture of a structural node age and an estimated birth,
      # so score it only where the birth is explicit in every retained sample
      if (!is.null(k_mat)) {
        pinned <- mean(rowSums(abs(k_mat - b_mat[, j]) < 1e-6 * pmax(1, abs(b_mat[, j]))) > 0)
        if (pinned > 0) next
      }
      first_age <- ranges$first_age[as.character(ranges$taxon) == taxa[j]]
      if (length(first_age) != 1) next
      cands[[length(cands) + 1]] <- if (nonext)
        list(cols = c(b = bcols[j], tau1 = tcols[j]), truth = c(b = true_b_all[j], tau1 = first_age), taxon = taxa[j])
      else
        list(cols = c(b = bcols[j], d = dcols[j], tau1 = tcols[j]),
             truth = c(b = true_b_all[j], d = true_d, tau1 = first_age), taxon = taxa[j])
    }
    list(cands = cands, n_taxa = n_taxa)
  }
  found <- par_map(ids, candidates)

  # one taxon per replicate: taxa within a replicate share a history, so their ranks are
  # not independent and pooling them narrows the band the ranks are judged against.
  # Drawn here in replicate order, so the same seed picks the same taxa as a serial read.
  set.seed(20260719)   # the per-replicate taxon draw, fixed so the figures reproduce
  picked <- list()
  for (i in seq_along(ids)) {
    cands <- found[[i]]$cands
    if (!length(cands)) next
    picked[[length(picked) + 1]] <- c(cands[[sample.int(length(cands), 1)]],
                                      rep = ids[i], n_taxa = found[[i]]$n_taxa)
  }
  if (!length(picked)) return(NULL)

  # the chosen taxon's three columns, and nothing else, from each log
  draws <- par_map(picked, function(p) {
    log <- file.path(mf$run_dir, "output", sprintf("rep_%d.log", p$rep))
    hdr <- log_header(log)
    ll <- read_columns(log, p$cols, hdr)
    m <- as.matrix(ll[post_burnin(nrow(ll)), hdr[p$cols], drop = FALSE])
    colnames(m) <- names(p$cols)
    posterior::as_draws_matrix(m)
  })

  cat("   latents: one", if (status_is_data) "extinct taxon" else "taxon", "from each of", length(draws), "replicates\n")
  res <- stats_for(draws, posterior::as_draws_matrix(do.call(rbind, lapply(picked, `[[`, "truth"))))
  # sim_id is a position in this list, not a replicate number, so the pair is recorded here
  # or there is no way to split a rank by anything about the record it came from
  res$scored_from <- data.frame(sim_id = seq_along(picked),
                                rep = vapply(picked, `[[`, 0L, "rep"),
                                taxon = vapply(picked, `[[`, "", "taxon"),
                                n_taxa = vapply(picked, `[[`, 0L, "n_taxa"))
  # youngest to oldest, the order the range is read in (rank_summary splits on this;
  # the plots take the same order through the variables= argument)
  res$stats$variable <- factor(res$stats$variable, levels = latent_order_for(mf))
  res
}

# SBC in one number per parameter: under calibration the rank of the true value among
# the posterior draws is uniform, so mean(rank/max) is 0.5 and z is standard normal.
#
# z only sees a location shift. A symmetric departure -- the U-shape that autocorrelated
# draws produce, where the rank piles up at both ends -- leaves the mean at 0.5 and is
# invisible to z, so also test the ranks for uniformity outright.
rank_summary <- function(res, what) {
  s <- res$stats
  d <- do.call(rbind, lapply(split(s, s$variable), function(x) {
    r <- x$rank / x$max_rank
    p <- suppressWarnings(chisq.test(table(cut(r, seq(0, 1, length.out = 21),
                                               include.lowest = TRUE)))$p.value)
    data.frame(variable = x$variable[1], n = nrow(x), mean_rank = mean(r),
               z = (mean(r) - 0.5) / sqrt(1 / 12 / nrow(x)), unif_p = p)
  }))
  d$flag <- ifelse(d$unif_p < 0.01 | abs(d$z) > 3, "  <-- NOT uniform", "")
  cat("  ", what, ": mean|rank-0.5| =", sprintf("%.4f", mean(abs(d$mean_rank - 0.5))), "\n")
  for (i in seq_len(nrow(d)))
    cat(sprintf("     %-10s rank %.4f  z %+6.2f  unif p %.3g%s\n",
                d$variable[i], d$mean_rank[i], d$z[i], d$unif_p[i], d$flag[i]))
  invisible(d)
}

# One 4x3 grid: the 9 rates then b, d, tau1, with a single legend under the whole thing.
# The two halves stay separate plot objects because their n differs (one rank per
# replicate for a rate, one per taxon for a latent) and so does their confidence band;
# only the bottom half carries the legend and the axis strip, so they read as one grid.
save_stacked <- function(path, top, bottom, height = 8) {
  # png only: the svgs are large and nothing reads them
  for (dev_name in c("png")) {
    png(paste0(path, ".png"), width = 1100, height = height * 110)
    grid::grid.newpage()
    n_row <- if (is.null(bottom)) 1 else 2
    heights <- if (is.null(bottom)) 1 else c(3, 1.5)
    grid::pushViewport(grid::viewport(layout = grid::grid.layout(n_row, 1,
                                       heights = grid::unit(heights, "null"))))
    print(top, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
    if (!is.null(bottom)) print(bottom, vp = grid::viewport(layout.pos.row = 2, layout.pos.col = 1))
    invisible(dev.off())
  }
}

# SBC carries the panel name as a character, so the facets come out alphabetical
order_facets <- function(p, lv) {
  for (v in c("group", "variable")) {
    if (!is.null(p$data) && v %in% names(p$data)) p$data[[v]] <- factor(p$data[[v]], levels = lv)
  }
  p
}

make_figures <- function(run) {
  run <- sub("/$", "", run)
  mf <- manifest_for(run)
  mf$run_dir <- run
  out <- file.path(run, "results")
  dir.create(out, recursive = TRUE, showWarnings = FALSE)

  logs <- mixedsort(list.files(file.path(run, "output"),
                               pattern = "^rep_.*log$", full.names = TRUE))
  if (!length(logs)) { cat(run, ": no complete logs, skipping\n"); return(invisible()) }
  plan <- rate_plan(mf, log_header(logs[[1]]))
  got <- par_map(logs, read_rates, plan, complete_lines_of(mf))
  status <- vapply(got, `[[`, "", "status")
  if (any(status == "nonfinite")) cat("  ", sum(status == "nonfinite"), "rep(s) dropped: non-finite posterior\n")
  logs <- logs[status == "ok"]
  draws <- setNames(lapply(got[status == "ok"], `[[`, "draws"), logs)
  if (!length(logs)) { cat(run, ": no complete logs, skipping\n"); return(invisible()) }
  ids <- as.integer(gsub("[^0-9]", "", basename(logs)))

  label <- sprintf("%s (batch %s, %s, reporting=%s)", basename(run), mf$batch,
                   mf$REPORTING, mf$INFER)
  cat(label, ":", length(logs), "complete reps\n")

  res <- rate_stats(plan, ids, draws)
  lat <- latent_stats(mf, ids)
  if (is.null(lat)) cat("   no tau1 traces: plotting the rates only\n")

  rank_summary(res, "skyline rates")
  if (!is.null(lat)) rank_summary(lat, "per-taxon latents")
  saveRDS(list(rates = res, latents = lat, manifest = mf),
          file.path(out, "sbc_results.rds"))

  # the run names the whole grid; the legend rides on the bottom half only
  upper <- function(p) p + labs(title = label) + theme(legend.position = "none")
  lower <- function(p) order_facets(p, latent_order_for(mf)) + theme(legend.position = "bottom")
  stack <- function(name, plot_fn) save_stacked(file.path(out, name),
             upper(plot_fn(res)), if (!is.null(lat)) lower(plot_fn(lat)))
  stack("output_ecdf_diff", plot_ecdf_diff)
  cat("   figure ->", file.path(out, "output_ecdf_diff.png"), "\n")
}

runs <- commandArgs(TRUE)
if (!length(runs)) stop("usage: simfbdr summarize <batch>/runs/<name> ...")
for (r in runs) make_figures(r)
