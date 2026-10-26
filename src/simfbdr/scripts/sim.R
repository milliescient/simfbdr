######################################################
##                    Chapter 3                     ##
##      Treeless skyline FBDR simulation study      ##
##                Simulation script                 ##
##            Bruno do Rosario Petrucci             ##
######################################################

###
# load packages

# paleobuddy 
library(paleobuddy)

# The fork's closed-form draw carries CRAN's version string, so a reinstall replaces it
# silently and shows up only as a run that is orders of magnitude slower.
if (!exists("rexp.var") || !any(grepl("stepfun", deparse(body(rexp.var))))) {
  warning("paleobuddy has no stepfun fast path, so this will be far slower than it should be",
          call. = FALSE, immediate. = TRUE)
}

# ape
library(ape)
library(TreeSim)
library(FossilSim)




###
# write minor auxiliary functions

# smarter dir.create function
smart_dir_create <- function(dir) {
  if (!dir.exists(dir)) dir.create(dir)
}

# equivalent to colSums for max
colMaxes <- function(df) {
  # make return
  res <- c()
  
  # iterate through columns
  for (i in 1:ncol(df)) {
    # add max to res
    res <- c(res, max(df[, i]))
  }
  
  # name res
  names(res) <- colnames(df)
  
  # return res
  return(res)
}

# same for min
colMins <- function(df) {
  # make return
  res <- c()
  
  # iterate through columns
  for (i in 1:ncol(df)) {
    # add max to res
    res <- c(res, min(df[, i]))
  }
  
  # name res
  names(res) <- colnames(df)
  
  # return res
  return(res)
}

# making time bins
make_bins <- function(age, unc) {
  # find range based on unc
  range <- switch(unc,
                  "mid" = {c(age / 6, age / 5)},
                  "low" = {c(age / 10, age / 8)},
                  "high" = {c(age / 4, age / 3)})
  
  # create bins vector
  bins <- c(0)
  
  # while bins doesn't have age
  while (!(age %in% bins)) {
    # if the highest bin is less than the max range, just add age
    if (max(bins) + range[2] > age) {
      bins <- c(bins, age)
    } else {
      # if not, draw a new bin
      bin <- runif(1, range[1], range[2])
      
      # if max(bins) + bin is higher than age, add age instead
      if (max(bins) + bin > age) {
        bins <- c(bins, age)
      } else {
        # add to bins
        bins <- c(bins, max(bins) + bin)
        
      }
    }
  }
  
  # return bins
  return(bins)
}

###
# SBC simulation constants
#
# Every generative parameter comes from the batch's config.yml, which simfbdr exports into
# the environment. A missing key is fatal: a silent default is how the origin prior drifted
# from the value the analysis assumed.
cfg <- function(k) {
  v <- Sys.getenv(k)
  if (v == "") stop("config key ", k, " is not set; simfbdr must export it")
  as.numeric(v)
}

# a key a batch may leave unset
cfg_opt <- function(k, default) {
  v <- Sys.getenv(k)
  if (v == "") default else as.numeric(v)
}

NINTERVALS    <- as.integer(cfg("NINTERVALS"))
INTERVAL_WIDTH <- cfg("INTERVAL_WIDTH")
# PyRate's rate prior: an independent gamma draw per interval. A speciation or extinction
# rate parameter left unset is PyRate's hyperprior, drawn per replicate and estimated.
GAMMA_L_SHAPE <- cfg("GAMMA_L_SHAPE"); GAMMA_L_RATE <- cfg_opt("GAMMA_L_RATE", NA)
GAMMA_M_SHAPE <- cfg_opt("GAMMA_M_SHAPE", NA); GAMMA_M_RATE <- cfg_opt("GAMMA_M_RATE", NA)
# Or one turnover shared by every interval, mu = turnover * lambda: Beta keeps it below 1, LogNormal
# lets a clade decline. Exactly one of mu's gamma and a turnover prior is stated.
TURN_A <- cfg_opt("TURN_A", NA); TURN_B <- cfg_opt("TURN_B", NA)
TURN_LMEAN <- cfg_opt("TURN_LMEAN", NA); TURN_LSD <- cfg_opt("TURN_LSD", NA)
TURNOVER <- if (!is.na(TURN_A) || !is.na(TURN_B)) "beta" else
            if (!is.na(TURN_LMEAN) || !is.na(TURN_LSD)) "lognormal" else "none"
if (TURNOVER == "beta" && (is.na(TURN_A) || is.na(TURN_B) || !is.na(TURN_LMEAN) || !is.na(TURN_LSD)))
  stop("rate_prior.turnover.beta needs alpha and beta, and no lognormal beside it")
if (TURNOVER == "lognormal" && (is.na(TURN_LMEAN) || is.na(TURN_LSD)))
  stop("rate_prior.turnover.lognormal needs log_mean and log_sd")
if ((TURNOVER == "none") == is.na(GAMMA_M_SHAPE))
  stop("state the extinction prior as rate_prior.mu or as rate_prior.turnover, one of the two")
GAMMA_P_SHAPE <- cfg("GAMMA_P_SHAPE"); GAMMA_P_RATE <- cfg("GAMMA_P_RATE")
GAMMA_RATE_MAX <- cfg_opt("GAMMA_RATE_MAX", 20)
AGE_MIN <- cfg("AGE_MIN"); AGE_MAX <- cfg("AGE_MAX")
BIN_WIDTH <- cfg("BIN_WIDTH")

# a key holding several numbers, comma separated; unset is none
cfg_list <- function(k) {
  v <- Sys.getenv(k)
  if (v == "") numeric(0) else as.numeric(strsplit(v, ",")[[1]])
}
# The same record reported to other grids as well, every width at every offset, so analyses that
# differ only in dating resolution are paired on one set of fossils. Width 0 is the exact ages.
# An offset is a fraction of the width: 0.5 moves a grid's edges by half a bin.
BIN_WIDTHS <- cfg_list("BIN_WIDTHS")
BIN_OFFSETS <- cfg_list("BIN_OFFSETS")
if (!length(BIN_OFFSETS)) BIN_OFFSETS <- 0
if (anyNA(BIN_WIDTHS) || any(BIN_WIDTHS < 0)) stop("BIN_WIDTHS must be numbers, none below 0")
if (anyNA(BIN_OFFSETS) || any(BIN_OFFSETS < 0 | BIN_OFFSETS >= 1))
  stop("BIN_OFFSETS are fractions of a bin, from 0 up to but not 1")
more_grids <- unique(do.call(rbind, c(list(data.frame(w = numeric(0), off = numeric(0))),
  lapply(BIN_WIDTHS, function(w) data.frame(w = w, off = if (w == 0) 0 else BIN_OFFSETS)))))
grid_label <- function(w, off) paste0("w", w, if (off != 0) paste0("_o", off))

# An age's cell on a grid, (lower, upper] as paleobuddy's are. The grid stops at PRESENT, and a
# bin straddling AGE_MAX stops there, which the analysis knows anyway, so an offset grid has no
# bin reaching past the origin prior. Bins wholly above it, which only a tree drawn backward from
# its survivors reaches, are left whole.
rebin <- function(age, w, off) {
  if (w == 0) return(list(min = age, max = age))
  hi <- w * (off + ceiling(age / w - off))
  list(min = pmax(hi - w, PRESENT), max = ifelse(hi - w < AGE_MAX, pmin(hi, AGE_MAX), hi))
}

# Off the grid, the top occurrence bin runs past AGE_MAX, so no origin the prior can draw
# covers it and tau_1 is squeezed into what is left.
if (abs(AGE_MAX %% BIN_WIDTH) > 1e-9) {
  warning(sprintf("AGE_MAX %g is off the BIN_WIDTH %g grid, so the top bin extends past it",
                  AGE_MAX, BIN_WIDTH), call. = FALSE, immediate. = TRUE)
}

# Ceiling on lineages ever born. Unset means no ceiling, which is right when the origin
# is shallow. Deeper origins need one: the oldest interval absorbs the extra depth, so
# diversity is exponential in it and a draw in the upper tail of lambda never terminates.
MAX_LINEAGES <- cfg_opt("MAX_LINEAGES", Inf)
# Ceiling on the work one bd.sim call may do, counted in the waiting times it draws. A count
# and not a number of seconds, so which draws are refused follows from the seed alone.
SIM_MAX_DRAWS <- cfg_opt("SIM_MAX_DRAWS", Inf)
if (Sys.getenv("SIM_TIMEOUT") != "")
  stop("SIM_TIMEOUT is a wall-clock limit, so a batch made under it does not reproduce from ",
       "its seeds. Set SIM_MAX_DRAWS instead: 100000 refuses about what 20 seconds did, 200000 ",
       "what 60 did.")
RHO <- cfg_opt("RHO", 1)                     # chance an extant lineage is seen at the present
# The record stops here rather than at zero. Occurrences younger than PRESENT are discarded
# and a lineage still alive there is a boundary survivor, which rho then reports or misses.
PRESENT <- cfg_opt("PRESENT", 0)
if (PRESENT > 0 && abs(PRESENT %% BIN_WIDTH) > 1e-9) {
  warning(sprintf("PRESENT %g is off the BIN_WIDTH %g grid, so a bin straddles it",
                  PRESENT, BIN_WIDTH), call. = FALSE, immediate. = TRUE)
}
# Fossil-sampled taxa a replicate must hold to be kept, both bounds inclusive. Rejection is on
# a function of the record and redraws theta, so it leaves the posterior alone. At rho < 1 the
# extant singletons are added after this test, so it counts the fossil record only.
MIN_TAXA <- as.integer(cfg_opt("MIN_TAXA", 1))
MAX_TAXA <- cfg_opt("MAX_TAXA", Inf)
# "record" makes MIN_TAXA and MAX_TAXA count every reported taxon, the survivors seen at the
# present included; unset counts the fossil-sampled ones only.
TAXA_COUNT <- Sys.getenv("TAXA_COUNT", "fossil")
if (!TAXA_COUNT %in% c("fossil", "record")) stop("TAXA_COUNT must be fossil or record")
# a floor on fossil-sampled taxa alone, for analyses that leave the present-only taxa out
MIN_FOSSIL_TAXA <- as.integer(cfg_opt("MIN_FOSSIL_TAXA", 0))

# Set, trees are drawn backward from exactly this many extant species and the origin is an
# output. AGE_MAX still bounds an offset grid's bins, and both bound a run's uniform origin prior. Rejection on extant count cannot reach the
# near-critical regime: at turnover 0.9 it needs over 480s a tree against 0.2s here.
N_EXTANT <- suppressWarnings(as.integer(cfg_opt("N_EXTANT", NA)))
BACKWARD <- !is.na(N_EXTANT)

# A psi mixture: each species draws one of MIX_K relative sampling rates at its birth,
# independently of its parent, from a discretized gamma of shape alpha with equal weights. Unset
# or 1 is none, and then no extra random number is drawn, so such a batch reproduces exactly.
MIX_K <- as.integer(cfg_opt("MIX_K", 1))
MIX_ALPHA_RATE <- cfg_opt("MIX_ALPHA_RATE", 1)
if (MIX_K < 1) stop("psi_mixture.components must be at least 1")

# Characters that change only at speciation, each of a new species' changing with probability
# epsilon to a uniformly drawn other state. Unset is none, again with no extra draws.
N_CHARS <- as.integer(cfg_opt("N_CHARS", 0))
CHAR_STATES <- cfg_list("CHAR_STATES")
if (N_CHARS > 0) {
  if (!length(CHAR_STATES)) CHAR_STATES <- 2
  if (length(CHAR_STATES) == 1) CHAR_STATES <- rep(CHAR_STATES, N_CHARS)
  if (length(CHAR_STATES) != N_CHARS || any(CHAR_STATES < 2) || any(CHAR_STATES > 10))
    stop("characters.states must be one number or one per character, each from 2 to 10")
}
EPS_A <- cfg_opt("EPS_A", 1); EPS_B <- cfg_opt("EPS_B", 9)

# Anagenesis: each lineage also splits into a new species at rate lambda_a, drawn per replicate
# from the speciation prior and constant in time. Off, nothing extra is drawn.
ANAGENESIS <- Sys.getenv("ANAGENESIS", "false") == "true"

# RevBayes's fnDiscretizeGamma(alpha, alpha, K): each category's mean, so the rates average one
discretized_gamma <- function(alpha, K) {
  q <- c(0, qgamma(seq_len(K - 1) / K, alpha, rate = alpha), Inf)
  K * diff(pgamma(q, alpha + 1, rate = alpha))
}

# Rate breakpoints (before the present) and fossil bins. Both are fixed rather than
# derived from age: the analysis estimates the origin, so anything it reads that
# scales with age hands it the answer. paleobuddy counts shift times forward from
# the start, so these are age - cov_breaks there.
# the intervals partition the process's own span, which starts at PRESENT rather than zero
cov_breaks <- PRESENT + INTERVAL_WIDTH * seq_len(NINTERVALS - 1)
# Occurrence bins are every multiple of BIN_WIDTH, without end. Each replicate builds them as far
# back as its tree reaches; the edges do not depend on how deep that is.
bins_to <- function(age) seq(0, BIN_WIDTH * (floor(age / BIN_WIDTH) + 1), BIN_WIDTH)

# Truth, named, in one place, so the names travel with the values.
truth_of <- function(th, lambda, mu, psi, age) {
  K <- length(lambda)
  v <- c(setNames(rev(lambda), paste0("lambda", seq_len(K))),
         setNames(rev(mu),     paste0("mu", seq_len(K))),
         setNames(rev(psi),    paste0("psi", seq_len(K))))
  v <- c(v, age = age)
  if (is.na(GAMMA_L_RATE)) v <- c(v, gamma_rate_lambda = th$gamma_rate_lambda)
  if (TURNOVER == "none" && is.na(GAMMA_M_RATE)) v <- c(v, gamma_rate_mu = th$gamma_rate_mu)
  if (TURNOVER != "none") v <- c(v, turnover = th$turnover)
  if (MIX_K > 1) v <- c(v, alpha = th$alpha)
  if (N_CHARS > 0) v <- c(v, epsilon = th$epsilon)
  if (ANAGENESIS) v <- c(v, lambda_a = th$lambda_a)
  v
}

# The generative draw, in one place. Both the initial draw and the rejection redraw call
# this, so they cannot diverge.
draw_theta <- function() {
  b_lambda <- if (is.na(GAMMA_L_RATE)) runif(1, 0, GAMMA_RATE_MAX) else GAMMA_L_RATE
  if (TURNOVER == "none") {
    b_mu <- if (is.na(GAMMA_M_RATE)) runif(1, 0, GAMMA_RATE_MAX) else GAMMA_M_RATE
    th <- list(age    = runif(1, AGE_MIN, AGE_MAX),
               lambda = rgamma(NINTERVALS, GAMMA_L_SHAPE, rate = b_lambda),
               mu     = rgamma(NINTERVALS, GAMMA_M_SHAPE, rate = b_mu),
               psi    = rgamma(NINTERVALS, GAMMA_P_SHAPE, rate = GAMMA_P_RATE),
               gamma_rate_lambda = b_lambda, gamma_rate_mu = b_mu)
  } else {
    # the turnover takes mu's place in the stream, so a batch without one draws as it always did
    age <- runif(1, AGE_MIN, AGE_MAX)
    lambda <- rgamma(NINTERVALS, GAMMA_L_SHAPE, rate = b_lambda)
    turnover <- if (TURNOVER == "beta") rbeta(1, TURN_A, TURN_B) else rlnorm(1, TURN_LMEAN, TURN_LSD)
    th <- list(age    = age,
               lambda = lambda,
               mu     = turnover * lambda,
               psi    = rgamma(NINTERVALS, GAMMA_P_SHAPE, rate = GAMMA_P_RATE),
               gamma_rate_lambda = b_lambda, turnover = turnover)
  }
  # drawn last, and only when asked for, so the draws before them are what they always were
  if (MIX_K > 1) th$alpha <- rexp(1, MIX_ALPHA_RATE)
  if (N_CHARS > 0) th$epsilon <- rbeta(1, EPS_A, EPS_B)
  if (ANAGENESIS) th$lambda_a <- rgamma(1, GAMMA_L_SHAPE, rate = b_lambda)
  th
}

# Anagenetic speciation relabels the same tree rather than reshaping it: split each species at
# Poisson(lambda_a) times in its lifespan and reassign occurrences by time. ORIG is the species
# each new one is a segment of, so what was drawn per species before the split follows it.
add_anagenesis <- function(sim, fossils, lambda_a) {
  n0 <- length(sim$TS)

  # segment boundaries per species, oldest first
  starts <- ends <- vector("list", n0)
  for (k in 1:n0) {
    te <- ifelse(is.na(sim$TE[k]), 0, sim$TE[k])
    ev <- sort(runif(rpois(1, lambda_a * (sim$TS[k] - te)), te, sim$TS[k]),
               decreasing = TRUE)
    starts[[k]] <- c(sim$TS[k], ev)
    ends[[k]] <- c(ev, te)
  }

  # new indices run in (species, segment) order, so a parent's index stays below its daughter's
  n_seg <- sapply(starts, length)
  offset <- cumsum(c(0, n_seg))[1:n0]

  # segment of species k alive at time t
  seg_at <- function(k, t) {
    s <- which(t <= starts[[k]] & t >= ends[[k]])
    if (length(s) == 0) n_seg[k] else s[1]
  }

  N <- sum(n_seg)
  TS <- TE <- rep(NA_real_, N)
  PAR <- rep(NA_integer_, N)
  EXTANT <- rep(FALSE, N)

  for (k in 1:n0) {
    for (s in 1:n_seg[k]) {
      i <- offset[k] + s
      last <- (s == n_seg[k]) && sim$EXTANT[k]
      TS[i] <- starts[[k]][s]
      TE[i] <- ifelse(last, NA_real_, ends[[k]][s])
      EXTANT[i] <- last
      # a later segment descends from the previous one; the first keeps the original
      # parent, joined at whichever of that parent's segments was then alive
      PAR[i] <- if (s > 1) offset[k] + s - 1 else
        if (is.na(sim$PAR[k])) NA_integer_ else
          offset[sim$PAR[k]] + seg_at(sim$PAR[k], sim$TS[k])
    }
  }

  if (nrow(fossils) > 0) {
    k <- as.integer(sub("^t", "", fossils$Species))
    i <- offset[k] + mapply(seg_at, k, fossils$SampT)
    fossils$Species <- paste0("t", i)
    fossils$Extant <- EXTANT[i]
    # retain_occs reads the first row as the oldest occurrence and the last as the youngest
    fossils <- fossils[order(i, -fossils$SampT), , drop = FALSE]
  }

  sim <- list(TE = TE, TS = TS, PAR = PAR, EXTANT = EXTANT)
  class(sim) <- "sim"

  list(SIM = sim, FOSSILS = fossils, ORIG = rep(seq_len(n0), n_seg))
}

# Each species' component, drawn at its birth and independent of its parent's. NULL without a
# mixture, which draws nothing.
draw_components <- function(n) if (MIX_K > 1) sample.int(MIX_K, n, replace = TRUE) else NULL

# Character states of every species, sampled or not, down the true tree. A species keeps its
# parent's states except that each changes with probability epsilon, to a uniformly drawn other
# state; the first species' states are uniform. A parent's index is always below its daughter's.
draw_characters <- function(par, epsilon) {
  n <- length(par)
  x <- matrix(0L, n, N_CHARS)
  for (i in seq_len(n)) {
    if (is.na(par[i])) {
      x[i, ] <- vapply(CHAR_STATES, function(k) sample.int(k, 1) - 1L, 0L)
      next
    }
    x[i, ] <- x[par[i], ]
    for (c in which(runif(N_CHARS) < epsilon)) {
      x[i, c] <- (x[i, c] + sample.int(CHAR_STATES[c] - 1, 1)) %% CHAR_STATES[c]
    }
  }
  x
}

# Expected lineages ever born, 1 + integral of lambda(t)E[N(t)]: E[N] grows as exp(lambda-mu)
# within an interval, so both integrals are closed form. Filtering on this before bd.sim is
# what keeps a runaway draw from hanging, and it is a function of theta alone, so a redraw
# leaves the accepted set a truncation of the prior rather than a conditioning on the data.
expected_lineages <- function(lambda, mu, age, shifts) {
  edges <- c(shifts, age)
  total <- 1
  n <- 1
  for (j in seq_along(lambda)) {
    dt <- edges[j + 1] - edges[j]
    r <- lambda[j] - mu[j]
    total <- total + lambda[j] * n * (if (abs(r) < 1e-12) dt else expm1(r * dt) / r)
    n <- n * exp(r * dt)
    if (!is.finite(total)) return(Inf)
  }
  total
}

# Expected size filters out the astronomical draws, but with rate shifts paleobuddy takes
# its general path, whose cost is not a function of tree size alone, so a rare draw still
# runs for minutes. Its cost is the waiting times it draws, over every tree bd.sim tries
# before one meets nFinal, so those are what is counted. NULL means the call went over
# SIM_MAX_DRAWS and the caller should redraw. The count is a function of the random stream,
# so the same seed refuses the same draws on any machine; sim_over reports how many.
sim_over <- 0
sim_draws_left <- Inf     # Inf outside a bounded call, so nothing else that samples is counted
sim_draws_most <- 0       # the most a call that finished used, which is what to set the ceiling by
over_budget <- structure(class = c("sim_over_budget", "error", "condition"),
                         list(message = "bd.sim went over SIM_MAX_DRAWS", call = NULL))
if (is.finite(SIM_MAX_DRAWS)) {
  invisible(suppressMessages(trace("rexp.var", where = asNamespace("paleobuddy"), print = FALSE,
                         tracer = function() {
                           sim_draws_left <<- sim_draws_left - 1
                           if (sim_draws_left < 0) stop(over_budget)
                         })))
}
# draws refused by the expected-clade filter. This one is a function of theta alone, so it
# truncates the prior and the ranks are non-uniform by exactly that much: the number belongs
# in the record of every batch, not just in a pilot.
sim_filtered <- 0
sim_declining <- 0
sim_drawn <- 0
bd_sim_bounded <- function(n0, lambda, mu, tMax, lShifts, mShifts, nFinal) {
  # paleobuddy refuses shifts alongside a constant rate, and the lone shift a single
  # interval leaves is at the process start, so it is not one
  if (length(lambda) < 2) lShifts <- NULL
  if (length(mu) < 2)     mShifts <- NULL
  call_it <- function() bd.sim(n0, lambda, mu, tMax, lShifts = lShifts,
                               mShifts = mShifts, nFinal = c(nFinal, Inf))
  if (!is.finite(SIM_MAX_DRAWS)) return(call_it())
  sim_draws_left <<- SIM_MAX_DRAWS
  on.exit(sim_draws_left <<- Inf)
  sim <- tryCatch(call_it(), sim_over_budget = function(e) {
    sim_over <<- sim_over + 1
    NULL
  })
  if (!is.null(sim)) sim_draws_most <<- max(sim_draws_most, SIM_MAX_DRAWS - sim_draws_left)
  sim
}

# TreeSim draws the tree backward from a fixed number of extant species, so the origin is an
# output rather than a draw. FossilSim then assigns budding species identity to the complete
# tree: beta = 0 is budding, lambda.a = 0 no anagenesis. Both packages index rates
# youngest-first, and the pipeline carries them oldest-first, so both arguments reverse.
backward_sim <- function(n_extant, lambda, mu, psi, breaks, alpha = NULL, lambda_a = 0) {
  tr <- TreeSim::sim.rateshift.taxa(n = n_extant, numbsim = 1,
                                    lambda = rev(lambda), mu = rev(mu),
                                    frac = rep(1, length(lambda)),
                                    times = c(0, breaks), complete = TRUE)[[1]]
  tx <- FossilSim::sim.taxonomy(tr, beta = 0, lambda.a = lambda_a)

  sp  <- sort(unique(tx$sp))
  TS  <- vapply(sp, function(s) max(tx$start[tx$sp == s]), 0)
  TE  <- vapply(sp, function(s) min(tx$end[tx$sp == s]), 0)
  par <- vapply(sp, function(s) tx$parent[tx$sp == s][1], 0)
  anag <- vapply(sp, function(s) tx$mode[tx$sp == s][1] == "a", FALSE)

  # species 1 has to be the origin lineage: the ancestor walk in summarize.R stops there and
  # the acceptance test asks whether it was sampled
  ord <- order(TS, decreasing = TRUE)
  idx <- integer(max(sp)); idx[sp[ord]] <- seq_along(ord)
  TS <- TS[ord]; TE <- TE[ord]; par <- par[ord]; anag <- anag[ord]

  # the root's parent is 0, and R drops a zero index, so idx[par] would come back short and
  # silently shift every parent after it
  PAR <- rep(NA_integer_, length(par))
  PAR[par != 0] <- idx[par[par != 0]]

  # a tip at the present carries accumulated floating-point error, not an extinction
  extant <- TE < 1e-8
  sim <- list(TS = TS, TE = ifelse(extant, NA_real_, TE), PAR = PAR, EXTANT = extant)

  # the grid must be monotone and reach past the tree whichever is deeper, the oldest break or
  # the origin, or the last interval is malformed
  top <- max(max(TS), if (length(breaks)) max(breaks) else 0) + 1
  comp <- draw_components(length(TS))
  # an anagenetic species continues its ancestor's lineage and keeps its sampling component, as
  # in the forward simulator. Oldest first, so an ancestor is settled before its descendant.
  if (!is.null(comp)) for (i in which(anag)) comp[i] <- comp[PAR[i]]
  if (is.null(comp)) {
    fos <- FossilSim::sim.fossils.intervals(taxonomy = tx, rates = rev(psi),
                                            interval.ages = c(0, breaks, top))
  } else {
    # each component's species sample at its own multiple of psi
    r <- discretized_gamma(alpha, MIX_K)
    fos <- do.call(rbind, lapply(seq_len(MIX_K), function(h) {
      keep <- idx[tx$sp] %in% which(comp == h)
      if (!any(keep)) return(NULL)
      FossilSim::sim.fossils.intervals(taxonomy = tx[keep, ], rates = rev(psi) * r[h],
                                       interval.ages = c(0, breaks, top))
    }))
  }

  # an occurrence cannot predate the tree it was sampled on, and a non-finite age would be
  # binned as the oldest edge with no upper one, which reads as NA downstream
  if ( nrow(fos) > 0 && ( any(!is.finite(fos$hmin)) || max(fos$hmin) > max(TS) + 1e-8 ) ) {
    stop(sprintf(paste("backward_sim: %d of %d occurrence ages are non-finite and %d exceed the",
                       "origin %.4f (max %.4f); interval.ages top %.4f, species %d"),
                 sum(!is.finite(fos$hmin)), nrow(fos),
                 sum(fos$hmin > max(TS) + 1e-8, na.rm = TRUE), max(TS),
                 suppressWarnings(max(fos$hmin[is.finite(fos$hmin)])), top, length(TS)))
  }
  cols <- c("Species", "Extant", "SampT", "MaxT", "MinT")
  if (nrow(fos) == 0) {
    f <- setNames(data.frame(character(0), logical(0), numeric(0), numeric(0), numeric(0),
                             stringsAsFactors = FALSE), cols)
    return(list(SIM = sim, FOSSILS = f, AGE = max(TS), COMP = comp))
  }
  k <- idx[fos$sp]
  # the origin is an output here, so the bins are built out to the tree that was drawn
  bins <- bins_to(max(TS))
  j <- findInterval(fos$hmin, bins)
  f <- data.frame(Species = paste0("t", k), Extant = extant[k], SampT = fos$hmin,
                  MaxT = bins[j + 1], MinT = bins[j], stringsAsFactors = FALSE)
  # retain_occs reads the first row as the oldest occurrence and the last as the youngest
  f <- f[order(k, -f$SampT), , drop = FALSE]
  list(SIM = sim, FOSSILS = f, AGE = max(TS), COMP = comp)
}

# Reporting models, by the sim's model number. An SBC rep emits both from one
# simulated record, so they share a tree, a fossil set and a timeline.
sbc_models <- c(complete = 1, firstlast = 2, presence = 3)

# Which of a species' occurrences get reported, under each model.
retain_occs <- function(occs, model) {
  # complete: the whole record
  if (model == 1) return(occs)

  # first/last (model 2): the oldest and youngest. Keep both even when they fall
  # in the same bin, so a bracketed pair is reported as two (possibly identical)
  # occurrences and the count >= 2 likelihood marginalizes the interior within the
  # bin. Only a genuine single fossil yields count == 1.
  if (model == 2) {
    if (nrow(occs) > 1) return(occs[c(1, nrow(occs)), ])
    return(occs[1, ])
  }

  # presence (model 3): each occupied bin once. The bins are only known once the columns are
  # named, so report_presence below does it.
  occs
}

# Each bin a taxon was sampled in, once. A bin left out then truly held nothing, which is what
# reporting="presence" assumes. On a grid of exact ages every occurrence is its own bin.
report_presence <- function(s) s[!duplicated(s[, c("taxon", "min_age", "max_age")]), , drop = FALSE]

###
# write auxiliary simulation functions

# simulate one rep
simulate_rep <- function(rates, age, shifts,
                         model, unc, extant_singletons,
                         sbc = FALSE, theta = NULL) {
  th <- theta
  ## set parameters
  # number of initial species
  n0 <- 1
  
  # total simulation time
  tMax <- age
  
  # lambda
  lambda <- rates[[1]]
  
  # mu
  mu <- rates[[2]]
  
  # psi
  psi <- rates[[3]]
  
  # shifts
  lShifts <- shifts[[1]]
  mShifts <- shifts[[2]]
  pShifts <- shifts[[3]]
  
  # nFinal based on sbc
  # SBC fix: do NOT impose a diversification floor on the tree
  # (the old value 5 conditions on >=5 total species, which biases the rates).
  nFinal <- ifelse(sbc, 1, 10)
  
  # get bins
  if (sbc) {
    # edges on the fixed grid, built out past the origin
    bins <- bins_to(age)
  } else {
    # if not, make bins
    bins <- make_bins(age, unc)
  }
  
  ##
  # run simulations
  
  # create conditions boolean
  cond <- FALSE

  # Every rejection path redraws theta jointly, so none of them can diverge. Deriving every
  # local from the one draw in one place is what keeps a redraw from leaving part of the old
  # one behind.
  use_theta <- function(t) {
    th <<- t
    age <<- t$age; tMax <<- t$age
    lambda <<- t$lambda; mu <<- t$mu; psi <<- t$psi
    lShifts <<- mShifts <<- pShifts <<- c(0, t$age - rev(cov_breaks))
    bins <<- bins_to(t$age)
  }
  redraw <- function() {
    sim_drawn <<- sim_drawn + 1
    use_theta(draw_theta())
  }

  # run until we fulfill conditions
  while (!cond) {
    if (BACKWARD) {
      # The extant count fixes the standing diversity, not the total: leaving N_EXTANT
      # survivors costs about N_EXTANT / (1 - mu/lambda) species, which diverges as turnover
      # approaches 1 and is unreachable above it, where the process declines. And a tree whose
      # oldest interval declines (mu > lambda) may never find a root, since that interval runs
      # back without end. Both are refused before simulating, as the forward path refuses its own
      # runaways, in one loop so that a draw one of them redraws faces the other again. Both
      # truncate the prior, so the counts belong in the batch record.
      total_species <- function() {
        r <- mean(mu) / mean(lambda)
        if (r >= 1) return(Inf)
        N_EXTANT / (1 - r)
      }
      repeat {
        if (sbc && mu[1] > lambda[1]) {
          sim_declining <<- sim_declining + 1
        } else if (!(sbc && is.finite(MAX_LINEAGES) && total_species() > MAX_LINEAGES)) {
          break
        }
        sim_filtered <<- sim_filtered + 1
        redraw()
      }

      bw <- backward_sim(N_EXTANT, lambda, mu, psi, cov_breaks, alpha = th$alpha,
                         lambda_a = if (ANAGENESIS) th$lambda_a else 0)
      sim <- bw$SIM
      fossils <- bw$FOSSILS
      comp <- bw$COMP
      # the origin is what the tree realized, and it is the truth written for it
      age <- bw$AGE
      tMax <- age
    } else {
    # a draw whose expected tree exceeds the ceiling is refused before it is simulated
    if (sbc && is.finite(MAX_LINEAGES)) {
      while (expected_lineages(lambda, mu, tMax, lShifts) > MAX_LINEAGES) {
        sim_filtered <<- sim_filtered + 1
        redraw()
      }
    }

    # run BD simulation - make sure we get 10+ species
    sim <- bd_sim_bounded(n0, lambda, mu, tMax, lShifts, mShifts, nFinal)

    # a draw that overshot its expectation, or its ceiling on work, is refused before it costs
    # a fossil sampling pass
    if (sbc && (is.null(sim) || length(sim$TS) > MAX_LINEAGES)) {
      redraw()
      next
    }

    # run fossil sampling, each component's species at its own multiple of psi
    sample_at <- function(rate, S = NULL)
      suppressMessages(sample.clade(sim, rate, tMax, S = S,
                                    rShifts = if (length(psi) < 2) NULL else pShifts,
                                    bins = bins, returnAll = TRUE))
    # kept beside sim, whose shape paleobuddy checks
    comp <- draw_components(length(sim$TS))
    if (is.null(comp)) {
      fossils <- sample_at(psi)
    } else {
      r <- discretized_gamma(th$alpha, MIX_K)
      fossils <- do.call(rbind, lapply(seq_len(MIX_K), function(h) {
        S <- which(comp == h)
        if (length(S)) sample_at(psi * r[h], S)
      }))
    }

    # split species at anagenetic events before any of them is reported. A segment keeps the
    # sampling component of the lineage it continues.
    if (ANAGENESIS) {
      ana <- add_anagenesis(sim, fossils, th$lambda_a)
      sim <- ana$SIM
      fossils <- ana$FOSSILS
      if (!is.null(comp)) comp <- comp[ana$ORIG]
    }
    }

    # drop what the boundary hides, and recompute who is alive at it
    if (PRESENT > 0) {
      fossils <- fossils[fossils$SampT > PRESENT, , drop = FALSE]
      alive <- is.na(sim$TE) | sim$TE < PRESENT
      if (nrow(fossils) > 0) {
        fossils$Extant <- alive[as.integer(sub("^t", "", fossils$Species))]
      }
      sim$EXTANT <- alive
    }
    # SBC fix: if no species were sampled, fall through and let the condition
    # below redraw the whole replicate, rather than resampling fossils on the
    # same tree (which can loop forever for tiny trees once the diversification
    # floor is removed). The species loop below uses seq_along to handle 0 safely.

    # an SBC rep reports the same record three ways; an accuracy rep just once
    models <- if (sbc) sbc_models else model
    specs <- rep(list(data.frame(matrix(nrow = 0, ncol = 5))), length(models))
    names(specs) <- names(models)
    ranges <- data.frame(matrix(nrow = 0, ncol = 3))

    # loop through species (seq_along handles the 0-sampled case safely)
    sp_list <- unique(fossils$Species)
    for (i in seq_along(sp_list)) {
      # which species is this
      sp <- sp_list[i]
      
      # get vector of occurrences for that species
      # the true time rides along as a fifth column, for the other grids
      occs <- fossils[fossils$Species == sp, ]
      
      # and get true occurrence times
      true_occs <- fossils[fossils$Species == sp, c("Species", "Extant", "SampT")]
      
      # check if it is extant
      ext_sp <- fossils$Extant[fossils$Species == sp][1]
      
      # report this species' record under each model
      for (m in seq_along(models)) {
        specs[[m]] <- rbind(specs[[m]], retain_occs(occs, models[m]))
      }

      # get true range
      range <- c(max(true_occs$SampT),
                 ifelse(sum(true_occs$Extant) > 0, PRESENT, min(true_occs$SampT)))
      
      # add to ranges
      ranges <- rbind(ranges, c(sp, range))
    }
    
    # reorder and rename columns
    for (m in seq_along(specs)) {
      specs[[m]] <- specs[[m]][, c(1, 5, 4, 2, 3)]
      colnames(specs[[m]]) <- c("taxon", "min_age", "max_age", "status", "age")

      # change status column to extant and extinct
      specs[[m]]$status <- c("extinct", "extant")[specs[[m]]$status + 1]
      if (names(specs)[m] == "presence") specs[[m]] <- report_presence(specs[[m]])
    }

    # An extant lineage is seen at the present with probability rho, and that observation is
    # the status flag rather than an occurrence. A survivor that is not seen has no present-day
    # observation at all: its record stops at its youngest fossil and reads as extinct. One
    # that left no fossil either is absent from the record entirely, below.
    rho_seen <- rep(TRUE, length(sim$TS))
    if (RHO < 1) rho_seen <- runif(length(sim$TS)) < RHO
    if (RHO < 1) {
      for (m in seq_along(specs)) {
        k <- as.integer(sub("^t", "", specs[[m]]$taxon))
        specs[[m]]$status[specs[[m]]$status == "extant" & !rho_seen[k]] <- "extinct"
      }
    }

    # the reporting models keep the same species, so any of them fixes the taxon set
    occurrences <- specs[[1]]
    
    # name columns for ranges
    colnames(ranges) <- c("taxon", "first_age", "last_age")
    
    # conditions (depending if it's an SBC sim or not)
    if (sbc) {
      # SBC fix: accept on >=1 sampled species, matching the analysis'
      # condition="sampling". The previous 5-50 window conditions on the number
      # of sampled species, which no dnFBDRMatrix condition can correct for and
      # which biases the estimated rates (selecting for high-diversification trees).
      # Rejection redraws theta below, so conditioning on a function of the reported record
      # leaves the posterior alone. occurrences holds only fossil-sampled taxa here; the extant
      # singletons are added after this loop.
      ntax <- length(unique(occurrences$taxon))
      if (TAXA_COUNT == "record")
        ntax <- length(union(unique(occurrences$taxon), sprintf("t%d", which(sim$EXTANT & rho_seen))))
      cond <- ntax >= MIN_TAXA && ntax <= MAX_TAXA &&
              length(unique(occurrences$taxon)) >= MIN_FOSSIL_TAXA

      # if cond is false, redraw everything
      if (!cond) redraw()
    }
    else {
      cond <- length(unique(occurrences$taxon)) > 10 &&
        length(unique(occurrences$taxon)) < 500 &&
        length(unique(occurrences$taxon)) / length(sim$TS) < 0.9 &&
        length(unique(occurrences$taxon)) / length(sim$TS) > 0.2
    }
  }

  # A survivor seen at the present that left no fossil still belongs in the record, as a
  # present-day (0,0) tip named t{k} to match its sim lineage index. A lone tip reports
  # identically under every model. Survivors rho missed are simply not here.
  if (sbc) {
    reported <- specs[[1]]$taxon
    for (k in which(sim$EXTANT & rho_seen)) {
      nm <- paste0("t", k)
      if (nm %in% reported) next
      for (m in seq_along(specs)) {
        specs[[m]] <- rbind(specs[[m]], data.frame(taxon = nm, min_age = PRESENT,
                            max_age = PRESENT, status = "extant", age = PRESENT,
                            stringsAsFactors = FALSE))
      }
      ranges <- rbind(ranges, data.frame(taxon = nm, first_age = PRESENT, last_age = PRESENT,
                                         stringsAsFactors = FALSE))
    }
  }

  # the characters of every species, drawn once the record is accepted so that a rejected
  # draw costs none, then kept for the reported ones
  chars <- NULL
  if (sbc && N_CHARS > 0) {
    x <- draw_characters(sim$PAR, th$epsilon)
    k <- as.integer(sub("^t", "", unique(specs[[1]]$taxon)))
    chars <- x[k, , drop = FALSE]
    rownames(chars) <- paste0("t", k)
  }

  # paleobuddy indexes intervals oldest-first because its shifts run forward from the origin.
  # Everything downstream indexes them youngest-first, so truth_of reverses, once.
  true_vals <- if (sbc) truth_of(th, lambda, mu, psi, age) else NULL

  # return sim, ranges and k
  return(list(SIM = sim, SPECIMENS = specs, RANGES = ranges,
              TV = true_vals, BINS = bins, CHARS = chars, COMP = comp))
}

# A standard-data NEXUS matrix, one row per taxon
write_nexus <- function(x, file) {
  symbols <- paste(seq_len(max(CHAR_STATES)) - 1, collapse = "")
  rows <- sprintf("  %s %s", rownames(x), apply(x, 1, paste, collapse = ""))
  writeLines(c("#NEXUS", "BEGIN DATA;",
               sprintf("  DIMENSIONS NTAX=%d NCHAR=%d;", nrow(x), ncol(x)),
               sprintf('  FORMAT DATATYPE=STANDARD SYMBOLS="%s" MISSING=? GAP=-;', symbols),
               "MATRIX", rows, ";", "END;"), file)
}

# simulate one set
simulate_set <- function(n_key, reps, rates, age, base_dir,
                         model = 2, unc = FALSE, extant_singletons = FALSE,
                         diff_bins = 0, sbc = FALSE,
                         first = 1) {
  # outside the SBC sets, rates are prepared as normal
  if (!sbc) {
    # create shifts list
    shifts <- vector("list", 3)
    
    # and number of stages per rate
    n_stages <- c()
    
    # iterate through rates to see how many shifts per rate
    for (i in 1:length(rates)) {
      # check number of stages for this rate
      n_stages <- c(n_stages, length(rates[[i]]))
      
      # make shifts vector
      if (n_stages[i] > 1) { 
        shifts[[i]] <- seq(0, age, age / n_stages[i])[-(n_stages[i] + 1)]
      }
    }
    
    # get bins to bin fossil data on
    bins <- seq(0, age, age / 3)
    if (max(n_stages) == 2) bins <- seq(0, age, age / 2)
    
    # if we're on the case where we need different bins
    if (diff_bins) {
      # set those bins
      bins <- seq(0, age, age / diff_bins)
    }
  }
  
  # create vectors for number of species, sampled, and percentage sampled
  n_sp <- c()
  n_sampled <- c()
  perc_sampled <- c()
  
  # create list for sim, fossil ranges, and k
  sims <- vector("list", reps)
  # An extension simulates replicates first..reps and leaves the earlier ones exactly as they
  # are. They are not redrawn: a batch made before SIM_MAX_DRAWS replaced the wall-clock limit
  # does not reproduce from its seeds, and the earlier replicates are the ones pages rest on.
  if (first > 1) {
    load(paste0(base_dir, "sim_list.RData"))
    if (length(sims) != first - 1)
      stop("sim_list.RData holds ", length(sims), " replicates, not the ", first - 1,
           " the manifest counts")
    length(sims) <- reps
  }

  # reuse existing seeds when present so a re-run reproduces the same trees and
  # rates (only the retention/output changes); otherwise draw and save fresh ones
  seeds_file <- paste0(base_dir, "seeds.RData")
  seeds <- numeric(0)
  if (file.exists(seeds_file)) load(seeds_file)
  if (!sbc && !length(seeds)) {
    # create seeds - reps*100 apart to ensure a lot of possible seeds
    seeds <- runif(reps, (n_key - 1) * reps * 100, n_key * reps * 100)
    save(seeds, file = seeds_file)
  }
  if (sbc) {
    # seeds past first - 1 belong to an extension that did not finish, and are drawn again
    if (first > 1) seeds <- head(seeds, first - 1)
    need <- reps - length(seeds)
    if (need > 0) {
      # Whole numbers, none repeated and none already used: set.seed truncates, so two seeds
      # in one unit interval are one replicate simulated twice.
      repeat {
        fresh <- sample.int(.Machine$integer.max, need)
        if (!any(fresh %in% as.integer(seeds))) break
      }
      seeds <- c(seeds, fresh)
      save(seeds, file = seeds_file)
    }
  }
  
  # one named row per replicate, bound at the end so the names become the columns
  true_vals <- vector("list", reps)
  
  # iterate through reps
  for (rep in first:reps) {
    # print some info
    print(paste0("key: ", n_key, " rep: ", rep, 
                 " seed: ", seeds[rep]))
    
    # set seed 
    set.seed(seeds[rep])
    
    # an SBC set draws its own values. AGE_MIN must exceed the oldest breakpoint so
    # that every replicate has lineages in all intervals.
    if (sbc) {
      th <- draw_theta()
      age <- th$age
      lambda <- th$lambda; mu <- th$mu; psi <- th$psi

      # make rates list
      rates <- list(lambda = lambda,
                    mu = mu,
                    psi = psi)
      
      # shifts list
      shifts <- rep(list(c(0, age - rev(cov_breaks))), 3)
    }
    
    # run sim
    sim_rep <- simulate_rep(rates, age, shifts,
                            model = model, unc = unc,
                            extant_singletons = extant_singletons,
                            sbc = sbc,
                            theta = if (sbc) th else NULL)
    
    # an SBC set shares one timeline and adds to true_vals
    if (sbc) {
      # write timeline. The bins no longer follow the rate intervals, so the
      # timeline is the breakpoints rather than the interior bin edges.
      smart_dir_create(paste0(base_dir, "/times"))
      write.table(t(cov_breaks), 
                  paste0(base_dir, "/times/times_", rep, ".tsv"),
                  col.names = FALSE, row.names = FALSE, 
                  quote = FALSE, sep = "\t")
      
      # add true values to true_vals
      true_vals[[rep]] <- sim_rep$TV
    }
    
    # get sim, occurrences and ranges
    sim <- sim_rep$SIM
    specs <- sim_rep$SPECIMENS
    occurrences <- specs[[1]]
    ranges <- sim_rep$RANGES

    # calculate numbers of interest
    n_sp <- c(n_sp, length(sim$TS))
    n_sampled <- c(n_sampled, length(unique(occurrences$taxon)))
    perc_sampled <- c(perc_sampled, tail(n_sampled, 1) / tail(n_sp, 1))
    
    # save occurrences. An SBC rep reports the same record under every model, so
    # each gets its own directory alongside the one shared timeline and true values.
    for (m in seq_along(specs)) {
      spec_dir <- if (sbc) paste0(base_dir, "/occurrences_", names(specs)[m]) else
                                paste0(base_dir, "/occurrences")
      smart_dir_create(spec_dir)
      write.table(specs[[m]][, 1:4], paste0(spec_dir, "/taxa_", rep, ".tsv"),
                  col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

      # and the same rows on each other grid, as occurrences_<model>_w<width>[_o<offset>]
      if (sbc) for (g in seq_len(nrow(more_grids))) {
        # presence on another grid thins the whole record on that grid, not this one's survivors
        s <- if (names(specs)[m] == "presence") specs[["complete"]] else specs[[m]]
        # a survivor seen at the present is a status, not an occurrence, and has no bin
        occ <- s$min_age != s$max_age
        cell <- rebin(s$age[occ], more_grids$w[g], more_grids$off[g])
        s$min_age[occ] <- cell$min
        s$max_age[occ] <- cell$max
        if (names(specs)[m] == "presence") s <- report_presence(s)
        grid_dir <- paste0(spec_dir, "_", grid_label(more_grids$w[g], more_grids$off[g]))
        smart_dir_create(grid_dir)
        write.table(s[, 1:4], paste0(grid_dir, "/taxa_", rep, ".tsv"),
                    col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
      }
    }
    
    # characters, as NEXUS for readDiscreteCharacterData
    if (!is.null(sim_rep$CHARS)) {
      smart_dir_create(paste0(base_dir, "/characters"))
      write_nexus(sim_rep$CHARS, paste0(base_dir, "/characters/characters_", rep, ".nex"))
    }

    # a species' psi mixture component is a latent of the analysis, so it is a truth too
    if (!is.null(sim_rep$COMP)) {
      smart_dir_create(paste0(base_dir, "/true_components"))
      k <- as.integer(sub("^t", "", ranges$taxon))
      write.table(data.frame(taxon = ranges$taxon, component = sim_rep$COMP[k]),
                  paste0(base_dir, "true_components/components_", rep, ".tsv"),
                  col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
    }

    # save true ranges
    smart_dir_create(paste0(base_dir, "/true_ranges"))
    write.table(ranges, paste0(base_dir, "true_ranges/ranges_", rep, ".tsv"),
                col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
    
    # append to list of simulations
    sims[[rep]] <- sim
  }
  
  # save sim list
  save(sims, file = paste0(base_dir, "sim_list.RData"))
  
  # outside an SBC set
  if (!sbc) {
    # write timeline
    write.table(t(bins[-c(1, length(bins))]), 
                paste0(base_dir, "times.tsv"),
                col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
    
    # get highest number of rates for this set
    max_len <- max(unlist(lapply(1:length(rates), 
                                 function(x) length(rates[[x]]))))
    
    # get make true rates by expanding rates that are not of len max_len
    true_rates <- lapply(1:length(rates), function(x) {
      if (length(rates[[x]]) < max_len) {
        rep(rates[[x]][1], max_len)
      } else {
        rates[[x]]
      }
    })
    
    # name them
    names(true_rates) <- c("lambda", "mu", "psi")
    
    # make it a data frame
    true_rates <- as.data.frame(true_rates)
    
    # write it to file
    write.table(true_rates, paste0(base_dir, "true_vals.tsv"),
                col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
  } else {
    # truth_of named every column when it built the row, so binding them is the whole job
    true_vals <- do.call(rbind, true_vals)
    tv_file <- paste0(base_dir, "/true_vals.tsv")
    # An extension appends its rows, so the earlier ones stay byte for byte what they were.
    # The columns are a function of the config, so a different header means a different one.
    if (first > 1) {
      had <- strsplit(readLines(tv_file, n = 1), "\t")[[1]]
      if (!identical(had, colnames(true_vals)))
        stop("true_vals.tsv has columns ", paste(had, collapse = " "),
             " but these replicates have ", paste(colnames(true_vals), collapse = " "))
    }

    # save true_vals
    write.table(true_vals, tv_file, append = first > 1,
                col.names = first == 1, row.names = FALSE, quote = FALSE, sep = "\t")
  }

  # make numbers into a data frame
  nums <- data.frame(n_sp = n_sp, n_sampled = n_sampled,
                     perc_sampled = perc_sampled)
  
  # return numbers of interest
  return(nums)
}


# create function for simulating an SBC set
simulate_batch <- function(reps, reps_dir, first = 1) {
  # run simulation set
  nums <- simulate_set("cov", reps, NULL, NULL, reps_dir,
                       2, FALSE, FALSE, 0, sbc = TRUE,
                       first = first)
  
  # make nums a data frame
  nums <- as.data.frame(nums)
  
  # name nums data frames
  colnames(nums) <- c("n_sp", "n_sampled", "perc_sampled")
  
  # return the data frames
  return(nums)
}

###
###
# run simulations
#
# One pass reports each simulated record under every reporting model, so this fills
# occurrences_{complete,firstlast,presence}/ alongside the one shared timeline,
# true_vals.tsv, true_ranges/ and sim_list.RData.
#
# seeds.RData is reused when present, so reruns of a batch stay comparable.
#
# Driven entirely by the batch's config.yml, which simfbdr exports; BATCHDIR and NREPS
# come from simfbdr too. NREPS is the count the batch ends with and REPS_HAVE the count it
# already holds, so an extension simulates only the replicates after REPS_HAVE.
# Sourcing with sim_functions_only set loads the functions only.
if (!exists("sim_functions_only")) {
  reps <- as.integer(cfg("NREPS"))
  first <- as.integer(cfg_opt("REPS_HAVE", 0)) + 1L
  reps_dir <- Sys.getenv("BATCHDIR")
  if (reps_dir == "") stop("BATCHDIR is not set; simfbdr must export it")
  if (!grepl("/$", reps_dir)) reps_dir <- paste0(reps_dir, "/")
  smart_dir_create(reps_dir)

  # checked before the hours of simulating, not after: rows past the count would mean an
  # earlier extension wrote its truth and then never reached the manifest
  if (first > 1) {
    had <- length(readLines(paste0(reps_dir, "true_vals.tsv"))) - 1
    if (had != first - 1)
      stop("true_vals.tsv holds ", had, " replicates, not the ", first - 1, " the manifest counts")
  }
  cat("simulating reps", first, "to", reps, "->", reps_dir, "\n")
  nums <- simulate_batch(reps, reps_dir, first = first)

  write.table(nums, paste0(reps_dir, "nums.tsv"), append = first > 1,
              row.names = FALSE, col.names = first == 1, quote = FALSE, sep = "\t")

  # what the draws actually used, for simfbdr to check against the config
  keys <- c("NINTERVALS", "INTERVAL_WIDTH", "AGE_MIN", "AGE_MAX", "BIN_WIDTH", "NREPS",
            "GAMMA_L_SHAPE", "GAMMA_L_RATE", "GAMMA_M_SHAPE", "GAMMA_M_RATE",
            "TURN_A", "TURN_B", "TURN_LMEAN", "TURN_LSD",
            "GAMMA_P_SHAPE", "GAMMA_P_RATE", "GAMMA_RATE_MAX")
  vals <- c(NINTERVALS, INTERVAL_WIDTH, AGE_MIN, AGE_MAX, BIN_WIDTH, reps,
            GAMMA_L_SHAPE, GAMMA_L_RATE, GAMMA_M_SHAPE, GAMMA_M_RATE,
            TURN_A, TURN_B, TURN_LMEAN, TURN_LSD,
            GAMMA_P_SHAPE, GAMMA_P_RATE, GAMMA_RATE_MAX)
  if (MIX_K > 1) {
    keys <- c(keys, "MIX_K", "MIX_ALPHA_RATE"); vals <- c(vals, MIX_K, MIX_ALPHA_RATE)
  }
  if (N_CHARS > 0) {
    keys <- c(keys, "N_CHARS", "EPS_A", "EPS_B"); vals <- c(vals, N_CHARS, EPS_A, EPS_B)
  }
  if (ANAGENESIS) {
    keys <- c(keys, "ANAGENESIS"); vals <- c(vals, "true")
  }
  if (is.finite(MAX_LINEAGES)) {
    keys <- c(keys, "MAX_LINEAGES"); vals <- c(vals, MAX_LINEAGES)
  }
  if (is.finite(SIM_MAX_DRAWS)) {
    keys <- c(keys, "SIM_MAX_DRAWS"); vals <- c(vals, SIM_MAX_DRAWS)
    cat(sprintf("draws refused for going over SIM_MAX_DRAWS: %d (the largest that finished used %d of %d)\n",
                sim_over, sim_draws_most, SIM_MAX_DRAWS))
  }
  if (BACKWARD) cat(sprintf("draws refused for an oldest interval that declines: %d\n", sim_declining))
  if (is.finite(MAX_LINEAGES) || BACKWARD) {
    tot <- sim_drawn + reps - first + 1
    cat(sprintf("draws refused by the prior filters: %d of %d (%.2f%%)\n",
                sim_filtered, tot, 100 * sim_filtered / max(tot, 1)))
    # the batch's count is over every replicate it holds, not only the ones just added
    filter_file <- paste0(reps_dir, "filtered.tsv")
    # a batch begun before this file was called that has its earlier counts under the old name
    was_file <- if (file.exists(filter_file)) filter_file else paste0(reps_dir, "screened.tsv")
    if (first > 1 && file.exists(was_file)) {
      was <- read.table(was_file, sep = "\t", col.names = c("k", "v"))
      sim_filtered <- sim_filtered + was$v[was$k %in% c("filtered", "screened")]
      tot <- tot + was$v[was$k %in% c("filtered_of", "screened_of")]
    }
    writeLines(sprintf("filtered\t%d\nfiltered_of\t%d\nfiltered_frac\t%s",
                       sim_filtered, tot, signif(sim_filtered / max(tot, 1), 4)),
               filter_file)
  }
  writeLines(sprintf("%s\t%s", keys, vals), paste0(reps_dir, "sim_params.tsv"))
  cat("done\n")
}
