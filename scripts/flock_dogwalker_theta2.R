# ============================================================================
# FLOCK illustration: two posterior modes from a correctly specified likelihood
# Base R only.  Run: source("flock_dogwalker_theta2.R")
#
# TRUE MODEL: y_i | theta ~ Normal(theta^2, 1); theta_true = +3
# FITTED MODEL: same likelihood; theta ~ Normal(2, 3^2).
# Because y depends on theta^2, the posterior has modes near -3 and +3.
# An asymmetric prior gives the modes UNEQUAL posterior mass.
#
# Eight independent POSTERIOR runs per method, with the SAME number of target
# density calls in every run.  FLOCK scouts/followers/pioneers and the adaptive
# mixture exploration are NOT treated as posterior draws.  Parallel tempering
# uses ONLY its cold (temperature 1) replica as posterior draws.
#
# Methods: RandomWalkMH, AdaptiveMetropolis, IndependenceMH,
#          ParallelTempering, AdaptiveMixtureMH, SAChES_inspired,
#          FLOCK_Novelty, FLOCK_Unexplored, DogWalker90,
#          DogWalker90_NoRepulsion, DogWalker65.
# SAChES_inspired is a teaching approximation, NOT the published SAChES.
# AdaptiveMixtureMH is a custom baseline, NOT a published AIMM implementation.
#
# Configuration can be overridden BEFORE sourcing:
# Sys.setenv(BUDGET="12000", SEED="20261025", OUTDIR="flock_theta2_outputs")
# ============================================================================

B <- as.integer(Sys.getenv("BUDGET", "12000"))   # TOTAL target evaluations per independent algorithm run
N_CHAINS <- 8L                                # Eight independent final runs
SEED <- as.integer(Sys.getenv("SEED", "20261025"))
OUTDIR <- Sys.getenv("OUTDIR", "flock_theta2_outputs")
stopifnot(is.finite(B), B >= 2000L, is.finite(SEED), nzchar(OUTDIR))
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
set.seed(SEED)

# ----------------------------------------------------------------------------
# 1. Generate data and compute the EXACT 1D posterior by grid integration.
# ----------------------------------------------------------------------------
theta_true <- 3
n <- 12L
y <- rnorm(n, mean = theta_true^2, sd = 1)
prior_mu <- 2
prior_sd <- 3

# The likelihood differs from the full product of normal densities by a
# theta-independent factor. That factor cancels in all MH acceptance ratios.
log_target <- function(theta) {
  -0.5 * n * (theta^2 - mean(y))^2 +
    dnorm(theta, mean = prior_mu, sd = prior_sd, log = TRUE)
}

prior_draw <- function() rnorm(1L, prior_mu, prior_sd)
log_prior <- function(theta) dnorm(theta, prior_mu, prior_sd, log = TRUE)

ref_grid <- seq(-9, 9, length.out = 12001L)
step <- ref_grid[2L] - ref_grid[1L]
ref_log <- log_target(ref_grid)
ref_raw <- exp(ref_log - max(ref_log))
trap <- rep(1, length(ref_grid)); trap[c(1L, length(trap))] <- 0.5
ref_pdf <- ref_raw / sum(ref_raw * trap * step)
ref_cdf <- cumsum(ref_pdf * trap * step)
ref_mean <- sum(ref_grid * ref_pdf * trap * step)
ref_prob_pos <- sum(ref_pdf[ref_grid > 0] * trap[ref_grid > 0] * step)
ref_quantile <- function(p) ref_grid[which(ref_cdf >= p)[1L]]
ref_ci <- vapply(c(0.025, 0.975), ref_quantile, numeric(1))
write.csv(data.frame(y = y), file.path(OUTDIR, "generated_data.csv"), row.names = FALSE)
write.csv(data.frame(theta = ref_grid, exact_posterior_density = ref_pdf),
          file.path(OUTDIR, "exact_posterior.csv"), row.names = FALSE)

cat("\nTRUE theta =", theta_true, "  observed ybar =", round(mean(y), 3), "\n")
cat("Exact posterior: mean =", round(ref_mean, 4),
    "; P(theta > 0 | y) =", round(ref_prob_pos, 4),
    "; 95% CrI = [", round(ref_ci[1], 3), ", ", round(ref_ci[2], 3), "]\n")
cat("Target evaluations per run =", B, "; independent runs per method =", N_CHAINS, "\n")
cat("DogWalker90: 90% exploration + 10% final MH; 10% of final draws discarded.\n\n")

# ----------------------------------------------------------------------------
# 2. Counting posterior evaluations and elementary MH kernels.
# ----------------------------------------------------------------------------
new_counter <- function() {
  e <- new.env(parent = emptyenv())
  e$count <- 0L
  e$eval <- function(x) {
    e$count <- e$count + length(x)
    log_target(x)
  }
  e$get <- function() e$count
  e
}

# Every call to this routine consumes nsteps log-target evaluations; the
# incoming current log target is cached from initialization/previous stage.
explore_rw <- function(f, x, lp, nsteps, step_sd,
                       global_prob = 0, large_step_sd = NULL) {
  stopifnot(nsteps >= 0L, step_sd > 0)
  h <- numeric(nsteps + 1L)
  h[1L] <- x
  accepted <- 0L
  if (nsteps > 0L) for (i in seq_len(nsteps)) {
    move <- if (runif(1) < global_prob) "prior" else "rw"
    if (move == "prior") {
      proposal <- prior_draw()
      lproposal <- f$eval(proposal)
      # Independence MH: q(theta') = prior(theta'); include q(x)/q(theta').
      log_alpha <- lproposal - lp + log_prior(x) - log_prior(proposal)
    } else {
      sd_here <- if (!is.null(large_step_sd) && runif(1) < 0.4)
        large_step_sd else step_sd
      proposal <- x + rnorm(1L, sd = sd_here)
      lproposal <- f$eval(proposal)
      log_alpha <- lproposal - lp  # Symmetric random-walk proposal
    }
    if (log(runif(1)) < log_alpha) {
      x <- proposal; lp <- lproposal; accepted <- accepted + 1L
    }
    h[i + 1L] <- x
  }
  list(x = x, lp = lp, h = h, accept = if (nsteps) accepted/nsteps else NA_real_)
}

retain <- function(h, fraction) {
  discard <- floor(length(h) * fraction)
  stopifnot(discard < length(h))
  as.numeric(h[seq.int(discard + 1L, length(h))])
}

# A mixture of *separately valid* MH kernels: local symmetric random walk and
# independence MH with learned global mixture q.  q is FROZEN throughout.
final_global_mh <- function(f, x, lp, nsteps, q,
                            local_probability = 0.20, local_sd = 0.15) {
  h <- numeric(nsteps + 1L); h[1L] <- x
  lqx <- q$logq(x)
  acc <- 0L
  if (nsteps > 0L) for (i in seq_len(nsteps)) {
    use_local <- runif(1) < local_probability
    if (use_local) {
      z <- x + rnorm(1L, sd = local_sd)
      lz <- f$eval(z)
      log_alpha <- lz - lp
      lqz <- q$logq(z)
    } else {
      z <- q$draw()
      lz <- f$eval(z)
      lqz <- q$logq(z)
      log_alpha <- lz - lp + lqx - lqz
    }
    if (log(runif(1)) < log_alpha) {
      x <- z; lp <- lz; lqx <- lqz; acc <- acc + 1L
    }
    h[i + 1L] <- x
  }
  list(h = h, accept = if (nsteps) acc/nsteps else NA_real_)
}

# Learned 1D Gaussian-mixture independence proposal.  Equal regional proposal
# weights DO NOT encode posterior mass; the MH correction handles that.
# A broad normal component guarantees q(x) > 0 for all real x.
make_global_proposal <- function(histories) {
  pool <- unlist(lapply(histories, function(h) {
    idx <- unique(as.integer(round(seq(1, length(h), length.out =
                                         min(120L, length(h))))))
    h[idx]
  }), use.names = FALSE)
  pool <- pool[is.finite(pool)]
  if (length(pool) < 4L) stop("Too few exploratory draws")

  if (sd(pool) < 0.10 || length(unique(pool)) < 4L) {
    centers <- mean(pool)
    sds <- max(0.12, sd(pool))
  } else {
    fit <- tryCatch(kmeans(matrix(pool, ncol = 1L), centers = 2L,
                           nstart = 5L, iter.max = 60L),
                    error = function(e) NULL)
    if (is.null(fit)) {
      centers <- mean(pool); sds <- max(0.12, sd(pool))
    } else {
      centers <- as.numeric(fit$centers)
      sds <- vapply(seq_along(centers), function(j) {
        z <- pool[fit$cluster == j]
        max(0.12, if (length(z) > 1L) sd(z) else 0.12)
      }, numeric(1))
    }
  }

  # A broad component derived from the known prior, not true mode locations.
  centers <- c(centers, prior_mu)
  sds <- c(sds, 2 * prior_sd)
  k <- length(centers) - 1L
  probs <- c(rep(0.80/k, k), 0.20)
  logq <- function(z) {
    parts <- log(probs) + dnorm(z, mean = centers, sd = sds, log = TRUE)
    a <- max(parts)
    a + log(sum(exp(parts - a)))
  }
  draw <- function() {
    j <- sample.int(length(probs), 1L, prob = probs)
    rnorm(1L, centers[j], sds[j])
  }
  list(draw = draw, logq = logq, centers = centers,
       sds = sds, proposal_weights = probs)
}

# ----------------------------------------------------------------------------
# 3. Single-run implementations.  Each uses exactly B target evaluations.
# ----------------------------------------------------------------------------
random_walk_mh <- function(B, x0) {
  f <- new_counter(); lp <- f$eval(x0)
  a <- explore_rw(f, x0, lp, B - 1L, step_sd = 0.35)
  stopifnot(f$get() == B)
  list(draws = retain(a$h, 0.25), full = a$h, evals = f$get(),
       accept = a$accept, stages = list())
}

adaptive_metropolis <- function(B, x0) {
  f <- new_counter(); x <- x0; lp <- f$eval(x)
  h <- numeric(B); h[1L] <- x
  burn <- floor(0.25 * B); scale <- 0.35; acc <- 0L
  for (t in 2:B) {
    if (t <= burn && t >= 101L && t %% 100L == 0L) {
      # 1D AM covariance estimated only from warmup, then frozen.
      scale <- max(0.09, min(1.5, 2.38 * sd(h[1L:(t - 1L)])))
    }
    z <- x + rnorm(1L, sd = scale)
    lz <- f$eval(z)
    if (log(runif(1)) < lz - lp) {x <- z; lp <- lz; acc <- acc + 1L}
    h[t] <- x
  }
  stopifnot(f$get() == B)
  list(draws = h[(burn + 1L):B], full = h,
       evals = f$get(), accept = acc/(B - 1L), stages = list())
}

independence_mh <- function(B, x0) {
  f <- new_counter(); x <- x0; lp <- f$eval(x)
  h <- numeric(B); h[1L] <- x; acc <- 0L
  for (t in 2:B) {
    z <- prior_draw(); lz <- f$eval(z)
    log_alpha <- lz - lp + log_prior(x) - log_prior(z)
    if (log(runif(1)) < log_alpha) {x <- z; lp <- lz; acc <- acc + 1L}
    h[t] <- x
  }
  stopifnot(f$get() == B)
  list(draws = retain(h, 0.25), full = h,
       evals = f$get(), accept = acc/(B - 1L), stages = list())
}

parallel_tempering <- function(B, x0) {
  f <- new_counter()
  temps <- c(1, 2.5, 6, 15)
  beta <- 1/temps
  nt <- length(temps)
  x <- c(x0, replicate(nt - 1L, prior_draw()))
  lp <- f$eval(x)
  n_sweeps <- floor((B - nt)/nt)
  cold <- numeric(n_sweeps + 1L)
  cold[1L] <- x[1L]
  cold_acc <- 0L; swap_acc <- 0L
  for (i in seq_len(n_sweeps)) {
    for (j in seq_len(nt)) {
      z <- x[j] + rnorm(1L, sd = 0.35 * sqrt(temps[j]))
      lz <- f$eval(z)
      if (log(runif(1)) < beta[j]*(lz - lp[j])) {
        x[j] <- z; lp[j] <- lz
        if (j == 1L) cold_acc <- cold_acc + 1L
      }
    }
    # Swap neighboring temperatures. NO new log-target evaluation needed.
    j <- sample.int(nt - 1L, 1L); k <- j + 1L
    log_alpha <- (beta[j] - beta[k]) * (lp[k] - lp[j])
    if (log(runif(1)) < log_alpha) {
      tmp <- x[j]; x[j] <- x[k]; x[k] <- tmp
      tmp <- lp[j]; lp[j] <- lp[k]; lp[k] <- tmp
      swap_acc <- swap_acc + 1L
    }
    cold[i + 1L] <- x[1L]
  }
  # Any remainder calls are spent on the cold chain, and retained.
  rest <- B - f$get()
  if (rest > 0L) {
    a <- explore_rw(f, x[1L], lp[1L], rest, step_sd = 0.35)
    cold <- c(cold, a$h[-1L])
  }
  stopifnot(f$get() == B)
  list(draws = retain(cold, 0.25), full = cold,
       evals = f$get(), accept = cold_acc/n_sweeps,
       swap_accept = swap_acc/n_sweeps, stages = list())
}

adaptive_mixture_mh <- function(B, x0) {
  f <- new_counter()
  exploration_budget <- floor(0.50 * B)
  n_workers <- 8L
  allot <- rep((exploration_budget - n_workers) %/% n_workers, n_workers)
  allot[1L] <- allot[1L] + exploration_budget - n_workers - sum(allot)
  histories <- vector("list", n_workers)
  endpoints <- numeric(n_workers)
  endpoint_lp <- numeric(n_workers)
  for (j in seq_len(n_workers)) {
    start <- if (j == 1L) x0 else prior_draw()
    lp <- f$eval(start)
    z <- explore_rw(f, start, lp, allot[j], step_sd = 0.38,
                    global_prob = 0.18)
    histories[[j]] <- z$h
    endpoints[j] <- z$x; endpoint_lp[j] <- z$lp
  }
  stopifnot(f$get() == exploration_budget)
  q <- make_global_proposal(histories)
  j <- which.max(endpoint_lp)
  fin <- final_global_mh(f, endpoints[j], endpoint_lp[j], B - f$get(), q)
  stopifnot(f$get() == B)
  list(draws = retain(fin$h, 0.10), full = fin$h,
       evals = f$get(), accept = fin$accept,
       stages = list(exploration = histories), proposal = q)
}

saches_inspired <- function(B, x0) {
  # Teaching approximation ONLY: DE ensemble exploration followed by a
  # covariance-tuned, LOCAL random-walk final sampler.  Not published SAChES.
  f <- new_counter(); K <- 8L
  exploration_budget <- floor(0.50*B)
  x <- c(x0, replicate(K - 1L, prior_draw()))
  lp <- f$eval(x)
  histories <- lapply(x, function(z) z)
  step_count <- exploration_budget - K
  for (i in seq_len(step_count)) {
    j <- (i - 1L) %% K + 1L
    others <- setdiff(seq_len(K), j)
    pair <- sample(others, 2L, replace = FALSE)
    z <- x[j] + 0.9 * (x[pair[1L]] - x[pair[2L]]) + rnorm(1L, sd = 0.18)
    lz <- f$eval(z)
    # Conditional on other walkers, signed DE differences are symmetric.
    if (log(runif(1)) < lz - lp[j]) {x[j] <- z; lp[j] <- lz}
    histories[[j]] <- c(histories[[j]], x[j])
  }
  stopifnot(f$get() == exploration_budget)
  local_sds <- vapply(histories, function(z) if(length(z)>2L) sd(z) else NA_real_, numeric(1))
  local_sds <- local_sds[is.finite(local_sds)]
  scale <- if (length(local_sds)) max(0.12, min(0.65, 2.38*median(local_sds))) else 0.35
  j <- which.max(lp)
  fin <- explore_rw(f, x[j], lp[j], B - f$get(), step_sd = scale)
  stopifnot(f$get() == B)
  list(draws = retain(fin$h, 0.25), full = fin$h,
       evals = f$get(), accept = fin$accept,
       stages = list(de_exploration = histories))
}

flock <- function(B, x0, variant = c("Novelty", "Unexplored")) {
  variant <- match.arg(variant)
  f <- new_counter()
  scout_budget <- floor(0.18 * B)
  mid_budget <- floor(0.32 * B)

  # STAGE 1: Two mixed-scale aggressive scouts; their draws inform proposals
  # but are NEVER part of the final posterior sample.
  sstarts <- c(x0, prior_draw())
  ns <- rep((scout_budget - 2L) %/% 2L, 2L)
  ns[1L] <- ns[1L] + scout_budget - 2L - sum(ns)
  scouts <- vector("list", 2L)
  sx <- slp <- numeric(2L)
  for (j in seq_len(2L)) {
    lp <- f$eval(sstarts[j])
    z <- explore_rw(f, sstarts[j], lp, ns[j], step_sd = 0.65,
                    global_prob = 0.20, large_step_sd = 2.8)
    scouts[[j]] <- z$h; sx[j] <- z$x; slp[j] <- z$lp
  }
  stopifnot(f$get() == scout_budget)

  # STAGE 2: Three followers and three novelty/unexplored pioneers.
  allot <- rep(mid_budget %/% 6L, 6L)
  allot[1L] <- allot[1L] + mid_budget - sum(allot)
  middle <- vector("list", 6L)
  ranks <- order(slp, decreasing = TRUE)
  for (j in 1:3) {
    scout_id <- ranks[(j - 1L) %% 2L + 1L]
    x <- sx[scout_id]; lp <- slp[scout_id]; steps <- allot[j]
    if (j == 3L) {
      x <- x + rnorm(1L, sd = 0.20)
      lp <- f$eval(x)
      steps <- steps - 1L
    }
    z <- explore_rw(f, x, lp, steps, step_sd = 0.22)
    middle[[j]] <- z$h
  }

  scout_pool <- unlist(lapply(scouts, function(h)
    h[unique(as.integer(round(seq(1, length(h), length.out =
                                  min(200L, length(h))))))]),
    use.names = FALSE)
  map_scale <- max(0.50, sd(scout_pool))
  for (j in 4:6) {
    # Explicitly account for all candidate scoring as target evaluations.
    nc <- 8L
    cand <- replicate(nc, prior_mu + rnorm(1L, sd = 2 * prior_sd))
    clp <- f$eval(cand)
    novelty <- vapply(cand, function(c) min(abs(c - scout_pool))/map_scale,
                      numeric(1))
    if (variant == "Novelty") {
      utility <- novelty - 0.30 * pmax(0, max(clp) - clp)
      selected <- which.max(utility)
    } else {
      eligible <- which(novelty >= median(novelty))
      selected <- eligible[sample.int(length(eligible), 1L)]
    }
    z <- explore_rw(f, cand[selected], clp[selected], allot[j] - nc,
                    step_sd = 0.45)
    middle[[j]] <- z$h
  }
  stopifnot(f$get() == scout_budget + mid_budget)

  # STAGE 3: Freeze the learned global proposal.  This is the ONLY stage
  # contributing posterior draws to FLOCK results.
  q <- make_global_proposal(c(scouts, middle))
  j <- which.max(slp)
  fin <- final_global_mh(f, sx[j], slp[j], B - f$get(), q)
  stopifnot(f$get() == B)
  list(draws = retain(fin$h, 0.10), full = fin$h,
       evals = f$get(), accept = fin$accept,
       stages = list(scouts = scouts, middle = middle), proposal = q)
}

# ----------------------------------------------------------------------------
# DOG WALKER: interacting exploratory walkers, then a single final MH sampler.
# In the exploration stage, 7 walkers deliberately REPEL each other.
# Their exploratory ensemble has a different (repulsive) joint distribution;
# THESE DRAWS ARE NEVER USED FOR POSTERIOR INFERENCE.
# In a coordinate update, conditioned on other walkers, the log joint target is
#    log pi(x_i) + strength * sum_{j != i} log{1-exp[-(x_i-x_j)^2/(2h^2)]}.
# This favors occupied-space separation, but does not guarantee mode coverage.
# The last chain uses a FROZEN learned mixture proposal, with exact MH ratio.
# ----------------------------------------------------------------------------
repulsion_log <- function(z, others, radius = 0.18) {
  if (!length(others)) return(0)
  u <- (z - others)/radius
  # -expm1(-a) = 1-exp(-a), stable when 'a' is small. Clamp to avoid -Inf.
  sum(log(pmax(-expm1(-0.5*u^2), 1e-12)))
}

dog_walker <- function(B, x0, explore_frac = 0.90, repel_strength = 0.35,
                       radius = 0.18, n_walkers = 7L) {
  stopifnot(B >= 2000L, explore_frac > 0, explore_frac < 1,
            repel_strength >= 0, radius > 0, n_walkers >= 2L)
  f <- new_counter()
  expl_budget <- floor(B * explore_frac)
  x <- c(x0, replicate(n_walkers - 1L, prior_draw()))
  lp <- f$eval(x)
  histories <- lapply(x, function(z) z)
  acc <- 0L
  close_before <- numeric()
  nsteps <- expl_budget - n_walkers

  for (i in seq_len(nsteps)) {
    j <- (i - 1L) %% n_walkers + 1L
    others <- x[-j]
    # Mixture of symmetric local and large random-walk proposals.
    jump_sd <- if (runif(1) < 0.30) 2.5 else 0.36
    z <- x[j] + rnorm(1L, sd = jump_sd)
    lz <- f$eval(z)
    loga <- lz - lp[j]
    if (repel_strength > 0) {
      loga <- loga + repel_strength *
        (repulsion_log(z, others, radius) -
         repulsion_log(x[j], others, radius))
    }
    if (log(runif(1)) < loga) {
      x[j] <- z; lp[j] <- lz; acc <- acc + 1L
    }
    histories[[j]] <- c(histories[[j]], x[j])
    if (i %% (n_walkers * 20L) == 0L) {
      close_before <- c(close_before, min(dist(matrix(x,ncol=1L))))
    }
  }
  stopifnot(f$get() == expl_budget)
  q <- make_global_proposal(histories)
  best <- which.max(lp)
  fin <- final_global_mh(f, x[best], lp[best], B - f$get(), q)
  stopifnot(f$get() == B)
  list(draws = retain(fin$h, 0.10), full = fin$h,
       evals = f$get(), accept = fin$accept,
       exploration_acceptance = acc/nsteps,
       stages = list(walkers = histories, nearest_walker_distance = close_before),
       proposal = q,
       explorer_fraction = explore_frac,
       repulsion_strength = repel_strength)
}

# ----------------------------------------------------------------------------
# 4. Run eight separate final posterior chains for EACH method.
#    Under a given method, the 8 runs are independent conditional on the data.
# ----------------------------------------------------------------------------
methods <- c("RandomWalkMH", "AdaptiveMetropolis", "IndependenceMH",
             "ParallelTempering", "AdaptiveMixtureMH", "SAChES_inspired",
             "FLOCK_Novelty", "FLOCK_Unexplored",
             "DogWalker90", "DogWalker90_NoRepulsion", "DogWalker65")

run_method <- function(method, B, x0) switch(method,
  RandomWalkMH = random_walk_mh(B, x0),
  AdaptiveMetropolis = adaptive_metropolis(B, x0),
  IndependenceMH = independence_mh(B, x0),
  ParallelTempering = parallel_tempering(B, x0),
  AdaptiveMixtureMH = adaptive_mixture_mh(B, x0),
  SAChES_inspired = saches_inspired(B, x0),
  FLOCK_Novelty = flock(B, x0, "Novelty"),
  FLOCK_Unexplored = flock(B, x0, "Unexplored"),
  DogWalker90 = dog_walker(B, x0, explore_frac = .90, repel_strength = .35),
  DogWalker90_NoRepulsion = dog_walker(B, x0, explore_frac = .90, repel_strength = 0),
  DogWalker65 = dog_walker(B, x0, explore_frac = .65, repel_strength = .35),
  stop("Unknown method: ", method))

# Deliberately dispersed INITIAL STATES to show local chains getting trapped.
# These are not posterior draws and their allocation does not reflect mode mass.
starts <- c(-5, -4, -2, -0.5, 0.5, 2, 4, 5) + rnorm(N_CHAINS, sd = 0.12)
results <- setNames(vector("list", length(methods)), methods)
for (m in seq_along(methods)) {
  name <- methods[m]
  cat("Running", name, "\n")
  runs <- vector("list", N_CHAINS)
  for (i in seq_len(N_CHAINS)) {
    set.seed(SEED + 10000L*m + i)
    t0 <- proc.time()[["elapsed"]]
    a <- run_method(name, B, starts[i])
    a$elapsed_seconds <- as.numeric(proc.time()[["elapsed"]] - t0)
    stopifnot(a$evals == B, length(a$draws) > 100L,
              all(is.finite(a$draws)))
    runs[[i]] <- a
    cat("  chain", i, ": retained", length(a$draws),
        "; P(theta>0) =", sprintf("%.3f", mean(a$draws > 0)), "\n")
  }
  results[[name]] <- runs
  # Crash-resistant method-level checkpoints; not a replacement for full runs.
  saveRDS(runs, file.path(OUTDIR, paste0(name, "_runs.rds")))
}
saveRDS(list(data = y, starts = starts, methods = results,
             truth = list(theta_true = theta_true, grid = ref_grid,
                          density = ref_pdf, mean = ref_mean,
                          p_positive = ref_prob_pos, ci = ref_ci)),
        file.path(OUTDIR, "all_runs.rds"))

# ----------------------------------------------------------------------------
# 5. Per-chain posterior uncertainty, mode occupancy and movement diagnostics.
# ----------------------------------------------------------------------------
mode_switches <- function(x) sum(diff(as.integer(x > 0)) != 0)
chain_rows <- do.call(rbind, lapply(methods, function(m) {
  do.call(rbind, lapply(seq_len(N_CHAINS), function(i) {
    a <- results[[m]][[i]]
    x <- a$draws
    data.frame(method = m, chain = i, evaluations = a$evals,
               elapsed_seconds = a$elapsed_seconds, retained_draws = length(x), mean = mean(x), sd = sd(x),
               lower95 = unname(quantile(x, 0.025)),
               upper95 = unname(quantile(x, 0.975)),
               prob_positive = mean(x > 0), mode_switches = mode_switches(x),
               acceptance = a$accept,
               abs_mean_error = abs(mean(x) - ref_mean),
               abs_prob_error = abs(mean(x > 0) - ref_prob_pos))
  }))
}))
write.csv(chain_rows, file.path(OUTDIR, "chain_posterior_summaries.csv"),
          row.names = FALSE)

# Classic R-hat (not rank-normalized); meant only to illustrate disagreement
# between eight independent final chains, not certify convergence.
basic_rhat <- function(chains) {
  k <- length(chains); nmin <- min(lengths(chains))
  x <- sapply(chains, function(z) tail(z, nmin))
  W <- mean(apply(x, 2L, var))
  if (!is.finite(W) || W < 1e-15) return(NA_real_)
  Bbetween <- nmin * var(colMeans(x))
  sqrt(((nmin - 1)/nmin * W + Bbetween/nmin)/W)
}

method_rows <- do.call(rbind, lapply(methods, function(m) {
  runs <- results[[m]]
  chains <- lapply(runs, `[[`, "draws")
  xx <- unlist(chains, use.names = FALSE)
  data.frame(method = m, total_target_evaluations = B*N_CHAINS,
             total_elapsed_seconds = sum(vapply(runs, `[[`, numeric(1), "elapsed_seconds")),
             total_retained_draws = length(xx),
             pooled_mean = mean(xx), pooled_sd = sd(xx),
             pooled_lower95 = unname(quantile(xx, 0.025)),
             pooled_upper95 = unname(quantile(xx, 0.975)),
             pooled_prob_positive = mean(xx > 0),
             mean_abs_error = abs(mean(xx) - ref_mean),
             prob_abs_error = abs(mean(xx > 0) - ref_prob_pos),
             mean_mode_switches = mean(vapply(chains, mode_switches, numeric(1))),
             classical_Rhat = basic_rhat(chains),
             mean_acceptance = mean(vapply(runs, `[[`, numeric(1), "accept")))
}))
write.csv(method_rows, file.path(OUTDIR, "method_posterior_summaries.csv"),
          row.names = FALSE)
write.csv(data.frame(theta_true = theta_true, exact_mean = ref_mean,
                     exact_p_positive = ref_prob_pos,
                     exact_lower95 = ref_ci[1], exact_upper95 = ref_ci[2]),
          file.path(OUTDIR, "reference_summary.csv"), row.names = FALSE)

# ----------------------------------------------------------------------------
# 6. Visuals.  For every method, each of 8 panels shows that run's FINAL
#    retained posterior density (blue) against the exact posterior (black).
#    Hot PT replicas and all exploratory trajectories are EXCLUDED.
# ----------------------------------------------------------------------------
plot_chain_density_pages <- function(file) {
  pdf(file, width = 13, height = 8.5, onefile = TRUE)
  on.exit(dev.off(), add = TRUE)
  for (m in methods) {
    par(mfrow = c(2, 4), mar = c(3.5, 3.5, 3, 1), oma = c(0, 0, 3, 0))
    for (i in seq_len(N_CHAINS)) {
      z <- results[[m]][[i]]$draws
      kde <- density(z, bw = 0.04, from = -5, to = 5, n = 1000)
      yrange <- c(0, 1.08 * max(ref_pdf, kde$y))
      plot(ref_grid, ref_pdf, type = "l", col = "black", lwd = 2,
           xlim = c(-5, 5), ylim = yrange,
           xlab = expression(theta), ylab = "Posterior density",
           main = paste0("Chain ", i, " | P(+)=", sprintf("%.2f", mean(z > 0))))
      lines(kde$x, kde$y, col = "#2468A0", lwd = 2)
      abline(v = theta_true, col = "grey55", lty = 3)
      abline(v = -theta_true, col = "grey70", lty = 3)
      if (i == 1L) legend("topright", c("Exact posterior", "Chain KDE"),
                          col = c("black", "#2468A0"), lwd = 2, cex = 0.7, bty = "n")
    }
    mtext(paste0(m, " | 8 independent retained posterior chains"),
          outer = TRUE, line = 0.9, cex = 1.3, font = 2)
  }
}

plot_trace_pages <- function(file) {
  pdf(file, width = 13, height = 8.5, onefile = TRUE)
  on.exit(dev.off(), add = TRUE)
  for (m in methods) {
    par(mfrow = c(2, 4), mar = c(3.5, 3.5, 2.5, 1), oma = c(0, 0, 3, 0))
    for (i in seq_len(N_CHAINS)) {
      z <- results[[m]][[i]]$draws
      plot(z, type = "l", col = "#2468A0", lwd = 0.7,
           ylim = c(-5, 5), xlab = "Retained iteration", ylab = expression(theta),
           main = paste0("Chain ", i, " | switches=", mode_switches(z)))
      abline(h = c(-3, 3), lty = 3, col = "grey50")
    }
    mtext(paste0(m, " | retained posterior traces (after burn-in)"),
          outer = TRUE, line = 0.9, cex = 1.3, font = 2)
  }
}

plot_pooled_comparison <- function(file) {
  pdf(file, width = 12, height = 8, onefile = TRUE)
  on.exit(dev.off(), add = TRUE)
  par(mfrow = c(2, 4), mar = c(4, 4, 2.5, 1))
  for (m in methods) {
    z <- unlist(lapply(results[[m]], `[[`, "draws"), use.names = FALSE)
    kde <- density(z, bw = .04, from = -5, to = 5, n = 1000)
    plot(ref_grid, ref_pdf, type = "l", col = "black", lwd = 2,
         xlim = c(-5, 5), ylim = c(0, 1.08 * max(ref_pdf, kde$y)),
         xlab = expression(theta), ylab = "Density", main = m)
    lines(kde$x, kde$y, col = "#2468A0", lwd = 2)
    mtext(sprintf("P(+)=%.2f; exact=%.2f", mean(z > 0), ref_prob_pos),
          side = 3, line = 0.2, cex = .75)
  }
  # Pooled overlays are only valid posterior estimates if chains adequately
  # sample the SAME target; inspect chain-to-chain heterogeneity first.
}

plot_flock_stages <- function(file) {
  pdf(file, width = 12, height = 8, onefile = TRUE)
  on.exit(dev.off(), add = TRUE)
  for (m in c("FLOCK_Novelty", "FLOCK_Unexplored")) {
    a <- results[[m]][[1L]]
    par(mfrow = c(3, 1), mar = c(4, 4, 3, 1))
    for (stage in c("scouts", "middle")) {
      histories <- a$stages[[stage]]
      plot(NA, xlim = c(1, max(lengths(histories))), ylim = c(-6, 6),
           xlab = "Exploratory step", ylab = expression(theta),
           main = paste(m, stage, "(NOT posterior draws)"))
      for (j in seq_along(histories))
        lines(histories[[j]], col = grDevices::rainbow(length(histories))[j])
      abline(h = c(-3, 3), lty = 3)
    }
    plot(a$draws, type = "l", col = "#2468A0", ylim = c(-6, 6),
         xlab = "Retained step", ylab = expression(theta),
         main = paste(m, "final posterior chain ONLY (after burn-in)"))
    abline(h = c(-3, 3), lty = 3)
  }
}

# Report which modes were encountered in exploration vs the retained final sample.
# Discovery means a visit to theta<0 / >0, not calibrated mass or convergence.
exploration_rows <- do.call(rbind, lapply(methods, function(m) {
  do.call(rbind, lapply(seq_len(N_CHAINS), function(i) {
    a <- results[[m]][[i]]
    # Ignore the nearest-walker-distance diagnostic in discovery counting.
    stage_hist <- if (m %in% c("DogWalker90", "DogWalker65", "DogWalker90_NoRepulsion"))
      unlist(a$stages$walkers, use.names = FALSE)
    else if (m == "AdaptiveMixtureMH") unlist(a$stages$exploration, use.names = FALSE)
    else if (m == "SAChES_inspired") unlist(a$stages$de_exploration, use.names = FALSE)
    else if (grepl("^FLOCK_", m)) unlist(c(a$stages$scouts,a$stages$middle),use.names=FALSE)
    else numeric()
    data.frame(method=m,chain=i,
               exploration_negative=if(length(stage_hist)) any(stage_hist<0) else NA,
               exploration_positive=if(length(stage_hist)) any(stage_hist>0) else NA,
               final_negative=any(a$draws<0),final_positive=any(a$draws>0))
  }))
}))
write.csv(exploration_rows,file.path(OUTDIR,"stage_mode_coverage.csv"),row.names=FALSE)

plot_dog_walker_stages <- function(file) {
  pdf(file, width=12,height=8,onefile=TRUE)
  on.exit(dev.off(),add=TRUE)
  for (m in c("DogWalker90", "DogWalker90_NoRepulsion", "DogWalker65")) {
    a<-results[[m]][[1L]]
    h<-a$stages$walkers
    par(mfrow=c(2,1),mar=c(4,4,3,1))
    plot(NA,xlim=c(1,max(lengths(h))),ylim=c(-6,6),
         xlab="Exploration updates per walker",ylab=expression(theta),
         main=paste(m,"| 7 interacting exploratory trajectories (NOT posterior draws)"))
    cols<-grDevices::rainbow(length(h))
    for(j in seq_along(h)) lines(h[[j]],col=cols[j],lwd=.7)
    abline(h=c(-3,3),lty=3)
    plot(a$draws,type="l",col="#2468A0",ylim=c(-6,6),
         xlab="Retained iteration",ylab=expression(theta),
         main=paste(m,"| final corrected MH chain (post burn-in)"))
    abline(h=c(-3,3),lty=3)
  }
}

plot_chain_density_pages(file.path(OUTDIR, "eight_chains_vs_exact.pdf"))
plot_trace_pages(file.path(OUTDIR, "eight_chain_traces.pdf"))
plot_pooled_comparison(file.path(OUTDIR, "pooled_vs_exact_CAUTION.pdf"))
plot_flock_stages(file.path(OUTDIR, "flock_exploration_vs_final.pdf"))
plot_dog_walker_stages(file.path(OUTDIR, "dog_walker_exploration_vs_final.pdf"))

cat("\n======= Method summaries (pooled estimates require adequate mixing) =======\n")
print(method_rows[, c("method", "pooled_mean", "pooled_prob_positive",
                      "mean_abs_error", "prob_abs_error", "classical_Rhat")],
      row.names = FALSE, digits = 3)
cat("\nExact mean:", round(ref_mean, 4),
    " | exact P(theta > 0):", round(ref_prob_pos, 4), "\n")
cat("\nSaved outputs in:", normalizePath(OUTDIR), "\n")
cat("IMPORTANT: burn-in alone does NOT guarantee mixing or convergence.\n")
cat("The repelled exploratory walker trajectories DO NOT target the original posterior.\n")