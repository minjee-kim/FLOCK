# FLOCK EXPERIMENT 5 -- version 1.0.0
# Standalone, base R only; no other FLOCK files need to be sourced.
#
# source("flock_experiment_5.R")                   # default: smoke, 3 x 5,000
# Sys.setenv(RUN_MODE="pilot"); source(...)       # 10 x 20,000
# Sys.setenv(RUN_MODE="full"); source(...)        # 100 x 30,000
# Overrides: REPS, BUDGET, OUTDIR, METHODS (comma-separated), BASE_SEED,
#            RESUME=TRUE|FALSE, SAVE_DRAWS=TRUE|FALSE, BRIDGE_CHECK=TRUE|FALSE.
# To load functions without running: Sys.setenv(FLOCK_AUTORUN="FALSE").
#
# Default methods: FLOCK_Novelty, FLOCK_Unexplored, FLOCK_Map,
# FLOCK_MapMerge, AdaptiveMixtureMH, ParallelTempering, SAChES_inspired,
# RWMH_4chains. No method is given mixture parameters or reference labels.
#
# IMPORTANT SCIENTIFIC QUALIFICATIONS
# * This is a research prototype, not a validated general-purpose sampler.
# * AdaptiveMixtureMH is a custom full-covariance-adaptation baseline, NOT AIMM.
# * SAChES_inspired implements DE exploration -> independent adaptive local
#   chains. It is NOT a reproduction or assessment of published SAChES.
#   Source: https://www.sandia.gov/research/publications/details/
#           saches-scalable-adaptive-chain-ensemble-sampling-2017-08-01/
# * Mixture components, map clusters, and density modes are different things.
#   Diagnostics below use component cores and prespecified component GROUPS.
#   They do NOT claim an oracle density-basin classification or merge precision.
# * n in a map counts stored observations, NOT independent information or mass.
# * Larger proposal steps alone need not make scouts efficient in 10 dimensions.
# * Final kernels freeze adaptation. Their MH corrections target the original
#   density, but do not guarantee adequate finite-run exploration.
# * All algorithms execute serially here. FLOCK middle jobs depend only on the
#   completed scout stage, not on other middle jobs; no live communication.
#
# FAIRNESS / COST
# * Every run uses EXACTLY B target calls, including starts and bridge checks.
# * FLOCK uses 20% scouting, 45% middle/preparation, 35% final transitions.
# * Map and MapMerge both reserve 2% of B INSIDE the 45% preparation budget.
#   Merge checks spend this reserve; unused reserve becomes discarded global-MH
#   settling steps. Therefore MapMerge does not receive free bridge evaluations.
# * Map and MapMerge use the same scouts, middle-job budgets, policy, seeds,
#   covariance regularization and final kernel. Merging and its consequences
#   (including paid check cost) are the only design difference.
# * FLOCK_Novelty / Unexplored retain Experiment 3's exploration settings.
#   Gaussian draw/density calculations use a consistent Cholesky convention.
# * External methods have their own fixed tuning settings, logged to disk.
#   Equal evaluation budget is NOT a claim of equal tuning or wall-clock cost.
# * Reference scoring is retrospective, outside sampling time and algorithm
#   budgets; scoring time is reported separately.
# * Checkpoints are per replicate AND method. A manifest prevents accidental
#   mixing of settings or script versions; increasing REPS is permitted.

FLOCK5_VERSION <- "1.0.0"

ff_assert <- function(ok, message) {
  if (length(ok) != 1L || is.na(ok) || !ok) stop(message, call. = FALSE)
}
ff_int <- function(name, default, minimum = 1L) {
  txt <- Sys.getenv(name, unset = as.character(default))
  val <- suppressWarnings(as.numeric(txt))
  ff_assert(length(val) == 1L && is.finite(val) && val == floor(val) &&
              val >= minimum && val <= .Machine$integer.max,
            paste("Invalid", name, ":", txt))
  as.integer(val)
}
ff_bool <- function(name, default = TRUE) {
  x <- toupper(Sys.getenv(name, unset = if (default) "TRUE" else "FALSE"))
  ff_assert(x %in% c("TRUE", "FALSE", "1", "0"), paste("Invalid", name))
  x %in% c("TRUE", "1")
}
ff_config <- function() {
  mode <- tolower(Sys.getenv("RUN_MODE", unset = "smoke"))
  specs <- list(smoke = c(3L, 5000L), pilot = c(10L, 20000L),
                full = c(100L, 30000L))
  ff_assert(mode %in% names(specs), "RUN_MODE must be smoke, pilot, or full.")
  methods <- c("FLOCK_Novelty", "FLOCK_Unexplored", "FLOCK_Map", "FLOCK_MapMerge",
               "AdaptiveMixtureMH", "ParallelTempering", "SAChES_inspired",
               "RWMH_4chains")
  selected <- Sys.getenv("METHODS", unset = paste(methods, collapse = ","))
  selected <- unique(trimws(strsplit(selected, ",", fixed = TRUE)[[1L]]))
  ff_assert(length(selected) > 0L && all(selected %in% methods),
            paste("Unknown method. Choose from", paste(methods, collapse = ", ")))
  budget <- ff_int("BUDGET", specs[[mode]][2L], 3000L)
  list(version = FLOCK5_VERSION, mode = mode,
       reps = ff_int("REPS", specs[[mode]][1L]), budget = budget,
       methods = selected, all_methods = methods,
       outdir = Sys.getenv("OUTDIR", unset = paste0("flock_exp5_", mode, "_B", budget)),
       seed = ff_int("BASE_SEED", 205100L),
       resume = ff_bool("RESUME", TRUE), save_draws = ff_bool("SAVE_DRAWS", FALSE),
       bridge = ff_bool("BRIDGE_CHECK", TRUE),
       d = 10L, scout_fraction = .20, middle_fraction = .45,
       map_reserve_fraction = .02, final_burn_fraction = .10,
       scout_sd = 5, follower_sd = 1.15, pioneer_sd = 2.3,
       initial_sd = 7, candidate_sd = 8, n_candidates = 12L,
       local_sd = .8, local_probability = .20,
       broad_sd = 15, broad_weight = .25, eigen_floor = .04,
       max_regions = 8L, points_per_job = 150L, merge_mahal = 1.5,
       merge_symkl = 1, merge_drop = 2, bridge_points = 9L,
       bootstrap = 2000L)
}
ff_seed <- function(cfg, rep, stream) {
  s <- (as.double(cfg$seed) + 100003 * rep + 7919 * stream) %% 2147483646
  set.seed(as.integer(s + 1))
}
ff_split <- function(total, k) {
  ff_assert(total >= k && k >= 1L, "Too few evaluations for the requested jobs.")
  ans <- rep(as.integer(total %/% k), k)
  rem <- as.integer(total %% k)
  if (rem > 0L) ans[seq_len(rem)] <- ans[seq_len(rem)] + 1L
  ans
}
ff_pick <- function(ids) {
  ff_assert(length(ids) > 0L, "Cannot select from an empty set.")
  # sample(ids, 1) is unsafe when ids is a single numeric value.
  ids[sample.int(length(ids), 1L)]
}
ff_lse <- function(x) {
  a <- max(x)
  if (!is.finite(a)) return(a)
  a + log(sum(exp(x - a)))
}
ff_spd <- function(S, floor = .04) {
  S <- (S + t(S)) / 2
  ee <- eigen(S, symmetric = TRUE)
  tcrossprod(sweep(ee$vectors, 2L, sqrt(pmax(ee$values, floor)), "*"))
}
ff_gaussian <- function(mu, S, floor = .04) {
  S <- ff_spd(S, floor)
  R <- chol(S)                     # t(R) %*% R == S
  list(mu = as.numeric(mu), S = S, R = R,
       constant = -length(mu) / 2 * log(2 * pi) - sum(log(diag(R))))
}
ff_gdraw <- function(g) as.numeric(g$mu + t(g$R) %*% rnorm(length(g$mu)))
ff_glog <- function(x, g) {
  z <- backsolve(g$R, x - g$mu, transpose = TRUE)
  g$constant - sum(z * z) / 2
}
ff_glog_batch <- function(X, g) {
  X <- as.matrix(X)
  z <- backsolve(g$R, t(sweep(X, 2L, g$mu, "-")), transpose = TRUE)
  g$constant - colSums(z * z) / 2
}
ff_target <- function() {
  d <- 10L
  w <- c(.30, .25, .20, .15, .08, .02)
  mu <- rbind(c(-5,-5,0,0,0,0,0,0,0,0),
              c(-4.1,-4.5,.3,0,0,0,0,0,0,0),
              c(5,5,0,0,0,0,0,0,0,0),
              c(-5,5,3,0,0,0,0,0,0,0),
              c(-4.1,5.7,3.2,0,0,0,0,0,0,0),
              c(6,-6,-3,0,0,0,0,0,0,0))
  scales <- c(1.2, 1.1, 1.4, 1.1, 1.0, .65)
  S <- lapply(seq_len(6L), function(j) {
    A <- diag(scales[j]^2, d)
    A[1,2] <- A[2,1] <- c(.45,.35,-.45,.30,.35,0)[j] * scales[j]^2
    A[3,4] <- A[4,3] <- .25 * scales[j]^2
    A
  })
  g <- lapply(seq_len(6L), function(j) ff_gaussian(mu[j, ], S[[j]], 1e-10))
  log_density <- local({ gg <- g; ww <- w
    function(x) ff_lse(log(ww) + vapply(gg, function(a) ff_glog(x, a), numeric(1)))
  })
  mean <- colSums(mu * w)
  covariance <- Reduce("+", lapply(seq_len(6L), function(j) {
    z <- mu[j, ] - mean
    w[j] * (S[[j]] + tcrossprod(z))
  }))
  probs <- c(.025, .50, .975)
  quantiles <- sapply(seq_len(d), function(k) {
    ss <- sqrt(vapply(S, function(a) a[k,k], numeric(1)))
    vapply(probs, function(p) {
      uniroot(function(t) sum(w * pnorm(t, mu[,k], ss)) - p,
              c(min(mu[,k] - 12 * ss), max(mu[,k] + 12 * ss)),
              tol = 1e-11)$root
    }, numeric(1))
  })
  list(d = d, w = w, mu = mu, S = S, g = g, log_density = log_density,
       mean = mean, covariance = covariance, quantiles = quantiles,
       quantile_probs = probs, groups = c(1L,1L,2L,3L,3L,4L),
       group_names = c("AB", "C", "DE", "F"), group_w = c(.55,.20,.23,.02),
       core_cutoff = qchisq(.95, d),
       tail = sum(w * pnorm(0, mu[,1], sqrt(vapply(S, function(a) a[1,1], numeric(1))),
                            lower.tail = FALSE)))
}
ff_counter <- function(log_density, budget, d) {
  e <- new.env(parent = emptyenv())
  e$n <- 0L; e$stage <- "initialization"
  e$points <- matrix(NA_real_, budget, d); e$stages <- character(budget)
  e$set_stage <- function(s) { e$stage <- s; invisible(NULL) }
  e$eval <- function(x) {
    ff_assert(length(x) == d && all(is.finite(x)), "Invalid parameter vector.")
    ff_assert(e$n < budget, "Target-evaluation budget would be exceeded.")
    e$n <- e$n + 1L
    e$points[e$n, ] <- x; e$stages[e$n] <- e$stage
    ans <- log_density(x)
    ff_assert(length(ans) == 1L && !is.na(ans) && ans != Inf, "Invalid log target.")
    ans
  }
  e$count <- function() e$n
  e$trace <- function() list(points = e$points[seq_len(e$n), , drop = FALSE],
                             stage = e$stages[seq_len(e$n)])
  e
}
ff_rw <- function(f, x, lp, n, R) {
  d <- length(x)
  h <- matrix(NA_real_, n + 1L, d); ll <- numeric(n + 1L)
  h[1L, ] <- x; ll[1L] <- lp; accepted <- 0L
  for (i in seq_len(n)) {
    y <- as.numeric(x + t(R) %*% rnorm(d)); ly <- f$eval(y)
    if (log(runif(1)) < ly - lp) { x <- y; lp <- ly; accepted <- accepted + 1L }
    h[i + 1L, ] <- x; ll[i + 1L] <- lp
  }
  list(x = x, lp = lp, h = h, ll = ll, accepted = accepted, proposals = n)
}
ff_am_walk <- function(f, x, lp, n, cfg, initial_sd = 1, adapt_until = n) {
  d <- length(x); R <- diag(initial_sd, d)
  h <- matrix(NA_real_, n + 1L, d); ll <- numeric(n + 1L)
  h[1L, ] <- x; ll[1L] <- lp; acc <- 0L
  for (i in seq_len(n)) {
    if (i > max(100L, 5L*d) && i <= adapt_until && i %% 100L == 0L) {
      # Use the latter half of history to reduce transient-scale inflation.
      z <- h[seq.int(max(1L, floor(i/2)), i), , drop = FALSE]
      S <- (2.38^2 / d) * cov(z) + diag(.01, d)
      R <- chol(ff_spd(S, cfg$eigen_floor))
    }
    # A fixed small symmetric proposal prevents a singular adaptive covariance
    # from being the only proposal. Each branch is symmetric conditional on past.
    y <- if (runif(1) < .05) x + rnorm(d, sd = .15) else
      as.numeric(x + t(R) %*% rnorm(d))
    ly <- f$eval(y)
    if (log(runif(1)) < ly - lp) { x <- y; lp <- ly; acc <- acc + 1L }
    h[i + 1L, ] <- x; ll[i + 1L] <- lp
  }
  list(x = x, lp = lp, h = h, ll = ll, accepted = acc, proposals = n)
}
ff_pool <- function(jobs, max_each = 150L) {
  chunks <- lapply(jobs, function(z) {
    ix <- unique(as.integer(round(seq(1L, nrow(z$h), length.out = min(max_each, nrow(z$h))))))
    list(X = z$h[ix, , drop = FALSE], lp = z$ll[ix])
  })
  list(X = do.call(rbind, lapply(chunks, `[[`, "X")),
       lp = unlist(lapply(chunks, `[[`, "lp"), use.names = FALSE))
}
ff_region <- function(X, lp, id, cfg, provenance = id) {
  n <- nrow(X); nu <- nrow(unique(X)); d <- ncol(X)
  S <- if (nu >= d + 3L) cov(X) + diag(.25, d) else diag(2, d)
  best <- which.max(lp)
  list(id = id, provenance = provenance, X = X, lp = lp, n = n, n_unique = nu,
       mu = colMeans(X), S = ff_spd(S, cfg$eigen_floor),
       peak = lp[best], anchor = as.numeric(X[best, ]), anchor_lp = lp[best])
}
ff_make_map <- function(jobs, cfg, rep, stream) {
  p <- ff_pool(jobs, cfg$points_per_job)
  U <- unique(p$X); k <- min(cfg$max_regions, nrow(U))
  ff_assert(k >= 1L, "Empty exploration pool.")
  ff_seed(cfg, rep, stream)
  if (k == 1L) lab <- rep.int(1L, nrow(p$X)) else {
    centers <- U[sample.int(nrow(U), k), , drop = FALSE]
    km <- tryCatch(kmeans(p$X, centers = centers, iter.max = 100L,
                           algorithm = "Lloyd"), error = function(e) NULL)
    if (is.null(km)) {
      warning("Map clustering failed; using one region.")
      lab <- rep.int(1L, nrow(p$X))
    } else lab <- km$cluster
  }
  ids <- sort(unique(lab))
  regions <- lapply(seq_along(ids), function(j) {
    keep <- lab == ids[j]
    ff_region(p$X[keep, , drop = FALSE], p$lp[keep], paste0("R", j), cfg)
  })
  list(regions = regions, pool = p, original_regions = length(regions))
}
ff_q_from_regions <- function(regions, cfg) {
  g <- lapply(regions, function(z) ff_gaussian(z$mu, z$S, cfg$eigen_floor))
  k <- length(g)
  g[[k + 1L]] <- ff_gaussian(rep(0, cfg$d), diag(cfg$broad_sd^2, cfg$d), 1e-8)
  w <- c(rep((1 - cfg$broad_weight) / k, k), cfg$broad_weight)
  logq <- function(x) ff_lse(log(w) + vapply(g, function(a) ff_glog(x, a), numeric(1)))
  draw <- function() ff_gdraw(g[[sample.int(length(w), 1L, prob = w)]])
  list(draw = draw, logq = logq, k = k, weights = w, gaussians = g)
}
ff_global <- function(f, q, x, lp, n, cfg) {
  h <- matrix(NA_real_, n + 1L, cfg$d); ll <- numeric(n + 1L)
  h[1L, ] <- x; ll[1L] <- lp; lq <- q$logq(x)
  acc <- glob_acc <- glob_try <- 0L
  for (i in seq_len(n)) {
    local <- runif(1) < cfg$local_probability
    y <- if (local) x + rnorm(cfg$d, sd = cfg$local_sd) else q$draw()
    ly <- f$eval(y); lqy <- q$logq(y)
    loga <- if (local) ly - lp else ly - lp + lq - lqy
    if (!local) glob_try <- glob_try + 1L
    if (log(runif(1)) < loga) {
      x <- y; lp <- ly; lq <- lqy; acc <- acc + 1L
      if (!local) glob_acc <- glob_acc + 1L
    }
    h[i + 1L, ] <- x; ll[i + 1L] <- lp
  }
  list(x = x, lp = lp, h = h, ll = ll, accepted = acc, proposals = n,
       global_accepted = glob_acc, global_proposals = glob_try)
}
ff_retain <- function(job, fraction = .10) {
  n <- nrow(job$h) - 1L
  first <- 2L + floor(n * fraction)  # initial state and warmup transitions excluded
  ff_assert(first <= nrow(job$h), "No retained samples after warmup.")
  job$h[seq.int(first, nrow(job$h)), , drop = FALSE]
}

# --- Mapping, allocation and conservative merge audit -------------------------
ff_pair <- function(a, b, cfg) {
  A <- ff_spd(a$S, cfg$eigen_floor); B <- ff_spd(b$S, cfg$eigen_floor)
  delta <- a$mu - b$mu
  P <- ff_spd((A+B)/2, cfg$eigen_floor)
  mah <- sqrt(max(0, as.numeric(crossprod(delta, solve(P, delta)))))
  kl <- .25 * (sum(diag(solve(A, B))) + sum(diag(solve(B, A))) - 2*cfg$d +
                 as.numeric(crossprod(delta, solve(A, delta) + solve(B, delta))))
  list(mahal = mah, symkl = max(0, kl))
}
ff_empty_audit <- function() data.frame(
  phase = character(), left = character(), right = character(),
  n_left = integer(), n_right = integer(), unique_left = integer(),
  unique_right = integer(), mahalanobis = double(), symmetric_kl = double(),
  bridge_drop = double(), bridge_evals = integer(), merged = logical(),
  reason = character(), stringsAsFactors = FALSE)
ff_merge_map <- function(map, f, cap, cfg, phase) {
  initial <- f$count(); audit <- ff_empty_audit(); snapshots <- list()
  repeat {
    regs <- map$regions; nr <- length(regs)
    if (nr < 2L) break
    pairs <- t(combn(seq_len(nr), 2L))
    metrics <- lapply(seq_len(nrow(pairs)), function(k)
      ff_pair(regs[[pairs[k,1L]]], regs[[pairs[k,2L]]], cfg))
    ord <- order(vapply(metrics, `[[`, numeric(1), "symkl"))
    changed <- FALSE
    for (k in ord) {
      i <- pairs[k,1L]; j <- pairs[k,2L]; a <- regs[[i]]; b <- regs[[j]]
      m <- metrics[[k]]; bridge_drop <- NA_real_; calls <- 0L
      enough <- min(a$n_unique, b$n_unique) >= cfg$d + 3L
      geometry <- m$mahal <= cfg$merge_mahal && m$symkl <= cfg$merge_symkl
      pass <- enough && geometry
      reason <- if (!enough) "insufficient_unique_points" else
        if (!geometry) "geometry_reject" else "geometry_accept"
      if (pass && cfg$bridge) {
        if (f$count() - initial + cfg$bridge_points > cap) {
          pass <- FALSE; reason <- "bridge_budget_exhausted"
        } else {
          f$set_stage(paste0(phase, "_bridge"))
          vals <- vapply(seq(0,1,length.out=cfg$bridge_points), function(t)
            f$eval((1-t)*a$mu + t*b$mu), numeric(1))
          calls <- length(vals)
          bridge_drop <- if (all(is.finite(vals)))
            max(0, min(vals[c(1L,length(vals))]) - min(vals)) else Inf
          pass <- is.finite(bridge_drop) && bridge_drop <= cfg$merge_drop
          reason <- if (pass) "geometry_and_bridge_accept" else "density_valley_reject"
        }
      }
      row <- data.frame(phase=phase, left=a$id, right=b$id, n_left=a$n, n_right=b$n,
                        unique_left=a$n_unique, unique_right=b$n_unique,
                        mahalanobis=m$mahal, symmetric_kl=m$symkl,
                        bridge_drop=bridge_drop, bridge_evals=calls,
                        merged=pass, reason=reason, stringsAsFactors=FALSE)
      audit <- rbind(audit, row)
      if (pass) {
        provenance <- c(a$provenance, b$provenance)
        merged <- ff_region(rbind(a$X,b$X), c(a$lp,b$lp),
                            paste(provenance,collapse="+"), cfg, provenance)
        snapshots[[length(snapshots)+1L]] <- list(phase=phase, left=a, right=b,
                                                 merged=merged)
        map$regions[[i]] <- merged; map$regions[[j]] <- NULL
        changed <- TRUE
        break
      }
    }
    if (!changed) break
  }
  ff_assert(f$count()-initial <= cap, "Merge check exceeded its reserve.")
  list(map=map, audit=audit, snapshots=snapshots, calls=f$count()-initial)
}
ff_assign <- function(map, number = 3L) {
  regs <- map$regions
  peaks <- vapply(regs, `[[`, numeric(1), "peak")
  nu <- vapply(regs, `[[`, integer(1), "n_unique")
  # An explicitly heuristic priority, NOT estimated posterior probability.
  # Reward poorly mapped but not wildly implausible discovered regions.
  score <- exp(pmax(-40, peaks-max(peaks))/10) / sqrt(1+nu)
  allocations <- integer(length(regs)); chosen <- integer(number)
  for (j in seq_len(number)) {
    chosen[j] <- which.max(score/(1+allocations))
    allocations[chosen[j]] <- allocations[chosen[j]]+1L
  }
  list(ids=chosen, score=score)
}
ff_coverage <- function(candidates, scout_history) {
  ix <- unique(as.integer(round(seq(1,nrow(scout_history),
                                    length.out=min(250L,nrow(scout_history))))))
  h <- scout_history[ix,,drop=FALSE]
  scale <- pmax(apply(h,2L,sd),2)
  vapply(seq_len(nrow(candidates)), function(i) {
    z <- sweep(sweep(h,2L,candidates[i,],"-"),2L,scale,"/")
    sqrt(min(rowSums(z*z)))
  }, numeric(1))
}
ff_allocation_row <- function(job, role, source, n, x, lp, R) {
  data.frame(job=job, role=role, source_region=source, target_calls=n,
             start_log_target=lp, proposal_rms_sd=sqrt(mean(colSums(R*R))),
             start=paste(signif(x,7),collapse=";"), stringsAsFactors=FALSE)
}
ff_finish <- function(draws, stages, f, accept, chain_id = rep(1L,nrow(draws)),
                      map_snapshots=list(), audit=ff_empty_audit(), merge_snapshots=list(),
                      allocation=data.frame(), q_components=NA_integer_, extra=list()) {
  list(draws=draws, chain_id=chain_id, stages=stages, evals=f$count(),
       evaluated=f$trace(), accept=accept, map_snapshots=map_snapshots,
       merge_audit=audit, merge_snapshots=merge_snapshots, allocation=allocation,
       proposal_components=q_components, extra=extra)
}
ff_flock <- function(method, f, init, cfg, rep_id) {
  mapped <- method %in% c("FLOCK_Map","FLOCK_MapMerge")
  merging <- method == "FLOCK_MapMerge"
  bs <- floor(cfg$budget*cfg$scout_fraction)
  bm <- floor(cfg$budget*cfg$middle_fraction)
  reserve <- if (mapped) floor(cfg$budget*cfg$map_reserve_fraction) else 0L
  scouts <- vector("list",2L); scouts_alloc <- ff_split(bs,2L)
  for (k in seq_len(2L)) {
    ff_seed(cfg,rep_id,100L+k); f$set_stage("scout")
    x <- init[k,]; lp <- f$eval(x)
    scouts[[k]] <- ff_rw(f,x,lp,scouts_alloc[k]-1L,diag(cfg$scout_sd,cfg$d))
  }
  ff_assert(f$count()==bs,"Scout budget mismatch.")
  ends <- do.call(rbind,lapply(scouts,`[[`,"x"))
  lps <- vapply(scouts,`[[`,numeric(1),"lp"); best <- order(lps,decreasing=TRUE)
  scout_history <- do.call(rbind,lapply(scouts,`[[`,"h"))
  audit <- ff_empty_audit(); merge_snapshots <- list(); snapshots <- list()
  merge_calls <- 0L
  if (mapped) {
    smap <- ff_make_map(scouts,cfg,rep_id,300L)
    snapshots$scout_before <- smap
    if (merging) {
      z <- ff_merge_map(smap,f,floor(reserve/2),cfg,"scout_map")
      smap <- z$map; audit <- rbind(audit,z$audit)
      merge_snapshots <- c(merge_snapshots,z$snapshots); merge_calls <- z$calls
    }
    snapshots$scout_after <- smap
    assignment <- ff_assign(smap,3L)
  }
  alloc <- ff_split(bm-reserve,6L)
  middle <- vector("list",6L); allocation <- vector("list",6L)
  # Every middle job receives the SAME frozen scout map. No job sees outputs
  # from any other middle job. Seeds are isolated by job to improve pairing.
  for (k in seq_len(6L)) {
    ff_seed(cfg,rep_id,400L+k); n <- alloc[k]
    role <- if (k<=3L) "follower" else "pioneer"
    f$set_stage(role)
    if (k<=3L) {
      if (mapped) {
        region <- smap$regions[[assignment$ids[k]]]
        x <- region$anchor; lp <- region$anchor_lp
        R <- chol(ff_spd((2.38^2/cfg$d)*region$S, cfg$eigen_floor))
        source <- region$id
      } else {
        j <- best[(k-1L)%%2L+1L]; x <- ends[j,]; lp <- lps[j]
        if (k>2L) { x <- x+rnorm(cfg$d,sd=.7); lp <- f$eval(x); n <- n-1L }
        R <- diag(cfg$follower_sd,cfg$d); source <- paste0("scout_endpoint_",j)
      }
    } else {
      nc <- min(cfg$n_candidates,max(2L,n%/%20L))
      candidates <- matrix(rnorm(nc*cfg$d,sd=cfg$candidate_sd),nc,cfg$d)
      f$set_stage("pioneer_screen")
      lc <- apply(candidates,1L,f$eval)
      novelty <- ff_coverage(candidates,scout_history)
      if (method=="FLOCK_Novelty") {
        chosen <- which.max(novelty-.15*pmax(0,max(lc)-lc))
      } else {
        # Relative coverage rule; NOT a claim these points lie in a truly
        # unexplored posterior region. Both map variants use this same policy.
        chosen <- ff_pick(which(novelty>=median(novelty)))
      }
      x <- candidates[chosen,]; lp <- lc[chosen]; n <- n-nc
      R <- diag(cfg$pioneer_sd,cfg$d); source <- "candidate_pool"
      f$set_stage("pioneer")
    }
    allocation[[k]] <- ff_allocation_row(k,role,source,alloc[k],x,lp,R)
    middle[[k]] <- ff_rw(f,x,lp,n,R)
  }
  alljobs <- c(scouts,middle)
  fmap <- ff_make_map(alljobs,cfg,rep_id,800L)
  if (mapped) snapshots$final_before <- fmap
  if (merging) {
    z <- ff_merge_map(fmap,f,reserve-floor(reserve/2),cfg,"final_map")
    fmap <- z$map; audit <- rbind(audit,z$audit)
    merge_snapshots <- c(merge_snapshots,z$snapshots); merge_calls <- merge_calls+z$calls
  }
  if (mapped) snapshots$final_after <- fmap
  q <- ff_q_from_regions(fmap$regions,cfg)
  x <- ends[best[1L],]; lp <- lps[best[1L]]
  stages <- list(scout=scout_history,
                 follower=do.call(rbind,lapply(middle[1:3],`[[`,"h")),
                 pioneer=do.call(rbind,lapply(middle[4:6],`[[`,"h")))
  settle <- reserve-merge_calls
  if (settle>0L) {
    ff_seed(cfg,rep_id,850L); f$set_stage("map_settle")
    z <- ff_global(f,q,x,lp,settle,cfg)
    x <- z$x; lp <- z$lp; stages$map_settle <- z$h
  }
  ff_assert(f$count()==bs+bm,"Exploration/preparation budget mismatch.")
  ff_seed(cfg,rep_id,900L); f$set_stage("final")
  fin <- ff_global(f,q,x,lp,cfg$budget-f$count(),cfg)
  draws <- ff_retain(fin,cfg$final_burn_fraction)
  stages$final_all <- fin$h; stages$final_retained <- draws
  ff_assert(f$count()==cfg$budget,"FLOCK total budget mismatch.")
  ff_finish(draws,stages,f,fin$accepted/fin$proposals,
            map_snapshots=snapshots,audit=audit,merge_snapshots=merge_snapshots,
            allocation=do.call(rbind,allocation),q_components=q$k,
            extra=list(merge_calls=merge_calls,settle_calls=settle,
                       scout_calls=bs,middle_preparation_calls=bm,
                       final_calls=fin$proposals,
                       global_accept=if(fin$global_proposals>0L)
                         fin$global_accepted/fin$global_proposals else NA_real_))
}

# --- External/custom comparators ---------------------------------------------
ff_adaptive_mixture <- function(f,init,cfg,rep_id) {
  ne <- floor(cfg$budget*cfg$scout_fraction)+floor(cfg$budget*cfg$middle_fraction)
  alloc <- ff_split(ne,8L); jobs <- vector("list",8L)
  for (k in seq_len(8L)) {
    ff_seed(cfg,rep_id,1000L+k); f$set_stage("adaptive_exploration")
    x <- init[k,]; lp <- f$eval(x)
    jobs[[k]] <- ff_am_walk(f,x,lp,alloc[k]-1L,cfg,initial_sd=1.5)
  }
  ff_assert(f$count()==ne,"Adaptive-mixture exploration budget mismatch.")
  map <- ff_make_map(jobs,cfg,rep_id,800L); q <- ff_q_from_regions(map$regions,cfg)
  best <- which.max(vapply(jobs,`[[`,numeric(1),"lp"))
  ff_seed(cfg,rep_id,900L); f$set_stage("final")
  z <- ff_global(f,q,jobs[[best]]$x,jobs[[best]]$lp,cfg$budget-f$count(),cfg)
  draws <- ff_retain(z,cfg$final_burn_fraction)
  ff_finish(draws,list(exploration=do.call(rbind,lapply(jobs,`[[`,"h")),
                       final_all=z$h,final_retained=draws),f,z$accepted/z$proposals,
            q_components=q$k,extra=list(final_calls=z$proposals))
}
ff_rwmh <- function(f,init,cfg,rep_id) {
  alloc <- ff_split(cfg$budget,4L); jobs <- retained <- vector("list",4L)
  acc <- proposals <- 0L
  for (k in seq_len(4L)) {
    ff_seed(cfg,rep_id,1100L+k); f$set_stage("rwmh")
    x <- init[k,]; lp <- f$eval(x)
    z <- ff_rw(f,x,lp,alloc[k]-1L,diag(2.38/sqrt(cfg$d),cfg$d))
    jobs[[k]] <- z; retained[[k]] <- ff_retain(z,.35)
    acc <- acc+z$accepted; proposals <- proposals+z$proposals
  }
  draws <- do.call(rbind,retained)
  ff_finish(draws,list(all_chains=do.call(rbind,lapply(jobs,`[[`,"h")),
                       final_retained=draws),f,acc/proposals,
            chain_id=rep(seq_len(4L),vapply(retained,nrow,integer(1))))
}
ff_pt <- function(f,init,cfg,rep_id) {
  ff_seed(cfg,rep_id,1200L)
  K <- 8L; betas <- exp(seq(log(1),log(.03),length.out=K))
  x <- init; f$set_stage("pt_initialization"); lp <- apply(x,1L,f$eval)
  updates <- cfg$budget-f$count(); sweeps <- ceiling(updates/K)
  histories <- lapply(seq_len(K),function(k) matrix(NA_real_,sweeps+1L,cfg$d))
  for(k in seq_len(K)) histories[[k]][1L,] <- x[k,]
  count_at <- integer(sweeps+1L); count_at[1L] <- f$count()
  warmup <- floor(.35*cfg$budget)
  logstep <- log(2.38/sqrt(cfg$d)/sqrt(betas))
  ba <- bn <- integer(K); windows <- integer(K)
  accept <- 0L; swap_accept <- swap_attempt <- integer(K-1L)
  for (s in seq_len(sweeps)) {
    for (k in seq_len(K)) {
      if(f$count()>=cfg$budget) break
      f$set_stage(if(f$count()<warmup) "pt_warmup" else "pt_sampling")
      y <- x[k,]+rnorm(cfg$d,sd=exp(logstep[k])); ly <- f$eval(y)
      ok <- log(runif(1)) < betas[k]*(ly-lp[k])
      if(ok) { x[k,] <- y; lp[k] <- ly; accept <- accept+1L }
      if(f$count()<=warmup) {
        bn[k] <- bn[k]+1L; ba[k] <- ba[k]+as.integer(ok)
        if(bn[k]>=50L) {
          windows[k] <- windows[k]+1L
          gain <- min(.5,1/sqrt(windows[k]))
          logstep[k] <- pmin(log(20),pmax(log(.03),
                                      logstep[k]+gain*(ba[k]/bn[k]-.234)))
          bn[k] <- 0L; ba[k] <- 0L
        }
      }
    }
    # Disjoint adjacent swaps, alternating parity; cached log targets cost zero.
    for (k in seq.int(if(s%%2L==1L) 1L else 2L,K-1L,by=2L)) {
      swap_attempt[k] <- swap_attempt[k]+1L
      loga <- (betas[k]-betas[k+1L])*(lp[k+1L]-lp[k])
      if(log(runif(1))<loga) {
        tmp <- x[k,]; x[k,] <- x[k+1L,]; x[k+1L,] <- tmp
        tt <- lp[k]; lp[k] <- lp[k+1L]; lp[k+1L] <- tt
        swap_accept[k] <- swap_accept[k]+1L
      }
    }
    for(k in seq_len(K)) histories[[k]][s+1L,] <- x[k,]
    count_at[s+1L] <- f$count()
  }
  keep <- which(count_at>warmup)
  draws <- histories[[1L]][keep,,drop=FALSE]
  ff_finish(draws,list(cold_all=histories[[1L]],
                       hot_all=do.call(rbind,histories[-1L]),final_retained=draws),
            f,accept/updates,extra=list(betas=betas,steps=exp(logstep),
                swap_accept=swap_accept,swap_attempt=swap_attempt))
}
ff_saches <- function(f,init,cfg,rep_id) {
  # Architecture-inspired control only. Published SAChES uses parallel chains;
  # this code is serial and does not reproduce its adaptation/stopping design.
  ff_seed(cfg,rep_id,1300L); K <- 8L; x <- init
  ne <- floor(cfg$budget*cfg$scout_fraction)+floor(cfg$budget*cfg$middle_fraction)
  f$set_stage("de_exploration"); lp <- apply(x,1L,f$eval)
  nupdates <- ne-K; capacities <- ff_split(nupdates,K)+1L
  h <- lapply(capacities,function(n) matrix(NA_real_,n,cfg$d)); pos <- rep(1L,K)
  for(k in seq_len(K)) h[[k]][1L,] <- x[k,]
  for(t in seq_len(nupdates)) {
    k <- (t-1L)%%K+1L
    pair <- sample(setdiff(seq_len(K),k),2L,replace=FALSE)
    gamma <- if(runif(1)<.1) 1 else 2.38/sqrt(2*cfg$d)
    y <- x[k,]+gamma*(x[pair[1L],]-x[pair[2L],])+rnorm(cfg$d,sd=.1)
    ly <- f$eval(y)
    if(log(runif(1))<ly-lp[k]) { x[k,] <- y; lp[k] <- ly }
    pos[k] <- pos[k]+1L; h[[k]][pos[k],] <- x[k,]
  }
  ff_assert(f$count()==ne,"DE exploration budget mismatch.")
  alloc <- ff_split(cfg$budget-ne,K); retained <- jobs <- vector("list",K)
  acc <- proposals <- 0L
  for(k in seq_len(K)) {
    ff_seed(cfg,rep_id,1350L+k); f$set_stage("adaptive_refinement")
    nw <- floor(.30*alloc[k])
    z <- ff_am_walk(f,x[k,],lp[k],alloc[k],cfg,initial_sd=2.38/sqrt(cfg$d),
                    adapt_until=nw)
    jobs[[k]] <- z; retained[[k]] <- ff_retain(z,.30)
    acc <- acc+z$accepted; proposals <- proposals+z$proposals
  }
  draws <- do.call(rbind,retained)
  ff_finish(draws,list(de_exploration=do.call(rbind,h),
                       refinement=do.call(rbind,lapply(jobs,`[[`,"h")),
                       final_retained=draws),f,acc/proposals,
            chain_id=rep(seq_len(K),vapply(retained,nrow,integer(1))))
}
ff_run_method <- function(method,log_density,init,cfg,rep_id) {
  f <- ff_counter(log_density,cfg$budget,cfg$d)
  before <- proc.time()
  ans <- if(grepl("^FLOCK_",method)) ff_flock(method,f,init,cfg,rep_id) else
    switch(method,AdaptiveMixtureMH=ff_adaptive_mixture(f,init,cfg,rep_id),
           ParallelTempering=ff_pt(f,init,cfg,rep_id),
           SAChES_inspired=ff_saches(f,init,cfg,rep_id),
           RWMH_4chains=ff_rwmh(f,init,cfg,rep_id),stop("Unknown method"))
  timing <- proc.time()-before
  ans$elapsed <- unname(timing["elapsed"])
  ans$cpu <- unname(timing["user.self"]+timing["sys.self"])
  ff_assert(ans$evals==cfg$budget,"Total budget mismatch.")
  ff_assert(nrow(ans$draws)>=2L && all(is.finite(ans$draws)),"Invalid final draws.")
  ans
}

# --- Retrospective reference diagnostics: NEVER passed to the samplers ---------
ff_oracle <- function(X,target) {
  X <- as.matrix(X); n <- nrow(X)
  if(n==0L) return(list(r=matrix(numeric(),0,6),core=matrix(logical(),0,6),
                       group_r=matrix(numeric(),0,4)))
  logg <- matrix(vapply(target$g,function(g) ff_glog_batch(X,g),numeric(n)),n,6L)
  lcomp <- sweep(logg,2L,log(target$w),"+")
  mx <- apply(lcomp,1L,max)
  r <- exp(lcomp-mx); r <- r/rowSums(r)
  constants <- vapply(target$g,`[[`,numeric(1),"constant")
  distance2 <- -2*sweep(logg,2L,constants,"-")
  core <- distance2<=target$core_cutoff
  group_r <- sapply(seq_len(4L),function(k)
    rowSums(r[,target$groups==k,drop=FALSE]))
  list(r=r,core=core,group_r=matrix(group_r,n,4L))
}
ff_discovery_block <- function(X,target,rep_id,method,source,stage) {
  z <- ff_oracle(X,target)
  component <- data.frame(rep=rep_id,method=method,source=source,stage=stage,
                          unit="component_core",region=LETTERS[1:6],
                          visited=colSums(z$core)>0L,count=colSums(z$core),
                          states=nrow(X),stringsAsFactors=FALSE)
  group_core <- matrix(vapply(seq_len(4L),function(k)
    rowSums(z$core[,target$groups==k,drop=FALSE])>0L, logical(nrow(X))), nrow(X), 4L)
  group <- data.frame(rep=rep_id,method=method,source=source,stage=stage,
                      unit="component_group_core",region=target$group_names,
                      visited=colSums(group_core)>0L,count=colSums(group_core),
                      states=nrow(X),stringsAsFactors=FALSE)
  rbind(component,group)
}
ff_score <- function(ans,target,cfg,rep_id,method) {
  before <- proc.time()["elapsed"]
  X <- ans$draws; o <- ff_oracle(X,target); wh <- colMeans(o$r)
  gh <- colMeans(o$group_r); mh <- colMeans(X); Sh <- cov(X)
  scales <- sqrt(diag(target$covariance))
  qh <- apply(X,2L,quantile,probs=target$quantile_probs,names=FALSE,type=8)
  groups <- max.col(o$group_r,ties.method="first")
  transitions <- sum(vapply(split(seq_len(nrow(X)),ans$chain_id),function(ii) {
    if(length(ii)<2L) return(0L)
    sum(diff(groups[ii])!=0L)
  },integer(1)))
  # Ratios use analytical reference scales, not information given to samplers.
  result <- data.frame(rep=rep_id,method=method,evals=ans$evals,
    retained=nrow(X),chains_retained=length(unique(ans$chain_id)),
    mean_rmse=sqrt(mean((mh-target$mean)^2)),
    standardized_mean_rmse=sqrt(mean(((mh-target$mean)/scales)^2)),
    covariance_relative_error=sqrt(sum((Sh-target$covariance)^2))/
                              sqrt(sum(target$covariance^2)),
    weight_l1=sum(abs(wh-target$w)),group_weight_l1=sum(abs(gh-target$group_w)),
    rare_est=wh[6L],rare_abs_error=abs(wh[6L]-target$w[6L]),
    rare_core_visited=as.integer(any(o$core[,6L])),
    quantile_standardized_rmse=sqrt(mean(sweep(qh-target$quantiles,2L,scales,"/")^2)),
    tail_abs_error=abs(mean(X[,1L]>0)-target$tail),
    macro_group_transitions=transitions,accept=ans$accept,
    elapsed_seconds=ans$elapsed,cpu_seconds=ans$cpu,
    proposal_components=ans$proposal_components,
    merges=sum(ans$merge_audit$merged),
    bridge_evaluations=sum(ans$merge_audit$bridge_evals),stringsAsFactors=FALSE)
  for(j in seq_len(6L)) result[[paste0("weight_",LETTERS[j])]] <- wh[j]
  for(j in seq_len(cfg$d)) result[[paste0("mean_theta",j)]] <- mh[j]
  discovery <- list(); i <- 0L
  for(s in names(ans$stages)) {
    i <- i+1L
    discovery[[i]] <- ff_discovery_block(ans$stages[[s]],target,rep_id,method,
                                         "chain_states",s)
  }
  for(s in unique(ans$evaluated$stage)) {
    i <- i+1L; keep <- ans$evaluated$stage==s
    discovery[[i]] <- ff_discovery_block(ans$evaluated$points[keep,,drop=FALSE],
                          target,rep_id,method,"evaluated_locations",s)
  }
  oe <- ff_oracle(ans$evaluated$points,target)
  first <- vapply(seq_len(6L),function(j) {
    hit <- which(oe$core[,j]); if(length(hit)) hit[1L] else NA_integer_
  },integer(1))
  first_evaluated <- data.frame(rep=rep_id,method=method,component=LETTERS[1:6],
                               first_evaluation=first,
                               stage=ans$evaluated$stage[first],stringsAsFactors=FALSE)
  counts <- table(ans$evaluated$stage)
  budgets <- data.frame(rep=rep_id,method=method,stage=names(counts),
                         calls=as.integer(counts),stringsAsFactors=FALSE)
  ff_assert(sum(budgets$calls)==cfg$budget,"Stage accounting does not sum to B.")
  map_rows <- list(); k <- 0L
  for(phase in names(ans$map_snapshots)) {
    map <- ans$map_snapshots[[phase]]
    for(region in map$regions) {
      k <- k+1L
      diagnostic_groups <- colMeans(ff_oracle(region$X,target)$group_r)
      map_rows[[k]] <- data.frame(rep=rep_id,method=method,phase=phase,
        region=region$id,n_stored=region$n,n_unique=region$n_unique,
        regions_in_map=length(map$regions),peak_log_target=region$peak,
        center=paste(signif(region$mu,7),collapse=";"),
        covariance_trace=sum(diag(region$S)),
        dominant_reference_group=target$group_names[which.max(diagnostic_groups)],
        reference_group_concentration=max(diagnostic_groups),stringsAsFactors=FALSE)
    }
  }
  audit <- ans$merge_audit
  if(nrow(audit)) { audit$rep <- rep_id; audit$method <- method }
  merge_diag <- list()
  for(k in seq_along(ans$merge_snapshots)) {
    z <- ans$merge_snapshots[[k]]
    a <- colMeans(ff_oracle(z$left$X,target)$group_r)
    b <- colMeans(ff_oracle(z$right$X,target)$group_r)
    merge_diag[[k]] <- data.frame(rep=rep_id,method=method,phase=z$phase,
       left=z$left$id,right=z$right$id,group_left=target$group_names[which.max(a)],
       group_right=target$group_names[which.max(b)],
       different_dominant_groups=which.max(a)!=which.max(b),
       group_responsibility_overlap=sum(pmin(a,b)),stringsAsFactors=FALSE)
  }
  allocation <- ans$allocation
  if(nrow(allocation)) { allocation$rep <- rep_id; allocation$method <- method }
  result$diagnostic_seconds <- unname(proc.time()["elapsed"]-before)
  list(result=result,discovery=do.call(rbind,discovery),first_evaluated=first_evaluated,
       budgets=budgets,maps=if(length(map_rows))do.call(rbind,map_rows) else data.frame(),
       merge_audit=audit,merge_diagnostics=if(length(merge_diag))
         do.call(rbind,merge_diag) else data.frame(),allocation=allocation,
       algorithm_details=ans$extra,
       draws=if(cfg$save_draws || rep_id==1L) ans$draws else NULL,
       chain_id=if(cfg$save_draws || rep_id==1L) ans$chain_id else NULL)
}
ff_bind <- function(bundles,field) {
  z <- lapply(bundles,`[[`,field)
  z <- Filter(function(x) is.data.frame(x) && nrow(x)>0L,z)
  if(length(z)) do.call(rbind,z) else data.frame()
}
ff_summary <- function(res) {
  metrics <- c("mean_rmse","standardized_mean_rmse","covariance_relative_error",
               "weight_l1","group_weight_l1","rare_abs_error",
               "quantile_standardized_rmse","tail_abs_error","elapsed_seconds")
  do.call(rbind,lapply(split(res,res$method),function(z) {
    a <- data.frame(method=z$method[1L],n=nrow(z),evals=min(z$evals),
                    rare_core_discovery_rate=mean(z$rare_core_visited),
                    mean_rare_est=mean(z$rare_est),mean_merges=mean(z$merges),
                    mean_bridge_evaluations=mean(z$bridge_evaluations),
                    mean_retained=mean(z$retained),stringsAsFactors=FALSE)
    for(metric in metrics) {
      a[[paste0("mean_",metric)]] <- mean(z[[metric]])
      a[[paste0("median_",metric)]] <- median(z[[metric]])
      a[[paste0("sd_",metric)]] <- if(nrow(z)>1L) sd(z[[metric]]) else NA_real_
    }
    a
  }))
}
ff_paired <- function(res,cfg) {
  present <- unique(res$method)
  pairs <- list()
  if("FLOCK_MapMerge" %in% present)
    for(other in setdiff(present,"FLOCK_MapMerge"))
      pairs[[length(pairs)+1L]] <- c("FLOCK_MapMerge",other)
  if("FLOCK_Map" %in% present)
    for(other in intersect(c("FLOCK_Novelty","FLOCK_Unexplored","AdaptiveMixtureMH"),present))
      pairs[[length(pairs)+1L]] <- c("FLOCK_Map",other)
  metrics <- c("mean_rmse","covariance_relative_error","weight_l1",
               "rare_abs_error","elapsed_seconds")
  rows <- list(); idx <- 0L
  for(pair in pairs) for(metric in metrics) {
    a <- res[res$method==pair[1L],c("rep",metric),drop=FALSE]
    b <- res[res$method==pair[2L],c("rep",metric),drop=FALSE]
    z <- merge(a,b,by="rep",suffixes=c("_a","_b"))
    if(!nrow(z)) next
    delta <- z[[paste0(metric,"_a")]]-z[[paste0(metric,"_b")]]
    ci <- c(NA_real_,NA_real_)
    if(length(delta)>=2L) {
      set.seed(1000L+sum(utf8ToInt(paste(pair,metric,collapse="|"))))
      boot <- replicate(cfg$bootstrap,mean(delta[sample.int(length(delta),length(delta),TRUE)]))
      ci <- unname(quantile(boot,c(.025,.975)))
    }
    idx <- idx+1L
    rows[[idx]] <- data.frame(focal=pair[1L],competitor=pair[2L],metric=metric,n=length(delta),
       mean_difference=mean(delta),median_difference=median(delta),
       bootstrap_low=ci[1L],bootstrap_high=ci[2L],
       fraction_focal_better=mean(delta<0),stringsAsFactors=FALSE)
  }
  if(length(rows)) do.call(rbind,rows) else data.frame()
}
ff_discovery_summary <- function(discovery) {
  if(!nrow(discovery)) return(data.frame())
  z <- discovery
  out <- aggregate(as.numeric(z$visited),
                   by=z[c("method","source","stage","unit","region")],FUN=mean)
  names(out)[ncol(out)] <- "discovery_rate"
  out
}
ff_pdf <- function(res,summary,paired,discovery_summary,maps,audit,cfg,file) {
  pdf(file,width=11.7,height=8.3,onefile=TRUE)
  on.exit(dev.off(),add=TRUE)
  par(mar=c(5,5,4,2)+.1)
  plot.new()
  text(.5,.88,"FLOCK Experiment 5",cex=2,font=2)
  text(.5,.80,"10-dimensional, six-component Gaussian mixture",cex=1.2)
  lines <- c(paste("Version",cfg$version,"| mode:",cfg$mode,
                    "| target evaluations per method:",cfg$budget),
             paste("Completed replicate-method pairs:",nrow(res)),
             "Component cores and component groups are diagnostics, not density modes.",
             "Map vs MapMerge is the primary merging ablation.",
             "Negative paired differences favor the named focal method.",
             "Bootstrap intervals are exploratory and unadjusted for multiple comparisons.",
             "SAChES-inspired and AdaptiveMixtureMH are custom implementations.",
             "Equal target-call budgets do not imply equal computational overhead.")
  text(.08,seq(.66,.18,length.out=length(lines)),lines,adj=0,cex=.95)
  pretty <- function(x) gsub("_"," ",x,fixed=TRUE)
  order_methods <- cfg$methods[cfg$methods %in% unique(res$method)]
  metrics <- c(mean_rmse="Posterior mean RMSE",
    covariance_relative_error="Relative covariance error",
    weight_l1="Component-weight L1 error",
    quantile_standardized_rmse="Marginal-quantile standardized RMSE",
    elapsed_seconds="Sampling wall-clock seconds (diagnostics excluded)",
    rare_est="Estimated rare-component probability")
  for(metric in names(metrics)) {
    par(mar=c(5,13,4,2)+.1)
    groups <- lapply(order_methods,function(m) res[res$method==m,metric])
    names(groups) <- pretty(order_methods)
    boxplot(groups,horizontal=TRUE,las=1,main=metrics[[metric]],
            xlab=if(metric=="rare_est")"Reference probability = 0.02" else "Lower is better",
            cex.axis=.85)
    if(metric=="rare_est") abline(v=.02,lty=2)
  }
  z <- discovery_summary[discovery_summary$source=="chain_states" &
        discovery_summary$stage=="final_retained" &
        discovery_summary$unit=="component_core",,drop=FALSE]
  if(nrow(z)) {
    mat <- matrix(0,length(order_methods),6L,
                  dimnames=list(pretty(order_methods),LETTERS[1:6]))
    for(i in seq_along(order_methods)) for(j in seq_len(6L)) {
      v <- z$discovery_rate[z$method==order_methods[i] & z$region==LETTERS[j]]
      if(length(v)) mat[i,j] <- v[1L]
    }
    par(mar=c(5,13,4,2)+.1)
    barplot(mat,horiz=TRUE,beside=TRUE,las=1,xlim=c(0,1.32),
            main="Retained-chain component-core visits",xlab="Fraction of runs",
            legend.text=rownames(mat),args.legend=list(x="topright",cex=.7,bty="n"))
  }
  if(nrow(paired)) {
    z <- paired[paired$metric=="mean_rmse",,drop=FALSE]
    if(nrow(z)) {
      par(mar=c(5,17,4,2)+.1)
      xx <- c(z$mean_difference,z$bootstrap_low,z$bootstrap_high,0)
      xx <- xx[is.finite(xx)]; span <- range(xx)
      if(diff(span)==0) span <- span+c(-.01,.01)
      plot(z$mean_difference,seq_len(nrow(z)),xlim=span,
           ylim=c(.5,nrow(z)+.5),yaxt="n",ylab="",xlab="Mean error difference (focal - competitor)",
           main="Paired posterior-mean RMSE differences",pch=19)
      axis(2,at=seq_len(nrow(z)),labels=paste(pretty(z$focal),"-",pretty(z$competitor)),
           las=1,cex.axis=.72)
      ok <- is.finite(z$bootstrap_low)&is.finite(z$bootstrap_high)
      segments(z$bootstrap_low[ok],which(ok),z$bootstrap_high[ok],which(ok))
      abline(v=0,lty=2)
    }
  }
  if(nrow(maps)) {
    z <- unique(maps[c("rep","method","phase","regions_in_map")])
    z <- z[z$phase %in% c("scout_before","scout_after","final_before","final_after"),]
    par(mar=c(8,5,4,2)+.1)
    boxplot(regions_in_map~interaction(method,phase,sep=" / "),data=z,las=2,
            cex.axis=.7,ylab="Map regions",main="Map size before and after merging")
  }
  if(nrow(audit)) {
    par(mar=c(5,5,4,2)+.1)
    plot(audit$mahalanobis,audit$symmetric_kl,pch=ifelse(audit$merged,19,1),
         xlab="Mahalanobis center distance",ylab="Symmetric Gaussian KL",
         main="Merge audit: filled points accepted, open points rejected")
    abline(v=cfg$merge_mahal,h=cfg$merge_symkl,lty=2)
  }
  invisible(file)
}

# --- Self-tests, provenance, persistence and driver ----------------------------
ff_self_tests <- function() {
  old_exists <- exists(".Random.seed",envir=.GlobalEnv,inherits=FALSE)
  if(old_exists) old <- get(".Random.seed",envir=.GlobalEnv)
  on.exit(if(old_exists) assign(".Random.seed",old,envir=.GlobalEnv) else
    if(exists(".Random.seed",envir=.GlobalEnv,inherits=FALSE))
      rm(".Random.seed",envir=.GlobalEnv),add=TRUE)
  set.seed(47)
  # In particular, guard against transposing a non-Cholesky covariance root.
  S <- matrix(c(2,.8,.8,1),2,2); mu <- c(.3,-.4)
  g <- ff_gaussian(mu,S,1e-10); x <- c(.9,.2)
  analytic <- -log(2*pi)-.5*as.numeric(determinant(S,logarithm=TRUE)$modulus)-
    .5*as.numeric(crossprod(x-mu,solve(S,x-mu)))
  ff_assert(abs(ff_glog(x,g)-analytic)<1e-10,"Gaussian density self-test failed.")
  ff_assert(max(abs(crossprod(g$R)-S))<1e-10,"Cholesky convention self-test failed.")
  X <- matrix(replicate(8000L,ff_gdraw(g)),nrow=2L)
  ff_assert(max(abs(cov(t(X))-S))<.12,"Gaussian draw covariance self-test failed.")
  # Verify the actual q mixture used in final MH, including broad support.
  cfg <- list(d=2L,eigen_floor=1e-10,broad_sd=5,broad_weight=.25)
  regs <- list(list(mu=mu,S=S),list(mu=c(-1,2),S=diag(c(.5,1))))
  q <- ff_q_from_regions(regs,cfg)
  Y <- t(replicate(12000L,q$draw()))
  qmean <- Reduce("+",Map(function(a,w) w*a$mu,q$gaussians,as.list(q$weights)))
  qcov <- Reduce("+",Map(function(a,w) w*(a$S+tcrossprod(a$mu-qmean)),
                          q$gaussians,as.list(q$weights)))
  ff_assert(max(abs(colMeans(Y)-qmean))<.14,"Mixture draw mean self-test failed.")
  ff_assert(max(abs(cov(Y)-qcov))/max(diag(qcov))<.12,
            "Mixture draw covariance self-test failed.")
  # Pointwise independence-MH flux identity, using the same evaluable q.
  a <- c(.2,.5); b <- c(-.7,1.1)
  lpa <- ff_glog(a,g); lpb <- ff_glog(b,g)
  lqa <- q$logq(a); lqb <- q$logq(b)
  forward <- lpa+lqb+min(0,lpb-lpa+lqa-lqb)
  reverse <- lpb+lqa+min(0,lpa-lpb+lqb-lqa)
  ff_assert(abs(forward-reverse)<1e-10,"MH flux self-test failed.")
  # Merge geometry must discriminate overlapping and distant synthetic groups.
  cfg$d <- 2L
  r1 <- list(mu=c(0,0),S=diag(2)); r2 <- list(mu=c(.2,.2),S=diag(2))
  r3 <- list(mu=c(8,8),S=diag(2))
  ff_assert(ff_pair(r1,r2,cfg)$mahal<1.5 && ff_pair(r1,r3,cfg)$mahal>1.5,
            "Merge geometry self-test failed.")
  ff_assert(identical(ff_pick(7L),7L),"Singleton candidate selection failed.")
  ff_assert(sum(ff_split(103L,8L))==103L,"Budget allocation self-test failed.")
  cat("Built-in numerical self-tests passed.\n")
  invisible(TRUE)
}
ff_save_rds <- function(object,path) {
  tmp <- paste0(path,".tmp")
  saveRDS(object,tmp,version=2)
  if(file.exists(path)) {
    file.copy(path,paste0(path,".bak"),overwrite=TRUE)
    unlink(path)
  }
  ff_assert(file.rename(tmp,path),paste("Could not write",path))
  invisible(path)
}
ff_script_path <- function() {
  frames <- sys.frames()
  files <- unlist(lapply(frames,function(f) if(!is.null(f$ofile)) f$ofile else NULL),
                  use.names=FALSE)
  if(length(files)) return(normalizePath(tail(files,1L),mustWork=FALSE))
  arg <- grep("^--file=",commandArgs(FALSE),value=TRUE)
  if(length(arg)) return(normalizePath(sub("^--file=","",arg[1L]),mustWork=FALSE))
  NULL
}
ff_write_csv <- function(x,file) {
  if(!ncol(x)) x <- data.frame(note="No rows were produced for this diagnostic.")
  write.csv(x,file,row.names=FALSE,na="")
}
ff_export <- function(bundles,cfg,plots=FALSE) {
  if(!length(bundles)) return(invisible(NULL))
  res <- ff_bind(bundles,"result"); summary <- ff_summary(res)
  discovery <- ff_bind(bundles,"discovery")
  ds <- ff_discovery_summary(discovery); maps <- ff_bind(bundles,"maps")
  audit <- ff_bind(bundles,"merge_audit")
  tab <- list(results=res,summary=summary,component_discovery=discovery,
              discovery_summary=ds,first_evaluated=ff_bind(bundles,"first_evaluated"),
              budget_audit=ff_bind(bundles,"budgets"),map_summary=maps,
              merge_audit=audit,merge_reference_diagnostics=ff_bind(bundles,"merge_diagnostics"),
              follower_assignments=ff_bind(bundles,"allocation"))
  if(plots) tab$paired <- ff_paired(res,cfg)
  for(name in names(tab)) ff_write_csv(tab[[name]],
      file.path(cfg$outdir,paste0("flock_experiment_5_",name,".csv")))
if(plots) ff_pdf(res,summary,tab$paired,ds,maps,audit,cfg,
                   file.path(cfg$outdir,"flock_experiment_5_diagnostics.pdf"))
  invisible(tab)
}
ff_run_experiment <- function() {
  cfg <- ff_config()
  ff_assert(cfg$seed+100003*cfg$reps<2^53,"Seed range too large.")
  dir.create(cfg$outdir,recursive=TRUE,showWarnings=FALSE)
  checkpoint_dir <- file.path(cfg$outdir,"checkpoints")
  dir.create(checkpoint_dir,showWarnings=FALSE)
  script <- ff_script_path()
  algorithm_config <- cfg[setdiff(names(cfg),c("reps","methods","mode","outdir",
                                              "resume","save_draws","bootstrap"))]
  signature <- list(config=algorithm_config,
      script_md5=if(!is.null(script) && file.exists(script)) unname(tools::md5sum(script)) else NULL,
      r_version=paste(R.version$major,R.version$minor,sep="."),rng=RNGkind())
  manifest_path <- file.path(cfg$outdir,"manifest.rds")
  if(file.exists(manifest_path)) {
    old <- readRDS(manifest_path)
    ff_assert(cfg$resume,
      "Output directory already has a run. Choose a new OUTDIR or use RESUME=TRUE.")
    ff_assert(identical(old$signature,signature),
      "Settings, R version, or script changed. Use a NEW OUTDIR; do not mix results.")
  } else {
    ff_assert(!length(list.files(checkpoint_dir,pattern="\\.rds$")),
              "Checkpoints exist without a manifest. Use a new output directory.")
    ff_save_rds(list(signature=signature,created=as.character(Sys.time()),
                    original_config=cfg),manifest_path)
  }
  cat("\nFLOCK Experiment 5 |",cfg$mode,"|",cfg$reps,"replicates | B =",cfg$budget,"\n")
  cat("Methods:",paste(cfg$methods,collapse=", "),"\n")
  cat("Output:",normalizePath(cfg$outdir),"\n")
  cat("Budget: target calls, not iterations. Reference diagnostics are post-hoc.\n")
  if(nzchar(Sys.getenv("REPS")) || nzchar(Sys.getenv("BUDGET")))
    cat("NOTE: environment REPS/BUDGET overrides are active; settings above are authoritative.\n")
  ff_self_tests()
  target <- ff_target()
  ff_save_rds(target[setdiff(names(target),"log_density")],
              file.path(cfg$outdir,"reference_target.rds"))
  capture.output(dput(cfg),file=file.path(cfg$outdir,"settings.txt"))
  capture.output(sessionInfo(),file=file.path(cfg$outdir,"session_info.txt"))
  if(!is.null(script) && file.exists(script))
    file.copy(script,file.path(cfg$outdir,"script_used.R"),overwrite=TRUE)
  description <- c(
    "FLOCK Experiment 5: interpretation notes",
    "",
    "Six Gaussian components with weights .30,.25,.20,.15,.08,.02 in 10 dimensions.",
    "Groups AB, C, DE, F are prespecified reference component groups, NOT proven density basins.",
    "Component core means squared Mahalanobis distance <= qchisq(.95,10).",
    "Cores can overlap. A visited core is not proof of accurate probability-mass recovery.",
    "Component weight estimates average exact responsibilities; their reference means are w.",
    "Sampling methods receive only the log density, common dispersed starts and tuning settings.",
    "Scoring labels, analytical moments and component locations are never used by samplers.",
    "",
    "Map stores provisional clusters of capped exploratory histories, including repeated states.",
    "n_unique is NOT ESS. Allocation uses density-penalized inverse-square-root unique count.",
    "Map vs MapMerge: same scouts, middle budgets, pioneer policy, seeds and proposal family.",
    "Both reserve 2% of B inside preparation. Bridge checks are charged; remainder is MH settling.",
    "All cumulative map construction is at stage boundaries; middle jobs never communicate.",
    "Merge uses Mahalanobis <=1.5, symmetrized KL <=1, sufficient unique points, and optional bridge.",
    "The bridge is only a straight-line check, not a connectivity theorem.",
    "Different dominant reference groups in merge diagnostics are flags, not true error labels.",
    "",
    "AdaptiveMixtureMH: custom dispersed AM explorers -> shared frozen-mixture final MH.",
    "SAChES_inspired: custom DE ensemble -> locally adaptive chains, frozen for retained sampling.",
    "Neither is a validated reproduction of AIMM or SAChES; no published-method superiority claim.",
    "ParallelTempering: 8 fixed inverse temperatures from 1 to .03, adjacent swaps, warmup scale tuning.",
    "RWMH_4chains: 4 common dispersed starts, step 2.38/sqrt(d), 35% per-chain burn-in.",
    "Acceptance rates refer to different phases across methods; they are not an accuracy ranking.",
    "",
    "Checkpoints save each completed replicate-method pair. Resume skips these pairs.",
    "First-replicate draws are saved by default; SAVE_DRAWS=TRUE retains every run.",
    "Paired bootstrap intervals resample replicate pairs and are exploratory, not multiplicity-adjusted.",
    "An easy map demo and passing self-tests do not validate finite-budget posterior recovery.",
    "Do not import the earlier experiment results into this comparison; target and implementation differ.")
  writeLines(description,file.path(cfg$outdir,"INTERPRETATION.txt"))
  bundles <- list(); idx <- 0L
  for(r in seq_len(cfg$reps)) {
    ff_seed(cfg,r,1L)
    init <- matrix(rnorm(8L*cfg$d,sd=cfg$initial_sd),8L,cfg$d)
    # Rotate method order by replicate to reduce systematic execution-order effects.
    mm <- cfg$methods
    offset <- (r-1L)%%length(mm)
    if(offset>0L) mm <- c(mm[(offset+1L):length(mm)],mm[seq_len(offset)])
    for(method in mm) {
      path <- file.path(checkpoint_dir,sprintf("rep_%04d_%s.rds",r,method))
      if(cfg$resume && file.exists(path)) {
        bundle <- readRDS(path)
        ff_assert(identical(bundle$signature,signature),"Checkpoint signature mismatch.")
        cat(sprintf("[%d/%d] %s: resumed\n",r,cfg$reps,method))
      } else {
        cat(sprintf("[%d/%d] %s ... ",r,cfg$reps,method)); flush.console()
        bundle <- tryCatch({
          ans <- ff_run_method(method,target$log_density,init,cfg,r)
          scored <- ff_score(ans,target,cfg,r,method)
          scored$signature <- signature
          scored
        },error=function(e) {
          error_file <- file.path(cfg$outdir,sprintf("ERROR_rep_%04d_%s.txt",r,method))
          writeLines(c(conditionMessage(e),"",capture.output(sessionInfo())),error_file)
          stop(paste("Failed:",method,"replicate",r,"-",conditionMessage(e),
                     "Completed checkpoints are preserved. Error file:",error_file),call.=FALSE)
        })
        ff_save_rds(bundle,path)
        cat(sprintf("RMSE %.4g | L1 %.4g | %d/%d calls | %.2fs sampling\n",
                    bundle$result$mean_rmse,bundle$result$weight_l1,
                    bundle$result$evals,cfg$budget,bundle$result$elapsed_seconds))
      }
      idx <- idx+1L; bundles[[idx]] <- bundle
    }
    # Update readable CSVs after every complete replicate; binary checkpoints
    # already protect each individual method if interruption occurs sooner.
    ff_export(bundles,cfg,plots=FALSE)
  }
  tables <- ff_export(bundles,cfg,plots=TRUE)
  cat("\nCompleted. Summary (lower errors are better):\n")
  print(tables$summary[,c("method","n","mean_mean_rmse","mean_weight_l1",
                          "mean_covariance_relative_error","rare_core_discovery_rate",
                          "mean_merges","mean_elapsed_seconds")],row.names=FALSE,digits=4)
  cat("\nOutputs saved in:",normalizePath(cfg$outdir),"\n")
  cat("Read INTERPRETATION.txt before interpreting rankings. No automatic winner declaration.\n")
  invisible(tables)
}

if (ff_bool("FLOCK_AUTORUN",TRUE)) ff_run_experiment()