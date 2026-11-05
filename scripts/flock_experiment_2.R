# FLOCK: paired ablation experiment (100 replications by default)
# Run: source("flock_experiment_2.R")
# Optional: Sys.setenv(REPS="20", BUDGET="20000") for a smoke test.
# Requires base R only. Does not require previous scripts to be sourced.
# Note: all variants use same scout seed within each replicate. The
# FLOCK-Random variant changes only candidate choice, not candidate count.
# The reference mixture parameters are used ONLY for diagnostics.

# ARPM-MCMC pilot: fixed-budget comparison against RWMH and Adaptive Metropolis
# Base R only. Run: Rscript arpm_pilot.R
# All methods use the same number of target-density evaluations per replicate.
# This is a heuristic pilot, not evidence of algorithmic superiority.

set.seed(20261024)
D <- 5L
weights <- c(.45,.30,.20,.05)
means <- rbind(c(-6,-6,0,0,0), c(6,6,0,0,0),
               c(-6,6,4,0,0), c(6,-6,-4,0,0))
scales <- c(.75,1.7,1.0,.35)
Sigmas <- lapply(scales, function(s) diag(s^2,D))
# Component 2 is broad, component 4 is narrow and rare.
chols <- lapply(Sigmas, chol)
logsumexp <- function(z) { a <- max(z); a + log(sum(exp(z-a))) }
log_target <- function(x) {
  lp <- numeric(length(weights))
  for (j in seq_along(weights)) {
    z <- backsolve(chols[[j]], x-means[j,], transpose=TRUE)
    lp[j] <- log(weights[j]) - D/2*log(2*pi) -
      sum(log(diag(chols[[j]]))) - sum(z*z)/2
  }
  logsumexp(lp)
}
# Known truth; regions are nearest mean under component-specific Mahalanobis
# distance. For mode-mass diagnostics, use a soft responsibility instead of
# assuming these hard regions have exact mixture weights.
responsibility <- function(x) {
  lp <- numeric(4)
  for (j in 1:4) {
    z <- backsolve(chols[[j]], x-means[j,], transpose=TRUE)
    lp[j] <- log(weights[j]) - D/2*log(2*pi) -
      sum(log(diag(chols[[j]]))) - sum(z*z)/2
  }
  exp(lp-logsumexp(lp))
}
truth_mean <- colSums(means * weights)
truth_second <- colSums((means^2 + matrix(scales^2,4,D)) * weights)
truth_var <- truth_second - truth_mean^2

# All target calls are counted by explicit wrappers.
make_counter <- function() {
  e <- new.env(parent=emptyenv()); e$n <- 0L
  e$eval <- function(x) { e$n <- e$n+1L; log_target(x) }
  e$count <- function() e$n
  e
}
rmv <- function(mu, L) as.numeric(mu + t(L)%*%rnorm(length(mu)))
# Cholesky of positive definite covariance, with regularization.
safe_chol <- function(S, floor=.03) {
  S <- (S+t(S))/2
  ee <- eigen(S, symmetric=TRUE)
  ee$vectors %*% diag(sqrt(pmax(ee$values,floor)),nrow(S))
}
# Proposal q(x): mixture of Gaussian local kernels and a broad Gaussian.
# q must be evaluable exactly for independence MH.
log_mvnorm <- function(x,mu,L) {
  z <- solve(L, x-mu)
  -D/2*log(2*pi)-sum(log(abs(diag(qr.R(qr(L))))))-sum(z*z)/2
}
make_proposal <- function(centers, covs, masses, broad_weight=.15) {
  centers <- rbind(centers, rep(0,D))
  covs[[length(covs)+1L]] <- diag(12^2,D)
  masses <- c((1-broad_weight)*masses/sum(masses), broad_weight)
  Ls <- lapply(covs, safe_chol)
  logq <- function(x) {
    terms <- vapply(seq_along(masses), function(j)
      log(masses[j])+log_mvnorm(x,centers[j,],Ls[[j]]), numeric(1))
    logsumexp(terms)
  }
  draw <- function() {
    j <- sample.int(length(masses),1,prob=masses)
    rmv(centers[j,],Ls[[j]])
  }
  list(logq=logq,draw=draw)
}
# RWMH baseline. B target evaluations, including initialization.
rwmh <- function(B, x0, step=1.8) {
  f <- make_counter(); x <- x0; lp <- f$eval(x)
  out <- matrix(NA_real_,B,D); out[1,] <- x; accepts <- 0L
  for (t in 2:B) {
    y <- x+step*rnorm(D); ly <- f$eval(y)
    if (log(runif(1)) < ly-lp) { x<-y; lp<-ly; accepts<-accepts+1L }
    out[t,] <- x
  }
  list(draws=out,evals=f$count(),accept=accepts/(B-1))
}
# Adaptive Metropolis: covariance updated from past chain states during warmup,
# then frozen. Adaptation is NOT used in retained estimation samples.
am <- function(B,x0,warmup,initial_step=1.8,update_every=100L) {
  f<-make_counter(); x<-x0; lp<-f$eval(x)
  out<-matrix(NA_real_,B,D); out[1,]<-x; accepts<-0L
  L<-diag(initial_step,D)
  for(t in 2:B) {
    if(t>max(100L,D*20L) && t<=warmup && t%%update_every==0L) {
      S<-cov(out[1:(t-1),,drop=FALSE]); L<-safe_chol((2.38^2/D)*S+.01*diag(D))
    }
    y<-rmv(x,L); ly<-f$eval(y)
    if(log(runif(1))<ly-lp) {x<-y;lp<-ly;accepts<-accepts+1L}
    out[t,]<-x
  }
  list(draws=out,evals=f$count(),accept=accepts/(B-1))
}

# SURVEY CORPS: sequential, role-specialized exploration + corrected final MH.
# Base R only. Does NOT use the true mixture parameters to choose proposals.
# Run: source('survey_corps.R') or Rscript survey_corps.R
# Environment overrides: BUDGET=20000 REPS=20
# NOTE: target density is a known Gaussian mixture solely for evaluation.

# Counted MH local exploration; every target evaluation is charged.
local_explore <- function(f, x, lp, n, sd, history) {
  acc <- 0L
  for (j in seq_len(n)) {
    y <- x + rnorm(D, sd=sd)
    ly <- f$eval(y)
    if (log(runif(1)) < ly-lp) {x <- y; lp <- ly; acc <- acc+1L}
    history <- rbind(history,x)
  }
  list(x=x, lp=lp, history=history, accept=acc)
}

# Coverage is relative to previously visited points, measured on a
# standardized scale. This is a search heuristic, NOT a target distribution.
coverage_distance <- function(candidates, visited, scale_vec) {
  # Cap historical pool for computational practicality; deterministic thinning.
  visited <- visited[unique(round(seq(1,nrow(visited),length.out=min(350L,nrow(visited))))),,drop=FALSE]
  v <- sweep(visited,2,scale_vec,'/')
  vapply(seq_len(nrow(candidates)),function(i) {
    z <- candidates[i,]/scale_vec
    sqrt(min(rowSums((v-matrix(z,nrow(v),D,byrow=TRUE))^2)))
  },numeric(1))
}

# Cluster visits to build a frozen proposal, with a heavy broad component.
# Pooling uses equal subsampling per role to avoid domination by followers.
learn_survey_proposal <- function(histories, max_clusters=8L) {
  pool <- do.call(rbind,lapply(histories,function(z) {
    ids <- unique(round(seq(1,nrow(z),length.out=min(180L,nrow(z)))))
    z[ids,,drop=FALSE]
  }))
  pool <- unique(pool)
  nc <- min(max_clusters,nrow(pool)-1L)
  stopifnot(nc>=1L)
  km <- kmeans(pool,centers=nc,nstart=5,iter.max=60)
  covs <- lapply(seq_len(nc),function(j) {
    z <- pool[km$cluster==j,,drop=FALSE]
    if(nrow(z)<D+2L) diag(2,D) else cov(z)+diag(.25,D)
  })
  make_proposal(km$centers,covs,rep(1,nc),broad_weight=.25)
}

survey_corps <- function(B,x0,scout_frac=.20,middle_frac=.45,
                         final_warmup_frac=.10, pioneer=TRUE,
                         n_scout=2L,n_follower=3L,n_pioneer=3L) {
  stopifnot(B>=200L,scout_frac>0,middle_frac>0,
            scout_frac+middle_frac<.85)
  f <- make_counter()
  scout_budget <- floor(B*scout_frac)
  middle_budget <- floor(B*middle_frac)
  # Scouts first: larger step sizes and widely dispersed starts.
  starts <- matrix(rnorm(n_scout*D,sd=7),n_scout,D)
  starts[1,] <- x0
  scout_hist <- vector('list',n_scout)
  scout_end <- matrix(NA_real_,n_scout,D)
  scout_lp <- numeric(n_scout)
  scout_steps <- rep(5,n_scout)
  scout_n <- rep((scout_budget-n_scout)%/%n_scout,n_scout)
  if (sum(scout_n)<scout_budget-n_scout)
    scout_n[seq_len(scout_budget-n_scout-sum(scout_n))] <-
      scout_n[seq_len(scout_budget-n_scout-sum(scout_n))]+1L
  for(k in seq_len(n_scout)) {
    lp <- f$eval(starts[k,]); hist <- matrix(starts[k,],1,D)
    z <- local_explore(f,starts[k,],lp,scout_n[k],scout_steps[k],hist)
    scout_hist[[k]] <- z$history
    scout_end[k,] <- z$x; scout_lp[k] <- z$lp
  }
  stopifnot(f$count()==scout_budget)
  scout_pool <- do.call(rbind,scout_hist)
  # Select distinct high-density scout discoveries for followers.
  # Ranking uses cached/retrospective target calls? Avoid extra target calls:
  # track log density at candidate points via counted evaluations instead.
  # We evaluate a modest fixed number of scout candidates and charge the budget.
  # Here followers start from scout endpoints (already evaluated), ordered by lp.
  # Remaining followers start from jittered scout endpoints (new counted calls).
  follow_hist <- list(); pioneer_hist <- list()
  n_mid <- n_follower + n_pioneer
  allocations <- rep(middle_budget %/% n_mid,n_mid)
  if(middle_budget %% n_mid > 0L)
    allocations[seq_len(middle_budget %% n_mid)] <- allocations[seq_len(middle_budget %% n_mid)]+1L
  best_order <- order(scout_lp,decreasing=TRUE)
  for(k in seq_len(n_follower)) {
    idx <- best_order[(k-1L) %% n_scout + 1L]
    x <- scout_end[idx,]; lp <- scout_lp[idx]
    # Follow existing high-density locations; local moves map covariance.
    if(k>n_scout) {
      x <- x + rnorm(D,sd=.7)
      lp <- f$eval(x); n <- allocations[k]-1L
    } else n <- allocations[k]
    z <- local_explore(f,x,lp,n,sd=1.15,matrix(x,1,D))
    follow_hist[[k]] <- z$history
  }
  # Pioneers start in places missed by scouts. Candidate generation is cheap;
  # candidate scoring requires target evaluations and is explicitly budgeted.
  # With pioneer=FALSE, the same slots use follower-style exploration.
  for(k in seq_len(n_pioneer)) {
    slot <- n_follower+k
    n <- allocations[slot]
    if(pioneer) {
      n_candidates <- min(12L,max(2L,n %/% 20L))
      cand <- matrix(rnorm(n_candidates*D,sd=8),n_candidates,D)
      lp_cand <- apply(cand,1,f$eval)
      scale_vec <- pmax(apply(scout_pool,2,sd),2)
      novelty <- coverage_distance(cand,scout_pool,scale_vec)
      # Density-aware novelty: reject wildly implausible candidates via
      # a capped density penalty. Heuristic only, never used for inference.
      utility <- novelty - .15*pmax(0,max(lp_cand)-lp_cand)
      j <- which.max(utility)
      x <- cand[j,]; lp <- lp_cand[j]
      n <- n-n_candidates
      z <- local_explore(f,x,lp,n,sd=2.3,matrix(x,1,D))
    } else {
      idx <- best_order[(k-1L) %% n_scout + 1L]
      x <- scout_end[idx,]; lp <- scout_lp[idx]
      z <- local_explore(f,x,lp,n,sd=1.15,matrix(x,1,D))
    }
    pioneer_hist[[k]] <- z$history
  }
  stopifnot(f$count()==scout_budget+middle_budget)
  q <- learn_survey_proposal(c(scout_hist,follow_hist,pioneer_hist))
  # Final stage: frozen mixture proposal and corrected independence MH.
  # Starting point already evaluated by a scout, so no uncounted calls.
  x <- scout_end[best_order[1],]; lp <- scout_lp[best_order[1]]
  lq <- q$logq(x)
  remain <- B-f$count()
  out <- matrix(NA_real_,remain+1L,D); out[1,] <- x
  acc <- 0L
  for(t in seq_len(remain)) {
    is_local <- runif(1)<.2
    if(is_local) {
      y <- x+rnorm(D,sd=.8); ly <- f$eval(y)
      loga <- ly-lp
      lqy <- q$logq(y)
    } else {
      y <- q$draw(); ly <- f$eval(y); lqy <- q$logq(y)
      loga <- ly-lp+lq-lqy
    }
    if(log(runif(1))<loga) {x<-y;lp<-ly;lq<-lqy;acc<-acc+1L}
    out[t+1L,] <- x
  }
  stopifnot(f$count()==B)
  burn <- floor(nrow(out)*final_warmup_frac)
  list(draws=out[(burn+1L):nrow(out),,drop=FALSE],evals=f$count(),
       accept=acc/remain,scout_evals=scout_budget,middle_evals=middle_budget,
       final_evals=remain, proposal=q)
}


# Control: candidate points scored exactly as in FLOCK, but selection is random.
flock_random <- function(B,x0,scout_frac=.20,middle_frac=.45,
                         final_warmup_frac=.10, pioneer=TRUE,
                         n_scout=2L,n_follower=3L,n_pioneer=3L) {
  stopifnot(B>=200L,scout_frac>0,middle_frac>0,
            scout_frac+middle_frac<.85)
  f <- make_counter()
  scout_budget <- floor(B*scout_frac)
  middle_budget <- floor(B*middle_frac)
  # Scouts first: larger step sizes and widely dispersed starts.
  starts <- matrix(rnorm(n_scout*D,sd=7),n_scout,D)
  starts[1,] <- x0
  scout_hist <- vector('list',n_scout)
  scout_end <- matrix(NA_real_,n_scout,D)
  scout_lp <- numeric(n_scout)
  scout_steps <- rep(5,n_scout)
  scout_n <- rep((scout_budget-n_scout)%/%n_scout,n_scout)
  if (sum(scout_n)<scout_budget-n_scout)
    scout_n[seq_len(scout_budget-n_scout-sum(scout_n))] <-
      scout_n[seq_len(scout_budget-n_scout-sum(scout_n))]+1L
  for(k in seq_len(n_scout)) {
    lp <- f$eval(starts[k,]); hist <- matrix(starts[k,],1,D)
    z <- local_explore(f,starts[k,],lp,scout_n[k],scout_steps[k],hist)
    scout_hist[[k]] <- z$history
    scout_end[k,] <- z$x; scout_lp[k] <- z$lp
  }
  stopifnot(f$count()==scout_budget)
  scout_pool <- do.call(rbind,scout_hist)
  # Select distinct high-density scout discoveries for followers.
  # Ranking uses cached/retrospective target calls? Avoid extra target calls:
  # track log density at candidate points via counted evaluations instead.
  # We evaluate a modest fixed number of scout candidates and charge the budget.
  # Here followers start from scout endpoints (already evaluated), ordered by lp.
  # Remaining followers start from jittered scout endpoints (new counted calls).
  follow_hist <- list(); pioneer_hist <- list()
  n_mid <- n_follower + n_pioneer
  allocations <- rep(middle_budget %/% n_mid,n_mid)
  if(middle_budget %% n_mid > 0L)
    allocations[seq_len(middle_budget %% n_mid)] <- allocations[seq_len(middle_budget %% n_mid)]+1L
  best_order <- order(scout_lp,decreasing=TRUE)
  for(k in seq_len(n_follower)) {
    idx <- best_order[(k-1L) %% n_scout + 1L]
    x <- scout_end[idx,]; lp <- scout_lp[idx]
    # Follow existing high-density locations; local moves map covariance.
    if(k>n_scout) {
      x <- x + rnorm(D,sd=.7)
      lp <- f$eval(x); n <- allocations[k]-1L
    } else n <- allocations[k]
    z <- local_explore(f,x,lp,n,sd=1.15,matrix(x,1,D))
    follow_hist[[k]] <- z$history
  }
  # Pioneers start in places missed by scouts. Candidate generation is cheap;
  # candidate scoring requires target evaluations and is explicitly budgeted.
  # With pioneer=FALSE, the same slots use follower-style exploration.
  for(k in seq_len(n_pioneer)) {
    slot <- n_follower+k
    n <- allocations[slot]
    if(pioneer) {
      n_candidates <- min(12L,max(2L,n %/% 20L))
      cand <- matrix(rnorm(n_candidates*D,sd=8),n_candidates,D)
      lp_cand <- apply(cand,1,f$eval)
      scale_vec <- pmax(apply(scout_pool,2,sd),2)
      novelty <- coverage_distance(cand,scout_pool,scale_vec)
      # Density-aware novelty: reject wildly implausible candidates via
      # a capped density penalty. Heuristic only, never used for inference.
      utility <- novelty - .15*pmax(0,max(lp_cand)-lp_cand)
      j <- sample.int(n_candidates,1L)
      x <- cand[j,]; lp <- lp_cand[j]
      n <- n-n_candidates
      z <- local_explore(f,x,lp,n,sd=2.3,matrix(x,1,D))
    } else {
      idx <- best_order[(k-1L) %% n_scout + 1L]
      x <- scout_end[idx,]; lp <- scout_lp[idx]
      z <- local_explore(f,x,lp,n,sd=1.15,matrix(x,1,D))
    }
    pioneer_hist[[k]] <- z$history
  }
  stopifnot(f$count()==scout_budget+middle_budget)
  q <- learn_survey_proposal(c(scout_hist,follow_hist,pioneer_hist))
  # Final stage: frozen mixture proposal and corrected independence MH.
  # Starting point already evaluated by a scout, so no uncounted calls.
  x <- scout_end[best_order[1],]; lp <- scout_lp[best_order[1]]
  lq <- q$logq(x)
  remain <- B-f$count()
  out <- matrix(NA_real_,remain+1L,D); out[1,] <- x
  acc <- 0L
  for(t in seq_len(remain)) {
    is_local <- runif(1)<.2
    if(is_local) {
      y <- x+rnorm(D,sd=.8); ly <- f$eval(y)
      loga <- ly-lp
      lqy <- q$logq(y)
    } else {
      y <- q$draw(); ly <- f$eval(y); lqy <- q$logq(y)
      loga <- ly-lp+lq-lqy
    }
    if(log(runif(1))<loga) {x<-y;lp<-ly;lq<-lqy;acc<-acc+1L}
    out[t+1L,] <- x
  }
  stopifnot(f$count()==B)
  burn <- floor(nrow(out)*final_warmup_frac)
  list(draws=out[(burn+1L):nrow(out),,drop=FALSE],evals=f$count(),
       accept=acc/remain,scout_evals=scout_budget,middle_evals=middle_budget,
       final_evals=remain, proposal=q)
}



B <- as.integer(Sys.getenv('BUDGET', '20000'))
R <- as.integer(Sys.getenv('REPS', '100'))
stopifnot(B >= 2000L, R >= 2L)

# Soft component responsibilities are evaluation-only diagnostics. This
# estimator is unbiased in expectation under exact posterior draws.
measure <- function(obj, replicate, method) {
  z <- obj$draws
  rr <- vapply(seq_len(nrow(z)), function(i) responsibility(z[i,]), numeric(4))
  w <- rowMeans(rr)
  mu <- colMeans(z)
  data.frame(rep=replicate, method=method, evals=obj$evals,
             accept=obj$accept, retained=nrow(z),
             mean_rmse=sqrt(mean((mu-truth_mean)^2)),
             weight_l1=sum(abs(w-weights)),
             rare_est=w[4], rare_abs_error=abs(w[4]-weights[4]),
             rare_detected=as.integer(w[4] > 0.005),
             stringsAsFactors=FALSE)
}

rows <- vector('list', R*3L)
methods <- c('FLOCK','FLOCK_Random','FLOCK_Followers')
for (r in seq_len(R)) {
  seed <- 900000L + r
  set.seed(seed)
  x0 <- rnorm(D,sd=7)
  # Common random numbers: each variant begins with the same RNG state.
  # Differences in downstream random-number consumption are expected.
  set.seed(seed+100000L)
  full <- survey_corps(B,x0,pioneer=TRUE)
  set.seed(seed+100000L)
  rand <- flock_random(B,x0,pioneer=TRUE)
  set.seed(seed+100000L)
  follow <- survey_corps(B,x0,pioneer=FALSE)
  objs <- list(full,rand,follow)
  for (k in 1:3) rows[[(r-1L)*3L+k]] <- measure(objs[[k]],r,methods[k])
  if (r %% 10L == 0L || r == R) cat('Finished',r,'of',R,'replicates\n')
}
res <- do.call(rbind,rows)
stopifnot(all(res$evals == B))
write.csv(res,'flock_experiment_2_results.csv',row.names=FALSE)

summarize <- function(df) {
  do.call(rbind,lapply(split(df,df$method),function(z) {
    data.frame(method=z$method[1], n=nrow(z),
               mean_rmse=mean(z$mean_rmse), median_rmse=median(z$mean_rmse),
               sd_rmse=sd(z$mean_rmse), mean_weight_l1=mean(z$weight_l1),
               mean_rare_est=mean(z$rare_est), rare_detection_rate=mean(z$rare_detected),
               mean_rare_error=mean(z$rare_abs_error),
               mean_accept=mean(z$accept))
  }))
}
summary <- summarize(res)
write.csv(summary,'flock_experiment_2_summary.csv',row.names=FALSE)
print(summary,row.names=FALSE,digits=4)

# Paired differences: negative means FLOCK improves the error metric.
# Normal approximation and percentile bootstrap confidence intervals are
# descriptive; bootstrap resamples whole replicate pairs.
paired <- function(metric,competitor,bootstrap=2000L) {
  a <- res[res$method=='FLOCK',c('rep',metric)]
  b <- res[res$method==competitor,c('rep',metric)]
  zz <- merge(a,b,by='rep',suffixes=c('_flock','_control'))
  d <- zz[[paste0(metric,'_flock')]]-zz[[paste0(metric,'_control')]]
  set.seed(12345L+length(d)+nchar(metric)+nchar(competitor))
  bs <- replicate(bootstrap,mean(sample(d,length(d),replace=TRUE)))
  data.frame(metric=metric, competitor=competitor, n=length(d),
             mean_difference=mean(d), median_difference=median(d),
             bootstrap_low=unname(quantile(bs,.025)),
             bootstrap_high=unname(quantile(bs,.975)),
             fraction_flock_better=mean(d<0))
}
paired_table <- do.call(rbind,lapply(c('FLOCK_Random','FLOCK_Followers'),function(other)
  do.call(rbind,lapply(c('mean_rmse','weight_l1','rare_abs_error'),paired,
                       competitor=other))))
write.csv(paired_table,'flock_experiment_2_paired.csv',row.names=FALSE)
cat('\nPaired comparisons (negative mean difference favors FLOCK):\n')
print(paired_table,row.names=FALSE,digits=4)

pdf('flock_experiment_2_diagnostics.pdf',width=12,height=8)
par(mfrow=c(2,2),mar=c(7,4,3,1))
boxplot(mean_rmse~method,data=res,las=2,main='Posterior mean RMSE',ylab='RMSE')
boxplot(weight_l1~method,data=res,las=2,main='Mixture weight L1 error',ylab='L1 error')
boxplot(rare_est~method,data=res,las=2,main='Rare component estimated weight',ylab='Weight')
abline(h=.05,col='red',lty=2)
a <- res[res$method=='FLOCK',c('rep','mean_rmse')]
b <- res[res$method=='FLOCK_Random',c('rep','mean_rmse')]
zz <- merge(a,b,by='rep',suffixes=c('_flock','_random'))
plot(zz$mean_rmse_random,zz$mean_rmse_flock,
     xlab='FLOCK-Random RMSE',ylab='FLOCK RMSE',
     main='Paired RMSE: below line favors FLOCK',pch=16,cex=.7)
abline(0,1,col='red',lty=2)
dev.off()
cat('\nSaved flock_experiment_2_results.csv, flock_experiment_2_summary.csv,\n',
    'flock_experiment_2_paired.csv, flock_experiment_2_diagnostics.pdf\n')