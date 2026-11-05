# FLOCK: role-specialized exploration followed by importance consolidation.
# Standalone base-R research prototype; no supplementary package required.
#
# IMPORTANT: final samples are weighted importance samples, NOT one MCMC chain.
# Aggressive exploratory trajectories are NOT ordinary posterior draws.
# All methods use the same target-density budget and the same final proposal
# fitting / importance-sampling algorithms. Only exploration is different.
#
# Usage:
#   Sys.setenv(RUN_MODE='smoke', TARGET='5d'); source('flock_importance_consolidation.R')
#   Sys.setenv(RUN_MODE='pilot', TARGET='5d'); source('flock_importance_consolidation.R')
#   Sys.setenv(RUN_MODE='full', TARGET='5d');  source('flock_importance_consolidation.R')
# Override REPS, BUDGET, OUTDIR. TARGET can be 2d, 5d, or 10d.
#
# The true Gaussian mixture is NEVER accessed by exploratory algorithms;
# the simulated-data truth is only used for retrospective evaluation.

fc_mode <- Sys.getenv('RUN_MODE', 'smoke')
fc_defaults <- switch(fc_mode,
                      smoke = c(2L, 3000L),
                      pilot = c(10L, 12000L),
                      full  = c(50L, 30000L),
                      stop('RUN_MODE must be smoke, pilot, or full'))
FC_REPS <- as.integer(Sys.getenv('REPS', fc_defaults[1]))
FC_BUDGET <- as.integer(Sys.getenv('BUDGET', fc_defaults[2]))
FC_CASE <- Sys.getenv('TARGET', '5d')
stopifnot(FC_CASE %in% c('2d', '5d', '10d'),
          is.finite(FC_REPS), FC_REPS >= 1L,
          is.finite(FC_BUDGET), FC_BUDGET >= 1000L)
FC_OUTDIR <- Sys.getenv('OUTDIR',
                        sprintf('flock_importance_%s_%s_B%d',
                                FC_CASE, fc_mode, FC_BUDGET))
dir.create(FC_OUTDIR, recursive=TRUE, showWarnings=FALSE)
dir.create(file.path(FC_OUTDIR, 'checkpoints'), recursive=TRUE,
           showWarnings=FALSE)
FC_METHODS <- c('FLOCK_Roles', 'Undirected_Matched', 'Undirected_RWMH')

fc_logsumexp <- function(z) {
  a <- max(z)
  if (!is.finite(a)) return(a)
  a + log(sum(exp(z - a)))
}
fc_spd <- function(S, floor=0.09) {
  S <- (S + t(S)) / 2
  ee <- eigen(S, symmetric=TRUE)
  ee$vectors %*% diag(pmax(ee$values, floor), length(ee$values)) %*% t(ee$vectors)
}
fc_log_gaussian <- function(x, center, R) {
  # R = chol(S), upper triangular, with crossprod(R)=S.
  # x is generated as center + t(R) %*% z, z~N(0,I).
  z <- forwardsolve(t(R), x - center)
  -length(x)*log(2*pi)/2 - sum(log(diag(R))) - sum(z*z)/2
}
fc_draw_gaussian <- function(center, R) {
  as.numeric(center + t(R) %*% rnorm(length(center)))
}

fc_target <- function(which_case) {
  d <- switch(which_case, '2d'=2L, '5d'=5L, '10d'=10L)
  weights <- c(.44, .30, .21, .05)
  centers <- matrix(0, nrow=4, ncol=d)
  centers[,1:2] <- rbind(c(-5,-5), c(5,5), c(-5,5), c(5,-5))
  if (d >= 3) centers[,3] <- c(0, 0, 1.7, -1.7)
  scales <- c(.95, 1.5, 1.0, .55)
  covs <- lapply(seq_len(4L), function(j) {
    S <- diag(rep(scales[j]^2, d))
    S[1,2] <- S[2,1] <- c(.35, -.40, .25, 0)[j]*scales[j]^2
    if (d >= 4L) S[3,4] <- S[4,3] <- .20 * scales[j]^2
    S
  })
  Rs <- lapply(covs, chol)
  logp <- function(x) {
    terms <- vapply(seq_len(4L), function(j) {
      log(weights[j]) + fc_log_gaussian(x, centers[j,], Rs[[j]])
    }, numeric(1))
    fc_logsumexp(terms)
  }
  responsibilities <- function(x) {
    terms <- vapply(seq_len(4L), function(j)
      log(weights[j]) + fc_log_gaussian(x,centers[j,],Rs[[j]]),numeric(1))
    exp(terms-fc_logsumexp(terms))
  }
  true_mean <- colSums(sweep(centers,1,weights,'*'))
  true_cov <- Reduce('+',lapply(seq_len(4L),function(j) {
    delta <- centers[j,]-true_mean
    weights[j]*(covs[[j]]+tcrossprod(delta))
  }))
  true_q1 <- vapply(c(.025,.975),function(p) {
    marginal_cdf <- function(x) sum(weights*vapply(seq_len(4L),function(j)
      pnorm(x,centers[j,1],sqrt(covs[[j]][1,1])),numeric(1)))
    uniroot(function(x) marginal_cdf(x)-p,
            interval=c(-25,25),tol=1e-10)$root
  },numeric(1))
  start_sd <- c(8,8,rep(2.5,d-2L))
  q0R <- diag(start_sd, d)
  list(d=d, logp=logp, resp=responsibilities, weights=weights,
       means=centers, covs=covs, start_sd=start_sd, q0R=q0R,
       true_mean=true_mean,true_cov=true_cov,true_q1=true_q1)
}

# Target evaluation counter: ALL starts, candidate screening, exploration
# moves, and fresh final importance proposals are charged to the budget.
fc_counter <- function(logp) {
  e <- new.env(parent=emptyenv())
  e$n <- 0L
  e$eval <- function(x) {
    e$n <- e$n + 1L
    logp(x)
  }
  e$count <- function() e$n
  e
}

fc_explorer <- function(f, x0, n, kind, d, q0R, local_sd=1.0) {
  # One counted initial evaluation, n-1 counted proposal evaluations.
  # Scout portfolio (expected): 40% multi-scale ordinary MH,
  # 40% tempered proposals, 20% unconditional broad restarts.
  # Tempered/nonstationary trajectories are DISCOVERY ONLY.
  stopifnot(n>=2L, kind %in% c('scout','follower','pioneer','rw'))
  states <- matrix(NA_real_, n, d)
  lps <- numeric(n)
  x <- as.numeric(x0)
  lp <- f$eval(x)
  states[1L,] <- x; lps[1L] <- lp
  accepted <- 0L; restarts <- 0L
  for (i in 2:n) {
    if (kind=='scout') {
      u <- runif(1)
      if (u<.20) {
        x <- fc_draw_gaussian(rep(0,d), q0R)
        lp <- f$eval(x)
        restarts <- restarts + 1L
        states[i,] <- x; lps[i] <- lp
        next
      } else if (u<.60) {
        beta <- 0.35
        sd_step <- if (runif(1)<.5) 1.7 else 3.8
      } else {
        beta <- 1.0
        sd_step <- if (runif(1)<.5) .8 else 2.2
      }
    } else if (kind=='pioneer') {
      if (runif(1)<.07) {
        x <- fc_draw_gaussian(rep(0,d), q0R)
        lp <- f$eval(x)
        restarts <- restarts + 1L
        states[i,] <- x; lps[i] <- lp
        next
      }
      beta <- if (runif(1)<.6) .6 else 1.0
      sd_step <- if (runif(1)<.5) 1.2 else 2.8
    } else if (kind=='follower') {
      beta <- 1
      sd_step <- .85
    } else {
      beta <- 1
      sd_step <- local_sd
    }
    y <- x + rnorm(d,sd=sd_step)
    ly <- f$eval(y)
    if (log(runif(1)) < beta*(ly-lp)) {
      x <- y; lp <- ly; accepted <- accepted + 1L
    }
    states[i,] <- x; lps[i] <- lp
  }
  list(states=states, lps=lps, accept=accepted/(n-1),
       restarts=restarts, kind=kind)
}

fc_subsample <- function(states, max_keep=140L) {
  idx <- unique(as.integer(round(seq(1,nrow(states),
                                     length.out=min(max_keep,nrow(states))))))
  states[idx,,drop=FALSE]
}

fc_follower_starts <- function(scout,k,init_sd) {
  # Highest-density DISTINCT points, rejecting nearly duplicate candidates.
  o <- order(scout$lps,decreasing=TRUE)
  picked <- integer()
  for (i in o) {
    z <- scout$states[i,]
    if (!length(picked)) {
      picked <- c(picked,i)
    } else {
      ds <- vapply(picked, function(j)
        sqrt(sum(((z-scout$states[j,])/2)^2)),numeric(1))
      if (min(ds)>1.2) picked <- c(picked,i)
    }
    if (length(picked)>=k) break
  }
  ans <- scout$states[picked,,drop=FALSE]
  while (nrow(ans)<k) {
    ans <- rbind(ans,scout$states[o[1],]+rnorm(ncol(ans),sd=.6))
  }
  ans
}

fc_novelty <- function(x, ref, local_scale) {
  z <- sweep(ref,2,x,'-')
  min(sqrt(rowSums(sweep(z,2,local_scale,'/')^2)))
}

fc_pioneer_start <- function(f, history, d, q0R, n_candidates=12L,
                             selection=c('novelty','random')) {
  selection <- match.arg(selection)
  # Screening is charged: each candidate needs one posterior evaluation.
  cand <- t(replicate(n_candidates,fc_draw_gaussian(rep(0,d),q0R)))
  cand <- matrix(cand,n_candidates,d)
  lpc <- apply(cand,1L,f$eval)
  ref <- fc_subsample(history,180L)
  scale <- pmax(apply(ref,2L,sd),1.0)
  novelty <- vapply(seq_len(n_candidates),function(i)
    fc_novelty(cand[i,],ref,scale),numeric(1))
  utility <- novelty-.12*pmax(0,max(lpc)-lpc)
  # If no candidate has plausible density relative to the visited map,
  # prefer the highest-density candidate rather than a useless extreme.
  best <- if(selection=='novelty') which.max(utility) else
    sample.int(n_candidates,1L)
  list(x=cand[best,],spent=n_candidates,candidate_lp=lpc[best],
       novelty=novelty[best])
}

fc_exploration <- function(which_method, f, target, B) {
  d <- target$d
  exploration_B <- floor(.65*B)
  histories <- vector('list',7L)
  regions_by_stage <- list()
  if (which_method=='Undirected_RWMH') {
    alloc <- rep(exploration_B%/%7L,7L)
    alloc[1] <- alloc[1]+exploration_B-sum(alloc)
    for(j in seq_len(7L)) {
      start <- fc_draw_gaussian(rep(0,d),target$q0R)
      histories[[j]] <- fc_explorer(f,start,alloc[j],'rw',d,target$q0R,
                                    local_sd=1.2)
    }
    stage <- c('ordinary',rep('ordinary',6L))
  } else {
    scout_B <- floor(.20*B)
    mid_B <- exploration_B-scout_B
    alloc <- rep(mid_B%/%6L,6L)
    alloc[1] <- alloc[1]+mid_B-sum(alloc)
    scout_start <- fc_draw_gaussian(rep(0,d),target$q0R)
    histories[[1]] <- fc_explorer(f,scout_start,scout_B,'scout',d,target$q0R)
    stage <- c('scout',rep('follower',3L),rep('pioneer',3L))
    follower_starts <- fc_follower_starts(histories[[1]],3L,target$start_sd)
    for(j in 1:3) {
      if (which_method=='FLOCK_Roles') {
        start <- follower_starts[j,]
      } else {
        start <- fc_draw_gaussian(rep(0,d),target$q0R)
      }
      histories[[j+1L]] <- fc_explorer(f,start,alloc[j],'follower',d,target$q0R)
    }
    ref <- do.call(rbind,lapply(histories[1:4],function(z)
      fc_subsample(z$states,120L)))
    for(j in 1:3) {
      slot <- 3L+j
      n <- alloc[slot]
      nc <- min(12L,max(2L,n%/%25L))
      start <- fc_pioneer_start(
        f,ref,d,target$q0R,nc,
        selection=if(which_method=='FLOCK_Roles') 'novelty' else 'random')
      # Both variants spend EXACTLY the same screening budget and run
      # exactly the same pioneer kernel after selecting a starting point.
      n_steps <- n-start$spent
      x <- start$x
      histories[[slot+1L]] <- fc_explorer(f,x,n_steps,'pioneer',d,target$q0R)
      ref <- rbind(ref,fc_subsample(histories[[slot+1L]]$states,120L))
    }
  }
  stopifnot(f$count()==exploration_B)
  list(histories=histories,stage=stage,evals=exploration_B)
}

fc_build_q <- function(histories, d, q0R, region_max=7L, defensive=.15) {
  # All methods use IDENTICAL regional proposal learner.
  pool <- do.call(rbind,lapply(histories,function(z)
    fc_subsample(z$states,150L)))
  pool <- unique(pool)
  if(nrow(pool)<2L) stop('Insufficient distinct exploration points')
  k <- min(region_max,max(1L, floor(nrow(pool)/75L)))
  k <- min(k,nrow(pool))
  if(k>=2L) {
    fit <- tryCatch(stats::kmeans(pool,centers=k,iter.max=60,nstart=3),
                    error=function(e) NULL)
  } else fit <- NULL
  if(is.null(fit)) {
    k <- 1L
    cluster <- rep(1L,nrow(pool))
    centers <- matrix(colMeans(pool),nrow=1L)
  } else {
    cluster <- fit$cluster
    centers <- matrix(fit$centers,ncol=d)
  }
  Rs <- lapply(seq_len(k),function(j) {
    z <- pool[cluster==j,,drop=FALSE]
    S <- if(nrow(z)>d+3) cov(z) else diag(1,d)
    S <- .85 * matrix(S,d,d) + .15 * diag(diag(matrix(S,d,d)),d)
    chol(fc_spd(S + .18*diag(d),floor=.15))
  })
  centers <- rbind(centers,rep(0,d))
  Rs[[k+1L]] <- q0R
  mass <- c(rep((1-defensive)/k,k),defensive)
  logq <- function(x) {
    fc_logsumexp(vapply(seq_along(mass), function(j)
      log(mass[j])+fc_log_gaussian(x,centers[j,],Rs[[j]]),numeric(1)))
  }
  draw <- function() {
    j <- sample.int(length(mass),1L,prob=mass)
    fc_draw_gaussian(centers[j,],Rs[[j]])
  }
  list(logq=logq,draw=draw,regions=k,centers=centers,pool=pool)
}

fc_final_importance <- function(f,q,n,target,store_samples=FALSE) {
  d <- target$d
  x <- matrix(NA_real_,n,d)
  logweights <- numeric(n)
  for(i in seq_len(n)) {
    x[i,] <- q$draw()
    logweights[i] <- f$eval(x[i,]) - q$logq(x[i,])
  }
  log_norm <- fc_logsumexp(logweights)
  if(!is.finite(log_norm)) stop('Nonfinite importance weights')
  w <- exp(logweights-log_norm)
  ess <- 1/sum(w*w)
  mu <- colSums(sweep(x,1,w,'*'))
  z <- sweep(x,2,mu,'-')
  S <- crossprod(z,sweep(z,1,w,'*'))
  responsibilities <- t(vapply(seq_len(n),function(i)
    target$resp(x[i,]),numeric(4)))
  estimated_component_mass <- colSums(sweep(responsibilities,1,w,'*'))
  # Detection is retrospective only. The sampler never sees true components.
  component_cores <- colSums(responsibilities>.8)
  result <- list(mu=mu,cov=S,component_mass=estimated_component_mass,
                 marginal_ci_1=fc_wquantile(x[,1],w,c(.025,.975)),
                 final_core_counts=component_cores,ess=ess,ess_fraction=ess/n,
                 max_weight=max(w),draws=n)
  if(store_samples) {
    result$samples <- x
    result$weights <- w
  }
  result
}

fc_wquantile <- function(x,w,p) {
  o <- order(x)
  cs <- cumsum(w[o])
  x[o[pmax(1L, vapply(p,function(a)
    which(cs>=a)[1L],integer(1)))]]
}

fc_method <- function(method,target, B, store_samples=FALSE) {
  f <- fc_counter(target$logp)
  elapsed <- system.time({
    exploration <- fc_exploration(method,f,target,B)
    q <- fc_build_q(exploration$histories,target$d,target$q0R)
    nfinal <- B-f$count()
    stopifnot(nfinal>=10L)
    ans <- fc_final_importance(f,q,nfinal,target,store_samples=store_samples)
  })['elapsed']
  stopifnot(f$count()==B)
  # Retrospective component discoveries: require >=3 accepted/history states
  # with responsibility >0.8. Counts can include repeated states, so report
  # as coarse evidence of visits, not independent discoveries.
  stage_discovery <- lapply(seq_along(exploration$histories),function(j) {
    z <- fc_subsample(exploration$histories[[j]]$states,200L)
    resp <- t(vapply(seq_len(nrow(z)),function(i)
      target$resp(z[i,]),numeric(4)))
    data.frame(chain=j, stage=exploration$stage[j],
               component=1:4, core_visits=colSums(resp>.8),
               stringsAsFactors=FALSE)
  })
  stage_discovery <- do.call(rbind,stage_discovery)
  m_err <- sqrt(mean((ans$mu-target$true_mean)^2))
  w_err <- sum(abs(ans$component_mass-target$weights))
  cov_error <- sqrt(sum((ans$cov-target$true_cov)^2))/
    sqrt(sum(target$true_cov^2))
  ci_error <- mean(abs(ans$marginal_ci_1-target$true_q1))
  # Consequential discovery: all components with true mass >= 0.10.
  important <- which(target$weights>=.10)
  missed_important <- sum(ans$component_mass[important]<.025)
  result <- list(
    metrics=data.frame(method=method,evals=B,explore_evals=exploration$evals,
      final_evals=nfinal,elapsed_seconds=as.numeric(elapsed),
      proposal_regions=q$regions,rmse=m_err,component_L1=w_err,
      covariance_relative_error=cov_error,ci_endpoint_error=ci_error,
      rare_mass=ans$component_mass[4],rare_abs_error=abs(ans$component_mass[4]-target$weights[4]),
      is_ess=ans$ess,ess_fraction=ans$ess_fraction,
      max_normalized_weight=ans$max_weight,
      missed_important=missed_important),
    discovery=stage_discovery,
    estimate=ans$mu,
    ci_theta1=ans$marginal_ci_1,
    component_mass=ans$component_mass,
    final_core_counts=ans$final_core_counts,
    q_centers=q$centers
  )
  if(store_samples) {
    result$final_samples <- ans$samples
    result$final_weights <- ans$weights
    result$exploration_samples <- lapply(exploration$histories,function(z)
      fc_subsample(z$states,400L))
    result$exploration_stages <- exploration$stage
  }
  result
}

fc_checkpoint <- function(path,obj) {
  temp <- paste0(path,'.tmp')
  saveRDS(obj,temp)
  if(!file.rename(temp,path)) stop('Could not rename checkpoint ',path)
}
fc_summary <- function(results) {
  res <- do.call(rbind,lapply(results,function(z) z$metrics))
  do.call(rbind,lapply(split(res,res$method),function(z) {
    data.frame(method=z$method[1],n=nrow(z),
       mean_RMSE=mean(z$rmse),median_RMSE=median(z$rmse),
       mean_L1=mean(z$component_L1),
       mean_covariance_error=mean(z$covariance_relative_error),
       mean_CI_endpoint_error=mean(z$ci_endpoint_error),
       mean_rare_error=mean(z$rare_abs_error),
       mean_IS_ESS=mean(z$is_ess),
       mean_ESS_fraction=mean(z$ess_fraction),
       missed_important_rate=mean(z$missed_important>0),
       mean_seconds=mean(z$elapsed_seconds))
  }))
}
fc_pairs <- function(results) {
  res <- do.call(rbind,lapply(results,function(z) z$metrics))
  out <- list()
  for(comp in c('Undirected_Matched','Undirected_RWMH')) {
    a <- res[res$method=='FLOCK_Roles',]
    b <- res[res$method==comp,]
    pair <- merge(a,b,by='rep',suffixes=c('_flock','_control'))
    for(metric in c('rmse','component_L1','covariance_relative_error',
                    'ci_endpoint_error','rare_abs_error','ess_fraction')) {
      d <- pair[[paste0(metric,'_flock')]]-pair[[paste0(metric,'_control')]]
      # For ESS fraction, positive favors FLOCK. For errors, negative favors FLOCK.
      if(length(d)>=2L) {
        bs <- replicate(500L,mean(sample(d,length(d),replace=TRUE)))
        ci <- unname(quantile(bs,c(.025,.975)))
      } else ci <- c(NA_real_,NA_real_)
      out[[length(out)+1L]] <- data.frame(
        competitor=comp,metric=metric,n=length(d),
        mean_difference=mean(d),ci_low=ci[1],ci_high=ci[2],
        fraction_flock_better=mean(if(metric=='ess_fraction') d>0 else d<0))
    }
  }
  do.call(rbind,out)
}

fc_make_plots <- function(results,target,outdir) {
  res <- do.call(rbind,lapply(results,function(z) z$metrics))
  pdf(file.path(outdir,'flock_importance_diagnostics.pdf'),width=11,height=8)
  par(mfrow=c(2,2),mar=c(5,4,3,1))
  boxplot(rmse~method,res,las=2,main='Posterior mean RMSE',ylab='RMSE')
  boxplot(component_L1~method,res,las=2,main='Component mass L1 error',ylab='L1 error')
  boxplot(ess_fraction~method,res,las=2,main='Importance ESS / final N',ylab='ESS fraction')
  boxplot(elapsed_seconds~method,res,las=2,main='Elapsed seconds',ylab='Seconds')
  methods <- FC_METHODS
  par(mfrow=c(length(methods),2),mar=c(4,4,3,1))
  for(m in methods) {
    key <- which(vapply(results,function(z)
      z$metrics$method==m && z$metrics$rep==1L,logical(1)))
    if(!length(key) || is.null(results[[key]]$final_samples)) next
    z <- results[[key]]
    xx <- z$final_samples
    ww <- z$final_weights
    # Weighted scatter: intentionally stratify a small subsample to draw
    # points approximately according to the weighted posterior.
    sel <- sample.int(nrow(xx),min(1000L,nrow(xx)),replace=TRUE,prob=ww)
    plot(xx[sel,1],xx[sel,2],pch=16,cex=.4,
         xlab=expression(theta[1]),ylab=expression(theta[2]),
         main=paste(m,'weighted posterior'))
    points(target$means[,1],target$means[,2],pch=4,col='red',cex=1.2,lwd=2)
    barplot(rbind(target$weights,z$component_mass),beside=TRUE,
            col=c('grey60','steelblue'),names.arg=1:4,
            main=paste(m,'component weights'),ylab='Probability',
            ylim=c(0,.6))
    legend('topright',c('True','Estimated'),fill=c('grey60','steelblue'),
           bty='n',cex=.7)
  }
  # Show how the exploratory computational budget was actually deployed.
  par(mfrow=c(1,3),mar=c(4,4,3,1))
  stage_colors <- c(scout='firebrick', follower='darkgreen',
                    pioneer='steelblue',ordinary='grey45')
  for(m in methods) {
    key <- which(vapply(results,function(z)
      z$metrics$method==m && z$metrics$rep==1L,logical(1)))
    if(!length(key) || is.null(results[[key]]$exploration_samples)) {
      plot.new();next
    }
    z <- results[[key]]
    xs <- do.call(rbind,z$exploration_samples)
    labs <- rep(z$exploration_stages,
                times=vapply(z$exploration_samples,nrow,integer(1)))
    plot(xs[,1],xs[,2],pch=16,cex=.35,
         col=unname(stage_colors[labs]),
         xlab=expression(theta[1]),ylab=expression(theta[2]),
         main=paste(m,'exploration'))
    points(target$means[,1],target$means[,2],pch=4,col='black',cex=1.3,lwd=2)
  }
  dev.off()
}

# Algebraic self-tests: Cholesky covariance convention and normalization.
fc_test_cov <- matrix(c(2,.5,.5,1),2,2)
stopifnot(max(abs(crossprod(chol(fc_test_cov))-fc_test_cov))<1e-12)
fc_test_x <- c(.4,-.9)
fc_test_log_ref <- -log(2*pi)-sum(log(diag(chol(fc_test_cov))))-
  as.numeric(crossprod(fc_test_x,solve(fc_test_cov,fc_test_x)))/2
stopifnot(abs(fc_log_gaussian(fc_test_x,c(0,0),chol(fc_test_cov))-
              fc_test_log_ref)<1e-10)
# Draw and density conventions both use t(chol(S)), avoiding the earlier
# factor-transposition error in legacy FLOCK scripts.
set.seed(123L)
fc_test_draws <- t(replicate(3000L,
  fc_draw_gaussian(c(0,0),chol(fc_test_cov))))
stopifnot(max(abs(cov(fc_test_draws)-fc_test_cov))<.20)
fc_test_target <- fc_target(FC_CASE)
stopifnot(abs(sum(fc_test_target$weights)-1)<1e-12,
          all(vapply(seq_len(4L),function(j)
            abs(sum(fc_test_target$resp(fc_test_target$means[j,]))-1)<1e-12,
            logical(1))))
cat('Numerical checks passed. Target:',FC_CASE,
    '| replicates:',FC_REPS,'| budget/method:',FC_BUDGET,'\n')

fc_results <- list()
for(r in seq_len(FC_REPS)) {
  for(method in FC_METHODS) {
    file <- file.path(FC_OUTDIR,'checkpoints',
                      sprintf('%s_rep%03d.rds',method,r))
    if(file.exists(file)) {
      z <- readRDS(file)
    } else {
      # Per-method common initial seed gives a paired experimental design.
      # Different random consumption after the first stage is expected.
      set.seed(910000L + r)
      z <- fc_method(method,fc_test_target,FC_BUDGET,
                     store_samples=(r==1L))
      z$metrics$rep <- r
      fc_checkpoint(file,z)
    }
    fc_results[[length(fc_results)+1L]] <- z
    m <- z$metrics
    cat(sprintf('[%d/%d] %-20s RMSE=%.4f L1=%.4f ESS=%.1f/%d (%.1f%%) time=%.1fs\n',
        r,FC_REPS,method,m$rmse,m$component_L1,m$is_ess,m$final_evals,
        100*m$ess_fraction,m$elapsed_seconds))
  }
  all_res <- do.call(rbind,lapply(fc_results,function(z) z$metrics))
  write.csv(all_res,file.path(FC_OUTDIR,'replicate_results.csv'),row.names=FALSE)
  write.csv(fc_summary(fc_results),file.path(FC_OUTDIR,'method_summary.csv'),row.names=FALSE)
}

set.seed(810000L)
pair_table <- fc_pairs(fc_results)
write.csv(pair_table,file.path(FC_OUTDIR,'paired_comparisons.csv'),row.names=FALSE)
discovery <- do.call(rbind,lapply(fc_results,function(z) {
  cbind(method=z$metrics$method,rep=z$metrics$rep,z$discovery)
}))
write.csv(discovery,file.path(FC_OUTDIR,'exploratory_component_coverage.csv'),row.names=FALSE)

for(z in fc_results) if(z$metrics$rep==1L && !is.null(z$final_samples)) {
  frame <- as.data.frame(z$final_samples)
  names(frame) <- paste0('theta',seq_len(ncol(frame)))
  frame$importance_weight <- z$final_weights
  write.csv(frame,file.path(FC_OUTDIR,paste0(z$metrics$method,
                                             '_weighted_posterior.csv')),
            row.names=FALSE)
}
fc_make_plots(fc_results,fc_test_target,FC_OUTDIR)
cat('\nSummary:\n')
print(fc_summary(fc_results),row.names=FALSE,digits=4)
cat('\nPaired differences: errors < 0 favor FLOCK; ESS differences > 0 favor FLOCK.\n')
print(pair_table,row.names=FALSE,digits=4)
cat('\nFiles saved in:',normalizePath(FC_OUTDIR),'\n')