# FLOCK: Minimal staged sampler versus Adaptive Metropolis
# Six Bayesian scenarios: (1) product ridge, (2) weak-signal funnel,
# (3) three-component mixture, (4) correlated predictors with spike/slab,
# (5) diagnostic latent class analysis with two imperfect tests,
# (6) high-dimensional (40 covariates) Bayesian regression with sparse
#     ambiguous coefficients and analytic mixture posterior reference.
# Only FLOCK final-chain post-warmup draws are posterior draws.
# All target evaluations including candidate selection are counted.
# No external packages required.
#
# Run example:
#   Sys.setenv(RUN_MODE="smoke"); source("flock_simple_five_scenarios.R")
#   Sys.setenv(RUN_MODE="pilot"); source("flock_simple_five_scenarios.R")
#   Sys.setenv(RUN_MODE="full"); source("flock_simple_five_scenarios.R")
# Optional overrides: REPS, BUDGET, OUTDIR, SEED, SCENARIOS
# References are analytically/numerically integrated from specified models;
# algorithms are NEVER given the reference values.
# This is a research prototype, not a validated replacement for convergence
# diagnostics and not an implementation of an existing published FLOCK method.

RUN_MODE <- Sys.getenv("RUN_MODE", "smoke")
settings <- switch(RUN_MODE,
  smoke=c(2L,4000L), pilot=c(12L,30000L), full=c(60L,60000L),
  stop("RUN_MODE must be smoke, pilot, or full"))
REPS <- as.integer(Sys.getenv("REPS",settings[1]))
BUDGET <- as.integer(Sys.getenv("BUDGET",settings[2]))
SEED <- as.integer(Sys.getenv("SEED","20261025"))
OUTDIR <- Sys.getenv("OUTDIR",sprintf("flock_simple_%s_B%d",RUN_MODE,BUDGET))
stopifnot(REPS>=1L,BUDGET>=3000L,SEED>=0L,nzchar(OUTDIR))
dir.create(OUTDIR,showWarnings=FALSE,recursive=TRUE)
dir.create(file.path(OUTDIR,"checkpoints"),showWarnings=FALSE,recursive=TRUE)
methods <- c("AdaptiveMetropolis","FLOCK_Novelty","FLOCK_OrthogonalSoft")
all_scenario_names <- c("ProductRidge","WeakSignalFunnel","ThreeComponentMixture",
                        "CorrelatedPredictors","DiagnosticLCA","HighDimCorrelated")
selection <- Sys.getenv("SCENARIOS","all")
scenario_names <- if(tolower(selection)=="all") all_scenario_names else
  trimws(strsplit(selection,",",fixed=TRUE)[[1L]])
if(!length(scenario_names) || any(!scenario_names %in% all_scenario_names) ||
   anyDuplicated(scenario_names)) stop("SCENARIOS must be 'all' or comma-separated scenario names")
version <- "v2_flock_softorth_six_scenarios_20261025"
manifest <- list(version=version,methods=methods,scenarios=scenario_names,
                 budget=BUDGET,reps=REPS,seed=SEED)
manifest_file <- file.path(OUTDIR,"manifest.rds")
if (file.exists(manifest_file)) {
 if(!identical(readRDS(manifest_file),manifest)) stop("Output manifest mismatch; set a fresh OUTDIR")
} else saveRDS(manifest,manifest_file)

logsum <- function(x) {
  a<-max(x)
  if(!is.finite(a)) return(-Inf)
  a+log(sum(exp(x-a)))
}
logmean <- function(x) {
  a<-max(x)
  if(!is.finite(a)) return(-Inf)
  a+log(mean(exp(x-a)))
}
spd <- function(S,min_eigen=.035) {
  S<-(S+t(S))/2
  ee<-eigen(S,symmetric=TRUE)
  ee$vectors %*% diag(pmax(ee$values,min_eigen),nrow(S)) %*% t(ee$vectors)
}
# R's chol(S) is upper triangular U with crossprod(U)=S.
# Sampling uses t(U) z; Gaussian log density solves t(U) z = x - mu.
draw_norm <- function(mu,U) as.numeric(mu+crossprod(U,rnorm(length(mu))))
log_norm <- function(x,mu,U) {
  d<-length(mu)
  z<-forwardsolve(t(U),x-mu)
  -d*.5*log(2*pi)-sum(log(diag(U)))-sum(z*z)/2
}
make_counter <- function(target) {
  e<-new.env(parent=emptyenv())
  e$n<-0L
  e$eval<-function(x) {
    e$n<-e$n+1L
    as.numeric(target$logpost(as.numeric(x)))
  }
  e$get<-function() e$n
  e
}
# Common mixed-move exploratory MH kernel. All proposals count as target calls.
# The independence-from-prior move has its exact q(x)/q(y) correction.
walk <- function(f,target,x,lp,n,U,prior_jump=.08,adapt=FALSE) {
  d<-target$d
  H<-matrix(NA_real_,n+1L,d); L<-numeric(n+1L)
  H[1L,]<-x;L[1L]<-lp;acc<-0L
  for(t in seq_len(n)) {
    if(adapt && t>max(70L,15L*d) && t%%75L==0L) {
      S<-cov(H[seq_len(t),,drop=FALSE])
      U<-chol(spd((2.38^2/d)*S+.035*diag(d)))
    }
    do_prior<-runif(1)<prior_jump
    if(do_prior) {
      y<-target$prior_draw()
      ly<-f$eval(y)
      loga<-ly-lp+target$logprior(x)-target$logprior(y)
    } else {
      y<-draw_norm(x,U)
      ly<-f$eval(y)
      loga<-ly-lp
    }
    if(is.finite(loga) && log(runif(1))<loga) {
      x<-y;lp<-ly;acc<-acc+1L
    }
    H[t+1L,]<-x;L[t+1L]<-lp
  }
  list(x=x,lp=lp,h=H,lps=L,accept=acc/max(1L,n),last_U=U)
}
select_pool <- function(histories,limit=140L) {
  do.call(rbind,lapply(histories,function(z) {
    i<-unique(as.integer(round(seq(1,nrow(z),length.out=min(nrow(z),limit)))))
    z[i,,drop=FALSE]
  }))
}
# Approximate posterior geometry; the learned mixture is a PROPOSAL ONLY.
# Prior fallback ensures broad support as long as the prior has full support.
build_proposal <- function(target,histories) {
  d<-target$d
  pool<-unique(select_pool(histories,limit=150L))
  if(nrow(pool)<5L) stop('Insufficient distinct exploratory points')
  K<-min(8L,max(1L,floor(nrow(pool)/45L)))
  km<-tryCatch(kmeans(pool,centers=K,nstart=2,iter.max=45),
               error=function(e) NULL)
  if(is.null(km)) {
    K<-1L
    groups<-rep(1L,nrow(pool))
    centers<-matrix(colMeans(pool),1L,d)
  } else {
    groups<-km$cluster
    centers<-matrix(km$centers,ncol=d)
  }
  UU<-lapply(seq_len(K),function(k) {
    z<-pool[groups==k,,drop=FALSE]
    S<-if(nrow(z)>=d+3L) cov(z) else diag(2,d)
    chol(spd(1.3*S+.10*diag(d)))
  })
  masses<-rep(.82/K,K)
  qdraw<-function() {
    if(runif(1)<.18) return(target$prior_draw())
    k<-sample.int(K,1L,prob=masses)
    draw_norm(centers[k,],UU[[k]])
  }
  qlog<-function(x) {
    parts<-vapply(seq_len(K),function(k)
      log(masses[k])+log_norm(x,centers[k,],UU[[k]]),numeric(1))
    logsum(c(parts,log(.18)+target$logprior(x)))
  }
  list(draw=qdraw,logq=qlog,K=K)
}
# Mixture of two separately corrected MH kernels: local RWMH and independence MH.
final_sampler <- function(f,target,x,lp,n,q) {
  d<-target$d
  H<-matrix(NA_real_,n+1L,d);H[1L,]<-x
  lqx<-q$logq(x);acc<-0L
  local_sd <- if(d>10L) 2.0/sqrt(d) else .65
  Lloc<-chol(diag(local_sd^2,d))
  for(i in seq_len(n)) {
    local<-runif(1)<.20
    if(local) {
      y<-draw_norm(x,Lloc);ly<-f$eval(y)
      loga<-ly-lp
      lqy<-q$logq(y)
    } else {
      y<-q$draw();ly<-f$eval(y);lqy<-q$logq(y)
      loga<-ly-lp+lqx-lqy
    }
    if(is.finite(loga) && log(runif(1))<loga) {
      x<-y;lp<-ly;lqx<-lqy;acc<-acc+1L
    }
    H[i+1L,]<-x
  }
  burn<-floor(.10*n)
  list(draws=H[(burn+2L):(n+1L),,drop=FALSE],accept=acc/n)
}

# ---- Scenario 1: product ridge (nonidentifiability and sign symmetry) ------
make_ridge <- function() {
  set.seed(101)
  sd_y<-.7
  y<-rnorm(6L,mean=3.7,sd=sd_y)
  n<-length(y);yb<-mean(y);sp<-3
  lp<-function(x) {
    if(length(x)!=2L) stop('Incorrect dimension')
    sum(dnorm(y,mean=x[1]*x[2],sd=sd_y,log=TRUE))+
      sum(dnorm(x,0,sp,log=TRUE))
  }
  prior<-function(x) sum(dnorm(x,0,sp,log=TRUE))
  # Exact marginal theta1 is Gaussian prior times integrated normal likelihood.
  grid<-seq(-13,13,length.out=2601L)
  logw<-dnorm(grid,0,sp,log=TRUE)+
    dnorm(yb,0,sqrt(sd_y^2/n+grid^2*sp^2),log=TRUE)
  wg<-exp(logw-max(logw));wg<-wg/sum(wg)
  v2<-1/(1/sp^2+n*grid^2/sd_y^2)
  m2<-v2*n*grid*yb/sd_y^2
  ref<-c(theta1_sq=sum(wg*grid^2),
         theta2_sq=sum(wg*(v2+m2^2)),
         product=sum(wg*grid*m2),
         p_theta1_pos=sum(wg[grid>0]))
  metrics<-function(X) c(theta1_sq=mean(X[,1]^2),
    theta2_sq=mean(X[,2]^2),product=mean(X[,1]*X[,2]),
    p_theta1_pos=mean(X[,1]>0))
  list(name='ProductRidge',d=2L,logpost=lp,logprior=prior,
       prior_draw=function() rnorm(2,0,sp),reference=ref,
       metric_scales=c(4,4,4,.5),metrics=metrics,
       modes=function(X) ifelse(X[,1]>0,2L,1L),
       mode_probs=c(.5,.5),mode_type='branch',
       show_coord=1L,ref_grid=grid,ref_pdf=wg/(grid[2]-grid[1]),
       description='Y|theta ~ N(theta1*theta2,0.7^2): curved sign-symmetric ridge')
}

# ---- Scenario 2: weak-signal hierarchical funnel ---------------------------
make_funnel <- function() {
  y<-c(-.8,.25,1.1,-.55,.65,-1.05)
  sd_y<-2.3
  n<-length(y);d<-n+1L
  lprior<-function(x) {
    t<-x[1L];tau<-exp(t)
    if(!is.finite(tau) || tau<=0) return(-Inf)
    dnorm(t,-.6,1,log=TRUE)+sum(dnorm(x[-1L],0,tau,log=TRUE))
  }
  lp<-function(x) lprior(x)+sum(dnorm(y,x[-1L],sd_y,log=TRUE))
  prior_draw<-function() {
    t<-rnorm(1L,-.6,1)
    c(t,rnorm(n,0,exp(t)))
  }
  grid<-seq(-6.5,4.5,length.out=3301L)
  logw<-dnorm(grid,-.6,1,log=TRUE)+
    vapply(grid,function(t)
      sum(dnorm(y,0,sqrt(exp(2*t)+sd_y^2),log=TRUE)),numeric(1))
  wg<-exp(logw-max(logw));wg<-wg/sum(wg)
  tau2<-exp(2*grid)
  v<-tau2*sd_y^2/(tau2+sd_y^2)
  m<-tau2/(tau2+sd_y^2)*y[1]
  ref<-c(mean_logtau=sum(wg*grid),p_tau_gt_1=sum(wg[grid>0]),
         mean_eta1=sum(wg*m),eta1_sq=sum(wg*(m^2+v)))
  metrics<-function(X) c(mean_logtau=mean(X[,1]),
      p_tau_gt_1=mean(X[,1]>0),mean_eta1=mean(X[,2]),
      eta1_sq=mean(X[,2]^2))
  list(name='WeakSignalFunnel',d=d,logpost=lp,logprior=lprior,
       prior_draw=prior_draw,reference=ref,
       metric_scales=c(1,.5,1,1.5),metrics=metrics,
       modes=function(X) ifelse(X[,1]>0,2L,1L),
       mode_probs=c(1-ref['p_tau_gt_1'],ref['p_tau_gt_1']),
       mode_type='tail_threshold_not_mode',
       show_coord=1L,ref_grid=grid,ref_pdf=wg/(grid[2]-grid[1]),
       description='6 weakly observed random effects; log(tau) posterior has funnel geometry')
}

# ---- Scenario 3: genuine Bayesian finite mixture with label switching ------
make_mixture <- function() {
  y<-c(-4.15,-3.75,-4.4,-.30,.15,.38,3.7,4.1,4.35)
  n<-length(y);sd_y<-1;sp<-4.5
  weights_fn<-function(x) {
    z<-c(x[4],x[5],0);zz<-exp(z-max(z));zz/sum(zz)
  }
  lprior<-function(x) {
    if(length(x)!=5L || any(!is.finite(x))) return(-Inf)
    w<-weights_fn(x)
    sum(dnorm(x[1:3],0,sp,log=TRUE))+log(2)+sum(log(w))
  }
  lp<-function(x) {
    lpr<-lprior(x)
    if(!is.finite(lpr)) return(-Inf)
    w<-weights_fn(x)
    denmat<-matrix(dnorm(as.vector(outer(y,x[1:3],'-')),0,sd_y),nrow=n,ncol=3L)
    lpr+sum(log(as.numeric(denmat%*%w)))
  }
  prior_draw<-function() {
    g<-rgamma(3L,shape=1);w<-g/sum(g)
    c(rnorm(3L,0,sp),log(w[1]/w[3]),log(w[2]/w[3]))
  }
  # EXACT posterior reference by summing 3^9 allocations. Conjugate
  # N(0,sp^2) component means and Dirichlet(1,1,1) mixture weights.
  # No reference information is available to the samplers.
  assigns<-as.matrix(do.call(expand.grid,rep(list(1:3),n)))
  N<-nrow(assigns)
  ll<-numeric(N)
  summaries<-matrix(NA_real_,N,3L)
  mu_post<-matrix(NA_real_,N,3L)
  v_post<-matrix(NA_real_,N,3L)
  for(i in seq_len(N)) {
    z<-assigns[i,]
    nk<-tabulate(z,nbins=3L)
    log_ev<-lgamma(3)-lgamma(n+3)+sum(lgamma(nk+1))
    mm<-vv<-numeric(3)
    for(k in 1:3) {
      yy<-y[z==k]
      vv[k]<-1/(1/sp^2+nk[k]/sd_y^2)
      mm[k]<-vv[k]*sum(yy)/sd_y^2
      log_ev<-log_ev-nk[k]/2*log(2*pi*sd_y^2)+
        .5*log(vv[k]/sp^2)-
        .5*(sum(yy^2)/sd_y^2-mm[k]^2/vv[k])
    }
    ak<-(1+nk)/(n+3)
    summaries[i,]<-c(sum(ak*mm),mean(mm^2+vv),sum(ak*(mm^2+vv)))
    mu_post[i,]<-mm;v_post[i,]<-vv
    ll[i]<-log_ev
  }
  weights<-exp(ll-logsum(ll))
  ref<-c(pred_mean=sum(weights*summaries[,1]),
         mean_mu_sq=sum(weights*summaries[,2]),
         weighted_mu_sq=sum(weights*summaries[,3]),
         p_label_123=1/6)
  metrics<-function(X) {
    v1<-X[,4];v2<-X[,5];v3<-rep(0,nrow(X))
    zmax<-pmax(v1,v2,v3)
    a<-exp(v1-zmax);b<-exp(v2-zmax);c<-exp(v3-zmax)
    denom<-a+b+c
    w1<-a/denom;w2<-b/denom;w3<-c/denom
    mu<-X[,1:3,drop=FALSE]
    c(pred_mean=mean(w1*mu[,1]+w2*mu[,2]+w3*mu[,3]),
      mean_mu_sq=mean(rowMeans(mu^2)),
      weighted_mu_sq=mean(w1*mu[,1]^2+w2*mu[,2]^2+w3*mu[,3]^2),
      p_label_123=mean(mu[,1]<mu[,2] & mu[,2]<mu[,3]))
  }
  label_fn<-function(X) {
    mu<-X[,1:3,drop=FALSE]
    key<-apply(mu,1L,function(z) {
      ii<-order(z);100L*ii[1]+10L*ii[2]+ii[3]
    })
    as.integer(match(key,c(123,132,213,231,312,321)))
  }
  gx<-seq(-9,9,length.out=361L)
  pgrid<-vapply(gx,function(v)
    sum(weights*dnorm(v,mu_post[,1],sqrt(v_post[,1]))),numeric(1))
  list(name='ThreeComponentMixture',d=5L,logpost=lp,
       logprior=lprior,prior_draw=prior_draw,
       reference=ref,metric_scales=c(2,8,8,1/6),metrics=metrics,
       modes=label_fn,mode_probs=rep(1/6,6L),mode_type='six_label_permutations',
       show_coord=1L,ref_grid=gx,ref_pdf=pgrid,
       description='3-component Gaussian mixture with exchangeable priors: 6 labelings')
}


# ---- Scenario 4: highly correlated predictors, competing coefficient modes -
# Proper continuous spike/slab Gaussian prior and normal likelihood.
# Four latent prior-component assignments have an EXACT Gaussian-mixture
# posterior. We integrate them analytically for reference diagnostics.
make_correlated <- function() {
  set.seed(2203)
  n <- 30L
  x1 <- rnorm(n)
  x2 <- .9995*x1 + sqrt(1-.9995^2)*rnorm(n)
  X <- cbind(x1,x2)
  sigma <- .65
  y <- as.numeric(2*x1 + rnorm(n,sd=sigma))
  spike <- .09; slab <- 2.5; slab_prob <- .40
  log_mix <- function(x) {
    vapply(x,function(z) logsum(c(log(1-slab_prob)+dnorm(z,0,spike,log=TRUE),
                                     log(slab_prob)+dnorm(z,0,slab,log=TRUE))),numeric(1))
  }
  lprior <- function(b) sum(log_mix(b))
  lp <- function(b) {
    if(length(b)!=2L || any(!is.finite(b))) return(-Inf)
    lprior(b) + sum(dnorm(y,as.numeric(X%*%b),sigma,log=TRUE))
  }
  prior_draw <- function() rnorm(2L,sd=ifelse(runif(2L)<slab_prob,slab,spike))
  combos <- as.matrix(expand.grid(z1=0:1,z2=0:1))
  m <- matrix(NA_real_,4,2); S <- vector("list",4); lw <- numeric(4)
  xy <- as.numeric(crossprod(X,y)/sigma^2)
  for(j in seq_len(4L)) {
    variances <- ifelse(combos[j,]==1,slab^2,spike^2)
    A <- diag(1/variances,2)+crossprod(X)/sigma^2
    C <- solve(A); mm <- as.numeric(C%*%xy)
    m[j,] <- mm; S[[j]] <- C
    logdetA <- as.numeric(determinant(A,logarithm=TRUE)$modulus)
    lw[j] <- sum(ifelse(combos[j,]==1,log(slab_prob),log(1-slab_prob))) -
      n/2*log(2*pi*sigma^2)-sum(y^2)/(2*sigma^2) +
      .5*sum(mm*(A%*%mm)) - .5*(sum(log(variances))+logdetA)
  }
  ww <- exp(lw-logsum(lw))
  cvec <- c(1,-1)
  delta_m <- as.numeric(m%*%cvec)
  delta_s <- vapply(S,function(s) sqrt(as.numeric(crossprod(cvec,s%*%cvec))),numeric(1))
  p_gt <- sum(ww*(1-pnorm(1,delta_m,delta_s)))
  p_orient <- sum(ww*(1-pnorm(0,delta_m,delta_s)))
  ref <- c(beta1=sum(ww*m[,1]),beta2=sum(ww*m[,2]),
           contrast=sum(ww*delta_m),p_contrast_gt_1=p_gt,
           p_beta1_gt_beta2=p_orient)
  metrics <- function(Z) {
    contrast <- Z[,1]-Z[,2]
    c(beta1=mean(Z[,1]),beta2=mean(Z[,2]),
      contrast=mean(contrast),p_contrast_gt_1=mean(contrast>1),
      p_beta1_gt_beta2=mean(contrast>0))
  }
  grid <- seq(-4,5,length.out=1201L)
  density1 <- vapply(grid,function(g)
    sum(ww*vapply(seq_len(4L),function(j)
      dnorm(g,m[j,1],sqrt(S[[j]][1,1])),numeric(1))),numeric(1))
  list(name="CorrelatedPredictors",d=2L,logpost=lp,
       logprior=lprior,prior_draw=prior_draw,reference=ref,
       metric_scales=c(1.5,1.5,2,.5,.5),metrics=metrics,
       modes=function(Z) ifelse(Z[,1]>Z[,2],2L,1L),
       mode_probs=c(1-p_orient,p_orient),mode_type="coefficient_explanation",
       show_coord=1L,ref_grid=grid,ref_pdf=density1,
       description="Normal regression; two predictors correlated ~0.9995; independent spike/slab priors",
       data=list(n=n,X=X,y=y,sigma=sigma),
       param_draws=function(Z) {
         Z<-as.matrix(Z)
         cbind(beta1=Z[,1],beta2=Z[,2],contrast=Z[,1]-Z[,2])
       },
       display=function(Z) Z[,1])
}

# ---- Scenario 6: 40-dimensional nearly collinear Bayesian regression ------
# A REAL likelihood-based posterior: y = X beta + epsilon, with a proper
# two-component spike/slab mixture prior on beta[1:2] and regular Gaussian
# priors on beta[3:40]. Integrating over four inclusion configurations
# gives an EXACT four-Gaussian mixture posterior reference in 40 dimensions.
# True data-generating coefficient contrasts and posterior references differ:
# we score all algorithms against the Bayesian posterior, not the true beta.
make_highdim <- function() {
  set.seed(4262)
  n <- 72L; d <- 40L; sigma <- .90
  root <- rnorm(n)
  X <- matrix(NA_real_,n,d)
  X[,1] <- as.numeric(scale(root))
  X[,2] <- as.numeric(scale(.9992*X[,1] + sqrt(1-.9992^2)*rnorm(n)))
  groups <- replicate(5L,rnorm(n))
  for(j in 3:d) {
    g <- ((j-3L) %% 5L)+1L
    X[,j] <- as.numeric(scale(.55*groups[,g] + .84*rnorm(n)))
  }
  beta_truth <- numeric(d)
  beta_truth[c(1L,5L,9L)] <- c(2.15,.55,-.55)
  y <- as.numeric(X %*% beta_truth + rnorm(n,sd=sigma))
  G <- crossprod(X)
  b <- as.numeric(crossprod(X,y))
  yy <- sum(y*y)
  spike <- .10; slab <- 2.6; slab_prob <- .45; rest_sd <- .55
  log_mix <- function(v) vapply(v,function(z) logsum(c(
    log(1-slab_prob)+dnorm(z,0,spike,log=TRUE),
    log(slab_prob)+dnorm(z,0,slab,log=TRUE))),numeric(1))
  lprior <- function(beta) {
    if(length(beta)!=d || any(!is.finite(beta))) return(-Inf)
    sum(log_mix(beta[1:2]))+
      sum(dnorm(beta[-c(1L,2L)],0,rest_sd,log=TRUE))
  }
  lp <- function(beta) {
    lp0 <- lprior(beta)
    if(!is.finite(lp0)) return(-Inf)
    quad <- yy-2*sum(beta*b)+as.numeric(crossprod(beta,G%*%beta))
    lp0 - n/2*log(2*pi*sigma^2)-quad/(2*sigma^2)
  }
  prior_draw <- function() {
    c(rnorm(2L,sd=ifelse(runif(2L)<slab_prob,slab,spike)),
      rnorm(d-2L,sd=rest_sd))
  }
  comb <- as.matrix(expand.grid(c(0L,1L),c(0L,1L)))
  means <- matrix(NA_real_,4L,d)
  S <- vector('list',4L)
  logmass <- numeric(4L)
  for(j in 1:4) {
    vv <- c(ifelse(comb[j,]==1L,slab^2,spike^2),rep(rest_sd^2,d-2L))
    prec <- diag(1/vv,d)+G/sigma^2
    S[[j]] <- solve(prec)
    mm <- as.numeric(S[[j]] %*% (b/sigma^2))
    means[j,] <- mm
    detlog <- as.numeric(determinant(prec,logarithm=TRUE)$modulus)
    logmass[j] <- sum(ifelse(comb[j,]==1L,log(slab_prob),log(1-slab_prob))) -
      .5*(sum(log(vv))+detlog)+
      .5*as.numeric(crossprod(mm,prec%*%mm))
  }
  mass <- exp(logmass-logsum(logmass))
  contrast_mean <- means[,1L]-means[,2L]
  contrast_sd <- vapply(S,function(Q)
    sqrt(Q[1L,1L]+Q[2L,2L]-2*Q[1L,2L]),numeric(1))
  p_positive <- sum(mass*pnorm(0,contrast_mean,contrast_sd,lower.tail=FALSE))
  p_contrast_gt1 <- sum(mass*pnorm(1,contrast_mean,contrast_sd,lower.tail=FALSE))
  ref <- c(beta1=sum(mass*means[,1]), beta2=sum(mass*means[,2]),
           contrast=sum(mass*contrast_mean),
           p_beta1_gt_beta2=p_positive,
           p_contrast_gt_1=p_contrast_gt1,
           beta5=sum(mass*means[,5]))
  metrics <- function(Z) {
    delta <- Z[,1L]-Z[,2L]
    c(beta1=mean(Z[,1L]), beta2=mean(Z[,2L]),
      contrast=mean(delta),p_beta1_gt_beta2=mean(delta>0),
      p_contrast_gt_1=mean(delta>1),beta5=mean(Z[,5L]))
  }
  gg <- seq(-3,5,length.out=1001L)
  reference_pdf <- vapply(gg,function(z) sum(mass*vapply(1:4,function(j)
    dnorm(z,means[j,1L],sqrt(S[[j]][1L,1L])),numeric(1))),numeric(1))
  list(name="HighDimCorrelated",d=d,logpost=lp,logprior=lprior,
       prior_draw=prior_draw,reference=ref,
       metric_scales=c(1.5,1.5,2,.5,.5,.7),metrics=metrics,
       modes=function(Z) ifelse(Z[,1L]>Z[,2L],2L,1L),
       mode_probs=c(1-p_positive,p_positive),
       mode_type="coefficient_explanation_not_density_mode",
       show_coord=1L,ref_grid=gg,ref_pdf=reference_pdf,
       description="40 predictors; two nearly identical columns; sparse Gaussian-mixture prior",
       data=list(n=n,d=d,X=X,y=y,beta_truth=beta_truth,
                 reference_component_weights=mass),
       display=function(Z) Z[,1L],
       param_draws=function(Z) {
         Z <- as.matrix(Z)
         colnames(Z) <- paste0("beta",seq_len(ncol(Z)))
         cbind(Z,contrast=Z[,1L]-Z[,2L])
       })
}

# ---- Scenario 5: two-test, one-population latent-class model ---------------
# The model has 5 unknown probabilities but the observed 2x2 table has 3
# degrees of freedom. Hence it is not likelihood-identifiable. Symmetric
# Beta(2,2) priors produce disease-label orientation symmetry.
# Reference: EXACT enumeration of latent diseased counts in each observed
# test-result cell; conjugate Beta integration. No oracle information is passed
# to the MCMC algorithms.
make_lca <- function() {
  # Counts for (test1,test2)=(11,10,01,00)
  counts <- c(7L,2L,2L,19L)
  N <- sum(counts)
  a <- 2; b <- 2
  # theta = logits(prevalence, sensitivity1, sensitivity2,
  #                       specificity1, specificity2)
  lprior <- function(z) {
    if(length(z)!=5L || any(!is.finite(z))) return(-Inf)
    # Beta(a,b) on probability p transformed to z=logit(p):
    # log p(z) = a*log(p) + b*log(1-p) - log B(a,b)
    logsig <- function(t) -pmax(-t,0)-log1p(exp(-abs(t)))
    sum(a*logsig(z)+b*logsig(-z)-lbeta(a,b))
  }
  cellprobs <- function(z) {
    pp <- plogis(z)
    pi <- pp[1]; se1 <- pp[2]; se2 <- pp[3]
    sp1 <- pp[4]; sp2 <- pp[5]
    c(pi*se1*se2+(1-pi)*(1-sp1)*(1-sp2),
      pi*se1*(1-se2)+(1-pi)*(1-sp1)*sp2,
      pi*(1-se1)*se2+(1-pi)*sp1*(1-sp2),
      pi*(1-se1)*(1-se2)+(1-pi)*sp1*sp2)
  }
  lp <- function(z) {
    lp0 <- lprior(z)
    if(!is.finite(lp0)) return(-Inf)
    p <- cellprobs(z)
    if(any(p<=0) || any(!is.finite(p))) return(-Inf)
    lp0+sum(counts*log(p))
  }
  prior_draw <- function() qlogis(rbeta(5L,a,b))
  # Enumerate latent disease counts k11,k10,k01,k00.
  ks <- as.matrix(do.call(expand.grid,lapply(counts,function(t) 0:t)))
  K <- nrow(ks)
  betaA <- matrix(NA_real_,K,5L)
  betaB <- matrix(NA_real_,K,5L)
  lw <- numeric(K)
  beta_ratio <- function(s,f) lbeta(a+s,b+f)-lbeta(a,b)
  for(i in seq_len(K)) {
    k <- as.integer(ks[i,]); h <- counts-k
    counts_success <- c(sum(k), k[1]+k[2], k[1]+k[3],
                        h[3]+h[4],h[2]+h[4])
    counts_failure <- c(N-sum(k), k[3]+k[4],k[2]+k[4],
                        h[1]+h[2],h[1]+h[3])
    betaA[i,] <- a+counts_success
    betaB[i,] <- b+counts_failure
    lw[i] <- sum(lchoose(counts,k)) +
      sum(vapply(seq_len(5L),function(j)
        beta_ratio(counts_success[j],counts_failure[j]),numeric(1)))
  }
  ww <- exp(lw-logsum(lw))
  means <- betaA/(betaA+betaB)
  Epi <- means[,1]; Ese1 <- means[,2]
  ref <- c(mean_prevalence=sum(ww*Epi),
           p_prevalence_gt_065=sum(ww*(1-pbeta(.65,betaA[,1],betaB[,1]))),
           mean_sensitivity1=sum(ww*Ese1),
           p_prevalence_lt_030=sum(ww*pbeta(.30,betaA[,1],betaB[,1])),
           p_diseased_and_test1_pos=sum(ww*Epi*Ese1))
  metrics <- function(Z) {
    P <- matrix(plogis(Z),nrow=nrow(Z),ncol=ncol(Z))
    c(mean_prevalence=mean(P[,1]),
      p_prevalence_gt_065=mean(P[,1]>.65),
      mean_sensitivity1=mean(P[,2]),
      p_prevalence_lt_030=mean(P[,1]<.30),
      p_diseased_and_test1_pos=mean(P[,1]*P[,2]))
  }
  grid <- seq(.002,.998,length.out=1201L)
  density_pi <- vapply(grid,function(g)
    sum(ww*dbeta(g,betaA[,1],betaB[,1])),numeric(1))
  symmetry <- function(z) {  # swap latent disease labels
    c(-z[1],-z[4],-z[5],-z[2],-z[3])
  }
  prior_check <- prior_draw()
  stopifnot(abs(lp(prior_check)-lp(symmetry(prior_check)))<1e-8,
            abs(sum(ww)-1)<1e-10,
            abs(ref["mean_prevalence"]-.5)<1e-8)
  list(name="DiagnosticLCA",d=5L,logpost=lp,logprior=lprior,
       prior_draw=prior_draw,reference=ref,
       metric_scales=c(.35,.5,.35,.5,.30),metrics=metrics,
       modes=function(Z) ifelse(plogis(Z[,1])>.5,2L,1L),
       mode_probs=c(.5,.5),mode_type="latent_disease_orientation",
       show_coord=1L,ref_grid=grid,ref_pdf=density_pi,
       description="Two imperfect tests, one population; unidentified prevalence and accuracy",
       data=list(counts=counts,reference_components=K),
       param_draws=function(Z) {
         P<-matrix(plogis(Z),nrow=nrow(Z),ncol=ncol(Z))
         colnames(P)<-c("prevalence","Se1","Se2","Sp1","Sp2")
         P
       },
       display=function(Z) plogis(Z[,1]))
}

# ---- Consistent original-parameter transformations for posterior plots -----
complete_target <- function(target) {
  if(is.null(target$param_draws)) {
    if(target$name=="ProductRidge") {
      target$param_draws <- function(Z) {
        Z<-as.matrix(Z);colnames(Z)<-c("theta1","theta2");Z
      }
    } else if(target$name=="WeakSignalFunnel") {
      target$param_draws <- function(Z) {
        Z<-as.matrix(Z)
        colnames(Z)<-c("log_tau",paste0("eta",seq_len(ncol(Z)-1L)))
        cbind(Z,tau=exp(Z[,1]))
      }
    } else if(target$name=="ThreeComponentMixture") {
      target$param_draws <- function(Z) {
        Z<-as.matrix(Z)
        logits<-cbind(Z[,4:5,drop=FALSE],0)
        mx<-apply(logits,1,max)
        e<-exp(logits-mx)
        W<-e/rowSums(e)
        cbind(mu1=Z[,1],mu2=Z[,2],mu3=Z[,3],
              w1=W[,1],w2=W[,2],w3=W[,3])
      }
    }
  }
  if(is.null(target$display)) target$display <- function(Z) Z[,target$show_coord]
  target
}

# ---- Sampler A: 8 independent adaptive random-walk MH chains ---------------
# Adapt during each chain's initial 35% of steps; freeze thereafter.
walk_am <- function(f,target,x,lp,n,burn_frac=.35) {
  d<-target$d
  H<-matrix(NA_real_,n+1L,d);H[1L,]<-x
  start_sd<-if(d>10L) 1.6/sqrt(d) else 1.4
  U<-chol(diag(start_sd^2,d));acc<-0L
  warm_steps<-floor(n*burn_frac)
  for(t in seq_len(n)) {
    if(t>max(60L,12L*d) && t<=warm_steps && t%%75L==0L) {
      S<-cov(H[seq_len(t),,drop=FALSE])
      U<-chol(spd((2.38^2/d)*S+.035*diag(d)))
    }
    y<-draw_norm(x,U);ly<-f$eval(y)
    if(is.finite(ly) && log(runif(1))<ly-lp) {
      x<-y;lp<-ly;acc<-acc+1L
    }
    H[t+1L,]<-x
  }
  list(draws=H[(warm_steps+2L):(n+1L),,drop=FALSE],
       trace=H,accept=acc/n,final_cov=crossprod(U))
}

adaptive_metropolis <- function(target,B,x0) {
  f<-make_counter(target);C<-8L
  allocations<-rep(B%/%C,C)
  allocations[1L]<-allocations[1L]+B-sum(allocations)
  retained<-vector("list",C)
  traces<-vector("list",C)
  a<-numeric(C);covariances<-vector("list",C)
  for(j in seq_len(C)) {
    x<-if(j==1L) x0 else target$prior_draw()
    lp<-f$eval(x)
    z<-walk_am(f,target,x,lp,allocations[j]-1L,burn_frac=.35)
    retained[[j]]<-z$draws
    traces[[j]]<-z$trace
    a[j]<-z$accept;covariances[[j]]<-z$final_cov
  }
  stopifnot(f$get()==B)
  list(draws=do.call(rbind,retained),retained_by_chain=retained,
       all_traces=traces,accept=mean(a),evals=f$get(),
       covariance=covariances,stages=NULL,proposal_K=NA_integer_)
}

# ---- FLOCK samplers: staged exploration, one FINAL posterior chain ---------
novelty_score <- function(candidates,visited,scales,max_visited=180L) {
  candidates<-as.matrix(candidates)
  if(nrow(visited)>max_visited) visited<-visited[
    unique(as.integer(round(seq(1,nrow(visited),length.out=max_visited)))),,drop=FALSE]
  vapply(seq_len(nrow(candidates)),function(j) {
    dx<-sweep(visited,2L,candidates[j,],"-")
    sqrt(min(rowSums(sweep(dx,2L,scales,"/")^2)))
  },numeric(1))
}

# Soft probabilities for all eigenvector directions. No direction is excluded.
soft_direction_weights <- function(lambda,alpha=.5,rho=.7,delta=.20) {
  stopifnot(all(lambda>=0),alpha>=0,rho>=0,rho<=1,delta>0)
  d<-length(lambda)
  raw<-(lambda+delta)^(-alpha)
  (1-rho)/d+rho*raw/sum(raw)
}

# Select one credible candidate stochastically, favoring novelty; no oracle.
select_credible <- function(logp,novelty,max_logdrop=10,pressure=.75) {
  eligible<-which(logp>=max(logp)-max_logdrop)
  if(!length(eligible)) eligible<-which.max(logp)
  if(length(eligible)==1L) return(eligible[1L])
  score<-pressure*novelty[eligible]
  score<-score-max(score)
  eligible[sample.int(length(eligible),size=1L,prob=exp(score))]
}

flock_staged <- function(target,B,x0,kind=c("Novelty","OrthogonalSoft")) {
  kind<-match.arg(kind)
  f<-make_counter(target);d<-target$d
  scout_budget<-floor(.20*B)
  middle_budget<-floor(.40*B)
  histories<-list();scouts<-list();followers<-list();pioneers<-list()
  scout_lps<-list();best_lp<--Inf;best_x<-x0
  nscout<-if(kind=="Novelty") 2L else 1L
  alloc_s<-rep(scout_budget%/%nscout,nscout)
  alloc_s[1L]<-alloc_s[1L]+scout_budget-sum(alloc_s)
  for(j in seq_len(nscout)) {
    x<-if(j==1L) x0 else target$prior_draw()
    lp<-f$eval(x)
    scout_sd<-if(d>10L) (if(kind=="Novelty") 5.0 else 3.0)/sqrt(d)
               else (if(kind=="Novelty") 2.2 else 1.4)
    z<-walk(f,target,x,lp,alloc_s[j]-1L,
            chol(diag(scout_sd^2,d)),
            prior_jump=if(kind=="Novelty") .12 else .08,
            adapt=(kind=="OrthogonalSoft"))
    scouts[[j]]<-z$h;scout_lps[[j]]<-z$lps
    histories[[length(histories)+1L]]<-z$h
    i<-which.max(z$lps)
    if(z$lps[i]>best_lp) {best_lp<-z$lps[i];best_x<-z$h[i,]}
  }
  stopifnot(f$get()==scout_budget)
  visited<-do.call(rbind,scouts)
  scales<-pmax(apply(visited,2,sd),.8)
  eig<-NULL;direction_prob<-NULL
  if(kind=="OrthogonalSoft") {
    # The scout's LOCAL covariance: closest states to its high-density point.
    # Whole-chain covariance could mix different modes and isn't local geometry.
    scaled<-sweep(visited,2,pmax(scales,.5),"/")
    dist2<-rowSums((scaled-matrix(best_x/pmax(scales,.5),nrow(visited),d,
                                     byrow=TRUE))^2)
    choose<-order(dist2)[seq_len(min(nrow(visited),max(d+12L,
                                              floor(.40*nrow(visited)))))]
    Vloc<-spd(if(length(choose)>=d+3L) cov(visited[choose,,drop=FALSE])
              else diag(1,d))
    eig<-eigen(Vloc,symmetric=TRUE)
    direction_prob<-soft_direction_weights(pmax(eig$values,.035))
  }
  mid_alloc<-rep(middle_budget%/%6L,6L)
  mid_alloc[1L]<-mid_alloc[1L]+middle_budget-sum(mid_alloc)
  # For interpretability, exactly three followers then three pioneers.
  for(j in 1:6) {
    n<-mid_alloc[j]
    if(kind=="Novelty") {
      if(j<=3L) {
        k<-(j-1L)%%nscout+1L
        ix<-which.max(scout_lps[[k]])
        x<-scouts[[k]][ix,];lp<-scout_lps[[k]][ix]
        follow_sd<-if(d>10L) 2.4/sqrt(d) else .95
        U<-chol(diag(follow_sd^2,d))
      } else {
        nc<-min(14L,max(4L,n%/%30L))
        cand<-t(replicate(nc,target$prior_draw()))
        lpc<-vapply(seq_len(nc),function(i) f$eval(cand[i,]),numeric(1))
        nov<-novelty_score(cand,visited,scales)
        utility<-nov-.15*pmax(0,max(lpc)-lpc)
        ix<-which.max(utility)
        x<-cand[ix,];lp<-lpc[ix];n<-n-nc
        pioneer_sd<-if(d>10L) 4.5/sqrt(d) else 1.9
        U<-chol(diag(pioneer_sd^2,d))
      }
    } else {
      if(j<=3L) {
        # Followers: investigate locations NOT adequately visited by the scout.
        nc<-min(14L,max(4L,n%/%30L))
        cand<-t(replicate(nc,target$prior_draw()))
        lpc<-vapply(seq_len(nc),function(i) f$eval(cand[i,]),numeric(1))
        nov<-novelty_score(cand,visited,scales)
        ix<-select_credible(lpc,nov,max_logdrop=10,pressure=.7)
        x<-cand[ix,];lp<-lpc[ix];n<-n-nc
        # Ordinary scout-covariance informed local proposal, regularized.
        U<-chol(spd((2.38^2/d)*eig$vectors %*%
                      diag(pmax(eig$values,.035),d) %*% t(eig$vectors)+
                      .08*diag(d)))
      } else {
        # Pioneers: probabilistic LOW-variance eigenvector bias, both signs.
        # Bias is used for INITIALIZATION, not imposed on all later transitions.
        nc<-min(14L,max(4L,n%/%30L))
        cand<-matrix(NA_real_,nc,d)
        radius<-max(.5,sqrt(stats::median(pmax(eig$values,.035))))
        for(i in seq_len(nc)) {
          axis<-sample.int(d,1L,prob=direction_prob)
          sign_dir<-sample(c(-1,1),1L)
          magnitude<-runif(1L,.7,3.0)*radius
          cand[i,]<-best_x+sign_dir*magnitude*eig$vectors[,axis]+rnorm(d,sd=.12*radius)
        }
        lpc<-vapply(seq_len(nc),function(i) f$eval(cand[i,]),numeric(1))
        nov<-novelty_score(cand,visited,scales)
        ix<-select_credible(lpc,nov,max_logdrop=10,pressure=.65)
        x<-cand[ix,];lp<-lpc[ix];n<-n-nc
        U<-chol(spd((2.38^2/d)*eig$vectors %*%
                      diag(pmax(eig$values,.035),d) %*% t(eig$vectors)+
                      .08*diag(d)))
      }
    }
    stopifnot(n>=10L)
    z<-walk(f,target,x,lp,n,U,prior_jump=0,adapt=FALSE)
    histories[[length(histories)+1L]]<-z$h
    if(j<=3L) followers[[j]]<-z$h else pioneers[[j-3L]]<-z$h
    i<-which.max(z$lps)
    if(z$lps[i]>best_lp) {best_lp<-z$lps[i];best_x<-z$h[i,]}
  }
  stopifnot(f$get()==scout_budget+middle_budget)
  q<-build_proposal(target,histories)
  final<-final_sampler(f,target,best_x,best_lp,B-f$get(),q)
  stopifnot(f$get()==B)
  list(draws=final$draws,retained_by_chain=list(final$draws),
       accept=final$accept,evals=f$get(),proposal_K=q$K,
       direction_weights=direction_prob,
       stages=list(scout=scouts,followers=followers,pioneers=pioneers),
       geometry=if(!is.null(eig)) list(values=eig$values,
                                      vectors=eig$vectors) else NULL)
}

# ---- Evaluation: compare estimates to an independently known reference -----
count_switches <- function(X,target) {
  g<-target$modes(X)
  if(length(g)<2L) return(0L)
  sum(g[-1L]!=g[-length(g)])
}

split_rhat_basic <- function(chain_draws,target) {
  # Classical split R-hat on the displayed scalar; diagnostic only.
  vals<-lapply(chain_draws,target$display)
  n<-min(vapply(vals,length,integer(1)))
  m<-floor(n/2L)
  if(m<10L) return(NA_real_)
  M<-do.call(cbind,lapply(vals,function(x)
    cbind(x[seq_len(m)],x[(m+1L):(2L*m)])))
  means<-colMeans(M)
  W<-mean(apply(M,2,var))
  if(!is.finite(W) || W<=0) return(NA_real_)
  Bv<-m*var(means)
  sqrt(((m-1)/m*W+Bv/m)/W)
}

summarize_sampler <- function(object,target,scenario,replicate,method,elapsed,
                              keep_draws=FALSE) {
  draws<-object$draws
  estimate<-target$metrics(draws)
  if(!identical(names(estimate),names(target$reference)))
    stop("Metric names do not agree for ",scenario)
  exact<-target$reference
  error<-estimate-exact
  relative<-error/target$metric_scales
  groups<-as.integer(target$modes(draws))
  counts<-tabulate(groups,nbins=length(target$mode_probs))
  probs<-counts/sum(counts)
  all_chains<-object$retained_by_chain
  nswitch<-sum(vapply(all_chains,count_switches,integer(1),target=target))
  group_l1<-sum(abs(probs-target$mode_probs))
  group_count<-sum(counts>0)
  explored_groups<-integer(0)
  if(!is.null(object$stages)) {
    for(stage in object$stages) for(x in stage) {
      if(is.matrix(x) && nrow(x)>0L)
        explored_groups<-union(explored_groups,as.integer(target$modes(x)))
    }
  }
  missed_from_exploration<-length(setdiff(explored_groups,which(counts>0)))
  rhat<-if(method=="AdaptiveMetropolis") split_rhat_basic(all_chains,target)
         else NA_real_
  method_summary<-data.frame(scenario=scenario,replicate=replicate,
      method=method,evals=object$evals,retained=nrow(draws),
      seconds=elapsed,accept=object$accept,
      norm_rmse=sqrt(mean(relative^2)),
      group_l1=group_l1,group_count=group_count,
      known_groups=length(target$mode_probs),mode_switches=nswitch,
      split_rhat=rhat,exploration_group_count=length(explored_groups),
      explored_but_missing=missed_from_exploration,
      proposal_components=if(is.null(object$proposal_K)) NA_integer_
                              else object$proposal_K)
  metric_table<-data.frame(scenario=scenario,replicate=replicate,method=method,
      metric=names(estimate),estimate=as.numeric(estimate),
      reference=as.numeric(exact),error=as.numeric(error),
      abs_error=abs(as.numeric(error)))
  P<-as.matrix(target$param_draws(draws))
  param_table<-do.call(rbind,lapply(seq_len(ncol(P)),function(j) {
    v<-P[,j]
    data.frame(scenario=scenario,replicate=replicate,method=method,
      parameter=colnames(P)[j],mean=mean(v),sd=sd(v),
      ci_lower=unname(quantile(v,.025)),ci_upper=unname(quantile(v,.975)))
  }))
  group_table<-data.frame(scenario=scenario,replicate=replicate,method=method,
      group=seq_along(probs),true_mass=target$mode_probs,
      estimated_mass=probs,n_visits=counts,
      mode_interpretation=target$mode_type)
  chain_table<-do.call(rbind,lapply(seq_along(all_chains),function(j) {
    XX<-all_chains[[j]]
    g<-target$modes(XX)
    data.frame(scenario=scenario,replicate=replicate,method=method,
      chain=j,retained=nrow(XX),mean_display=mean(target$display(XX)),
      sd_display=sd(target$display(XX)),group_count=length(unique(g)),
      switches=count_switches(XX,target))
  }))
  keep<-if(keep_draws) list(internal=draws,original=P,
                             retained_chains=all_chains,stages=object$stages,
                             geometry=object$geometry) else NULL
  list(summary=method_summary,metrics=metric_table,parameters=param_table,
       groups=group_table,chains=chain_table,keep=keep)
}

# ---- Posterior-reference numerical tests and scenario setup ----------------
set.seed(SEED)
S_test<-matrix(c(1.5,.55,.55,.9),2,2)
U_test<-chol(S_test)
Z_test<-replicate(10000L,draw_norm(c(0,0),U_test))
stopifnot(max(abs(cov(t(Z_test))-S_test))<.10,
          abs(log_norm(c(0,0),c(0,0),U_test)+log(2*pi)+
                .5*log(det(S_test)))<1e-9)
stopifnot(abs(sum(soft_direction_weights(c(2,.1,.01)))-1)<1e-12)
factories<-list(ProductRidge=make_ridge,
                WeakSignalFunnel=make_funnel,
                ThreeComponentMixture=make_mixture,
                CorrelatedPredictors=make_correlated,
                DiagnosticLCA=make_lca,
                HighDimCorrelated=make_highdim)
scenarios<-lapply(factories[scenario_names],function(fun)complete_target(fun()))
for(s in scenarios) {
  stopifnot(all(is.finite(s$reference)),length(s$metric_scales)==length(s$reference),
            all(s$metric_scales>0),length(s$reference)==length(s$metrics(
              matrix(s$prior_draw(),nrow=1L))),
            abs(sum(s$mode_probs)-1)<1e-6)
  x<-s$prior_draw()
  stopifnot(length(x)==s$d,is.finite(s$logpost(x)),
            all(is.finite(s$param_draws(matrix(x,nrow=1L)))))
}
if("ProductRidge" %in% names(scenarios)) stopifnot(
  abs(scenarios$ProductRidge$logpost(c(1,4))-
      scenarios$ProductRidge$logpost(c(-1,-4)))<1e-10)
if("ThreeComponentMixture" %in% names(scenarios)) {
  mi<-scenarios$ThreeComponentMixture
  zz1<-c(-4,0,4,log(.2/.5),log(.3/.5))
  zz2<-c(0,-4,4,log(.3/.5),log(.2/.5))
  stopifnot(abs(mi$logpost(zz1)-mi$logpost(zz2))<1e-9)
}
cat("Model-reference and Gaussian-proposal self-tests passed.\n")

reference_table<-do.call(rbind,lapply(scenarios,function(s)
  data.frame(scenario=s$name,metric=names(s$reference),
             reference=as.numeric(s$reference))))
write.csv(reference_table,file.path(OUTDIR,"reference_metrics.csv"),row.names=FALSE)
reference_marginal<-do.call(rbind,lapply(scenarios,function(s) {
  w<-s$ref_pdf
  dx<-diff(s$ref_grid)[1L]
  cdf<-cumsum(w)*dx
  cdf<-cdf/max(cdf)
  data.frame(scenario=s$name,
    quantile025=s$ref_grid[which(cdf>=.025)[1L]],
    quantile975=s$ref_grid[which(cdf>=.975)[1L]],
    displayed_parameter=if(s$name=="DiagnosticLCA") "prevalence"
      else if(s$name=="WeakSignalFunnel") "log_tau"
      else if(s$name %in% c("CorrelatedPredictors","HighDimCorrelated")) "beta1"
      else if(s$name=="ProductRidge") "theta1" else "mu1")
}))
write.csv(reference_marginal,file.path(OUTDIR,"reference_marginal_intervals.csv"),
          row.names=FALSE)

# ---- Checkpointing ----------------------------------------------------------
checkpoint_path <- function(s,r,m) file.path(OUTDIR,"checkpoints",
  sprintf("%s_rep%03d_%s.rds",s,r,m))
read_all_results <- function() {
  files<-list.files(file.path(OUTDIR,"checkpoints"),pattern="[.]rds$",
                    full.names=TRUE)
  lapply(files,readRDS)
}
write_outputs <- function(completed) {
  if(!length(completed)) return(invisible(NULL))
  for(slot in c("summary","metrics","parameters","groups","chains")) {
    df<-do.call(rbind,lapply(completed,function(x)x[[slot]]))
    write.csv(df,file.path(OUTDIR,paste0(slot,".csv")),row.names=FALSE)
  }
  sm<-do.call(rbind,lapply(completed,function(x)x$summary))
  spl<-split(sm,interaction(sm$scenario,sm$method,drop=TRUE))
  agg<-do.call(rbind,lapply(spl,function(v)
    data.frame(scenario=v$scenario[1],method=v$method[1],replicates=nrow(v),
               mean_norm_rmse=mean(v$norm_rmse),median_norm_rmse=median(v$norm_rmse),
               mean_group_l1=mean(v$group_l1),
               mean_group_coverage=mean(v$group_count/v$known_groups),
               mean_switches=mean(v$mode_switches),
               mean_seconds=mean(v$seconds))))
  agg<-agg[order(agg$scenario,agg$mean_norm_rmse),,drop=FALSE]
  write.csv(agg,file.path(OUTDIR,"method_summary.csv"),row.names=FALSE)
  invisible(agg)
}

for(sidx in seq_along(scenarios)) {
  scn<-scenarios[[sidx]]
  for(r in seq_len(REPS)) {
    set.seed(SEED+100000L*sidx+997L*r)
    x0<-scn$prior_draw()  # identical first start across methods
    for(midx in seq_along(methods)) {
      method<-methods[midx]
      dest<-checkpoint_path(scn$name,r,method)
      if(file.exists(dest)) next
set.seed(SEED+100000L*sidx+997L*r+17L*midx)
      elapsed<-system.time({
        fit<-switch(method,
          AdaptiveMetropolis=adaptive_metropolis(scn,BUDGET,x0),
          FLOCK_Novelty=flock_staged(scn,BUDGET,x0,"Novelty"),
          FLOCK_OrthogonalSoft=flock_staged(scn,BUDGET,x0,"OrthogonalSoft"))
      })[["elapsed"]]
      stopifnot(fit$evals==BUDGET)
      out<-summarize_sampler(fit,scn,scn$name,r,method,elapsed,
                             keep_draws=(r==1L))
      saveRDS(out,dest)
      # Convenient one-chain posterior output for FLOCK, first replicate only.
      if(r==1L) {
        prefix<-file.path(OUTDIR,paste0(scn$name,"_",method))
        if(method!="AdaptiveMetropolis") {
          write.csv(as.data.frame(out$keep$original),
                    paste0(prefix,"_FINAL_posterior.csv"),row.names=FALSE)
        } else {
          long_chains<-do.call(rbind,lapply(seq_along(fit$retained_by_chain),
            function(cid) {
              data.frame(chain=cid,as.data.frame(
                scn$param_draws(fit$retained_by_chain[[cid]])))
            }))
          write.csv(long_chains,paste0(prefix,"_eight_chains.csv"),
                    row.names=FALSE)
        }
        if(method=="FLOCK_OrthogonalSoft") {
          geom<-fit$geometry
          df_geom<-data.frame(eigenvector=seq_along(geom$values),
              eigenvalue=geom$values,
              pioneer_direction_probability=fit$direction_weights,
              geom$vectors)
          colnames(df_geom)[4:ncol(df_geom)]<-
            paste0("component",seq_len(scn$d))
          write.csv(df_geom,paste0(prefix,"_scout_geometry.csv"),
                    row.names=FALSE)
        }
      }
      cat(sprintf("[%s %d/%d] %-23s error=%.4f groups=%d/%d elapsed=%.1fs\n",
            scn$name,r,REPS,method,out$summary$norm_rmse,
            out$summary$group_count,out$summary$known_groups,elapsed))
    }
    write_outputs(read_all_results())
  }
}

completed<-read_all_results()
agg<-write_outputs(completed)

# ---- Paired comparisons over independent runs ------------------------------
summary_df<-do.call(rbind,lapply(completed,function(x)x$summary))
paired<-list()
for(scn in scenario_names) for(other in c("FLOCK_Novelty","FLOCK_OrthogonalSoft")) {
  a<-summary_df[summary_df$scenario==scn &
                  summary_df$method=="AdaptiveMetropolis",]
  b<-summary_df[summary_df$scenario==scn & summary_df$method==other,]
  z<-merge(a[,c("replicate","norm_rmse","group_l1")],
           b[,c("replicate","norm_rmse","group_l1")],by="replicate",
           suffixes=c("_am","_flock"))
  if(!nrow(z)) next
  for(metric in c("norm_rmse","group_l1")) {
    delta<-z[[paste0(metric,"_flock")]]-z[[paste0(metric,"_am")]]
    interval<-c(NA_real_,NA_real_)
    if(length(delta)>=8L) {
      bootstrap<-replicate(2000L,mean(sample(delta,length(delta),replace=TRUE)))
      interval<-unname(quantile(bootstrap,c(.025,.975)))
    }
    paired[[length(paired)+1L]]<-data.frame(scenario=scn,method=other,
      metric=metric,n=length(delta),mean_difference=mean(delta),
      fraction_flock_better=mean(delta<0),
      ci_low=interval[1L],ci_high=interval[2L])
  }
}
if(length(paired)) write.csv(do.call(rbind,paired),
    file.path(OUTDIR,"paired_comparisons.csv"),row.names=FALSE)

# ---- Visualizations: exact marginal density versus retained posterior -------
# NOTE: 8 AM chains are pooled ONLY for visual diagnostics and must pass
# convergence diagnostics before the pooled posterior is trusted.
pdf(file.path(OUTDIR,"posterior_illustrations.pdf"),width=12,height=9)
for(sidx in seq_along(scenarios)) {
  s<-scenarios[[sidx]]
  par(mfrow=c(3,2),mar=c(4,4,3.2,1))
  for(method in methods) {
    chk<-readRDS(checkpoint_path(s$name,1L,method))
    fit<-chk$keep
    dd<-fit$internal
    scalar<-as.numeric(s$display(dd))
    xr<-range(s$ref_grid)
    yr<-range(s$ref_pdf)
    plot(s$ref_grid,s$ref_pdf,type="l",lwd=2,col="black",
         xlab=reference_marginal$displayed_parameter[sidx],
         ylab="Density",xlim=xr,ylim=c(0,max(yr)*1.35),
         main=paste(s$name,method))
    if(length(unique(scalar))>2L) {
      kde<-density(scalar,n=1024)
      lines(kde,col="steelblue",lwd=2)
    }
    legend("topright",c("Reference","Posterior estimate"),
           col=c("black","steelblue"),lty=1,bty="n",cex=.75)
    trace<-if(method=="AdaptiveMetropolis") fit$retained_chains[[1L]] else dd
    plot(s$display(trace),type="l",col="grey35",lwd=.7,
         xlab="Retained draw index",ylab="Displayed parameter",
         main=if(method=="AdaptiveMetropolis")
              "AM: first of 8 chains (others omitted)" else
              "FLOCK: FINAL chain only (after burn-in)")
  }
}
dev.off()

# ---- Visualizations: scout, follower, pioneer discoveries in 2D projections --
pdf(file.path(OUTDIR,"flock_stage_exploration.pdf"),width=11,height=5)
for(s in scenarios) {
  par(mfrow=c(1,2),mar=c(4,4,3,1))
  for(method in c("FLOCK_Novelty","FLOCK_OrthogonalSoft")) {
    k<-readRDS(checkpoint_path(s$name,1L,method))$keep
    pool<-k$stages
    phases<-list(scout=pool$scout,followers=pool$followers,pioneers=pool$pioneers)
    phasepoints<-lapply(phases,function(L) {
      x<-do.call(rbind,L)
      if(nrow(x)>250L) x<-x[unique(as.integer(round(seq(1,nrow(x),length.out=250L)))),,drop=FALSE]
      x[,1:2,drop=FALSE]
    })
    xx<-do.call(rbind,phasepoints)
    finals<-k$internal[,1:2,drop=FALSE]
    if(nrow(finals)>250L) finals<-finals[unique(as.integer(round(
        seq(1,nrow(finals),length.out=250L)))),,drop=FALSE]
    bounds<-rbind(xx,finals)
    plot(NA,xlim=range(bounds[,1]),ylim=range(bounds[,2]),
         xlab="Internal parameter 1",ylab="Internal parameter 2",
         main=paste(s$name,method))
    cols<-c(scout="orange",followers="darkgreen",pioneers="purple")
    for(p in names(phasepoints)) points(phasepoints[[p]],col=cols[p],
                                         pch=16,cex=.28)
    points(finals,col="steelblue",pch=16,cex=.28)
    legend("topright",c("Scout","Followers","Pioneers","Final"),
           col=c(cols,"steelblue"),pch=16,cex=.65,bty="n")
  }
}
dev.off()

cat("\nFinished. Results saved in: ",normalizePath(OUTDIR),"\n",sep="")
print(agg,row.names=FALSE,digits=3)
cat("FLOCK posterior draws: ONLY final post-warmup chain.\n")
cat("AM pooled chains may be invalid if modes have not mixed.\n")
cat("These five targets are illustrative; all methods require convergence validation.\n")