# FLOCK / Dog Walker x Adaptive Metropolis: controlled benchmark
# Standalone, base R only. This is an experimental research prototype.
# Two exact-reference Bayesian targets replicated from flock_softorth_six_scenarios.R:
#   DiagnosticLCA (5D): two imperfect tests, nonidentified prevalence/accuracy
#   HighDimCorrelated (40D): Bayesian correlated-predictor regression
# Same log posterior and prior across every method; common target-call budget.
#
# Methods:
#   IndependentAM: eight independent AM chains (retained chains pooled)
#   StagedAM:      one AM scout + six undirected AM explorers + one corrected final sampler
#   FLOCK_AM:      one AM scout + 3 novelty followers + 3 soft orthogonal pioneers
#                  + one corrected final AM/local-global sampler
#   FLOCK_RWMH:    same formation but fixed local random-walk MH (ablation)
#   DogWalker65:   seven mutually repelling AM walkers, 65% exploration
#   DogWalker90:   seven mutually repelling AM walkers, 90% exploration
#   DogWalker65_NoRepulsion and DogWalker90_NoRepulsion: ablations
#
# IMPORTANT: exploratory histories are NOT final posterior draws.
# IndependentAM has 8 separately burned-in retained chains.
# Staged/FLOCK have 1 final corrected MH chain; final warmup excluded.
# Optional FLOCK_RWMH is fixed local random-walk MH *plus* a corrected
# independence-MH global component: it is NOT purely local RWMH.
# Dog walkers adapt local covariance during exploration; their repulsive
# histories are NOT posterior draws. They share FLOCK's final MH sampler.
# A learned proposal is frozen before final sampling. Local covariance adapts
# ONLY in final warmup. Warmup excluded and global independence proposals
# use q(x)/q(y) correction. Kernel mixing probabilities are fixed.
# The scout's covariance does NOT certify global posterior coverage.
#
# Examples:
#   Sys.setenv(RUN_MODE="smoke", SCENARIOS="DiagnosticLCA"); source("flock_am_formation_comparison.R")
#   Sys.setenv(RUN_MODE="pilot", SCENARIOS="all"); source("flock_am_formation_comparison.R")
#   Sys.setenv(RUN_MODE="full", REPS="60", BUDGET="60000", SCENARIOS="all"); source("flock_am_formation_comparison.R")

RUN_MODE <- Sys.getenv("RUN_MODE", "smoke")
settings <- switch(RUN_MODE,
   smoke=c(2L,4000L),pilot=c(12L,30000L),full=c(60L,60000L),
   stop("RUN_MODE must be smoke, pilot, or full"))
REPS <- as.integer(Sys.getenv("REPS",settings[1]))
BUDGET <- as.integer(Sys.getenv("BUDGET",settings[2]))
SEED <- as.integer(Sys.getenv("SEED","20261025"))
OUTDIR <- Sys.getenv("OUTDIR",sprintf("flock_am_dogwalker_%s_B%d",RUN_MODE,BUDGET))
stopifnot(REPS>=1L,BUDGET>=3000L,SEED>=0L,nzchar(OUTDIR))
dir.create(OUTDIR,recursive=TRUE,showWarnings=FALSE)
dir.create(file.path(OUTDIR,"checkpoints"),recursive=TRUE,showWarnings=FALSE)
all_methods <- c("IndependentAM","StagedAM","FLOCK_AM","FLOCK_RWMH",
                 "DogWalker65","DogWalker65_NoRepulsion",
                 "DogWalker90","DogWalker90_NoRepulsion")
METHODS <- trimws(strsplit(Sys.getenv("METHODS",paste(all_methods,collapse=",")),",",fixed=TRUE)[[1L]])
stopifnot(length(METHODS)>0L,all(METHODS%in%all_methods),!anyDuplicated(METHODS))
selection <- Sys.getenv("SCENARIOS","all")
all_scenarios <- c("DiagnosticLCA","HighDimCorrelated")
SCENARIOS <- if(tolower(selection)=="all") all_scenarios else
 trimws(strsplit(selection,",",fixed=TRUE)[[1L]])
stopifnot(length(SCENARIOS)>0L,all(SCENARIOS%in%all_scenarios),!anyDuplicated(SCENARIOS))
VERSION <- "flock_am_dogwalker_v1_20261025"
manifest <- list(version=VERSION,reps=REPS,budget=BUDGET,seed=SEED,
                 scenarios=SCENARIOS,methods=METHODS)
manifest_path <- file.path(OUTDIR,"manifest.rds")
if(file.exists(manifest_path)) {
  if(!identical(readRDS(manifest_path),manifest))
    stop("Manifest/settings mismatch; choose a new OUTDIR to avoid mixing results")
} else saveRDS(manifest,manifest_path)

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



# ----- Generic Adaptive Metropolis / Random-Walk Metropolis kernel ----------
# Every proposal (even rejected) charges exactly one log-posterior evaluation.
# Optional prior draws form a second independence MH kernel; its correction is
# log prior(x)-log prior(y). Every method uses the same 2% prior jump frequency.
# When q is provided, an additional global independence MH kernel is used.
# During exploratory/final warmup, the covariance and scalar tuning adapt.
# After warmup they are frozen; only post-warmup states are retained.

split_calls <- function(total,k) {
  z <- rep(as.integer(total%/%k),k)
  if(total%%k) z[seq_len(total%%k)] <- z[seq_len(total%%k)]+1L
  stopifnot(sum(z)==total,all(z>=1L))
  z
}

am_kernel <- function(f,target,x,lp,n,adapt_n=0L,adaptive=TRUE,
                      global_q=NULL,p_global=0,prior_prob=.02,
                      initial_U=NULL) {
  d<-target$d
  stopifnot(n>=1L,adapt_n>=0L,adapt_n<=n,p_global>=0,
            prior_prob>=0,p_global+prior_prob<1)
  H<-matrix(NA_real_,n+1L,d);LP<-numeric(n+1L)
  H[1,]<-x;LP[1]<-lp
  sd_start <- if(d>15L) .48/sqrt(d) else .32
  U <- if(is.null(initial_U)) diag(sd_start,d) else initial_U
  lscale<-0
  accepted<-0L; n_local<-0L; acc_local<-0L
  n_prior<-0L; acc_prior<-0L; n_global<-0L; acc_global<-0L
  recent_local<-0L;recent_accepted<-0L
  logqx<-if(!is.null(global_q)) global_q$logq(x) else NA_real_
  for(t in seq_len(n)) {
    choice<-runif(1)
    if(!is.null(global_q) && choice<p_global) {
      typ<-"global";n_global<-n_global+1L
      y<-global_q$draw();ly<-f$eval(y);logqy<-global_q$logq(y)
      la<-ly-lp+logqx-logqy
    } else if(choice<p_global+prior_prob) {
      typ<-"prior";n_prior<-n_prior+1L
      y<-target$prior_draw();ly<-f$eval(y)
      la<-ly-lp+target$logprior(x)-target$logprior(y)
    } else {
      typ<-"local";n_local<-n_local+1L;recent_local<-recent_local+1L
      y<-draw_norm(x,U);ly<-f$eval(y)
      la<-ly-lp
    }
    ok<-is.finite(la) && (log(runif(1))<la)
    if(ok) {
      x<-y;lp<-ly;accepted<-accepted+1L
      if(typ=="local") {acc_local<-acc_local+1L;recent_accepted<-recent_accepted+1L}
      if(typ=="prior") acc_prior<-acc_prior+1L
      if(typ=="global") acc_global<-acc_global+1L
      if(!is.null(global_q)) {
        logqx<-if(typ=="global") logqy else global_q$logq(x)
      }
    }
    H[t+1L,]<-x;LP[t+1L]<-lp
    # Diminishing-ish scalar AM adaptation during warmup ONLY.
    # Uses a diagonal shrinkage and eigenvalue floor for short 40D chains.
    if(adaptive && t<=adapt_n && t>=max(100L,3L*d) && t%%100L==0L) {
      S<-stats::cov(H[seq_len(t+1L),,drop=FALSE])
      if(!all(is.finite(S))) S<-diag(sd_start^2,d)
      diagS<-pmax(diag(S),(.12/sqrt(d))^2)
      Sshrink<-.80*S+.20*diag(diagS,d)
      C<-(2.38^2/d)*Sshrink+diag((.08/sqrt(d))^2,d)
      if(recent_local>5L) {
        recent_rate<-recent_accepted/recent_local
        lscale<-max(-1.2,min(1.2,lscale+.65*(recent_rate-.234)))
      }
      U<-chol(spd(C,1e-7))*exp(lscale)
      recent_local<-0L;recent_accepted<-0L
    }
  }
  list(h=H,lps=LP,x=x,lp=lp,U=U,
       accept=accepted/n,local_accept=if(n_local) acc_local/n_local else NA_real_,
       global_accept=if(n_global) acc_global/n_global else NA_real_,
       prior_accept=if(n_prior) acc_prior/n_prior else NA_real_,
       local_trials=n_local,global_trials=n_global,
       distinct=nrow(unique(H)),warmup=adapt_n)
}

# ----- Standard eight-chain Adaptive Metropolis ---------------------------
independent_am <- function(target,B,x0) {
  f<-make_counter(target);a<-split_calls(B,8L)
  chains<-vector("list",8L);acc<-numeric(8L)
  for(j in seq_len(8L)) {
    x<-if(j==1L) x0 else target$prior_draw()
    lp<-f$eval(x)
    n<-a[j]-1L;warm<-floor(.35*n)
    z<-am_kernel(f,target,x,lp,n,adapt_n=warm,adaptive=TRUE)
    chains[[j]]<-z$h[(warm+2L):(n+1L),,drop=FALSE]
    acc[j]<-z$accept
  }
  stopifnot(f$get()==B)
  list(draws=do.call(rbind,chains),chains=chains,
       accept=mean(acc),local_accept=mean(acc),global_accept=NA_real_,
       evals=f$get(),exploration=NULL,final=NULL,
       learned_components=NA_integer_,model="IndependentAM")
}

# ----- Proposal constructed from exploratory histories ----------------------
# Identical proposal fitting and final MH used in StagedAM and FLOCK_AM.
# It learns a fixed multivariate Gaussian mixture; this mixture is NOT the
# posterior and its component weights are NOT posterior mode probabilities.
# Prior fallback q0 gives full support wherever the prior has support.
exploration_pool <- function(histories,per_chain=130L) {
  do.call(rbind,lapply(histories,function(H) {
    idx<-unique(as.integer(round(seq(1L,nrow(H),length.out=min(nrow(H),per_chain)))))
    H[idx,,drop=FALSE]
  }))
}

make_global_proposal <- function(target,histories) {
  d<-target$d
  pool<-unique(exploration_pool(histories,140L))
  if(nrow(pool)<2L) {
    # The final sampler must remain defined even if exploration failed.
    # This is a diagnostic failure case, not a reason to discard the run.
    pool<-rbind(pool,pool+rep(.01,d))
  }
  K<-min(5L,max(1L,floor(nrow(pool)/(max(30L,3L*d)))))
  # Standardization is for clustering ONLY, to avoid one parameter scale
  # dominating k-means. Fitting/drawing use original parameter coordinates.
  center<-colMeans(pool)
  scale<-pmax(apply(pool,2,stats::sd),.25)
  Z<-sweep(sweep(pool,2,center,"-"),2,scale,"/")
  km<-tryCatch(stats::kmeans(Z,centers=K,nstart=2,iter.max=45),
               error=function(e) NULL)
  if(is.null(km)) {K<-1L;groups<-rep(1L,nrow(pool))}
  else groups<-km$cluster
  means<-lapply(seq_len(K),function(j) colMeans(pool[groups==j,,drop=FALSE]))
  U<-lapply(seq_len(K),function(j) {
    X<-pool[groups==j,,drop=FALSE]
    C<-if(nrow(X)>=d+5L) stats::cov(X) else diag(.12^2,d)
    diagC<-pmax(diag(C),.06^2)
    C<-.80*C+.20*diag(diagC,d)+diag(.035^2,d)
    chol(spd(1.4*C,1e-6))
  })
  prior_weight<-.15
  logq<-function(x) {
    vals<-vapply(seq_len(K),function(j)
      log((1-prior_weight)/K)+log_norm(x,means[[j]],U[[j]]),numeric(1))
    logsum(c(vals,log(prior_weight)+target$logprior(x)))
  }
  draw<-function() {
    if(runif(1)<prior_weight) return(target$prior_draw())
    j<-sample.int(K,1L)
    draw_norm(means[[j]],U[[j]])
  }
  list(draw=draw,logq=logq,K=K)
}

# ----- Scout geometry and candidate choices -------------------------------
score_novelty <- function(candidates,visited) {
  visited<-as.matrix(visited);candidates<-as.matrix(candidates)
  ids<-unique(as.integer(round(seq(1L,nrow(visited),
                     length.out=min(160L,nrow(visited))))))
  visited<-visited[ids,,drop=FALSE]
  scl<-pmax(apply(visited,2,stats::sd),.25)
  vapply(seq_len(nrow(candidates)),function(j) {
    diff<-sweep(visited,2,candidates[j,],"-")
    sqrt(min(rowSums(sweep(diff,2,scl,"/")^2)))
  },numeric(1))
}
soft_direction_prob <- function(eigenvalues,alpha=.5,rho=.7,delta=.15) {
  d<-length(eigenvalues)
  ww<-(pmax(eigenvalues,0)+delta)^(-alpha)
  (1-rho)/d+rho*ww/sum(ww)
}
select_plausible <- function(logp,novelty,drop=12,weight=.7) {
  valid<-which(logp >= max(logp)-drop)
  if(length(valid)==1L) return(valid)
  score<-weight*novelty[valid]
  score<-score-max(score)
  sample(valid,1L,prob=exp(score))
}

# Formation is identical for StagedAM and FLOCK_AM except middle starts.
# The local transition kernel and final global proposal are IDENTICAL.
formation_sampler <- function(target,B,x0,kind=c("StagedAM","FLOCK_AM","FLOCK_RWMH")) {
  kind<-match.arg(kind);d<-target$d
  adaptive<-kind!="FLOCK_RWMH"
  f<-make_counter(target)
  scout_budget<-floor(.15*B)
  middle_budget<-floor(.30*B)
  stopifnot(scout_budget>=200L,middle_budget>=100L)
  histories<-list();lphist<-list();starts<-matrix(NA_real_,6L,d)
  # One initial scout, adapted only while its trajectory is exploratory.
  lp0<-f$eval(x0)
  sc<-am_kernel(f,target,x0,lp0,scout_budget-1L,
                adapt_n=if(adaptive) floor(.65*(scout_budget-1L)) else 0L,
                adaptive=adaptive)
  histories[[1L]]<-sc$h;lphist[[1L]]<-sc$lps
  i_best<-which.max(sc$lps);best_x<-sc$h[i_best,];best_lp<-sc$lps[i_best]
  stopifnot(f$get()==scout_budget)
  # Scout local geometry is learned only from nearby/strongly visited states.
  H<-sc$h
  # Build a LOCAL scout geometry around its best discovered state.
  # This prevents widely separated exploratory states from being treated as
  # one elliptical posterior mode by the orthogonal-initialization rule.
  scout_scl<-pmax(apply(H,2,stats::sd),.20)
  standardized<-sweep(sweep(H,2,best_x,"-"),2,scout_scl,"/")
  delta2<-rowSums(standardized^2)
  neighbor_count<-min(nrow(H),max(d+15L,floor(.45*nrow(H))))
  Xlocal<-H[order(delta2)[seq_len(neighbor_count)],,drop=FALSE]
  S<-if(nrow(Xlocal)>d+5L) stats::cov(Xlocal) else diag(1,d)
  S<-spd(.85*S+.15*diag(pmax(diag(S),.1^2),d),1e-6)
  e<-eigen(S,symmetric=TRUE)
  pdir<-soft_direction_prob(e$values)
  middle_alloc<-split_calls(middle_budget,6L)
  for(j in seq_len(6L)) {
    allot<-middle_alloc[j]
    n_candidates<-min(12L,max(3L,floor(allot/30L)))
    candidates<-matrix(NA_real_,n_candidates,d)
    if(kind=="StagedAM") {
      # Baseline: no covariance-guided exploration. First half local scout
      # starts, second half broad prior starts; highest-density candidate.
      for(k in seq_len(n_candidates)) {
        candidates[k,]<-if(j<=3L) best_x+rnorm(d,sd=if(d>15L) .20 else .55)
                        else target$prior_draw()
      }
    } else if(j<=3L) {
      # FLOCK followers: half dispersed prior candidates, half plausible
      # scout-neighborhood candidates, to avoid wasting all 40D starts in
      # remote prior regions. Select using scout coverage plus density gate.
      for(k in seq_len(n_candidates)) {
        candidates[k,]<-if(k<=floor(n_candidates/2L)) target$prior_draw()
          else best_x+rnorm(d,sd=if(d>15L) .20 else .55)
      }
    } else {
      # FLOCK pioneers: start near high-density scout point and probe soft
      # orthogonal directions chosen probabilistically from the eigensystem.
      radius<-max(.25,sqrt(stats::median(pmax(e$values,.02))))
      for(k in seq_len(n_candidates)) {
        axis<-sample.int(d,1L,prob=pdir)
        sgn<-sample(c(-1,1),1L)
        mag<-runif(1,.5,2.5)*radius
        candidates[k,]<-if(runif(1)<.15) target$prior_draw() else
          best_x+sgn*mag*e$vectors[,axis]+rnorm(d,sd=.12*radius)
      }
    }
    lp_cand<-vapply(seq_len(n_candidates),function(i)
      f$eval(candidates[i,]),numeric(1))
    if(kind=="StagedAM") chosen<-which.max(lp_cand)
    else {
      novelty<-score_novelty(candidates,H)
      chosen<-select_plausible(lp_cand,novelty,drop=12,weight=.7)
    }
    starts[j,]<-candidates[chosen,]
    n_transitions<-allot-n_candidates
    stopifnot(n_transitions>=1L)
    z<-am_kernel(f,target,starts[j,],lp_cand[chosen],n_transitions,
                 adapt_n=if(adaptive) floor(.65*n_transitions) else 0L,
                 adaptive=adaptive)
    histories[[j+1L]]<-z$h;lphist[[j+1L]]<-z$lps
    ib<-which.max(z$lps)
    if(z$lps[ib]>best_lp) {best_lp<-z$lps[ib];best_x<-z$h[ib,]}
  }
  stopifnot(f$get()==scout_budget+middle_budget)
  q<-make_global_proposal(target,histories)
  final_n<-B-f$get()
  # One final chain: 20% burn-in for local AM adaptation, then frozen.
  warm<-floor(.20*final_n)
  zfin<-am_kernel(f,target,best_x,best_lp,final_n,
                  adapt_n=if(adaptive) warm else 0L,
                  adaptive=adaptive,
                  global_q=q,p_global=if(d>15L) .12 else .23,
                  prior_prob=.02)
  draws<-zfin$h[(warm+2L):(final_n+1L),,drop=FALSE]
  stopifnot(nrow(draws)==final_n-warm,f$get()==B)
  list(draws=draws,chains=list(draws),evals=f$get(),
       accept=zfin$accept,local_accept=zfin$local_accept,
       global_accept=zfin$global_accept,global_trials=zfin$global_trials,
       local_trials=zfin$local_trials,final_distinct=nrow(unique(draws)),
       learned_components=q$K,
       exploration=list(scout=sc$h,middle=histories[-1L],starts=starts),
       direction_prob=pdir,scout_eigenvalues=e$values,
       final_trace=zfin$h,final_warmup=warm,method=kind)
}

# ----- DOG WALKER: seven interacting exploratory AM chains ------------------
# Coordinate-wise updates target a REPULSIVE JOINT exploration density:
#   prod_j posterior(x_j) * prod_{i<j} [1-exp(-0.5*r_ij^2)]^lambda,
# where r_ij^2 = sum_k ((x_ik-x_jk)/scale_k)^2 / (d * radius^2).
# This is an exploration heuristic. The walker histories are NEVER posterior
# draws for the original model. Exact conditional MH corrections include the
# change in repulsion and, for independent prior jumps, the proposal ratio.
# Local walker proposals adapt their covariance ONLY during exploration.
# Their final sampler is the SAME local/global corrected MH kernel as FLOCK_AM.
# All expensive log-posterior calls go through the same counted wrapper.

repulsion_log <- function(x, other_states, scales, radius=0.9) {
  if(nrow(other_states)<1L) return(0)
  d<-length(x)
  z<-sweep(sweep(other_states,2L,x,"-"),2L,scales,"/")
  dsq<-rowSums(z*z)/(d*radius^2)
  # Stable log(1-exp(-0.5*dsq)); same floor in numerator/denominator.
  sum(log(pmax(-expm1(-.5*dsq),1e-12)))
}

min_scaled_walker_distance <- function(X,scales) {
  if(nrow(X)<2L) return(NA_real_)
  S<-sweep(X,2L,scales,"/")
  as.numeric(min(stats::dist(S)))/sqrt(ncol(X))
}

dogwalker_sampler <- function(target,B,x0,explore_fraction=.65,
                              repulsion_strength=.35,n_walkers=7L,
                              radius=.90) {
  d<-target$d
  stopifnot(B>=3000L,n_walkers>=2L,repulsion_strength>=0,
            explore_fraction>0,explore_fraction<.95,radius>0)
  f<-make_counter(target)
  exploration_budget<-floor(explore_fraction*B)
  nsteps<-exploration_budget-n_walkers
  stopifnot(nsteps>n_walkers)
  X<-matrix(NA_real_,n_walkers,d)
  X[1L,]<-x0
  for(j in 2:n_walkers) X[j,]<-target$prior_draw()
  logp<-vapply(seq_len(n_walkers),function(j) f$eval(X[j,]),numeric(1))

  # Coordinate scales are constructed from the known prior, then FROZEN.
  # They are NOT derived from target truths or exact reference quantities.
  prior_calibration<-do.call(rbind,lapply(seq_len(48L),function(i) target$prior_draw()))
  scales<-pmax(apply(rbind(X,prior_calibration),2L,stats::sd),.15)
  per_walker<-tabulate((seq_len(nsteps)-1L)%%n_walkers+1L,nbins=n_walkers)
  histories<-lapply(seq_len(n_walkers),function(j) matrix(NA_real_,per_walker[j]+1L,d))
  for(j in seq_len(n_walkers)) histories[[j]][1L,]<-X[j,]

  sd_start<-if(d>15L) .48/sqrt(d) else .32
  Us<-lapply(seq_len(n_walkers),function(j) diag(sd_start,d))
  counts<-integer(n_walkers); accepted<-integer(n_walkers)
  local_trials<-integer(n_walkers);local_accepted<-integer(n_walkers)
  prior_trials<-integer(n_walkers);prior_accepted<-integer(n_walkers)
  recent_trials<-integer(n_walkers);recent_accept<-integer(n_walkers)
  log_scale<-numeric(n_walkers)
  adaptation_limit<-floor(.65*max(per_walker))
  for(t in seq_len(nsteps)) {
    j<-(t-1L)%%n_walkers+1L
    here<-X[j,]
    others<-X[-j,,drop=FALSE]
    current_rep<-if(repulsion_strength>0) repulsion_log(here,others,scales,radius) else 0
    kind<-runif(1L)
    if(kind<.02) {
      # Prior-independence kernel (correct for q(y)=prior(y)).
      which_kernel<-"prior"
      y<-target$prior_draw()
      log_correction<-target$logprior(here)-target$logprior(y)
      prior_trials[j]<-prior_trials[j]+1L
    } else if(kind<.30) {
      # Large but still symmetric AM-covariance-aligned random walk.
      which_kernel<-"large"
      y<-draw_norm(here,2.8*Us[[j]])
      log_correction<-0
    } else {
      which_kernel<-"local"
      y<-draw_norm(here,Us[[j]])
      log_correction<-0
      local_trials[j]<-local_trials[j]+1L
      recent_trials[j]<-recent_trials[j]+1L
    }
    ly<-f$eval(y)
    proposed_rep<-if(repulsion_strength>0) repulsion_log(y,others,scales,radius) else 0
    la<-ly-logp[j]+repulsion_strength*(proposed_rep-current_rep)+log_correction
    ok<-is.finite(la) && log(runif(1L))<la
    if(ok) {
      X[j,]<-y;logp[j]<-ly;accepted[j]<-accepted[j]+1L
      if(which_kernel=="local") {
        local_accepted[j]<-local_accepted[j]+1L
        recent_accept[j]<-recent_accept[j]+1L
      }
      if(which_kernel=="prior") prior_accepted[j]<-prior_accepted[j]+1L
    }
    counts[j]<-counts[j]+1L
    histories[[j]][counts[j]+1L,]<-X[j,]
    if(counts[j]<=adaptation_limit && counts[j]>=max(100L,3L*d) &&
       counts[j]%%100L==0L) {
      H<-histories[[j]][seq_len(counts[j]+1L),,drop=FALSE]
      C<-stats::cov(H)
      if(all(is.finite(C))) {
        diagC<-pmax(diag(C),(.12/sqrt(d))^2)
        Cshrink<-.80*C+.20*diag(diagC,d)
        Cprop<-(2.38^2/d)*Cshrink+diag((.08/sqrt(d))^2,d)
        if(recent_trials[j]>5L) {
          rate<-recent_accept[j]/recent_trials[j]
          log_scale[j]<-max(-1.2,min(1.2,log_scale[j]+.65*(rate-.234)))
        }
        Us[[j]]<-chol(spd(Cprop,1e-7))*exp(log_scale[j])
      }
      recent_trials[j]<-0L;recent_accept[j]<-0L
    }
  }
  stopifnot(f$get()==exploration_budget,all(counts==per_walker))
  q<-make_global_proposal(target,histories)
  best<-which.max(logp)
  final_n<-B-f$get()
  warm<-floor(.20*final_n)
  final<-am_kernel(f,target,X[best,],logp[best],final_n,
                   adapt_n=warm,adaptive=TRUE,
                   global_q=q,p_global=if(d>15L) .12 else .23,
                   prior_prob=.02)
  draws<-final$h[(warm+2L):(final_n+1L),,drop=FALSE]
  stopifnot(f$get()==B,nrow(draws)==final_n-warm)
  list(draws=draws,chains=list(draws),evals=f$get(),
       accept=final$accept,local_accept=final$local_accept,
       global_accept=final$global_accept,
       global_trials=final$global_trials,local_trials=final$local_trials,
       final_distinct=nrow(unique(draws)),learned_components=q$K,
       exploration=list(walkers=histories,starts=X,scales=scales),
       walker_accept=sum(accepted)/nsteps,
       walker_local_accept=if(sum(local_trials)) sum(local_accepted)/sum(local_trials) else NA_real_,
       walker_prior_accept=if(sum(prior_trials)) sum(prior_accepted)/sum(prior_trials) else NA_real_,
       walker_min_separation=min_scaled_walker_distance(X,scales),
       repulsion_strength=repulsion_strength,
       explore_fraction=explore_fraction,
       final_trace=final$h,final_warmup=warm,
       method=paste0("DogWalker",round(100*explore_fraction),
                     if(repulsion_strength==0) "_NoRepulsion" else ""))
}

# ----- Diagnostics: consequential posterior quantities vs exact reference ----
mode_switches <- function(Z,target) {
  g<-as.integer(target$modes(Z))
  if(length(g)<2L) return(0L)
  sum(diff(g)!=0L)
}
ess_simple <- function(v,maxlag=150L) {
  # Conservative truncated autocorrelation estimate for the displayed scalar.
  v<-as.numeric(v); n<-length(v)
  if(n<50L || stats::sd(v)<1e-13) return(NA_real_)
  ac<-as.numeric(stats::acf(v,plot=FALSE,lag.max=min(maxlag,n%/%3L))$acf)[-1L]
  if(!length(ac)) return(n)
  # Initial positive sequence of adjacent pair autocorrelation sums.
  k<-floor(length(ac)/2L)
  pairs<-vapply(seq_len(k),function(j) ac[2L*j-1L]+ac[2L*j],numeric(1))
  valid<-which(pairs<=0)
  keep<-if(length(valid)) min(valid)-1L else k
  tau<-1+2*sum(ac[seq_len(2L*keep)])
  min(n,max(1,n/max(tau,1)))
}

summarize_fit <- function(fit,target,scenario,replicate,method,seconds) {
  est<-target$metrics(fit$draws)
  ref<-target$reference
  if(!identical(names(est),names(ref))) stop("Metric/reference mismatch")
  err<-est-ref
  norm_rmse<-sqrt(mean((err/target$metric_scales)^2))
  gids<-as.integer(target$modes(fit$draws))
  pp<-tabulate(gids,nbins=length(target$mode_probs))/length(gids)
  per_chain<-lapply(seq_along(fit$chains),function(j) {
    z<-fit$chains[[j]]
    data.frame(chain=j,posterior_mean=mean(target$display(z)),
               ess_display=ess_simple(target$display(z)),
               mode_switches=mode_switches(z,target),
               n_retained=nrow(z),stringsAsFactors=FALSE)
  })
  chains<-do.call(rbind,per_chain)
  chains$scenario<-scenario;chains$replicate<-replicate;chains$method<-method
  # Explicit stage coverage: earlier exploratory draws remain diagnostic ONLY.
  stage_rows<-list()
  if(!is.null(fit$exploration)) {
    if(!is.null(fit$exploration$walkers)) {
      stage_histories<-list(walkers=do.call(rbind,fit$exploration$walkers),
                            final=fit$draws)
    } else {
      stage_histories<-list(scout=fit$exploration$scout,
                            middle=do.call(rbind,fit$exploration$middle),
                            final=fit$draws)
    }
    for(label in names(stage_histories)) {
      zz<-stage_histories[[label]]
      gp<-tabulate(as.integer(target$modes(zz)),nbins=length(target$mode_probs))
      stage_rows[[label]]<-data.frame(scenario=scenario,replicate=replicate,
         method=method,stage=label,group=seq_along(gp),
         count=gp,visited=as.integer(gp>0L),reference_mass=target$mode_probs)
    }
  } else {
    stage_rows$final<-data.frame(scenario=scenario,replicate=replicate,
       method=method,stage="final_pooled",group=seq_along(pp),
       count=tabulate(gids,nbins=length(pp)),visited=as.integer(pp>0),
       reference_mass=target$mode_probs)
  }
  stage_df<-do.call(rbind,stage_rows)
  explorer_ids<-unique(stage_df$group[stage_df$stage %in% c("scout","middle","walkers") &
                                         stage_df$visited==1L])
  missing_ids<-setdiff(explorer_ids,which(pp>0))
  sumrow<-data.frame(scenario=scenario,replicate=replicate,method=method,
      evals=fit$evals,seconds=seconds,retained=nrow(fit$draws),
      normalized_rmse=norm_rmse,
      group_l1=sum(abs(pp-target$mode_probs)),
      group_count=sum(pp>0),
      explored_but_missing=if(!is.null(fit$exploration)) length(missing_ids) else NA_integer_,
      transitions=sum(chains$mode_switches),
      mean_display_ess=sum(chains$ess_display,na.rm=TRUE),
      local_accept=fit$local_accept,global_accept=fit$global_accept,
      final_distinct=if(is.null(fit$final_distinct)) NA_integer_ else fit$final_distinct,
      learned_components=if(is.null(fit$learned_components)) NA_integer_
                         else fit$learned_components,
      walker_accept=if(is.null(fit$walker_accept)) NA_real_ else fit$walker_accept,
      walker_min_separation=if(is.null(fit$walker_min_separation)) NA_real_
                            else fit$walker_min_separation,
      explorer_fraction=if(is.null(fit$explore_fraction)) NA_real_ else fit$explore_fraction,
      repulsion_strength=if(is.null(fit$repulsion_strength)) NA_real_ else fit$repulsion_strength)
  met<-data.frame(scenario=scenario,replicate=replicate,method=method,
                  metric=names(est),estimate=as.numeric(est),
                  reference=as.numeric(ref),error=as.numeric(err),
                  abs_error=abs(as.numeric(err)))
  gr<-data.frame(scenario=scenario,replicate=replicate,method=method,
                 group=seq_along(pp),reference_mass=target$mode_probs,
                 estimated_mass=pp,mode_type=target$mode_type)
  Z<-as.matrix(target$param_draws(fit$draws))
  params<-do.call(rbind,lapply(seq_len(ncol(Z)),function(j) {
    v<-Z[,j]
    data.frame(scenario=scenario,replicate=replicate,method=method,
      parameter=colnames(Z)[j],mean=mean(v),sd=stats::sd(v),
      lower95=unname(stats::quantile(v,.025)),upper95=unname(stats::quantile(v,.975)))
  }))
  list(summary=sumrow,metrics=met,groups=gr,chains=chains,
       parameters=params,stage_coverage=stage_df)
}

# ----- Model setup + diagnostics verification ------------------------------
set.seed(SEED)
Tcov<-matrix(c(1,.65,.65,1.7),2,2)
Tchol<-chol(Tcov)
Tdraws<-replicate(8000L,draw_norm(c(0,0),Tchol))
stopifnot(max(abs(stats::cov(t(Tdraws))-Tcov))<.10,
          abs(log_norm(c(0,0),c(0,0),Tchol)+log(2*pi)+.5*log(det(Tcov)))<1e-10,
          abs(sum(soft_direction_prob(c(2,.1,.01)))-1)<1e-10)
# Dog Walker potential must penalize collisions while remaining finite.
check_others<-rbind(c(0,0),c(1,0))
stopifnot(is.finite(repulsion_log(c(0,0),check_others,c(1,1))),
          repulsion_log(c(0,0),check_others,c(1,1)) <
          repulsion_log(c(4,4),check_others,c(1,1)))
scenario_factories<-list(DiagnosticLCA=make_lca,HighDimCorrelated=make_highdim)
scenarios<-lapply(scenario_factories[SCENARIOS],function(make) complete_target(make()))
for(target in scenarios) {
 stopifnot(length(target$reference)==length(target$metric_scales),
           abs(sum(target$mode_probs)-1)<1e-8)
 z<-target$prior_draw()
 stopifnot(length(z)==target$d,is.finite(target$logpost(z)))
}
cat("Reference-model and proposal-Cholesky checks passed.\n")
write.csv(do.call(rbind,lapply(scenarios,function(t)
  data.frame(scenario=t$name,metric=names(t$reference),
             reference=as.numeric(t$reference)))),
  file.path(OUTDIR,"reference_metrics.csv"),row.names=FALSE)

checkpoint_path<-function(s,r,m) file.path(OUTDIR,"checkpoints",
                         sprintf("%s_rep%03d_%s.rds",s,r,m))
atomic_rds<-function(obj,path) {
 tmp<-paste0(path,".tmp")
 saveRDS(obj,tmp)
 if(!file.rename(tmp,path)) stop("Unable to atomically create checkpoint: ",path)
 invisible(path)
}
all_checkpoints<-function() {
 fs<-list.files(file.path(OUTDIR,"checkpoints"),pattern="[.]rds$",full.names=TRUE)
 lapply(fs,readRDS)
}
export_tables<-function() {
 recs<-all_checkpoints()
 if(!length(recs)) return(invisible(NULL))
 for(field in c("summary","metrics","groups","chains","parameters","stage_coverage")) {
   z<-do.call(rbind,lapply(recs,function(r) r[[field]]))
   write.csv(z,file.path(OUTDIR,paste0(field,".csv")),row.names=FALSE)
 }
 S<-do.call(rbind,lapply(recs,function(r) r$summary))
 by<-split(S,interaction(S$scenario,S$method,drop=TRUE))
 sm<-do.call(rbind,lapply(by,function(z)
   data.frame(scenario=z$scenario[1],method=z$method[1],reps=nrow(z),
              mean_rmse=mean(z$normalized_rmse),
              median_rmse=stats::median(z$normalized_rmse),
              mean_group_l1=mean(z$group_l1),
              mean_seconds=mean(z$seconds),
              mean_transitions=mean(z$transitions),
              mean_final_global_accept=if(all(is.na(z$global_accept))) NA_real_
                 else mean(z$global_accept,na.rm=TRUE),
              mean_walker_accept=if(all(is.na(z$walker_accept))) NA_real_
                 else mean(z$walker_accept,na.rm=TRUE),
              mean_walker_min_separation=if(all(is.na(z$walker_min_separation))) NA_real_
                 else mean(z$walker_min_separation,na.rm=TRUE),
              mean_final_distinct=if(all(is.na(z$final_distinct))) NA_real_
                else mean(z$final_distinct,na.rm=TRUE))))
 write.csv(sm,file.path(OUTDIR,"method_summary.csv"),row.names=FALSE)
 invisible(sm)
}

# ----- Main paired repeated-benchmark loop ---------------------------------
cat(sprintf("FLOCK formation benchmark | %s | %d replicates | budget %d\n",
            RUN_MODE,REPS,BUDGET))
cat("Scenarios:",paste(SCENARIOS,collapse=", "),"\nMethods:",paste(METHODS,collapse=", "),"\n")
for(sidx in seq_along(scenarios)) {
 target<-scenarios[[sidx]]
 for(r in seq_len(REPS)) {
   set.seed(SEED+100003L*sidx+997L*r)
   x0<-target$prior_draw()
   for(method in METHODS) {
     ck<-checkpoint_path(target$name,r,method)
     if(file.exists(ck)) next
     # Same pre-run seed and shared x0 for comparable setups. Methods diverge
     # as they consume different random numbers, as expected.
     set.seed(SEED+100003L*sidx+997L*r+4242L)
     elapsed<-system.time({
       fit<-switch(method,
         IndependentAM=independent_am(target,BUDGET,x0),
         StagedAM=formation_sampler(target,BUDGET,x0,"StagedAM"),
         FLOCK_AM=formation_sampler(target,BUDGET,x0,"FLOCK_AM"),
         FLOCK_RWMH=formation_sampler(target,BUDGET,x0,"FLOCK_RWMH"),
         DogWalker65=dogwalker_sampler(target,BUDGET,x0,.65,.35),
         DogWalker65_NoRepulsion=dogwalker_sampler(target,BUDGET,x0,.65,0),
         DogWalker90=dogwalker_sampler(target,BUDGET,x0,.90,.35),
         DogWalker90_NoRepulsion=dogwalker_sampler(target,BUDGET,x0,.90,0))
     })[["elapsed"]]
     stopifnot(fit$evals==BUDGET)
     out<-summarize_fit(fit,target,target$name,r,method,elapsed)
     atomic_rds(out,ck)
     if(r==1L) {
       prefix<-file.path(OUTDIR,paste0(target$name,"_",method))
       if(method=="IndependentAM") {
         draws<-do.call(rbind,lapply(seq_along(fit$chains),function(i)
            data.frame(chain=i,target$param_draws(fit$chains[[i]]))))
         write.csv(draws,paste0(prefix,"_eight_retained_chains.csv"),row.names=FALSE)
       } else {
         write.csv(as.data.frame(target$param_draws(fit$draws)),
                   paste0(prefix,"_FINAL_posterior.csv"),row.names=FALSE)
         if(!is.null(fit$scout_eigenvalues)) {
           write.csv(data.frame(eigenaxis=seq_along(fit$scout_eigenvalues),
                     eigenvalue=fit$scout_eigenvalues,
                     soft_orthogonal_probability=fit$direction_prob),
                     paste0(prefix,"_scout_eigenvalues.csv"),row.names=FALSE)
         }
         if(!is.null(fit$exploration$walkers)) {
           write.csv(data.frame(walker=seq_len(nrow(fit$exploration$starts)),
                     fit$exploration$starts),
                     paste0(prefix,"_walkers_final_positions.csv"),row.names=FALSE)
         }
       }
     }
     cat(sprintf("[%s %d/%d] %-15s error=%.4f | group L1=%.4f | evals=%d | %.1fs\n",
       target$name,r,REPS,method,out$summary$normalized_rmse,
       out$summary$group_l1,fit$evals,elapsed))
   }
   export_tables()
 }
}

# ----- Paired differences in normalized error; negative favors the method ---
completed<-all_checkpoints()
summary_df<-do.call(rbind,lapply(completed,function(z) z$summary))
paired_rows<-list()
for(s in SCENARIOS) for(other in setdiff(METHODS,"IndependentAM")) {
  a<-summary_df[summary_df$scenario==s & summary_df$method==other,
                c("replicate","normalized_rmse")]
  b<-summary_df[summary_df$scenario==s & summary_df$method=="IndependentAM",
                c("replicate","normalized_rmse")]
  p<-merge(a,b,by="replicate",suffixes=c("_method","_independent"))
  if(nrow(p)) {
   diff<-p$normalized_rmse_method-p$normalized_rmse_independent
   set.seed(SEED+length(diff)+nchar(other)+nchar(s))
   bootstrap<-replicate(1200L,mean(sample(diff,length(diff),replace=TRUE)))
   paired_rows[[length(paired_rows)+1L]]<-data.frame(scenario=s,method=other,
     n=nrow(p),mean_difference=mean(diff),
     bootstrap_low=unname(stats::quantile(bootstrap,.025)),
     bootstrap_high=unname(stats::quantile(bootstrap,.975)),
     fraction_better=mean(diff<0),median_difference=stats::median(diff))
  }
}
if(length(paired_rows)) write.csv(do.call(rbind,paired_rows),
                      file.path(OUTDIR,"paired_vs_independent_am.csv"),row.names=FALSE)
# Mechanistic ablation: SAME AM kernels/final sampler, but undirected
# versus deliberately structured middle-stage initialization.
if(all(c("StagedAM","FLOCK_AM") %in% METHODS)) {
  comparisons<-list()
  for(s in SCENARIOS) {
    aa<-summary_df[summary_df$scenario==s & summary_df$method=="FLOCK_AM",
                   c("replicate","normalized_rmse")]
    bb<-summary_df[summary_df$scenario==s & summary_df$method=="StagedAM",
                   c("replicate","normalized_rmse")]
    zz<-merge(aa,bb,by="replicate",suffixes=c("_flock","_staged"))
    if(nrow(zz)>0L) {
      diff<-zz$normalized_rmse_flock-zz$normalized_rmse_staged
      set.seed(SEED+length(diff)+nchar(s)+457L)
      bs<-replicate(1200L,mean(sample(diff,length(diff),replace=TRUE)))
      comparisons[[length(comparisons)+1L]]<-data.frame(
         scenario=s,n=nrow(zz),mean_difference=mean(diff),
         bootstrap_low=unname(stats::quantile(bs,.025)),
         bootstrap_high=unname(stats::quantile(bs,.975)),
         median_difference=stats::median(diff),fraction_flock_better=mean(diff<0))
}
  }
  if(length(comparisons)) write.csv(do.call(rbind,comparisons),
        file.path(OUTDIR,"paired_flock_am_vs_staged_am.csv"),row.names=FALSE)
}

# Paired Dog Walker ablations: negative difference favors repulsion.
repulsion_comparisons<-list()
for(s in SCENARIOS) for(pct in c("65","90")) {
  repelled<-paste0("DogWalker",pct)
  control<-paste0("DogWalker",pct,"_NoRepulsion")
  if(!all(c(repelled,control)%in%METHODS)) next
  for(metric in c("normalized_rmse","group_l1")) {
    aa<-summary_df[summary_df$scenario==s & summary_df$method==repelled,
                   c("replicate",metric)]
    bb<-summary_df[summary_df$scenario==s & summary_df$method==control,
                   c("replicate",metric)]
    xy<-merge(aa,bb,by="replicate",suffixes=c("_repelled","_control"))
    if(!nrow(xy)) next
    difference<-xy[[paste0(metric,"_repelled")]]-xy[[paste0(metric,"_control")]]
    set.seed(SEED+101L*nrow(xy)+nchar(s)+nchar(metric)+as.integer(pct))
    bs<-replicate(1200L,mean(sample(difference,length(difference),replace=TRUE)))
    repulsion_comparisons[[length(repulsion_comparisons)+1L]]<-data.frame(
      scenario=s,exploration_fraction=as.integer(pct)/100,metric=metric,
      n=nrow(xy),mean_difference=mean(difference),
      bootstrap_low=unname(stats::quantile(bs,.025)),
      bootstrap_high=unname(stats::quantile(bs,.975)),
      fraction_repelled_better=mean(difference<0))
  }
}
if(length(repulsion_comparisons)) write.csv(do.call(rbind,repulsion_comparisons),
  file.path(OUTDIR,"paired_dogwalker_repulsion.csv"),row.names=FALSE)

# ----- PDF diagnostics from first replicate retained posterior samples -----
pdf(file.path(OUTDIR,"posterior_diagnostics.pdf"),width=13,height=9)
for(s in SCENARIOS) {
  target<-scenarios[[s]]
  files<-lapply(METHODS,function(m) {
    f<-file.path(OUTDIR,paste0(s,"_",m,
         if(m=="IndependentAM") "_eight_retained_chains.csv" else "_FINAL_posterior.csv"))
    if(file.exists(f)) read.csv(f,check.names=FALSE) else NULL
  })
  names(files)<-METHODS
  if(any(vapply(files,is.null,logical(1)))) next
  par(mfrow=c(2,2),mar=c(4,4,3,1))
  # Reference scalar is prevalence for LCA, beta1 for 40D.
  par(mar=c(4,4,3,1))
  plot(target$ref_grid,target$ref_pdf,type="l",lwd=2,
       xlab=if(s=="DiagnosticLCA") "prevalence" else "beta1",
       ylab="Posterior density",main=paste0(s," | exact vs MCMC"))
  palette<-grDevices::rainbow(length(METHODS),s=.75,v=.85)
  for(j in seq_along(METHODS)) {
    Z<-files[[j]]
    v<-if(s=="DiagnosticLCA") Z$prevalence else Z$beta1
    if(length(unique(v))<3L) {
      abline(v=mean(v),col=palette[j],lty=3)
    } else {
      de<-stats::density(v,from=min(target$ref_grid),to=max(target$ref_grid),n=1024)
      lines(de$x,de$y,col=palette[j],lwd=1.5)
    }
  }
  legend("topright",legend=c("Exact",METHODS),
         col=c("black",palette[seq_along(METHODS)]),lty=1,cex=.72,bty="n")
  illustrative_method<-if("FLOCK_AM"%in%METHODS) "FLOCK_AM" else METHODS[1L]
  Z<-files[[illustrative_method]]
  val<-if(s=="DiagnosticLCA") Z$prevalence else Z$beta1
  plot(val,type="l",main=paste0(s," | final sampled trajectory"),
       xlab="Retained draw index",ylab=if(s=="DiagnosticLCA") "prevalence" else "beta1")
  if(s=="DiagnosticLCA") {
    plot(Z$prevalence,Z$Se1,pch=16,cex=.25,
         xlab="Prevalence",ylab="Sensitivity 1",main="Parameter relationship")
  } else {
    plot(Z$beta1,Z$beta2,pch=16,cex=.25,
         xlab="beta1",ylab="beta2",main="Competing coefficient explanations")
  }
  sub<-summary_df[summary_df$scenario==s&summary_df$replicate==1L,]
  barplot(sub$normalized_rmse,names.arg=sub$method,las=2,
          main="Posterior-functional error",ylab="Normalized RMSE")
}
dev.off()
cat("Saved summaries, checkpoints, posterior CSVs, and diagnostics under:",OUTDIR,"\n")