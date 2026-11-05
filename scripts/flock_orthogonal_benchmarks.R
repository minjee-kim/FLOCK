# ============================================================================
# FLOCK-Novelty vs FLOCK-Orthogonal vs Adaptive Metropolis
# Three different difficult Bayesian posteriors. Base R only.
#
# Run:
#   Sys.setenv(RUN_MODE="smoke"); source("flock_orthogonal_benchmarks.R")
#   Sys.setenv(RUN_MODE="pilot"); source("flock_orthogonal_benchmarks.R")
#   Sys.setenv(RUN_MODE="full");  source("flock_orthogonal_benchmarks.R")
# Optional: REPS, BUDGET, OUTDIR, SEED. Defaults: smoke (2,3000),
# pilot (10,20000), full (100,60000) for reps and calls per method per trial.
#
# The FLOCK strategies use exploration for proposals, NOT posterior estimates.
# All posterior summary estimates for FLOCK use ONLY final MH post-burn-in.
# AM uses eight independent locally adaptive MH chains, with frozen proposals
# after burn-in, as a conventional multi-start comparison. Combining chains
# that have not mixed is NOT a valid way to repair nonconvergence.
#
# The benchmark compares full algorithms, NOT an isolated orthogonal direction
# effect: Novelty has two scouts, Orthogonal has one adaptive scout. This is a
# pilot experiment, not published-method reproduction or a theorem.
# ============================================================================

RUN_MODE <- Sys.getenv('RUN_MODE','smoke')
settings <- switch(RUN_MODE,
  smoke=c(2L,3000L), pilot=c(10L,20000L), full=c(100L,60000L),
  stop('RUN_MODE must be smoke, pilot or full'))
REPS <- as.integer(Sys.getenv('REPS',settings[1]))
BUDGET <- as.integer(Sys.getenv('BUDGET',settings[2]))
SEED <- as.integer(Sys.getenv('SEED','20261025'))
OUTDIR <- Sys.getenv('OUTDIR',sprintf('flock_orthogonal_%s_B%d',RUN_MODE,BUDGET))
stopifnot(REPS>=1L,BUDGET>=2500L,SEED>=0L,nzchar(OUTDIR))
dir.create(OUTDIR,recursive=TRUE,showWarnings=FALSE)
dir.create(file.path(OUTDIR,'checkpoints'),recursive=TRUE,showWarnings=FALSE)
methods <- c('AdaptiveMetropolis','FLOCK_Novelty','FLOCK_Orthogonal')
scenarios <- c('ProductRidge','WeakSignalFunnel','ThreeComponentMixture')
version <- 'v1_orthogonal_three_targets_2026_10_25'
manifest <- list(version=version,reps=REPS,budget=BUDGET,seed=SEED,
                 methods=methods,scenarios=scenarios)
manifest_path <- file.path(OUTDIR,'run_manifest.rds')
if(file.exists(manifest_path)) {
  if(!identical(readRDS(manifest_path),manifest))
    stop('Incompatible output manifest; use a different OUTDIR.')
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
  Lloc<-chol(diag(.65^2,d))
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

# ---- Sampler A: eight independent Adaptive Metropolis chains ---------------
# Actual AM kernel: re-estimate covariance during warmup, then freeze.
walk_am <- function(f,target,x,lp,n,burn_frac=.35) {
  d<-target$d
  h<-matrix(NA_real_,n+1L,d);h[1L,]<-x
  L<-chol(diag(1.25^2,d));acc<-0L
  adapt_until<-floor(n*burn_frac)
  for(t in seq_len(n)) {
    if(t>max(50L,10L*d) && t<=adapt_until && t%%75L==0L) {
      S<-cov(h[seq_len(t),,drop=FALSE])
      L<-chol(spd((2.38^2/d)*S+.035*diag(d)))
    }
    y<-draw_norm(x,L);ly<-f$eval(y)
    if(log(runif(1))<ly-lp) {x<-y;lp<-ly;acc<-acc+1L}
    h[t+1L,]<-x
  }
  burn<-floor((n+1L)*burn_frac)
  list(draws=h[(burn+1L):(n+1L),,drop=FALSE],
       trace=h,accept=acc/n)
}
# Use a frozen proposal after adapting only in the warmup window.
am_sampler <- function(target,B,starting_point) {
  f<-make_counter(target);d<-target$d;C<-8L
  pieces<-rep(B%/%C,C);pieces[1L]<-pieces[1L]+B-sum(pieces)
  kept<-vector('list',C);traces<-vector('list',C)
  a<-numeric(C)
  for(j in seq_len(C)) {
    x<-if(j==1L) starting_point else target$prior_draw()
    lp<-f$eval(x)
    z<-walk_am(f,target,x,lp,pieces[j]-1L,burn_frac=.35)
    kept[[j]]<-z$draws;traces[[j]]<-z$trace;a[j]<-z$accept
  }
  stopifnot(f$get()==B)
  list(draws=do.call(rbind,kept),evals=f$get(),accept=mean(a),
       first_trace=traces[[1]],phase=list(chains=traces),
       retained_by_chain=kept)
}

# ---- Samplers B and C: staged global proposals ------------------------------
novelty_score <- function(candidates,history,scale) {
  if(nrow(history)>220L) history<-history[
    unique(as.integer(round(seq(1,nrow(history),length.out=220L)))),,drop=FALSE]
  vapply(seq_len(nrow(candidates)),function(i) {
    dz<-sweep(history,2,candidates[i,],'-')
    sqrt(min(rowSums(sweep(dz,2,scale,'/')^2)))
  },numeric(1))
}
flock_sampler <- function(target,B,starting_point,kind=c('Novelty','Orthogonal')) {
  kind<-match.arg(kind)
  f<-make_counter(target);d<-target$d
  scout_budget<-floor(.20*B)
  mid_budget<-floor(.40*B)
  histories<-list();scout_hist<-list();scout_logp<-list()
  ns<-if(kind=='Novelty') 2L else 1L
  sb<-rep(scout_budget%/%ns,ns)
  sb[1L]<-sb[1L]+scout_budget-sum(sb)
  best_lp<--Inf;best_x<-starting_point
  best_scout_lp<--Inf;best_scout_x<-starting_point
  for(j in seq_len(ns)) {
    x<-if(j==1L) starting_point else target$prior_draw()
    lp<-f$eval(x)
    z<-walk(f,target,x,lp,sb[j]-1L,
            chol(diag((if(kind=='Novelty')2.2 else 1.4)^2,d)),
            prior_jump=if(kind=='Novelty') .12 else .08,
            adapt=(kind=='Orthogonal'))
    scout_hist[[j]]<-z$h;scout_logp[[j]]<-z$lps
    histories[[length(histories)+1L]]<-z$h
    k<-which.max(z$lps)
    if(z$lps[k]>best_scout_lp) {
      best_scout_lp<-z$lps[k];best_scout_x<-z$h[k,]
    }
    if(z$lps[k]>best_lp) {best_lp<-z$lps[k];best_x<-z$h[k,]}
  }
  stopifnot(f$get()==scout_budget)
  scouting<-do.call(rbind,scout_hist)
  scan_scale<-pmax(apply(scouting,2,sd),.8)
  pooled_cov<-spd(if(nrow(scouting)>d+3) cov(scouting) else diag(2,d))
  if(kind=='Orthogonal') {
    eig<-eigen(pooled_cov,symmetric=TRUE)
    V<-eig$vectors
    eigval<-pmax(eig$values,.035)
    # Align follower movements to large-variance scout directions.
    Sf<-spd((2.38^2/d)*pooled_cov+.08*diag(d))
    # Invert but regularize and cap eigenvalue ratios for pioneers.
    inv<-1/(eigval+.20)
    inv<-pmin(inv,5*median(inv))
    inv<-inv/mean(inv)
    scale_mean<-median(sqrt(eigval))
    Sp<-spd(V%*%diag(inv*max(.55,scale_mean)^2,d)%*%t(V)+.04*diag(d))
  }
  alloc<-rep(mid_budget%/%6L,6L)
  alloc[1L]<-alloc[1L]+mid_budget-sum(alloc)
  mid_traces<-list()
  for(j in 1:6) {
    n<-alloc[j]
    if(kind=='Novelty') {
      if(j<=3L) {
        # Local followers start from discovered scout states.
        jj<-(j-1L)%%ns+1L
        ix<-which.max(scout_logp[[jj]])
        x<-scout_hist[[jj]][ix,]
        lp<-scout_logp[[jj]][ix]
        U<-chol(diag(.95^2,d));pj<-0
      } else {
        nc<-min(14L,max(4L,n%/%30L))
        cand<-t(replicate(nc,target$prior_draw()))
        lpc<-vapply(seq_len(nc),function(i)f$eval(cand[i,]),numeric(1))
        nov<-novelty_score(cand,scouting,scan_scale)
        score<-nov-.15*pmax(0,max(lpc)-lpc)
        ix<-which.max(score)
        x<-cand[ix,];lp<-lpc[ix];n<-n-nc
        U<-chol(diag(1.9^2,d));pj<-0
      }
    } else {
      if(j<=3L) {
        # Orthogonal FLOCK FOLLOWERS: scout-underexplored, plausible starts.
        nc<-min(14L,max(4L,n%/%30L))
        cand<-t(replicate(nc,target$prior_draw()))
        lpc<-vapply(seq_len(nc),function(i)f$eval(cand[i,]),numeric(1))
        nov<-novelty_score(cand,scouting,scan_scale)
        eligible<-which(lpc>=max(lpc)-10)
        score<-nov[eligible]
        ix<-eligible[which.max(score)]
        x<-cand[ix,];lp<-lpc[ix];n<-n-nc
        U<-chol(Sf);pj<-0
      } else {
        # Orthogonal FLOCK PIONEERS: start near a high-density scout point
        # and probe smaller-eigenvalue directions, both signs, with screening.
        nc<-min(14L,max(4L,n%/%30L))
        cand<-matrix(NA_real_,nc,d)
        for(i in seq_len(nc)) {
          axis<-d-((i-1L)%%min(3L,d))
          sign_dir<-if(i%%2L) 1 else -1
          magnitude<-c(.7,1.3,2.2,3.0)[((i-1L)%/%2L)%%4L+1L]
          cand[i,]<-best_scout_x+sign_dir*magnitude*
            max(.7,sqrt(median(eigval)))*V[,axis]
        }
        lpc<-vapply(seq_len(nc),function(i)f$eval(cand[i,]),numeric(1))
        valid<-which(lpc>=max(lpc)-9)
        # Prefer credible outward displacement; no posterior modes are known.
        radial<-sqrt(rowSums((cand-matrix(best_scout_x,nc,d,byrow=TRUE))^2))
        ix<-valid[which.max(radial[valid])]
        x<-cand[ix,];lp<-lpc[ix];n<-n-nc
        U<-chol(Sp);pj<-0
      }
    }
    stopifnot(n>=10L)
    z<-walk(f,target,x,lp,n,U,prior_jump=pj,adapt=FALSE)
    histories[[length(histories)+1L]]<-z$h
    mid_traces[[j]]<-z$h
    if(max(z$lps)>best_lp) {
      k<-which.max(z$lps);best_lp<-z$lps[k];best_x<-z$h[k,]
    }
  }
  stopifnot(f$get()==scout_budget+mid_budget)
  q<-build_proposal(target,histories)
  fin<-final_sampler(f,target,best_x,best_lp,B-f$get(),q)
  stopifnot(f$get()==B)
  list(draws=fin$draws,evals=f$get(),accept=fin$accept,
       first_trace=fin$draws,phase=list(scout=scout_hist,middle=mid_traces),
       proposal_K=q$K,
       scout_eigenvalues=if(kind=='Orthogonal') eig$values else numeric())
}

# ---- Tests and reference calculations --------------------------------------
# Validate that Gaussian proposal draws match the density's covariance system.
set.seed(9)
Ss<-matrix(c(1.5,.55,.55,.9),2,2)
U<-chol(Ss)
z<-replicate(12000L,draw_norm(c(0,0),U))
stopifnot(max(abs(cov(t(z))-Ss))<.09,
          abs(log_norm(c(0,0),c(0,0),U)-
              (-log(2*pi)-log(det(Ss))/2))<1e-9)
scenarios_obj<-list(ProductRidge=make_ridge(),
                    WeakSignalFunnel=make_funnel(),
                    ThreeComponentMixture=make_mixture())
for(scn in scenarios_obj) {
  stopifnot(scn$d>=2L,all(is.finite(scn$reference)),
            length(scn$metric_scales)==length(scn$reference),
            all(scn$metric_scales>0),
            abs(sum(scn$mode_probs)-1)<1e-6)
  x<-scn$prior_draw()
  stopifnot(length(x)==scn$d,is.finite(scn$logpost(x)))
}
# Invariance checks for the genuine symmetries of scenarios 1 and 3.
stopifnot(abs(scenarios_obj$ProductRidge$logpost(c(1,4)) -
              scenarios_obj$ProductRidge$logpost(c(-1,-4)))<1e-10)
mi<-scenarios_obj$ThreeComponentMixture
a<-c(-4,0,4,log(.2/.5),log(.3/.5))
b<-c(0,-4,4,log(.3/.5),log(.2/.5))
stopifnot(abs(mi$logpost(a)-mi$logpost(b))<1e-10)
ref_tab<-do.call(rbind,lapply(scenarios_obj,function(sc)
  data.frame(scenario=sc$name,metric=names(sc$reference),
             exact_value=as.numeric(sc$reference),
             scale=sc$metric_scales,
             stringsAsFactors=FALSE)))
write.csv(ref_tab,file.path(OUTDIR,'reference_values.csv'),row.names=FALSE)
ref_ci<-do.call(rbind,lapply(scenarios_obj,function(sc) {
  dx<-sc$ref_grid[2]-sc$ref_grid[1]
  pp<-sc$ref_pdf/sum(sc$ref_pdf*dx)
  cdf<-cumsum(pp*dx)
  data.frame(scenario=sc$name,coord=sc$show_coord,
             q025=sc$ref_grid[which(cdf>=.025)[1]],
             q975=sc$ref_grid[which(cdf>=.975)[1]])
}))
write.csv(ref_ci,file.path(OUTDIR,'reference_marginal_ci.csv'),row.names=FALSE)

# ---- Diagnostics, resumable benchmark --------------------------------------
make_result <- function(scn,obj,rep,method,elapsed) {
  X<-obj$draws
  stopifnot(ncol(X)==scn$d,nrow(X)>50L,all(is.finite(X)),
            obj$evals==BUDGET)
  est<-scn$metrics(X)
  truth<-scn$reference
  stopifnot(identical(names(est),names(truth)))
  err<-(est-truth)/scn$metric_scales
  coord<-X[,scn$show_coord]
  ci<-unname(quantile(coord,c(.025,.975)))
  mm<-as.integer(scn$modes(X))
  labels<-seq_along(scn$mode_probs)
  prop<-as.numeric(table(factor(mm,levels=labels)))/nrow(X)
  if(!is.null(obj$retained_by_chain)) {
    switches<-sum(vapply(obj$retained_by_chain,function(x) {
      r<-as.integer(scn$modes(x))
      sum(diff(r)!=0L)
    },numeric(1)))
  } else switches<-sum(diff(mm)!=0L)
  summary<-data.frame(
    scenario=scn$name,replicate=rep,method=method,
    evaluations=obj$evals,retained=nrow(X),
    seconds=elapsed,accept_rate=obj$accept,
    marginal_mean=mean(coord),marginal_q025=ci[1],marginal_q975=ci[2],
    standardized_rmse=sqrt(mean(err^2)),
    label_or_tail_l1=sum(abs(prop-scn$mode_probs)),
    transitions=switches,
    visited_groups=sum(prop>0),
    total_groups=length(labels),
    group_type=scn$mode_type,
    proposal_components=if(is.null(obj$proposal_K)) NA_integer_ else obj$proposal_K,
    stringsAsFactors=FALSE)
  estimates<-data.frame(scenario=scn$name,replicate=rep,method=method,
                        metric=names(est),exact=as.numeric(truth),
                        estimate=as.numeric(est),
                        error=as.numeric(est-truth),
                        standardized_error=as.numeric(err),
                        stringsAsFactors=FALSE)
  occupancy<-data.frame(scenario=scn$name,replicate=rep,method=method,
                        group=labels,estimated_prob=prop,
                        reference_prob=scn$mode_probs,
                        group_type=scn$mode_type)
  traces <- NULL
  if(rep==1L) {
    ids <- unique(as.integer(round(seq(1,nrow(X),
                      length.out=min(12000L,nrow(X))))))
    T0 <- obj$first_trace
    tids <- unique(as.integer(round(seq(1,nrow(T0),
                       length.out=min(12000L,nrow(T0))))))
    traces <- list(draws=X[ids,,drop=FALSE],trace=T0[tids,,drop=FALSE],
                   scout_eigen=obj$scout_eigenvalues)
  }
  list(summary=summary,estimates=estimates,occupancy=occupancy,example=traces)
}

atomic_save <- function(obj,path) {
  t<-tempfile(pattern='checkpoint_',tmpdir=dirname(path),fileext='.rds')
  saveRDS(obj,t)
  if(!file.rename(t,path)) stop('Cannot save checkpoint: ',path)
}
results<-list();iout<-0L
for(s in seq_along(scenarios_obj)) {
  scn<-scenarios_obj[[s]]
  cat('\n=== Scenario:',scn$name,'===\n',scn$description,'\n')
  cat('Reference:',paste(names(scn$reference),
                        signif(scn$reference,4),collapse='; '),'\n')
  for(rep in seq_len(REPS)) {
    # Same initial prior draw for the first chain of each method.
    set.seed(SEED+100000L*s+rep)
    x0<-scn$prior_draw()
    for(j in seq_along(methods)) {
      method<-methods[j]
      ck<-file.path(OUTDIR,'checkpoints',sprintf('%s_rep%04d_%s.rds',
                                                scn$name,rep,method))
      if(file.exists(ck)) {
        res<-readRDS(ck)
      } else {
        set.seed(SEED+100000L*s+100L*rep+j)
        timing<-system.time({
          output<-switch(method,
            AdaptiveMetropolis=am_sampler(scn,BUDGET,x0),
            FLOCK_Novelty=flock_sampler(scn,BUDGET,x0,'Novelty'),
            FLOCK_Orthogonal=flock_sampler(scn,BUDGET,x0,'Orthogonal'),
            stop('Unknown method'))
        })
        res<-make_result(scn,output,rep,method,as.numeric(timing['elapsed']))
        atomic_save(res,ck)
        rm(output)
      }
      iout<-iout+1L
      results[[iout]]<-res
      cat(sprintf('  rep %d/%d %-19s error=%.4f groupL1=%.4f (%d/%d evals)\n',
           rep,REPS,method,res$summary$standardized_rmse,
           res$summary$label_or_tail_l1,res$summary$evaluations,BUDGET))
    }
  }
}
all_results<-do.call(rbind,lapply(results,`[[`,'summary'))
all_estimates<-do.call(rbind,lapply(results,`[[`,'estimates'))
all_occupancy<-do.call(rbind,lapply(results,`[[`,'occupancy'))
write.csv(all_results,file.path(OUTDIR,'replicate_results.csv'),row.names=FALSE)
write.csv(all_estimates,file.path(OUTDIR,'metric_estimates.csv'),row.names=FALSE)
write.csv(all_occupancy,file.path(OUTDIR,'group_occupancy.csv'),row.names=FALSE)

key<-interaction(all_results$scenario,all_results$method,drop=TRUE)
summary_tbl<-do.call(rbind,lapply(split(all_results,key),function(z)
 data.frame(scenario=z$scenario[1],method=z$method[1],n=nrow(z),
            mean_error=mean(z$standardized_rmse),
            median_error=median(z$standardized_rmse),
            sd_error=if(nrow(z)>1L)sd(z$standardized_rmse) else NA_real_,
            mean_group_l1=mean(z$label_or_tail_l1),
            median_seconds=median(z$seconds),
            mean_transitions=mean(z$transitions),
            mean_visited=mean(z$visited_groups))))
summary_tbl<-summary_tbl[order(summary_tbl$scenario,summary_tbl$mean_error),]
rownames(summary_tbl)<-NULL
write.csv(summary_tbl,file.path(OUTDIR,'method_summary.csv'),row.names=FALSE)

paired_table<-list()
for(nm in scenarios) {
  data<-all_results[all_results$scenario==nm,]
  for(pair in list(c('FLOCK_Orthogonal','FLOCK_Novelty'),
                   c('FLOCK_Orthogonal','AdaptiveMetropolis'),
                   c('FLOCK_Novelty','AdaptiveMetropolis'))) {
    a<-data[data$method==pair[1],c('replicate','standardized_rmse')]
    b<-data[data$method==pair[2],c('replicate','standardized_rmse')]
    both<-merge(a,b,by='replicate',suffixes=c('_a','_b'))
    diff<-both$standardized_rmse_a-both$standardized_rmse_b
    if(length(diff)>=2L) {
      set.seed(SEED+length(diff)+nchar(nm)+nchar(pair[1]))
      boot<-replicate(1500L,mean(sample(diff,length(diff),replace=TRUE)))
      ci<-unname(quantile(boot,c(.025,.975)))
    } else ci<-c(NA_real_,NA_real_)
    paired_table[[length(paired_table)+1L]]<-data.frame(
      scenario=nm,method_a=pair[1],method_b=pair[2],n=length(diff),
      mean_error_difference=mean(diff),ci_low=ci[1],ci_high=ci[2],
      fraction_a_lower=mean(diff<0))
  }
}
paired<-do.call(rbind,paired_table)
write.csv(paired,file.path(OUTDIR,'paired_comparisons.csv'),row.names=FALSE)

# ---- Visualization: marginal posterior density and example traces ----------
pdf(file.path(OUTDIR,'posterior_diagnostics.pdf'),width=12,height=8)
for(nm in scenarios) {
  scn<-scenarios_obj[[nm]]
  saved<-lapply(methods,function(method) {
    p<-file.path(OUTDIR,'checkpoints',sprintf('%s_rep%04d_%s.rds',nm,1L,method))
    readRDS(p)$example
  })
  par(mfrow=c(2,3),mar=c(4,4,3,1))
  for(j in seq_along(methods)) {
    X<-saved[[j]]$draws[,scn$show_coord]
    # Density overlays are illustrative; functional errors remain primary.
    plot(scn$ref_grid,scn$ref_pdf,type='l',lwd=2,col='black',
         xlab=if(nm=='WeakSignalFunnel') 'log(tau)' else 'parameter 1',
         ylab='Density',main=paste(nm,methods[j],sep=': '),
         ylim=c(0,max(scn$ref_pdf)*1.3))
    if(length(unique(X))>3L) lines(density(X),col='steelblue',lwd=2)
    legend('topright',c('Reference','Estimated'),
           col=c('black','steelblue'),lty=1,bty='n',cex=.75)
  }
  # Trace of one retained chain; AM example is its first independent chain.
  for(j in seq_along(methods)) {
    xx<-saved[[j]]$trace[,scn$show_coord]
    plot(xx,type='l',col='grey30',main=paste(methods[j],'trace'),
         xlab='Iteration',ylab='Parameter value',cex.main=.85)
  }
}
# Errors and mode/tail occupancy across replicates.
par(mfrow=c(2,2),mar=c(7,4,3,1))
boxplot(standardized_rmse~interaction(scenario,method),data=all_results,
        las=2,main='Posterior functional error',ylab='Normalized RMSE')
boxplot(label_or_tail_l1~interaction(scenario,method),data=all_results,
        las=2,main='Branch/label/tail error',ylab='L1 error')
boxplot(seconds~method,data=all_results,las=2,main='Wall-clock cost',ylab='Seconds')
plot.new();text(.5,.7,'All methods: same target evaluation budget',cex=1.05)
text(.5,.5,paste('B =',BUDGET,'| R =',REPS),cex=1)
text(.5,.3,'FLOCK posterior estimates: final chain ONLY',cex=.9)
dev.off()

cat('\n=== Aggregate summary ===\n')
print(summary_tbl,row.names=FALSE,digits=4)
cat('\nPaired differences: negative means method_a has lower error.\n')
print(paired,row.names=FALSE,digits=4)
cat('\nSaved output to: ',normalizePath(OUTDIR),'\n',sep='')
cat('Notes: AM draws are pooled from 8 independent post-warmup chains; this\n')
cat('does not correct nonmixing. FLOCK draws come from one final corrected MH.\n')