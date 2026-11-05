# FLOCK vs strong MCMC competitors: multimodal, unidentified, weak-data targets
# Version v1.0 — base R, optional BayesianTools::DREAMzs / DEzs.
# Three original posterior families + an optional scientifically consequential case.
# No fitted model is misspecified here. All methods target the SAME posterior.
# Only final post-burn-in FLOCK draws are used for posterior inference.
# Run: Sys.setenv(RUN_MODE="smoke"); source("flock_hard_multimodal_benchmark.R")
# Then: Sys.setenv(RUN_MODE="pilot",HARD="TRUE"); source(...)

RUN_MODE <- Sys.getenv("RUN_MODE", "smoke")
base_settings <- switch(RUN_MODE,
  smoke=c(2L,3500L), pilot=c(10L,20000L), full=c(40L,50000L),
  stop("RUN_MODE must be smoke / pilot / full"))
REPS <- as.integer(Sys.getenv("REPS", base_settings[1]))
BUDGET <- as.integer(Sys.getenv("BUDGET", base_settings[2]))
HARD <- tolower(Sys.getenv("HARD","TRUE")) %in% c("true","1","yes")
SEED <- as.integer(Sys.getenv("SEED","20261025"))
INCLUDE_DREAM <- tolower(Sys.getenv("INCLUDE_DREAM","auto"))
HAS_BT <- requireNamespace("BayesianTools",quietly=TRUE)
if (INCLUDE_DREAM %in% c("true","1","yes") && !HAS_BT)
  stop("INCLUDE_DREAM=TRUE requires install.packages('BayesianTools')")
WITH_BT <- INCLUDE_DREAM %in% c("true","1","yes") ||
  (INCLUDE_DREAM=="auto" && HAS_BT)
all_scenarios <- c("ProductRidge","WeakSignalFunnel","ThreeComponentMixture",
                   "ConsequentialModes")
scn_option <- Sys.getenv("SCENARIOS","ProductRidge,WeakSignalFunnel,ThreeComponentMixture")
SCENARIOS <- if (tolower(scn_option)=="all") all_scenarios else
  trimws(strsplit(scn_option,",",fixed=TRUE)[[1]])
if(any(!SCENARIOS %in% all_scenarios) || anyDuplicated(SCENARIOS))
  stop("Invalid SCENARIOS: ",scn_option)
METHODS <- c("AdaptiveMH","FLOCK_Novelty","FLOCK_OrthogonalHard",
             "FLOCK_OrthogonalSoft","UndirectedMixtureMH","ParallelTempering")
if(WITH_BT) METHODS <- c(METHODS,"DREAMzs_Published","DEzs_Published")
OUTDIR <- Sys.getenv("OUTDIR",sprintf("flock_hard_%s_B%d_%s",RUN_MODE,BUDGET,
                                            if(HARD) "hard" else "original"))
stopifnot(REPS>=1L,BUDGET>=2500L,SEED>=0L,nzchar(OUTDIR))
dir.create(OUTDIR,showWarnings=FALSE,recursive=TRUE)
dir.create(file.path(OUTDIR,"checkpoints"),showWarnings=FALSE,recursive=TRUE)
manifest <- list(version="1.0",scenarios=SCENARIOS,methods=METHODS,
                 budget=BUDGET,reps=REPS,seed=SEED,hard=HARD)
manifest_path <- file.path(OUTDIR,"manifest.rds")
if(file.exists(manifest_path)) {
  if(!identical(readRDS(manifest_path),manifest))
    stop("Existing output manifest differs. Use a new OUTDIR.")
} else saveRDS(manifest,manifest_path)

logsum <- function(x) {
  a<-max(x)
  if(!is.finite(a)) return(-Inf)
  a+log(sum(exp(x-a)))
}
spd <- function(S,floor=.04) {
  S<-as.matrix(S);S<-(S+t(S))/2
  e<-eigen(S,symmetric=TRUE)
  e$vectors%*%diag(pmax(e$values,floor),length(e$values))%*%t(e$vectors)
}
# chol() is upper U, crossprod(U) = Sigma.
# Correct draw and density are consistent: x = mu + t(U) z.
draw_norm <- function(mu,U) as.numeric(mu+crossprod(U,rnorm(length(mu))))
log_norm <- function(x,mu,U) {
  z<-forwardsolve(t(U),x-mu)
  -length(mu)*log(2*pi)/2-sum(log(diag(U)))-sum(z*z)/2
}
make_counter <- function(target) {
  e<-new.env(parent=emptyenv());e$n<-0L
  e$eval<-function(x) {
    e$n<-e$n+1L
    as.numeric(target$logpost(as.numeric(x)))
  }
  e$get<-function() e$n
  e
}
make_ridge <- function(hard=TRUE) {
  set.seed(101)
  sd_y<-if(hard) .45 else .7
  y<-rnorm(if(hard) 8L else 6L,mean=if(hard) 5.8 else 3.7,sd=sd_y)
  n<-length(y);yb<-mean(y);sp<-if(hard) 4 else 3
  lp<-function(x) {
    if(length(x)!=2L) stop('Incorrect dimension')
    sum(dnorm(y,mean=x[1]*x[2],sd=sd_y,log=TRUE))+
      sum(dnorm(x,0,sp,log=TRUE))
  }
  prior<-function(x) sum(dnorm(x,0,sp,log=TRUE))
  # Exact marginal theta1 is Gaussian prior times integrated normal likelihood.
  grid<-seq(if(hard) -18 else -13,if(hard) 18 else 13,length.out=3601L)
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
make_funnel <- function(hard=TRUE) {
  y<-if(hard) c(-.8,.25,1.1,-.55,.65,-1.05,.48,-.39) else c(-.8,.25,1.1,-.55,.65,-1.05)
  sd_y<-if(hard) 3.2 else 2.3
  mean_t<-if(hard) -.85 else -.6
  sd_t<-if(hard) 1.35 else 1
  n<-length(y);d<-n+1L
  lprior<-function(x) {
    t<-x[1L];tau<-exp(t)
    if(!is.finite(tau) || tau<=0) return(-Inf)
    dnorm(t,mean_t,sd_t,log=TRUE)+sum(dnorm(x[-1L],0,tau,log=TRUE))
  }
  lp<-function(x) lprior(x)+sum(dnorm(y,x[-1L],sd_y,log=TRUE))
  prior_draw<-function() {
    t<-rnorm(1L,mean_t,sd_t)
    c(t,rnorm(n,0,exp(t)))
  }
  grid<-seq(if(hard) -9 else -6.5,if(hard) 6 else 4.5,length.out=4201L)
  logw<-dnorm(grid,mean_t,sd_t,log=TRUE)+
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
make_mixture <- function(hard=TRUE) {
  y<-if(hard) c(-8.4,-7.9,-7.6,-.5,.10,.47,7.6,8.05,8.39) else
    c(-4.15,-3.75,-4.4,-.30,.15,.38,3.7,4.1,4.35)
  n<-length(y);sd_y<-1;sp<-if(hard) 8 else 4.5
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
  gx<-seq(if(hard) -15 else -9,if(hard) 15 else 9,length.out=601L)
  pgrid<-vapply(gx,function(v)
    sum(weights*dnorm(v,mu_post[,1],sqrt(v_post[,1]))),numeric(1))
  list(name='ThreeComponentMixture',d=5L,logpost=lp,
       logprior=lprior,prior_draw=prior_draw,
       reference=ref,metric_scales=c(2,8,8,1/6),metrics=metrics,
       modes=label_fn,mode_probs=rep(1/6,6L),mode_type='six_label_permutations',
       show_coord=1L,ref_grid=gx,ref_pdf=pgrid,
       description='3-component Gaussian mixture with exchangeable priors: 6 labelings')
}

# ---- Optional fourth target: same training fit, different future conclusions
make_consequential <- function(hard=TRUE) {
  sd_y<-if(hard) .6 else 1
  y<-c(9.15,8.8,9.25,8.65,9.3,8.9)
  z<-c(.10,-.20,.3)
  sd_z<-1.2
  mu0<-1.15;prior_sd<-3
  lp0 <- function(x) sum(dnorm(x,c(mu0,0),c(prior_sd,2),log=TRUE))
  ll <- function(x) sum(dnorm(y,x[1]^2,sd_y,log=TRUE))+
                       sum(dnorm(z,x[2],sd_z,log=TRUE))
  lpost <- function(x) lp0(x)+ll(x)
  grid<-seq(-9,9,length.out=4501L)
  lw<-dnorm(grid,mu0,prior_sd,log=TRUE)+
    vapply(grid,function(v)sum(dnorm(y,v^2,sd_y,log=TRUE)),numeric(1))
  wg<-exp(lw-max(lw));wg<-wg/sum(wg)
  v2<-1/(1/4+length(z)/sd_z^2)
  m2<-v2*sum(z)/sd_z^2
  fut <- (grid-1)^2+m2
  ref<-c(p_theta_positive=sum(wg[grid>0]),
         future_mean=sum(wg*fut),
         p_future_gt_12=sum(wg*pnorm((grid-1)^2+m2-12,sd=sqrt(v2+1))),
         mean_theta=sum(wg*grid))
  metrics <- function(X) c(p_theta_positive=mean(X[,1]>0),
        future_mean=mean((X[,1]-1)^2+X[,2]),
        p_future_gt_12=mean(pnorm((X[,1]-1)^2+X[,2]-12)),
        mean_theta=mean(X[,1]))
  list(name="ConsequentialModes",d=2L,logpost=lpost,logprior=lp0,
       prior_draw=function() c(rnorm(1,mu0,prior_sd),rnorm(1,0,2)),
       reference=ref,metric_scales=c(.5,12,.5,3),metrics=metrics,
       modes=function(X)ifelse(X[,1]>0,2L,1L),
       mode_probs=c(1-ref[1],ref[1]),mode_type="consequential_prediction",
       show_coord=1L,ref_grid=grid,ref_pdf=wg/(grid[2]-grid[1]),
       description="Same training fit under theta and -theta; asymmetric prior; differing future predictions")
}

# ---- Core: local MH, optional fixed Gaussian adaptation during warmup ----
walk <- function(f,target,x,lp,n,U,prior_jump=.05,adapt=FALSE,
                 adapt_until=n,update_every=80L) {
  d<-target$d
  h<-matrix(NA_real_,n+1L,d);h[1,]<-x
  lps<-numeric(n+1L);lps[1]<-lp
  accepted<-0L;prior_moves<-0L
  for(t in seq_len(n)) {
    if(adapt && t>=max(100L,12L*d) && t<=adapt_until && t%%update_every==0L) {
      S<-cov(h[seq_len(t),,drop=FALSE])
      U<-chol(spd((2.38^2/d)*S+.035*diag(d)))
    }
    isprior<-runif(1)<prior_jump
    if(isprior) {
      y<-target$prior_draw()
      ly<-f$eval(y)
      loga<-ly-lp+target$logprior(x)-target$logprior(y)
      prior_moves<-prior_moves+1L
    } else {
      y<-draw_norm(x,U)
      ly<-f$eval(y)
      loga<-ly-lp
    }
    if(is.finite(loga) && log(runif(1))<loga) {
      x<-y;lp<-ly;accepted<-accepted+1L
    }
    h[t+1L,]<-x;lps[t+1L]<-lp
  }
  list(x=x,lp=lp,h=h,lps=lps,U=U,accept=accepted/max(n,1L),
       prior_attempts=prior_moves)
}
# Eight independent Adaptive Metropolis chains; no pooling prior to burn-in.
independent_am <- function(target,B,x0) {
  f<-make_counter(target);d<-target$d
  allocations<-rep(B%/%8L,8L);allocations[1]<-allocations[1]+B-sum(allocations)
  kept<-vector("list",8L);a<-numeric(8L)
  for(j in 1:8) {
    x<-if(j==1) x0 else target$prior_draw()
    lp<-f$eval(x)
    n<-allocations[j]-1L;burn<-floor(.35*n)
    z<-walk(f,target,x,lp,n,chol(diag((if(d>10) 2/sqrt(d) else 1.3)^2,d)),
            prior_jump=0,adapt=TRUE,adapt_until=burn)
    kept[[j]]<-z$h[(burn+2L):(n+1L),,drop=FALSE]
    a[j]<-z$accept
  }
  stopifnot(f$get()==B)
  list(draws=do.call(rbind,kept),chains=kept,evals=f$get(),accept=mean(a),
       stage_coverage=NA_integer_,global_accept=NA_real_,local_accept=mean(a))
}
# Pool historical states equally by exploratory chain, NOT by residence time.
thin_pool <- function(histories,limit=150L) {
  do.call(rbind,lapply(histories,function(H) {
    k<-unique(as.integer(round(seq(1,nrow(H),length.out=min(limit,nrow(H))))))
    H[k,,drop=FALSE]
  }))
}
novelty_score <- function(C,H,scale) {
  H<-as.matrix(H);C<-as.matrix(C)
  if(nrow(H)>230L) {
    ix<-unique(as.integer(round(seq(1,nrow(H),length.out=230L))))
    H<-H[ix,,drop=FALSE]
  }
  vapply(seq_len(nrow(C)),function(j) {
    d<-sweep(H,2,C[j,],"-")
    sqrt(min(rowSums(sweep(d,2,scale,"/")^2)))
  },numeric(1))
}
# A globally supported proposal. The prior component is evaluated exactly.
fit_global_q <- function(target,histories) {
  d<-target$d
  pool<-unique(thin_pool(histories,limit=140L))
  if(nrow(pool)<4L) stop("Insufficient exploratory states")
  K<-min(8L,max(1L,floor(nrow(pool)/60L)))
  km<-tryCatch(kmeans(pool,centers=K,nstart=3L,iter.max=45L),
               error=function(e) NULL)
  if(is.null(km)) {
    K<-1L; cl<-rep(1L,nrow(pool))
    centers<-matrix(colMeans(pool),1L,d)
  } else {
    cl<-km$cluster
    centers<-matrix(km$centers,ncol=d)
  }
  UU<-lapply(seq_len(K),function(k) {
    z<-pool[cl==k,,drop=FALSE]
    S<-if(nrow(z)>=d+3L) cov(z) else diag(1.5,d)
    chol(spd(1.3*S+.1*diag(d)))
  })
  mixprob<-.82/K
  draw<-function() {
    if(runif(1)<.18) return(target$prior_draw())
    k<-sample.int(K,1L)
    draw_norm(centers[k,],UU[[k]])
  }
  logq<-function(x) {
    z<-vapply(seq_len(K),function(k)
      log(mixprob)+log_norm(x,centers[k,],UU[[k]]),numeric(1))
    logsum(c(z,log(.18)+target$logprior(x)))
  }
  list(draw=draw,logq=logq,K=K)
}
# Frozen mixture of valid local RWMH and independence MH transition kernels.
final_mh <- function(f,target,x,lp,n,q) {
  d<-target$d
  h<-matrix(NA_real_,n+1L,d);h[1,]<-x
  locU<-chol(diag((if(d>10)2/sqrt(d) else .70)^2,d))
  lqx<-q$logq(x);accept_local<-accept_global<-0L
  attempt_local<-attempt_global<-0L
  for(t in seq_len(n)) {
    islocal<-runif(1)<.20
    if(islocal) {
      attempt_local<-attempt_local+1L
      y<-draw_norm(x,locU);ly<-f$eval(y)
      lqy<-q$logq(y)
      loga<-ly-lp
    } else {
      attempt_global<-attempt_global+1L
      y<-q$draw();ly<-f$eval(y);lqy<-q$logq(y)
      loga<-ly-lp+lqx-lqy
    }
    if(is.finite(loga) && log(runif(1))<loga) {
      x<-y;lp<-ly;lqx<-lqy
      if(islocal) accept_local<-accept_local+1L else accept_global<-accept_global+1L
    }
    h[t+1L,]<-x
  }
  burn<-floor(.10*n)
  list(draws=h[(burn+2L):(n+1L),,drop=FALSE],
       local_accept=accept_local/max(1L,attempt_local),
       global_accept=accept_global/max(1L,attempt_global),
       accept=(accept_local+accept_global)/n)
}
soft_weights <- function(lambda,alpha=.5,rho=.7,delta=.2) {
  raw<-(pmax(lambda,0)+delta)^(-alpha)
  (1-rho)/length(lambda)+rho*raw/sum(raw)
}
# Four variants share the exact same final transition kernels and budget split.
staged_sampler <- function(target,B,x0,kind=c("Novelty","OrthogonalHard",
                                           "OrthogonalSoft","Undirected")) {
  kind<-match.arg(kind);d<-target$d
  f<-make_counter(target)
  scout_B<-floor(.20*B);mid_B<-floor(.40*B)
  nscout<-if(kind%in%c("OrthogonalHard","OrthogonalSoft"))1L else 2L
  sb<-rep(scout_B%/%nscout,nscout)
  sb[1]<-sb[1]+scout_B-sum(sb)
  scoutH<-scoutLP<-vector("list",nscout);histories<-list()
  best_lp<--Inf;best_x<-x0
  for(j in seq_len(nscout)) {
    x<-if(j==1L) x0 else target$prior_draw()
    lp<-f$eval(x)
    s<-if(d>10)3/sqrt(d) else if(kind=="Novelty")2.2 else 1.5
    z<-walk(f,target,x,lp,sb[j]-1L,chol(diag(s*s,d)),
            prior_jump=.10,adapt=kind%in%c("OrthogonalHard","OrthogonalSoft"),
            adapt_until=floor((sb[j]-1L)*.8))
    scoutH[[j]]<-z$h;scoutLP[[j]]<-z$lps
    histories[[length(histories)+1L]]<-z$h
    k<-which.max(z$lps)
    if(z$lps[k]>best_lp) {best_lp<-z$lps[k];best_x<-z$h[k,]}
  }
  stopifnot(f$get()==scout_B)
  scoutpool<-do.call(rbind,scoutH)
  scale<-pmax(apply(scoutpool,2,sd),.8)
  if(kind%in%c("OrthogonalHard","OrthogonalSoft")) {
    # Local geometry close to scout's best point; not necessarily global posterior.
    delta<-sweep(scoutpool,2,best_x,"-")
    scores<-rowSums(sweep(delta,2,scale,"/")^2)
    ids<-order(scores)[seq_len(min(nrow(scoutpool),max(d+12L,
                                                floor(.45*nrow(scoutpool)))))]
    S<-spd(if(length(ids)>=d+3L)cov(scoutpool[ids,,drop=FALSE]) else diag(1,d))
    eig<-eigen(S,symmetric=TRUE)
    probs<-soft_weights(eig$values)
    Uf<-chol(spd((2.38^2/d)*S+.10*diag(d)))
  }
  alloc<-rep(mid_B%/%6L,6L);alloc[1]<-alloc[1]+mid_B-sum(alloc)
  for(j in 1:6) {
    n<-alloc[j]
    if(kind=="Undirected") {
      x<-target$prior_draw();lp<-f$eval(x);n<-n-1L
      U<-chol(diag((if(d>10)2/sqrt(d) else 1.1)^2,d))
    } else if(kind=="Novelty" && j<=3) {
      jj<-(j-1L)%%nscout+1L;k<-which.max(scoutLP[[jj]])
      x<-scoutH[[jj]][k,];lp<-scoutLP[[jj]][k]
      U<-chol(diag((if(d>10)2/sqrt(d) else .95)^2,d))
    } else {
      nc<-min(14L,max(4L,floor(n/35L)))
      cand<-matrix(NA_real_,nc,d)
      if(kind%in%c("OrthogonalHard","OrthogonalSoft") && j>3L) {
        for(i in seq_len(nc)) {
          axis<-if(kind=="OrthogonalHard") {
            d-((i-1L)%%min(3L,d))
          } else sample.int(d,1L,prob=probs)
          sgn<-if(i%%2L==0L)1 else -1
          r<-c(.8,1.5,2.5,3.5)[((i-1L)%%4L)+1L]
          cand[i,]<-best_x+sgn*r*max(.7,sqrt(median(eig$values)))*eig$vectors[,axis]
          if(kind=="OrthogonalSoft")cand[i,]<-cand[i,]+rnorm(d,sd=.10)
        }
      } else {
        for(i in seq_len(nc))cand[i,]<-target$prior_draw()
      }
      lp_cand<-vapply(seq_len(nc),function(i)f$eval(cand[i,]),numeric(1))
      nv<-novelty_score(cand,scoutpool,scale)
      eligible<-which(lp_cand>=max(lp_cand)-10)
      if(kind=="Novelty"){
        score<-nv-.15*pmax(0,max(lp_cand)-lp_cand)
        chosen<-which.max(score)
      } else if(kind=="OrthogonalHard") {
        chosen<-eligible[which.max(nv[eligible])]
      } else if(kind=="OrthogonalSoft") {
        score<-.7*nv[eligible]-.3*(max(lp_cand)-lp_cand[eligible])
        score<-score-max(score)
        chosen<-eligible[sample.int(length(eligible),1L,prob=exp(score))]
      } else {
        chosen<-eligible[sample.int(length(eligible),1L)]
      }
      x<-cand[chosen,];lp<-lp_cand[chosen];n<-n-nc
      if(kind%in%c("OrthogonalHard","OrthogonalSoft")&&j<=3L) U<-Uf else
        U<-chol(diag((if(d>10)2.4/sqrt(d) else 1.7)^2,d))
    }
    z<-walk(f,target,x,lp,n,U,prior_jump=.025,adapt=FALSE)
    histories[[length(histories)+1L]]<-z$h
    k<-which.max(z$lps)
    if(z$lps[k]>best_lp) {best_lp<-z$lps[k];best_x<-z$h[k,]}
  }
  stopifnot(f$get()==scout_B+mid_B)
  q<-fit_global_q(target,histories)
  stage_labels<-target$modes(thin_pool(histories,limit=120L))
  final<-final_mh(f,target,best_x,best_lp,B-f$get(),q)
  stopifnot(f$get()==B)
  seen_stage<-unique(stage_labels)
  seen_final<-unique(target$modes(final$draws))
  list(draws=final$draws,evals=f$get(),accept=final$accept,
       local_accept=final$local_accept,global_accept=final$global_accept,
       stage_coverage=length(seen_stage),
       discovered_but_missing=sum(!seen_stage %in% seen_final),
       proposal_K=q$K)
}

# ---- Parallel tempering: likelihood tempering with valid adjacent swaps ----
parallel_tempering <- function(target,B,x0) {
  f<-make_counter(target);d<-target$d;C<-8L
  # Spaced inverse temperatures: 1 (cold) to 0.08 (hot).
  beta<-exp(seq(log(1),log(.08),length.out=C))
  starts<-matrix(NA_real_,C,d); lp0<-ll<-numeric(C)
  for(j in seq_len(C)) {
    starts[j,]<-if(j==1L)x0 else target$prior_draw()
    lfull<-f$eval(starts[j,])
    lp0[j]<-target$logprior(starts[j,])
    ll[j]<-lfull-lp0[j]
  }
  rounds<-(B-C)%/%C
  leftover<-B-C-rounds*C
  stopifnot(rounds>0L)
  cold<-matrix(NA_real_,rounds+leftover+1L,d);cold[1,]<-starts[1,]
  swaps<-attempt_swaps<-0L
  cold_accept<-0L
  for(t in seq_len(rounds)) {
    for(j in seq_len(C)) {
      # Proposals scaled to the likelihood temperature; symmetric RWMH.
      sd0<-if(d>10L) 1.5/sqrt(d) else .75
      y<-starts[j,]+rnorm(d,sd=sd0/sqrt(beta[j]))
      full_new<-f$eval(y)
      pr_new<-target$logprior(y)
      ll_new<-full_new-pr_new
      loga<-pr_new-lp0[j]+beta[j]*(ll_new-ll[j])
      if(is.finite(loga) && log(runif(1))<loga) {
        starts[j,]<-y;lp0[j]<-pr_new;ll[j]<-ll_new
        if(j==1L)cold_accept<-cold_accept+1L
      }
    }
    # Alternating neighbor pairs on every round; detailed balance for swaps.
    pairs<-if(t%%2L) seq(1L,C-1L,by=2L) else seq(2L,C-1L,by=2L)
    for(j in pairs) {
      attempt_swaps<-attempt_swaps+1L
      loga<-(beta[j]-beta[j+1L])*(ll[j+1L]-ll[j])
      if(log(runif(1))<loga) {
        tmp<-starts[j,];starts[j,]<-starts[j+1L,];starts[j+1L,]<-tmp
        tmp<-ll[j];ll[j]<-ll[j+1L];ll[j+1L]<-tmp
        tmp<-lp0[j];lp0[j]<-lp0[j+1L];lp0[j+1L]<-tmp
        swaps<-swaps+1L
      }
    }
    cold[t+1L,]<-starts[1,]
  }
  if(leftover>0L) for(k in seq_len(leftover)) {
    y<-starts[1,]+rnorm(d,sd=if(d>10L)1.5/sqrt(d) else .75)
    full_new<-f$eval(y);pr_new<-target$logprior(y)
    ll_new<-full_new-pr_new
    loga<-pr_new-lp0[1]+ll_new-ll[1]
    if(is.finite(loga)&&log(runif(1))<loga) {
      starts[1,]<-y;lp0[1]<-pr_new;ll[1]<-ll_new
      cold_accept<-cold_accept+1L
    }
    cold[rounds+k+1L,]<-starts[1,]
  }
  stopifnot(f$get()==B)
  burn<-floor(.35*nrow(cold))
  list(draws=cold[(burn+1L):nrow(cold),,drop=FALSE],
       evals=B,accept=cold_accept/(rounds+leftover),
       local_accept=cold_accept/(rounds+leftover),
       global_accept=swaps/max(attempt_swaps,1L),stage_coverage=NA_integer_)
}

# ---- Published BayesianTools implementations, optional dependency ----------
# These are the actual CRAN implementations, not hand-built substitutes.
# Their target calls can exceed or fall below B due to initialization/internal
# operations. We REPORT the observed number and flag unmatched budgets.
bayesian_tools_method <- function(target,B,which=c("DREAMzs","DEzs")) {
  which<-match.arg(which)
  if(!requireNamespace("BayesianTools",quietly=TRUE))stop("Install BayesianTools")
  f<-make_counter(target)
  sampler_prior<-function(n=1L) {
    n<-as.integer(n)
    if(n<=1L) return(target$prior_draw())
    t(replicate(n,target$prior_draw()))
  }
  prior<-BayesianTools::createPrior(
    density=function(x) target$logprior(x),sampler=sampler_prior)
  setup<-BayesianTools::createBayesianSetup(
    likelihood=function(x) f$eval(x)-target$logprior(x),
    prior=prior,names=paste0("theta",seq_len(target$d)))
  itr<-as.integer(B-10L)
  fit<-BayesianTools::runMCMC(
    bayesianSetup=setup,sampler=which,
    settings=list(iterations=itr,burnin=as.integer(floor(.35*itr)),
                  message=FALSE,consoleUpdates=1000000L))
  X<-as.matrix(BayesianTools::getSample(fit,start=1L,thin=1L))
  if(ncol(X)!=target$d) stop("BayesianTools extraction returned incorrect number of parameters")
  list(draws=X,evals=f$get(),accept=NA_real_,
       local_accept=NA_real_,global_accept=NA_real_,stage_coverage=NA_integer_)
}

# ---- Diagnostics (reference known but never used by algorithm) ------------
mode_stats <- function(target,Z) {
  lab<-target$modes(Z)
  K<-length(target$mode_probs)
  counts<-tabulate(lab,nbins=K)
  freq<-counts/sum(counts)
  c(mode_l1=sum(abs(freq-target$mode_probs)),
    modes_visited=sum(counts>0),
    mode_switches=sum(diff(lab)!=0))
}
measure <- function(target,obj,scenario,method,rep,elapsed) {
  Z<-as.matrix(obj$draws)
  if(ncol(Z)!=target$d||nrow(Z)<10L || any(!is.finite(Z)))
    stop("Invalid posterior draws for ",scenario," / ",method)
  est<-target$metrics(Z)
  ref<-target$reference
  if(!identical(names(est),names(ref)))
    stop("Metric-name mismatch ",scenario," / ",method)
  err<-(est-ref)/target$metric_scales
  m<-mode_stats(target,Z)
  data.frame(scenario=scenario,method=method,rep=rep,
             evals=obj$evals,budget=BUDGET,
             comparable_budget=abs(obj$evals-BUDGET)<=.05*BUDGET,
             retained=nrow(Z),elapsed_seconds=elapsed,
             normalized_rmse=sqrt(mean(err^2)),
             max_abs_scaled_error=max(abs(err)),
             mode_l1=as.numeric(m["mode_l1"]),
             modes_visited=as.integer(m["modes_visited"]),
             mode_switches=as.integer(m["mode_switches"]),
             mode_type=target$mode_type,
             stage_coverage=if(is.null(obj$stage_coverage))NA_integer_ else obj$stage_coverage,
             discovered_but_missing=if(is.null(obj$discovered_but_missing))NA_integer_
                                    else obj$discovered_but_missing,
             accept=if(is.null(obj$accept))NA_real_ else obj$accept,
             local_accept=if(is.null(obj$local_accept))NA_real_ else obj$local_accept,
             global_accept=if(is.null(obj$global_accept))NA_real_ else obj$global_accept,
             stringsAsFactors=FALSE)
}
checkpoint_path <- function(scn,method,r) {
  file.path(OUTDIR,"checkpoints",sprintf("%s__%s__%03d.rds",scn,method,r))
}
write_summaries <- function() {
  files<-list.files(file.path(OUTDIR,"checkpoints"),pattern="[.]rds$",full.names=TRUE)
  if(!length(files))return(invisible(NULL))
  objects<-lapply(files,function(z)tryCatch(readRDS(z),error=function(e)NULL))
  objects<-Filter(Negate(is.null),objects)
  if(!length(objects))return(invisible(NULL))
  res<-do.call(rbind,lapply(objects,function(z)z$row))
  res<-res[order(res$scenario,res$method,res$rep),,drop=FALSE]
  write.csv(res,file.path(OUTDIR,"replicate_results.csv"),row.names=FALSE)
  key<-split(res,paste(res$scenario,res$method,sep="__"))
  summary<-do.call(rbind,lapply(key,function(z) data.frame(
    scenario=z$scenario[1],method=z$method[1],n=nrow(z),
    mean_rmse=mean(z$normalized_rmse),median_rmse=median(z$normalized_rmse),
    mean_mode_l1=mean(z$mode_l1),
    failure_rate=mean(z$normalized_rmse>.25),
    mean_mode_visited=mean(z$modes_visited),
    median_elapsed=median(z$elapsed_seconds),
    mean_evals=mean(z$evals),budget_matched=all(z$comparable_budget))))
  rownames(summary)<-NULL
  write.csv(summary,file.path(OUTDIR,"method_summary.csv"),row.names=FALSE)
  invisible(res)
}
# Paired differences (FLOCK - comparator): negative favors FLOCK.
write_paired <- function(res) {
  flocks<-c("FLOCK_Novelty","FLOCK_OrthogonalHard","FLOCK_OrthogonalSoft")
  comps<-setdiff(METHODS,flocks)
  # Do not make paired fixed-budget claims for unmatched external methods.
  res<-res[res$comparable_budget,,drop=FALSE]
  rows<-list()
  for(scn in SCENARIOS)for(fm in flocks)for(cm in comps) {
    a<-res[res$scenario==scn & res$method==fm,c("rep","normalized_rmse")]
    b<-res[res$scenario==scn & res$method==cm,c("rep","normalized_rmse")]
    ab<-merge(a,b,by="rep",suffixes=c("_flock","_other"))
    if(nrow(ab)<2)next
    d<-ab$normalized_rmse_flock-ab$normalized_rmse_other
    # Descriptive paired percentile bootstrap; not evidence at tiny R.
    set.seed(100000+length(d)+nchar(fm)+nchar(cm)+nchar(scn))
    bs<-replicate(1500L,mean(sample(d,length(d),replace=TRUE)))
    rows[[length(rows)+1L]]<-data.frame(scenario=scn,flock=fm,competitor=cm,
      n=length(d),mean_difference=mean(d),
      lower=unname(quantile(bs,.025)),upper=unname(quantile(bs,.975)),
      fraction_flock_better=mean(d<0))
  }
  if(length(rows))write.csv(do.call(rbind,rows),
     file.path(OUTDIR,"paired_comparisons.csv"),row.names=FALSE)
}

# ---- Numerical self-tests before running the samplers ----------------------
set.seed(77)
S<-matrix(c(1.5,.55,.55,.9),2,2)
U<-chol(S)
xx<-replicate(14000L,draw_norm(c(0,0),U))
stopifnot(max(abs(cov(t(xx))-S))<.075)
stopifnot(abs(log_norm(c(0,0),c(0,0),U)-
              (-log(2*pi)-log(det(S))/2))<1e-9)
for(lambda in list(c(1.9,.1),c(8,.2,.04))) {
  pr<-soft_weights(lambda)
  stopifnot(all(pr>0),abs(sum(pr)-1)<1e-12)
}
cat("Passed proposal covariance and soft-orthogonal self-tests.\n")

# ---- Build and validate posterior references once -------------------------
construct_target <- function(s) {
  switch(s,ProductRidge=make_ridge(HARD),
         WeakSignalFunnel=make_funnel(HARD),
         ThreeComponentMixture=make_mixture(HARD),
         ConsequentialModes=make_consequential(HARD))
}
targets <- lapply(SCENARIOS,construct_target);names(targets)<-SCENARIOS
for(targ in targets) {
  stopifnot(all(is.finite(targ$reference)),
            abs(sum(targ$mode_probs)-1)<1e-8,
            length(targ$reference)==length(targ$metric_scales))
}
ref<-do.call(rbind,lapply(targets,function(targ)data.frame(
  scenario=targ$name,metric=names(targ$reference),
  reference=as.numeric(targ$reference),scale=targ$metric_scales)))
write.csv(ref,file.path(OUTDIR,"reference_values.csv"),row.names=FALSE)

run_one <- function(targ,method,B,x0) {
  switch(method,
    AdaptiveMH=independent_am(targ,B,x0),
    FLOCK_Novelty=staged_sampler(targ,B,x0,"Novelty"),
    FLOCK_OrthogonalHard=staged_sampler(targ,B,x0,"OrthogonalHard"),
    FLOCK_OrthogonalSoft=staged_sampler(targ,B,x0,"OrthogonalSoft"),
    UndirectedMixtureMH=staged_sampler(targ,B,x0,"Undirected"),
    ParallelTempering=parallel_tempering(targ,B,x0),
    DREAMzs_Published=bayesian_tools_method(targ,B,"DREAMzs"),
    DEzs_Published=bayesian_tools_method(targ,B,"DEzs"),
    stop("Unknown method: ",method))
}

cat(sprintf("FLOCK hard multimodal benchmark | %s | HARD=%s | B=%d | REPS=%d\n",
            RUN_MODE,HARD,BUDGET,REPS))
cat("Methods: ",paste(METHODS,collapse=", "),"\n",sep="")
cat("Scenarios: ",paste(SCENARIOS,collapse=", "),"\n",sep="")
if(!WITH_BT)cat("BayesianTools unavailable/disabled: published DREAMzs and DEzs skipped.\n")
for(scn in SCENARIOS) {
  target<-targets[[scn]]
  for(r in seq_len(REPS)) {
    # Keep initial x0 the same across methods within replicate.
    set.seed(SEED+as.integer(match(scn,SCENARIOS))*100000L+r)
    x0<-target$prior_draw()
    for(method in METHODS) {
      path<-checkpoint_path(scn,method,r)
      if(file.exists(path)) next
      # Same initial random-number state across methods (except as algorithms
      # consume different numbers of random draws after they begin).
      set.seed(SEED+as.integer(match(scn,SCENARIOS))*1000000L+r*100L)
      t0<-proc.time()[[3]]
      obj<-tryCatch(run_one(target,method,BUDGET,x0),error=function(e)e)
      elapsed<-proc.time()[[3]]-t0
      if(inherits(obj,"error")) {
        cat(sprintf("[%s %d/%d] %s FAILED: %s\n",
                    scn,r,REPS,method,conditionMessage(obj)))
        next
      }
      row<-measure(target,obj,scn,method,r,elapsed)
      # Save a single representative posterior for visual inspection.
      if(r==1L) {
        colnames(obj$draws)<-paste0("theta",seq_len(target$d))
        write.csv(obj$draws,file.path(OUTDIR,paste0(scn,"__",method,
                              "__retained_posterior.csv")),row.names=FALSE)
      }
      saveRDS(list(row=row),path)
      cat(sprintf("[%s %d/%d] %-23s RMSE %.4f | region L1 %.4f | %d calls | %.2fs\n",
                  scn,r,REPS,method,row$normalized_rmse,row$mode_l1,
                  row$evals,elapsed))
    }
    write_summaries()
  }
}
res<-write_summaries()
if(!is.null(res)) {
  write_paired(res)
  print(read.csv(file.path(OUTDIR,"method_summary.csv")),row.names=FALSE)
}
# Reference-overlay PDF for retained samples from first replicate.
pdf_path<-file.path(OUTDIR,"posterior_reference_overlays.pdf")
pdf(pdf_path,width=12,height=7)
for(scn in SCENARIOS) {
  target<-targets[[scn]]
  par(mfrow=c(2,4),mar=c(4,4,3,1))
  for(method in METHODS) {
    fn<-file.path(OUTDIR,paste0(scn,"__",method,"__retained_posterior.csv"))
    if(!file.exists(fn)) {plot.new();title(main=paste(scn,method));next}
    X<-as.matrix(read.csv(fn))
    display<-if(scn=="ConsequentialModes")X[,1] else X[,target$show_coord]
    dd<-density(display,n=1024L)
    plot(target$ref_grid,target$ref_pdf,type="l",lwd=2,col="black",
         xlab="Parameter",ylab="Density",main=method,
         ylim=c(0,max(dd$y,target$ref_pdf)*1.07))
    lines(dd,col="steelblue",lwd=2)
    legend("topright",legend=c("Reference","Retained"),
           col=c("black","steelblue"),lty=1,bty="n",cex=.65)
  }
  par(mfrow=c(1,1));mtext(scn,outer=TRUE,line=-1.5)
}
dev.off()
cat("Done. Results in ",normalizePath(OUTDIR),"\n",sep="")
cat("Key files: method_summary.csv, paired_comparisons.csv, reference_values.csv,\n")
cat("replicate_results.csv, posterior_reference_overlays.pdf, and checkpoints/.\n")