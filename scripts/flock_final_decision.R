# FLOCK: pre-registered final decision benchmark (experimental prototype)
# Standalone BASE R only. R version >= 4.0 recommended.
# See flock_final_decision_README.md before interpreting outcomes.
#
# Scenarios: KnownMixture5D, DiagnosticLCA, HighDimCorrelated (40 coefficients)
# Methods: IndependentAM, IndependentAM_Shifts, StagedAM,
#          FLOCK_AM, FLOCK_AM_LocalOnly
#
# All methods use exactly B calls to the identical unnormalized log posterior.
# All locally adaptive proposals use THE SAME Adaptive Metropolis transition code.
# Regional shifts are state-independent symmetric proposals with MH acceptance
# min(1, exp(logpi(y)-logpi(x))); the shift pairs are frozen before sampling.
# All exploratory/warmup samples excluded from reported posterior estimates.
# The final output of staged/FLOCK algorithms is one post-warmup chain only.
# The eight-chain baselines pool post-warmup draws, and record per-chain diagnostics.
# This is NOT a claim that one final chain proves convergence.
#
# RUN_MODE=smoke (2 reps x 6000), pilot (10 x 30000), full (60 x 60000).
# Optional REPS, BUDGET, SCENARIOS, METHODS, OUTDIR, SEED.
# e.g. Sys.setenv(RUN_MODE='smoke'); source('flock_final_decision.R')

FF_VERSION <- 'flock_decision_fixedkernel_v1_20261025'
FF_MODE <- Sys.getenv('RUN_MODE','smoke')
FF_DEFAULTS <- switch(FF_MODE,smoke=c(2L,6000L),pilot=c(10L,30000L),
                     full=c(60L,60000L),stop('RUN_MODE must be smoke/pilot/full'))
FF_REPS <- as.integer(Sys.getenv('REPS',FF_DEFAULTS[1]))
FF_B <- as.integer(Sys.getenv('BUDGET',FF_DEFAULTS[2]))
FF_SEED <- as.integer(Sys.getenv('SEED','20261025'))
FF_ALLOWED_SCENARIOS <- c('KnownMixture5D','DiagnosticLCA','HighDimCorrelated')
FF_ALLOWED_METHODS <- c('IndependentAM','IndependentAM_Shifts',
                         'StagedAM','FLOCK_AM','FLOCK_AM_LocalOnly')
ff_parse_selection <- function(arg,valid) {
  z <- Sys.getenv(arg,'all')
  v <- if(tolower(z)=='all') valid else trimws(strsplit(z,',',fixed=TRUE)[[1L]])
  stopifnot(length(v)>0L,!anyDuplicated(v),all(v%in%valid))
  v
}
FF_SCENARIOS <- ff_parse_selection('SCENARIOS',FF_ALLOWED_SCENARIOS)
FF_METHODS <- ff_parse_selection('METHODS',FF_ALLOWED_METHODS)
FF_OUT <- Sys.getenv('OUTDIR',paste0('flock_final_',FF_MODE,'_B',FF_B))
stopifnot(FF_REPS>=1L,FF_B>=6000L,is.finite(FF_SEED),nzchar(FF_OUT))
dir.create(FF_OUT,showWarnings=FALSE,recursive=TRUE)
dir.create(file.path(FF_OUT,'checkpoints'),showWarnings=FALSE)
ff_manifest <- list(version=FF_VERSION,R=FF_REPS,B=FF_B,seed=FF_SEED,
                    scenarios=FF_SCENARIOS,methods=FF_METHODS)
ff_manifest_file <- file.path(FF_OUT,'manifest.rds')
if(file.exists(ff_manifest_file)) {
  if(!identical(readRDS(ff_manifest_file),ff_manifest))
    stop('Output folder has different settings/code version; change OUTDIR')
} else saveRDS(ff_manifest,ff_manifest_file)

ff_logsum <- function(z) {
  m<-max(z)
  if(!is.finite(m)) return(-Inf)
  m+log(sum(exp(z-m)))
}
ff_spd <- function(S,floor=1e-8,cap=Inf) {
  S<-(S+t(S))/2
  ee<-eigen(S,symmetric=TRUE)
  ee$vectors %*% diag(pmin(pmax(ee$values,floor),cap),length(ee$values)) %*% t(ee$vectors)
}
ff_multinorm_log <- function(x,mu,U) {
  z<-forwardsolve(t(U),x-mu)
  -.5*length(mu)*log(2*pi)-sum(log(diag(U)))-sum(z*z)/2
}
ff_multinorm_draw <- function(mu,U) as.numeric(mu+crossprod(U,rnorm(length(mu))))
ff_counter <- function(target) {
  env<-new.env(parent=emptyenv())
  env$n<-0L
  env$eval<-function(x) {
    env$n<-env$n+1L
    z<-target$logpost(as.numeric(x))
    if(length(z)!=1L || is.na(z)) stop('logpost returned NA/non-scalar')
    as.numeric(z)
  }
  env$count<-function() env$n
  env
}
ff_split_budget <- function(n,k) {
  a<-rep(as.integer(n%/%k),k)
  if(n%%k>0L) a[seq_len(n%%k)]<-a[seq_len(n%%k)]+1L
  a
}

# ---------- Exactly known reference: 5D Gaussian mixture ----------
ff_target_mixture <- function() {
  d<-5L
  weights<-c(.45,.30,.20,.05)
  mu<-rbind(c(-5,-5,0,0,0),c(5,5,0,0,0),
            c(-5,5,3,0,0),c(5,-5,-3,0,0))
  s<-c(.85,1.5,1.0,.55)
  S<-lapply(seq_len(4L),function(i) {
    A<-diag(s[i]^2,d)
    A[1,2]<-A[2,1]<-.35*s[i]^2
    A
  })
  U<-lapply(S,chol)
  lp_parts<-function(x) vapply(seq_len(4L),function(j)
    log(weights[j])+ff_multinorm_log(x,mu[j,],U[[j]]),numeric(1))
  logpost<-function(x) ff_logsum(lp_parts(x))
  resp<-function(x) {z<-lp_parts(x); exp(z-ff_logsum(z))}
  mean_true<-colSums(weights*mu)
  pplus<-sum(weights*pnorm(0,mu[,1],sqrt(vapply(S,function(A) A[1,1],numeric(1))),lower.tail=FALSE))
  ref<-c(mean_theta1=mean_true[1],mean_theta2=mean_true[2],
         p_theta1_gt0=pplus,component3_mass=weights[3],
         rare_component_mass=weights[4])
  cdf<-function(z) sum(weights*pnorm(z,mu[,1],sqrt(vapply(S,function(A) A[1,1],numeric(1)))))
  quant<-vapply(c(.025,.975),function(t) uniroot(function(z)cdf(z)-t,c(-30,30))$root,numeric(1))
  list(name='KnownMixture5D',d=d,logpost=logpost,
       prior_draw=function() rnorm(d,sd=6),base_step=.55,
       ref=ref,scales=c(5,5,.35,.3,.12),
       metric=function(Z) {
         rr<-vapply(seq_len(nrow(Z)),function(i) resp(Z[i,]),numeric(4))
         c(mean_theta1=mean(Z[,1]),mean_theta2=mean(Z[,2]),
           p_theta1_gt0=mean(Z[,1]>0),
           component3_mass=mean(rr[3,]),rare_component_mass=mean(rr[4,]))
       },
       display=function(Z) Z[,1],
       region=function(Z) ifelse(Z[,1]>0,2L,1L),
       region_truth=c(1-pplus,pplus),
       important='rare_component_mass',fail_tol=.035,
       interval=quant,
       meaning='Known target; component memberships are diagnostic only, not given to samplers')
}

# ---------- Exactly known reference: two-test, one-population LCA ----------
ff_target_lca <- function() {
  counts<-c(7L,2L,2L,19L); N<-sum(counts)
  aa<-2;bb<-2
  logsig<-function(x) -pmax(-x,0)-log1p(exp(-abs(x)))
  prior<-function(z) {
    if(length(z)!=5L || any(!is.finite(z))) return(-Inf)
    sum(aa*logsig(z)+bb*logsig(-z)-lbeta(aa,bb))
  }
  probs<-function(z) {
    p<-plogis(z); pi<-p[1];se1<-p[2];se2<-p[3];sp1<-p[4];sp2<-p[5]
    c(pi*se1*se2+(1-pi)*(1-sp1)*(1-sp2),
      pi*se1*(1-se2)+(1-pi)*(1-sp1)*sp2,
      pi*(1-se1)*se2+(1-pi)*sp1*(1-sp2),
      pi*(1-se1)*(1-se2)+(1-pi)*sp1*sp2)
  }
  logpost<-function(z) {
    ll<-prior(z)
    if(!is.finite(ll)) return(-Inf)
    q<-probs(z)
    if(any(q<=0) || any(!is.finite(q))) return(-Inf)
    ll+sum(counts*log(q))
  }
  ks<-as.matrix(expand.grid(lapply(counts,function(n)0:n)))
  kk<-nrow(ks)
  BA<-matrix(0,kk,5L);BB<-matrix(0,kk,5L);lw<-numeric(kk)
  for(i in seq_len(kk)) {
    k<-as.integer(ks[i,]);h<-counts-k
    hits<-c(sum(k), k[1]+k[2],k[1]+k[3],h[3]+h[4],h[2]+h[4])
    misses<-c(N-sum(k),k[3]+k[4],k[2]+k[4],h[1]+h[2],h[1]+h[3])
    BA[i,]<-aa+hits;BB[i,]<-bb+misses
    lw[i]<-sum(lchoose(counts,k))+sum(lbeta(BA[i,],BB[i,])-lbeta(aa,bb))
  }
  wt<-exp(lw-ff_logsum(lw))
  ref<-c(mean_prevalence=sum(wt*BA[,1]/(BA[,1]+BB[,1])),
         p_prevalence_gt065=sum(wt*pbeta(.65,BA[,1],BB[,1],lower.tail=FALSE)),
         mean_Se1=sum(wt*BA[,2]/(BA[,2]+BB[,2])),
         p_prevalence_lt030=sum(wt*pbeta(.30,BA[,1],BB[,1])),
         orientation_prob=.5)
  cdf<-function(x) sum(wt*pbeta(x,BA[,1],BB[,1]))
  qtl<-vapply(c(.025,.975),function(t) uniroot(function(x)cdf(x)-t,c(0,1))$root,numeric(1))
  stopifnot(abs(ref['mean_prevalence']-.5)<1e-7,abs(sum(wt)-1)<1e-9)
  list(name='DiagnosticLCA',d=5L,logpost=logpost,
       prior_draw=function() qlogis(rbeta(5L,aa,bb)),base_step=.55,
       ref=ref,scales=c(.4,.4,.4,.4,.5),
       metric=function(Z) {
         pp<-plogis(Z)
         c(mean_prevalence=mean(pp[,1]),
           p_prevalence_gt065=mean(pp[,1]>.65),
           mean_Se1=mean(pp[,2]),
           p_prevalence_lt030=mean(pp[,1]<.30),
           orientation_prob=mean(pp[,1]>.5))
       },
       display=function(Z) plogis(Z[,1]),
       region=function(Z) ifelse(plogis(Z[,1])>.5,2L,1L),
       region_truth=c(.5,.5),
       important='p_prevalence_gt065',fail_tol=.10,
       interval=qtl,
       meaning='Exact Beta-mixture posterior; orientations reflect disease-label symmetry')
}

# ---------- Exactly known reference: 40D spike/slab prior regression ----------
ff_target_highdim <- function() {
  set.seed(4262)
  n<-72L;d<-40L;sigma<-.90
  root<-rnorm(n)
  X<-matrix(0,n,d)
  X[,1]<-as.numeric(scale(root))
  X[,2]<-as.numeric(scale(.9992*X[,1]+sqrt(1-.9992^2)*rnorm(n)))
  groups<-replicate(5L,rnorm(n))
  for(j in 3:d) {
    g<-((j-3L)%%5L)+1L
    X[,j]<-as.numeric(scale(.55*groups[,g]+.84*rnorm(n)))
  }
  beta<-numeric(d);beta[c(1L,5L,9L)]<-c(2.15,.55,-.55)
  y<-as.numeric(X%*%beta+rnorm(n,sd=sigma))
  G<-crossprod(X);bb<-as.numeric(crossprod(X,y));yy<-sum(y*y)
  spike<-.10;slab<-2.6;slab_prob<-.45;rest_sd<-.55
  logp_mix<-function(v) vapply(v,function(z) ff_logsum(c(
    log(1-slab_prob)+dnorm(z,0,spike,log=TRUE),
    log(slab_prob)+dnorm(z,0,slab,log=TRUE))),numeric(1))
  prior<-function(z) {
    if(length(z)!=d || any(!is.finite(z))) return(-Inf)
    sum(logp_mix(z[1:2]))+sum(dnorm(z[-c(1L,2L)],0,rest_sd,log=TRUE))
  }
  logpost<-function(z) {
    pr<-prior(z)
    if(!is.finite(pr)) return(-Inf)
    quad<-yy-2*sum(z*bb)+as.numeric(crossprod(z,G%*%z))
    pr-n/2*log(2*pi*sigma^2)-quad/(2*sigma^2)
  }
  comb<-as.matrix(expand.grid(c(0L,1L),c(0L,1L)))
  mu<-matrix(0,4L,d);SS<-vector('list',4L);lm<-numeric(4L)
  for(j in 1:4) {
    vv<-c(ifelse(comb[j,]==1L,slab^2,spike^2),rep(rest_sd^2,d-2L))
    Q<-diag(1/vv,d)+G/sigma^2
    SS[[j]]<-solve(Q)
    mu[j,]<-as.numeric(SS[[j]]%*%(bb/sigma^2))
    ld<-as.numeric(determinant(Q,logarithm=TRUE)$modulus)
    lm[j]<-sum(ifelse(comb[j,]==1L,log(slab_prob),log(1-slab_prob)))-
      .5*(sum(log(vv))+ld)+.5*as.numeric(crossprod(mu[j,],Q%*%mu[j,]))
  }
  wt<-exp(lm-ff_logsum(lm))
  cm<-mu[,1]-mu[,2]
  cs<-vapply(SS,function(S) sqrt(S[1,1]+S[2,2]-2*S[1,2]),numeric(1))
  ref<-c(mean_beta1=sum(wt*mu[,1]),mean_beta2=sum(wt*mu[,2]),
         mean_contrast=sum(wt*cm),
         p_beta1_gt_beta2=sum(wt*pnorm(0,cm,cs,lower.tail=FALSE)),
         p_contrast_gt1=sum(wt*pnorm(1,cm,cs,lower.tail=FALSE)),
         mean_beta5=sum(wt*mu[,5]))
  cdf<-function(x) sum(wt*pnorm(x,mu[,1],sqrt(vapply(SS,function(S)S[1,1],numeric(1)))))
  quant<-vapply(c(.025,.975),function(t)uniroot(function(x)cdf(x)-t,c(-20,20))$root,numeric(1))
  list(name='HighDimCorrelated',d=d,logpost=logpost,
       prior_draw=function()c(rnorm(2L,sd=ifelse(runif(2L)<slab_prob,slab,spike)),
                               rnorm(d-2L,sd=rest_sd)),base_step=.040,
       ref=ref,scales=c(1.5,1.5,2,.5,.5,.7),
       metric=function(Z) {
         contrast<-Z[,1]-Z[,2]
         c(mean_beta1=mean(Z[,1]),mean_beta2=mean(Z[,2]),
           mean_contrast=mean(contrast),
           p_beta1_gt_beta2=mean(contrast>0),
           p_contrast_gt1=mean(contrast>1),mean_beta5=mean(Z[,5]))
       },
       display=function(Z) Z[,1],
       region=function(Z) ifelse(Z[,1]-Z[,2]>0,2L,1L),
       region_truth=c(1-ref['p_beta1_gt_beta2'],ref['p_beta1_gt_beta2']),
       important='p_beta1_gt_beta2',fail_tol=.20,
       interval=quant,
       meaning='Exact 4-normal posterior mixture; sign of contrast is not necessarily a density mode')
}

ff_targets<-list(KnownMixture5D=ff_target_mixture,
                 DiagnosticLCA=ff_target_lca,
                 HighDimCorrelated=ff_target_highdim)

# ----------- Adaptive Metropolis kernel shared by EVERY method -----------
# Online covariance uses Welford's recurrence, including repeated rejected states.
# Local AM adapts during warmup only; then proposal is frozen.
ff_am_segment <- function(f,target,x,lp=NULL,evals,
                          adapt_frac=.35,initial_cov=NULL,
                          shifts=NULL,jump_prob=0,record_history=TRUE) {
  d<-target$d
  stopifnot(evals>=1L,adapt_frac>=0,adapt_frac<=1,jump_prob>=0,jump_prob<1)
  if(is.null(lp)) {
    lp<-f$eval(x);steps<-evals-1L
  } else steps<-evals
  h<-matrix(NA_real_,steps+1L,d)
  lps<-numeric(steps+1L)
  h[1,]<-x;lps[1]<-lp
  base_cov<-diag(target$base_step^2,d)
  Q<-if(is.null(initial_cov)) base_cov else ff_spd(initial_cov,1e-9)
  U<-chol(Q)
  warmup<-as.integer(floor(steps*adapt_frac))
  updates<-0L;nlocal<-0L;alocal<-0L;njump<-0L;ajump<-0L
  m<-as.numeric(x);M2<-matrix(0,d,d); nn<-1L
  for(t in seq_len(steps)) {
    use_jump<-length(shifts)>0L && runif(1)<jump_prob
    if(use_jump) {
      j<-sample.int(length(shifts),1L)
      sign<-sample(c(-1,1),1L)
      # Symmetric +/- shift mixture with state-INDEPENDENT small Gaussian noise.
      delta<-shifts[[j]]
      y<-x+sign*delta+rnorm(d,sd=target$base_step*.15/sqrt(d))
      njump<-njump+1L
    } else {
      y<-ff_multinorm_draw(x,U)
      nlocal<-nlocal+1L
    }
    ly<-f$eval(y)
    if(is.finite(ly) && log(runif(1)) < ly-lp) {
      x<-y;lp<-ly
      if(use_jump) ajump<-ajump+1L else alocal<-alocal+1L
    }
    h[t+1L,]<-x;lps[t+1L]<-lp
    if(t<=warmup) {
      nn<-nn+1L
      dd<-x-m
      m<-m+dd/nn
      M2<-M2+tcrossprod(dd,x-m)
      if(t>=max(100L,2L*d) && t%%100L==0L) {
        S<-M2/(nn-1L)
        # Shrinkage/regularization prevents singularity and overexpansion.
        Q<-(2.38^2/d)*S+.08*base_cov
        Q<-ff_spd(Q,floor=(target$base_step^2)*.02,
                  cap=max(10,target$base_step^2*100*d))
        U<-chol(Q);updates<-updates+1L
      }
    }
  }
  ret<-if(warmup<steps) h[(warmup+2L):(steps+1L),,drop=FALSE] else matrix(numeric(),0,d)
  list(x=x,lp=lp,h=h,logp=lps,retained=ret,cov=Q,
       accept_local=if(nlocal)alocal/nlocal else NA_real_,
       accept_jump=if(njump)ajump/njump else NA_real_,
       nlocal=nlocal,njump=njump,accept_jumps=ajump,
       adapt_updates=updates,retained_n=nrow(ret),warmup=warmup)
}

# ---------- A symmetric region-shift MH kernel ----------
# Construct tentative centers using *only* exploratory visited states.
# No true mixture centers, model labels, or true posterior values are consulted.
ff_make_shifts <- function(chains,d,max_k=4L,max_dims=6L) {
  pool<-do.call(rbind,lapply(chains,function(z) {
    h<-z$h
    ids<-unique(as.integer(round(seq(1,nrow(h),length.out=min(175L,nrow(h))))))
    h[ids,,drop=FALSE]
  }))
  pool<-unique(pool)
  if(nrow(pool)<30L) return(list())
  scale<-pmax(apply(pool,2,sd),.05)
  Z<-sweep(pool,2,scale,'/')
  k<-min(max_k,max(2L,as.integer(floor(nrow(pool)/80L))))
  k<-min(k,nrow(unique(Z)))
  if(k<2L) return(list())
  fit<-tryCatch(stats::kmeans(Z,centers=k,nstart=2L,iter.max=35L),
                error=function(e) NULL)
  if(is.null(fit)) return(list())
  centers<-sweep(fit$centers,2,scale,'*')
  shifts<-list()
  for(i in seq_len(k-1L)) for(j in (i+1L):k) {
    diff<-as.numeric(centers[j,]-centers[i,])
    # Non-oracle sparse shift in the most separated standardized coordinates.
    if(d>max_dims) {
      keep<-order(abs(diff)/scale,decreasing=TRUE)[seq_len(max_dims)]
      diff[-keep]<-0
    }
    if(sqrt(sum((diff/scale)^2))>1.25 && any(abs(diff)>1e-10))
      shifts[[length(shifts)+1L]]<-diff
  }
  shifts
}

# The +/- shift proposal is symmetric, including the Gaussian jitter.
# This numerical self-check is a guard against regression of the earlier bug.
ff_check_symmetry <- function() {
  set.seed(711)
  U<-chol(matrix(c(1,.4,.4,1),2,2))
  z<-replicate(4000L,ff_multinorm_draw(c(0,0),U))
  stopifnot(max(abs(cov(t(z))-crossprod(U)))<.075)
  ss<-list(c(2,-1),c(-.5,3))
  lmix<-function(z) ff_logsum(unlist(lapply(ss,function(s)
    c(sum(dnorm(z-s,0,.15,log=TRUE)),sum(dnorm(z+s,0,.15,log=TRUE))))))
  stopifnot(abs(lmix(c(1.8,-1.2))-lmix(c(-1.8,1.2)))<1e-10)
  TRUE
}
stopifnot(ff_check_symmetry())

# Use accumulated explorer states to initialize the local covariance of the
# final chain. The covariance is estimated in the neighborhood of the selected
# starting state, rather than pooled across distant modes.
ff_local_covariance <- function(histories,center,target) {
  d<-target$d
  pool<-do.call(rbind,lapply(histories,function(a) {
    h<-a$h
    ix<-unique(as.integer(round(seq(1,nrow(h),length.out=min(175L,nrow(h))))))
    h[ix,,drop=FALSE]
  }))
  if(nrow(pool)<2L*d+10L) return(diag(target$base_step^2,d))
  scales<-pmax(apply(pool,2,sd),target$base_step*.5)
  delta<-sweep(pool,2,center,'-')
  distance<-rowSums(sweep(delta,2,scales,'/')^2)
  nn<-min(nrow(pool),max(2L*d+30L,floor(nrow(pool)*.20)))
  ids<-order(distance)[seq_len(nn)]
  S<-cov(pool[ids,,drop=FALSE])
  ff_spd((2.38^2/d)*S+.08*diag(target$base_step^2,d),
         floor=target$base_step^2*.02,
         cap=max(10,target$base_step^2*100*d))
}

# ---------- Fixed scout; same candidate banks for Staged and FLOCK ----------
ff_start_from_pool <- function(h,lp,k=4L) {
  # highest-density visited point and several near-high states
  ord<-order(lp,decreasing=TRUE)
  as.matrix(h[ord[seq_len(min(k,length(ord)))],,drop=FALSE])
}
ff_middle_candidates <- function(target,scout,j,seed) {
  set.seed(seed)
  anchor<-ff_start_from_pool(scout$h,scout$logp,4L)
  guided<-t(vapply(seq_len(4L),function(i)
    as.numeric(anchor[((i+j-2L)%%nrow(anchor))+1L,]+
                 rnorm(target$d,sd=target$base_step*1.8)),numeric(target$d)))
  prior<-t(replicate(4L,target$prior_draw()))
  rbind(guided,prior)
}
ff_pick_start <- function(target,strategy,j,scout,candidates,lp) {
  if(strategy=='staged') return(sample.int(nrow(candidates),1L))
  if(j<=3L) return(which.max(lp[1:4]))
  # Novelty pioneers score relative to *scout* history, not true modes.
  h<-scout$h
  ids<-unique(as.integer(round(seq(1,nrow(h),length.out=min(140L,nrow(h))))))
  H<-h[ids,,drop=FALSE]
  scale<-pmax(apply(H,2,sd),target$base_step*2)
  novelty<-vapply(seq_len(nrow(candidates)),function(k) {
    z<-sweep(H,2,candidates[k,],'-')
    sqrt(min(rowSums(sweep(z,2,scale,'/')^2)))
  },numeric(1))
  plausible_penalty<-pmax(0,max(lp)-lp)
  utility<-novelty-.08*plausible_penalty
  which.max(utility)
}

ff_run_staged <- function(target,B,strategy=c('staged','flock'),
                          final_jumps=TRUE,seed=1L) {
  strategy<-match.arg(strategy)
  f<-ff_counter(target); t0<-proc.time()[['elapsed']]
  bs<-as.integer(floor(B*.15))
  bm<-as.integer(floor(B*.30))
  alloc<-ff_split_budget(bm,6L)
  set.seed(seed+10L)
  x0<-target$prior_draw()
  scout<-ff_am_segment(f,target,x0,evals=bs,adapt_frac=.70)
  explorers<-vector('list',6L)
  for(j in seq_len(6L)) {
    cand<-ff_middle_candidates(target,scout,j,seed+1000L+j)
    lp_cand<-apply(cand,1,f$eval)
    set.seed(seed+2000L+j)
    index<-ff_pick_start(target,strategy,j,scout,cand,lp_cand)
    # Every method gets 8 counted candidate evaluations per middle chain.
    explorers[[j]]<-ff_am_segment(f,target,cand[index,],lp_cand[index],
                                  evals=alloc[j]-nrow(cand),adapt_frac=.70)
  }
  stopifnot(f$count()==bs+bm)
  histories<-c(list(scout),explorers)
  set.seed(seed+3000L)
  shifts<-ff_make_shifts(histories,target$d)
  # Choose same initialization rule for staged and FLOCK: best visited state.
  lpall<-unlist(lapply(histories,function(a) a$logp))
  imax<-which.max(lpall)
  lens<-vapply(histories,function(a) length(a$logp),integer(1))
  owner<-which(imax<=cumsum(lens))[1L]
  prev<-if(owner==1L)0L else sum(lens[seq_len(owner-1L)])
  xfinal<-histories[[owner]]$h[imax-prev,]
  lpfinal<-histories[[owner]]$logp[imax-prev]
  informed_cov<-ff_local_covariance(histories,xfinal,target)
  final<-ff_am_segment(f,target,xfinal,lpfinal,evals=B-f$count(),
                       adapt_frac=.18,initial_cov=informed_cov,
                       shifts=if(final_jumps)shifts else list(),
                       jump_prob=if(final_jumps).12 else 0)
  stopifnot(f$count()==B)
  list(draws=final$retained,evals=f$count(),
       local_accept=final$accept_local,jump_accept=final$accept_jump,
       accepted_jumps=final$accept_jumps,jump_attempts=final$njump,
       shifts=length(shifts),cov_updates=sum(vapply(histories,function(a)a$adapt_updates,numeric(1)))+
         final$adapt_updates,
       scout_regions=length(unique(target$region(scout$h))),
       explored_regions=length(unique(target$region(do.call(rbind,lapply(histories,function(a)a$h))))),
       final_regions=length(unique(target$region(final$retained))),
       seconds=proc.time()[['elapsed']]-t0,
       stage_info=list(scout_evals=bs,middle_evals=bm,final_evals=B-bs-bm))
}

# ---------- Independent AM with optional ensemble-shared symmetric jumps ----------
ff_run_independent <- function(target,B,use_shifts=FALSE,seed=1L) {
  f<-ff_counter(target);t0<-proc.time()[['elapsed']]
  alloc<-ff_split_budget(B,8L)
  chains<-vector('list',8L)
  if(!use_shifts) {
    for(i in 1:8) {
      set.seed(seed+4000L+i)
      chains[[i]]<-ff_am_segment(f,target,target$prior_draw(),
                                 evals=alloc[i],adapt_frac=.35)
    }
    shifts<-list();ret<-lapply(chains,function(a)a$retained)
    naccept<-0L;nattempt<-0L
  } else {
    warm<-pmax(2L,as.integer(floor(alloc*.45)))
    pre<-vector('list',8L)
    for(i in 1:8) {
      set.seed(seed+4000L+i)
      pre[[i]]<-ff_am_segment(f,target,target$prior_draw(),
                              evals=warm[i],adapt_frac=.95)
    }
    set.seed(seed+4100L)
    shifts<-ff_make_shifts(pre,target$d)
    for(i in 1:8) {
      set.seed(seed+5000L+i)
      chains[[i]]<-ff_am_segment(f,target,pre[[i]]$x,pre[[i]]$lp,
          evals=alloc[i]-warm[i],adapt_frac=0,
          initial_cov=pre[[i]]$cov,shifts=shifts,jump_prob=.12)
    }
    ret<-lapply(chains,function(a)a$retained)
    naccept<-sum(vapply(chains,function(a)a$accept_jumps,integer(1)))
    nattempt<-sum(vapply(chains,function(a)a$njump,integer(1)))
  }
  stopifnot(f$count()==B)
  allret<-do.call(rbind,ret)
  list(draws=allret,evals=f$count(),
       local_accept=mean(vapply(chains,function(a)a$accept_local,numeric(1)),na.rm=TRUE),
       jump_accept=if(nattempt) naccept/nattempt else NA_real_,
       accepted_jumps=naccept,jump_attempts=nattempt,shifts=length(shifts),
       cov_updates=sum(vapply(chains,function(a)a$adapt_updates,integer(1)))+
         if(use_shifts)sum(vapply(pre,function(a)a$adapt_updates,integer(1))) else 0L,
       scout_regions=NA_integer_,explored_regions=NA_integer_,
       final_regions=length(unique(target$region(allret))),
       seconds=proc.time()[['elapsed']]-t0,
       stage_info=NULL,per_chain=ret)
}

# ---------- Post-hoc evaluation; target reference NEVER given to samplers ----------
ff_evaluate <- function(target,res,scenario,method,rep_id) {
  Z<-res$draws
  est<-target$metric(Z)
  stopifnot(identical(names(est),names(target$ref)),all(is.finite(est)))
  dif<-(est-target$ref)/target$scales
  primary<-target$important
  region<-target$region(Z)
  if(!is.null(res$per_chain)) {
    switches<-sum(vapply(res$per_chain,function(a)sum(diff(target$region(a))!=0),integer(1)))
  } else switches<-sum(diff(region)!=0)
  qq<-as.numeric(stats::quantile(target$display(Z),c(.025,.975)))
  data.frame(scenario=scenario,method=method,replicate=rep_id,
             evals=res$evals,seconds=res$seconds,
             normalized_rmse=sqrt(mean(dif^2)),
             important_estimate=unname(est[primary]),
             important_reference=unname(target$ref[primary]),
             important_error=abs(unname(est[primary]-target$ref[primary])),
             important_failure=as.integer(abs(est[primary]-target$ref[primary])>target$fail_tol),
             interval_l1=sum(abs(qq-target$interval)),
             region_switches=switches,
             region_1_fraction=mean(region==1L),
             final_covered_both=as.integer(length(unique(region))>=2L),
             n_retained=nrow(Z),local_accept=res$local_accept,
             jump_accept=res$jump_accept,
             accepted_jumps=res$accepted_jumps,jump_attempts=res$jump_attempts,
             shifts=res$shifts,cov_updates=res$cov_updates,
             explored_regions=res$explored_regions,
             final_regions=res$final_regions,
             lost_known_regions=if(is.na(res$explored_regions)) NA_integer_ else
               as.integer(res$explored_regions>res$final_regions),
             no_accepted_jumps=as.integer(res$jump_attempts>0L && res$accepted_jumps==0L),
             stringsAsFactors=FALSE)
}
ff_metric_rows<-function(target,res,scenario,method,rep_id) {
  v<-target$metric(res$draws)
  data.frame(scenario=scenario,method=method,replicate=rep_id,
             metric=names(v),estimate=as.numeric(v),
             reference=as.numeric(target$ref),
             error=as.numeric(v-target$ref))
}
ff_write_csv<-function(df,path) {
  tmp<-paste0(path,'.tmp')
  utils::write.csv(df,tmp,row.names=FALSE)
  if(file.exists(path)) file.remove(path)
  if(!file.rename(tmp,path)) stop('Failed to move CSV into place: ',path)
}
ff_bootstrap<-function(a,b,seed,boots=1000L) {
  d<-a-b
  set.seed(seed)
  sims<-replicate(boots,mean(sample(d,length(d),replace=TRUE)))
  c(mean_diff=mean(d),lo=unname(quantile(sims,.025)),hi=unname(quantile(sims,.975)),
    frac_first_better=mean(d<0))
}
ff_rebuild<-function(targets) {
  paths<-list.files(file.path(FF_OUT,'checkpoints'),pattern='\\.rds$',full.names=TRUE)
  if(!length(paths)) return(invisible(NULL))
  objs<-lapply(paths,readRDS)
  summary_rows<-do.call(rbind,lapply(objs,function(x)x$summary))
  metric_rows<-do.call(rbind,lapply(objs,function(x)x$metrics))
  summary_rows<-summary_rows[order(summary_rows$scenario,summary_rows$replicate,
                                   summary_rows$method),]
  ff_write_csv(summary_rows,file.path(FF_OUT,'replicate_results.csv'))
  ff_write_csv(metric_rows,file.path(FF_OUT,'metric_estimates.csv'))
  sums<-do.call(rbind,lapply(split(summary_rows,list(summary_rows$scenario,
                                 summary_rows$method),drop=TRUE),function(z) {
    data.frame(scenario=z$scenario[1],method=z$method[1],n=nrow(z),
               mean_rmse=mean(z$normalized_rmse),median_rmse=median(z$normalized_rmse),
               mean_important_error=mean(z$important_error),
               failure_rate=mean(z$important_failure),
               mean_interval_l1=mean(z$interval_l1),
               mean_seconds=mean(z$seconds),
               mean_accepted_jumps=mean(z$accepted_jumps),
               mean_retained=mean(z$n_retained))
  }))
  ff_write_csv(sums,file.path(FF_OUT,'method_summary.csv'))
  # Paired bootstrap comparisons only over replicate IDs present in both methods.
  comparisons<-list()
  pairs<-list(c('FLOCK_AM','StagedAM'),
              c('FLOCK_AM','IndependentAM'),
              c('FLOCK_AM','IndependentAM_Shifts'),
              c('FLOCK_AM','FLOCK_AM_LocalOnly'),
              c('IndependentAM_Shifts','IndependentAM'))
  for(sc in unique(summary_rows$scenario)) for(pair in pairs) {
    a<-summary_rows[summary_rows$scenario==sc & summary_rows$method==pair[1],]
    b<-summary_rows[summary_rows$scenario==sc & summary_rows$method==pair[2],]
    mm<-merge(a,b,by='replicate',suffixes=c('_A','_B'))
    if(nrow(mm)<2L) next
    for(metric in c('normalized_rmse','important_error','interval_l1')) {
      result<-ff_bootstrap(mm[[paste0(metric,'_A')]],
                           mm[[paste0(metric,'_B')]],
                           FF_SEED+nrow(mm)+nchar(sc)+nchar(metric))
      comparisons[[length(comparisons)+1L]]<-data.frame(
        scenario=sc,metric=metric,method_A=pair[1],method_B=pair[2],
        n=nrow(mm),t(result),row.names=NULL)
    }
  }
  if(length(comparisons)) {
    ff_write_csv(do.call(rbind,comparisons),file.path(FF_OUT,'paired_bootstrap.csv'))
  }
  invisible(sums)
}

# Avoid generating figures from assumed 'truth' in a way that leaks to samplers.
ff_plot_results<-function() {
  file<-file.path(FF_OUT,'replicate_results.csv')
  if(!file.exists(file)) return(invisible(NULL))
  dat<-read.csv(file)
  pdf(file.path(FF_OUT,'diagnostics.pdf'),width=12,height=8)
  for(sc in unique(dat$scenario)) {
    x<-dat[dat$scenario==sc,,drop=FALSE]
    par(mfrow=c(2,2),mar=c(8,4,3,1))
    boxplot(normalized_rmse~method,data=x,las=2,main=paste(sc,'RMSE'),ylab='Normalized RMSE')
    boxplot(important_error~method,data=x,las=2,main='Consequential event error',ylab='Absolute error')
    boxplot(seconds~method,data=x,las=2,main='Wall time',ylab='Seconds')
    boxplot(region_switches~method,data=x,las=2,main='Within-chain region transitions',
            ylab='Transitions (independent methods sum per-chain)')
  }
  dev.off()
}

# ---------- Execute with checkpoint per scenario x replicate x method ----------
ff_available<-lapply(FF_SCENARIOS,function(n)ff_targets[[n]]())
names(ff_available)<-FF_SCENARIOS
for(n in FF_SCENARIOS) {
  tar<-ff_available[[n]]
  stopifnot(is.finite(tar$logpost(tar$prior_draw())),
            all(is.finite(tar$ref)),all(tar$scales>0),
            all(is.finite(tar$interval)),length(tar$region_truth)==2L)
}
ref_file<-do.call(rbind,lapply(ff_available,function(t)
  data.frame(scenario=t$name,metric=names(t$ref),reference=as.numeric(t$ref))))
ff_write_csv(ref_file,file.path(FF_OUT,'reference_metrics.csv'))
cat('\nFLOCK final decision experiment | ',FF_MODE,' | ',FF_REPS,' replicates | B=',FF_B,'\n',sep='')
cat('Scenarios:',paste(FF_SCENARIOS,collapse=', '),'\n')
cat('Methods:',paste(FF_METHODS,collapse=', '),'\n')
for(s in FF_SCENARIOS) {
  tar<-ff_available[[s]]
  cat('\nScenario',s,'| dimension',tar$d,'| reference:',tar$meaning,'\n')
  for(r in seq_len(FF_REPS)) {
    for(method in FF_METHODS) {
      checkpoint<-file.path(FF_OUT,'checkpoints',paste(s,sprintf('%03d',r),method,sep='__'))
      checkpoint<-paste0(checkpoint,'.rds')
      if(file.exists(checkpoint)) next
      # Deterministic scenario-specific seed, paired across strategies.
      si<-match(s,FF_ALLOWED_SCENARIOS)
      seed<-FF_SEED+si*1000000L+r*10000L
      ans<-switch(method,
        IndependentAM=ff_run_independent(tar,FF_B,FALSE,seed),
        IndependentAM_Shifts=ff_run_independent(tar,FF_B,TRUE,seed),
        StagedAM=ff_run_staged(tar,FF_B,'staged',TRUE,seed),
        FLOCK_AM=ff_run_staged(tar,FF_B,'flock',TRUE,seed),
        FLOCK_AM_LocalOnly=ff_run_staged(tar,FF_B,'flock',FALSE,seed))
      stopifnot(ans$evals==FF_B,all(is.finite(ans$draws)),nrow(ans$draws)>20L)
      report<-list(summary=ff_evaluate(tar,ans,s,method,r),
                   metrics=ff_metric_rows(tar,ans,s,method,r))
      # Preserve posterior draws only for replicate 1 to limit disk usage.
      if(r==1L) {
        report$draws<-ans$draws
        report$per_chain<-ans$per_chain
      }
      temp<-paste0(checkpoint,'.tmp')
      saveRDS(report,temp)
      if(!file.rename(temp,checkpoint)) stop('Cannot write checkpoint ',checkpoint)
      cat(sprintf('[%s %d/%d] %-22s RMSE %.4f | event err %.4f | jumps %d/%d | %.1fs\n',
          s,r,FF_REPS,method,report$summary$normalized_rmse,
          report$summary$important_error,ans$accepted_jumps,ans$jump_attempts,ans$seconds))
    }
    ff_rebuild(ff_available)
  }
}
ff_rebuild(ff_available)
ff_plot_results()
cat('\nFinished. Output directory:',normalizePath(FF_OUT),'\n')
cat('Read method_summary.csv and paired_bootstrap.csv with the decision criteria in README.\n')