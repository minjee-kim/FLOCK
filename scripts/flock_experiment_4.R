# FLOCK Experiment 4: Bayesian mixture posterior with label switching
# Base R only. Rscript flock_experiment_4.R
# Smoke test: Sys.setenv(REPS='3', BUDGET='12000'); source('flock_experiment_4.R')
# Reference truth computed by exact enumeration of latent assignments.
set.seed(20261025)
D <- 3L
# 12 observations, two unknown component means and unknown mixing proportion.
y <- c(-4.4,-4.0,-3.8,-3.6,-3.4,-3.0,3.0,3.3,3.5,3.8,4.0,4.3)
nobs <- length(y); obs_sd <- .85; prior_sd <- 5; a0 <- 2
logsumexp <- function(z) { m<-max(z);m+log(sum(exp(z-m))) }
# x=(mu1,mu2,logit(p)); symmetric priors imply two separated label-switched modes.
log_target <- function(x) {
 p<-plogis(x[3]);if(p<=0||p>=1)return(-Inf)
 z<-dnorm(x[1],0,prior_sd,log=TRUE)+dnorm(x[2],0,prior_sd,log=TRUE)+
    dbeta(p,a0,a0,log=TRUE)+log(p)+log1p(-p)
 l1<-log(p)+dnorm(y,x[1],obs_sd,log=TRUE)
 l2<-log1p(-p)+dnorm(y,x[2],obs_sd,log=TRUE)
 z+sum(vapply(seq_along(y),function(i)logsumexp(c(l1[i],l2[i])),numeric(1)))
}
# Analytic reference: enumerate 2^12 latent allocations and integrate
# normal means and beta mixing weight. Used ONLY for scoring, never for proposals.
log_marg <- function(z) {
 k<-length(z);if(k==0)return(0)
 v<-1/(1/prior_sd^2+k/obs_sd^2);m<-v*sum(z)/obs_sd^2
 -k/2*log(2*pi*obs_sd^2)+.5*log(v/prior_sd^2)-sum(z*z)/(2*obs_sd^2)+m*m/(2*v)
}
N<-2^nobs; tab<-matrix(NA_real_,N,6); lw<-numeric(N)
for(i in 0:(N-1L)) {
 b<-as.logical(as.integer(intToBits(i))[seq_len(nobs)]);n1<-sum(b);n2<-nobs-n1
 v1<-1/(1/prior_sd^2+n1/obs_sd^2);v2<-1/(1/prior_sd^2+n2/obs_sd^2)
 tab[i+1L,]<-c(v1*sum(y[b])/obs_sd^2,v2*sum(y[!b])/obs_sd^2,
                sqrt(v1),sqrt(v2),a0+n1,a0+n2)
 lw[i+1L]<-lbeta(a0+n1,a0+n2)-lbeta(a0,a0)+log_marg(y[b])+log_marg(y[!b])
}
ww<-exp(lw-logsumexp(lw))
set.seed(987654)
ii<-sample.int(N,100000L,replace=TRUE,prob=ww)
ref<-cbind(rnorm(length(ii),tab[ii,1],tab[ii,3]),
           rnorm(length(ii),tab[ii,2],tab[ii,4]),
           qlogis(rbeta(length(ii),tab[ii,5],tab[ii,6])))
features<-function(z) {
 p<-plogis(z[,3]);c(low_mean=mean(pmin(z[,1],z[,2])),
 high_mean=mean(pmax(z[,1],z[,2])),
 separation=mean(abs(z[,1]-z[,2])),
 mixing_uncertainty=mean(p*(1-p)),
 label_probability=mean(z[,1]<z[,2]),
 predictive_tail=mean(p*pnorm(-5,z[,1],obs_sd)+(1-p)*pnorm(-5,z[,2],obs_sd)))
}
truth<-features(ref)
cat('Exact-mixture reference functionals (oracle MC approximation):\n');print(truth)
# A mode means one of the two label-switching configurations.
mode_label <- function(x) if(x[1]<x[2]) 1L else 2L
# All target calls are counted by explicit wrappers.
make_counter <- function() {
  e <- new.env(parent=emptyenv()); e$n <- 0L
  e$eval <- function(x) { e$n <- e$n+1L; log_target(x) }
  e$count <- function() e$n
  e
}
rmv <- function(mu, L) as.numeric(mu + L%*%rnorm(length(mu)))
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
  covs[[length(covs)+1L]] <- diag(c(9^2,9^2,5^2),D)
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
local_walk <- function(f,x,lp,n,step) {
 h <- matrix(NA_real_,n+1L,D);h[1,] <- x
 for(i in seq_len(n)) {
  y <- x+rnorm(D,sd=step);ly <- f$eval(y)
  if(log(runif(1))<ly-lp){x<-y;lp<-ly}
  h[i+1L,] <- x
 }
 list(x=x,lp=lp,h=h)
}
coverage <- function(cand,history) {
 ids <- unique(round(seq(1,nrow(history),length.out=min(250L,nrow(history)))))
 h <- history[ids,,drop=FALSE];scale <- pmax(apply(h,2,sd),2)
 vapply(seq_len(nrow(cand)),function(i) {
  sqrt(min(rowSums(sweep(sweep(h,2,cand[i,],'-'),2,scale,'/')^2)))
 },numeric(1))
}
fit_q <- function(histories) {
 pool <- do.call(rbind,lapply(histories,function(h){
  ids<-unique(round(seq(1,nrow(h),length.out=min(150L,nrow(h)))))
  h[ids,,drop=FALSE]
 }))
 pool <- unique(pool);n <- min(8L,nrow(pool)-1L)
 km <- tryCatch(kmeans(pool,centers=n,nstart=5,iter.max=60),error=function(e) kmeans(pool,centers=3,nstart=1))
 n <- nrow(km$centers)
 covs <- lapply(seq_len(n),function(j){
  z<-pool[km$cluster==j,,drop=FALSE]
  if(nrow(z)<D+3L) diag(2,D) else cov(z)+diag(.25,D)
 })
 make_proposal(km$centers,covs,rep(1,n),broad_weight=.25)
}
final_chain <- function(f,q,x,lp,n) {
 h<-matrix(NA_real_,n+1L,D);h[1,]<-x;lq<-q$logq(x);accept<-0L
 for(i in seq_len(n)) {
  if(runif(1)<.2) {
   y<-x+rnorm(D,sd=.8);ly<-f$eval(y);lqy<-q$logq(y);loga<-ly-lp
  } else {
   y<-q$draw();ly<-f$eval(y);lqy<-q$logq(y);loga<-ly-lp+lq-lqy
  }
  if(log(runif(1))<loga){x<-y;lp<-ly;lq<-lqy;accept<-accept+1L}
  h[i+1L,]<-x
 }
 burn<-floor(.1*nrow(h))
 list(draws=h[(burn+1L):nrow(h),,drop=FALSE],accept=accept/n)
}
choose_candidate <- function(policy,cand,lp,scouts,previous) {
 n<-nrow(cand)
 if(policy=='random') return(sample.int(n,1))
 nov<-coverage(cand,scouts)
 if(policy=='unexplored') { eligible<-which(nov>=median(nov));return(eligible[sample.int(length(eligible),1L)]) }
 if(policy=='density') {
  eligible<-which(lp>=max(lp)-8)
  return(eligible[which.max(nov[eligible])])
 }
 if(policy=='novelty') return(which.max(nov-.15*pmax(0,max(lp)-lp)))
 if(policy=='diversity') {
  if(nrow(previous)==0L) return(sample.int(n,1))
  scale<-pmax(apply(scouts,2,sd),2)
  d<-vapply(seq_len(n),function(i){
   min(sqrt(rowSums(sweep(sweep(previous,2,cand[i,],'-'),2,scale,'/')^2)))
  },numeric(1))
  return(which.max(d))
 }
 stop('Unknown policy')
}
flock <- function(B,x0,policy) {
 f<-make_counter();bs<-floor(.20*B);bm<-floor(.45*B)
 starts<-matrix(rnorm(2*D,sd=4),2,D);starts[1,]<-x0
 scouts<-vector('list',2);ends<-matrix(NA_real_,2,D);lps<-numeric(2)
 ns<-rep((bs-2L)%/%2L,2);ns[1]<-ns[1]+(bs-2L-sum(ns))
 for(k in 1:2) {
  lp<-f$eval(starts[k,]);z<-local_walk(f,starts[k,],lp,ns[k],2.5)
  scouts[[k]]<-z$h;ends[k,]<-z$x;lps[k]<-z$lp
 }
 scout_pool<-do.call(rbind,scouts);ord<-order(lps,decreasing=TRUE)
 alloc<-rep(bm%/%6L,6L);if(bm%%6L>0) alloc[seq_len(bm%%6L)]<-alloc[seq_len(bm%%6L)]+1L
 followers<-list();pioneers<-list();previous<-matrix(numeric(0),0,D)
 for(k in 1:3) {
  j<-ord[(k-1L)%%2L+1L];x<-ends[j,];lp<-lps[j];n<-alloc[k]
  if(k>2){x<-x+rnorm(D,sd=.7);lp<-f$eval(x);n<-n-1L}
  followers[[k]]<-local_walk(f,x,lp,n,.6)$h
 }
 for(k in 1:3) {
  n<-alloc[k+3L]
  if(policy=='followers') {
   j<-ord[(k-1L)%%2L+1L];z<-local_walk(f,ends[j,],lps[j],n,.6)
  } else {
   nc<-min(12L,max(2L,n%/%20L))
   cand<-matrix(rnorm(nc*D,sd=4),nc,D);lp<-apply(cand,1,f$eval)
   j<-choose_candidate(policy,cand,lp,scout_pool,previous)
   previous<-rbind(previous,cand[j,])
   z<-local_walk(f,cand[j,],lp[j],n-nc,2.3)
  }
  pioneers[[k]]<-z$h
 }
 stopifnot(f$count()==bs+bm)
 q<-fit_q(c(scouts,followers,pioneers))
 j<-ord[1];fin<-final_chain(f,q,ends[j,],lps[j],B-f$count())
 stopifnot(f$count()==B)
 list(draws=fin$draws,accept=fin$accept,evals=f$count(),
      stages=list(scout=scouts,follower=followers,pioneer=pioneers,final=list(fin$draws)))
}
# Competitive matched global-mixture baseline: eight dispersed explorers,
# exploratory covariance adaptation, same mixture fit and final MH kernel.
# Not a substitute for a published, fully tuned AIMM implementation.
adaptive_mixture_mh <- function(B,x0) {
 f<-make_counter();explore<-floor(.65*B);nc<-8L
 starts<-matrix(rnorm(nc*D,sd=4),nc,D);starts[1,]<-x0
 n<-rep((explore-nc)%/%nc,nc)
 if((explore-nc)%%nc>0)n[seq_len((explore-nc)%%nc)]<-n[seq_len((explore-nc)%%nc)]+1L
 h<-vector('list',nc);ends<-matrix(NA_real_,nc,D);lp_end<-numeric(nc)
 for(k in seq_len(nc)) {
  lp<-f$eval(starts[k,]);n1<-floor(n[k]/2)
  z1<-local_walk(f,starts[k,],lp,n1,1.3)
  S<-if(nrow(unique(z1$h))>D+3L) cov(z1$h) else diag(4,D)
  step<-sqrt(max(.5,min(16,mean(diag(S)))))*2.38/sqrt(D)
  z2<-local_walk(f,z1$x,z1$lp,n[k]-n1,step)
  h[[k]]<-rbind(z1$h,z2$h[-1,,drop=FALSE]);ends[k,]<-z2$x;lp_end[k]<-z2$lp
 }
 stopifnot(f$count()==explore)
 q<-fit_q(h);j<-which.max(lp_end)
 fin<-final_chain(f,q,ends[j,],lp_end[j],B-f$count())
 stopifnot(f$count()==B)
 list(draws=fin$draws,accept=fin$accept,evals=f$count(),
      stages=list(exploration=h,final=list(fin$draws)))
}

# SAChES-inspired: DE ensemble stage then covariance-adaptive MH stage.
# IMPORTANT: not a reproduction of the published SAChES/DREAM/DRAM method.
# Unlike FLOCK's frozen mixture final kernel, this final stage uses an AM kernel.
saches_inspired <- function(B,x0) {
 f<-make_counter();train<-floor(.65*B);nc<-8L
 starts<-matrix(rnorm(nc*D,sd=4),nc,D);starts[1,]<-x0
 x<-starts;lp<-apply(x,1,f$eval);h<-lapply(1:nc,function(k)matrix(x[k,],1,D))
 for(t in seq_len(train-nc)) {
  k<-(t-1L)%%nc+1L;ab<-sample(setdiff(1:nc,k),2)
  prop<-x[k,]+(2.38/sqrt(2*D))*(x[ab[1],]-x[ab[2],])+rnorm(D,sd=.08)
  lprop<-f$eval(prop)
  if(log(runif(1))<lprop-lp[k]){x[k,]<-prop;lp[k]<-lprop}
  h[[k]]<-rbind(h[[k]],x[k,])
 }
 # Second phase: a single frozen covariance learned from ensemble (valid MH).
 pool<-do.call(rbind,h);L<-safe_chol((2.38^2/D)*cov(pool)+diag(.04,D))
 j<-which.max(lp);xx<-x[j,];ll<-lp[j];nf<-B-f$count()
 out<-matrix(NA_real_,nf+1L,D);out[1,]<-xx;acc<-0L
 for(t in seq_len(nf)) {
  yy<-rmv(xx,L);ly<-f$eval(yy)
  if(log(runif(1))<ly-ll){xx<-yy;ll<-ly;acc<-acc+1L}
  out[t+1L,]<-xx
 }
 burn<-floor(nrow(out)*.1);draws<-out[(burn+1L):nrow(out),,drop=FALSE]
 stopifnot(f$count()==B)
 list(draws=draws,accept=acc/nf,evals=B,
      stages=list(exploration=h,final=list(draws)))
}
# Independent RWMH, four chains. Final 90% per chain scored; total budget counted.
rwmh4 <- function(B,x0) {
 f<-make_counter();nc<-4L;ns<-rep(B%/%nc,nc)
 if(B%%nc>0)ns[seq_len(B%%nc)]<-ns[seq_len(B%%nc)]+1L
 h<-list();fin<-list();acc<-numeric(nc)
 for(k in seq_len(nc)) {
  x<-if(k==1)x0 else rnorm(D,sd=4);lp<-f$eval(x)
  z<-local_walk(f,x,lp,ns[k]-1L,.7)
  h[[k]]<-z$h;burn<-floor(nrow(z$h)*.1)
  fin[[k]]<-z$h[(burn+1L):nrow(z$h),,drop=FALSE]
 }
 stopifnot(f$count()==B)
 list(draws=do.call(rbind,fin),accept=NA_real_,evals=B,
      stages=list(exploration=h,final=fin))
}
# Label-invariant functionals are primary; raw parameter means are misleading
# in symmetric label-switching posteriors.
measure <- function(z,r,method) {
 est<-features(z$draws);err<-est-truth
 data.frame(rep=r,method=method,evals=z$evals,accept=z$accept,
  functional_rmse=sqrt(mean(err^2)),
  location_rmse=sqrt(mean(err[c('low_mean','high_mean','separation')]^2)),
  label_error=abs(est['label_probability']-.5),
  tail_error=abs(err['predictive_tail']),
  n_retained=nrow(z$draws))
}
stage_diagnostics <- function(z,r,method) {
 do.call(rbind,lapply(names(z$stages),function(stage) {
  h<-do.call(rbind,z$stages[[stage]]);lab<-apply(h,1,mode_label)
  data.frame(rep=r,method=method,stage=stage,mode=1:2,
   visited=vapply(1:2,function(j)any(lab==j),logical(1)),
   occupancy=vapply(1:2,function(j)mean(lab==j),numeric(1)))
 }))
}
B<-as.integer(Sys.getenv('BUDGET','20000'))
R<-as.integer(Sys.getenv('REPS','100'))
stopifnot(B>=3000,R>=2)
methods<-c('FLOCK_Followers','FLOCK_Novelty','FLOCK_Unexplored',
           'AdaptiveMixtureMH','SAChES_inspired','RWMH_4chains')
rows<-vector('list',R*length(methods));stages<-vector('list',length(rows));ix<-0L
for(r in seq_len(R)) {
 seed<-810000L+r;set.seed(seed);x0<-rnorm(D,sd=4)
 for(k in seq_along(methods)) {
  set.seed(seed+100000L);method<-methods[k]
  z<-switch(method,
    FLOCK_Followers=flock(B,x0,'followers'),
    FLOCK_Novelty=flock(B,x0,'novelty'),
    FLOCK_Unexplored=flock(B,x0,'unexplored'),
    AdaptiveMixtureMH=adaptive_mixture_mh(B,x0),
    SAChES_inspired=saches_inspired(B,x0),
    RWMH_4chains=rwmh4(B,x0))
  ix<-ix+1L;rows[[ix]]<-measure(z,r,method)
  stages[[ix]]<-stage_diagnostics(z,r,method)
 }
 if(r%%5L==0L||r==R)cat('Finished',r,'of',R,'replicates\n')
}
res<-do.call(rbind,rows);stage_res<-do.call(rbind,stages)
stopifnot(all(res$evals==B))
write.csv(res,'flock_experiment_4_results.csv',row.names=FALSE)
write.csv(stage_res,'flock_experiment_4_mode_discovery.csv',row.names=FALSE)
summ<-do.call(rbind,lapply(split(res,res$method),function(z)
 data.frame(method=z$method[1],n=nrow(z),mean_rmse=mean(z$functional_rmse),
 median_rmse=median(z$functional_rmse),mean_location_rmse=mean(z$location_rmse),
 mean_label_error=mean(z$label_error),mean_tail_error=mean(z$tail_error))))
write.csv(summ,'flock_experiment_4_summary.csv',row.names=FALSE)
paired<-do.call(rbind,lapply(setdiff(methods,'FLOCK_Unexplored'),function(other){
 do.call(rbind,lapply(c('functional_rmse','location_rmse','label_error'),function(metric){
  a<-res[res$method=='FLOCK_Unexplored',c('rep',metric)]
  b<-res[res$method==other,c('rep',metric)]
  zz<-merge(a,b,by='rep',suffixes=c('_flock','_other'))
  d<-zz[[paste0(metric,'_flock')]]-zz[[paste0(metric,'_other')]]
  set.seed(4567+nchar(other)+nchar(metric))
  boot<-replicate(1000,mean(sample(d,length(d),replace=TRUE)))
  data.frame(competitor=other,metric=metric,mean_diff=mean(d),
   ci_low=unname(quantile(boot,.025)),ci_high=unname(quantile(boot,.975)),
   fraction_flock_better=mean(d<0))
 }))
}))
write.csv(paired,'flock_experiment_4_paired.csv',row.names=FALSE)
disc<-aggregate(visited~method+stage+mode,stage_res,mean)
write.csv(disc,'flock_experiment_4_discovery_summary.csv',row.names=FALSE)
pdf('flock_experiment_4_diagnostics.pdf',width=12,height=8)
par(mfrow=c(2,2),mar=c(9,4,3,1))
boxplot(functional_rmse~method,res,las=2,main='Functional RMSE',ylab='RMSE')
boxplot(location_rmse~method,res,las=2,main='Label-invariant location RMSE',ylab='RMSE')
boxplot(label_error~method,res,las=2,main='Label probability absolute error',ylab='Error')
boxplot(tail_error~method,res,las=2,main='Posterior predictive tail error',ylab='Error')
dev.off()
print(summ,row.names=FALSE,digits=4)
cat('Saved Experiment 4 results, paired comparisons, mode discovery, and diagnostics.\n')