# FLOCK experiment 3: targeted vs random pioneers and adaptive mixture MH
# Run: source('flock_experiment_3.R')
# Set Sys.setenv(REPS='5', BUDGET='20000') for smoke test.
# Base R only. Ground truth used exclusively in diagnostics.
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
 starts<-matrix(rnorm(2*D,sd=7),2,D);starts[1,]<-x0
 scouts<-vector('list',2);ends<-matrix(NA_real_,2,D);lps<-numeric(2)
 ns<-rep((bs-2L)%/%2L,2);ns[1]<-ns[1]+(bs-2L-sum(ns))
 for(k in 1:2) {
  lp<-f$eval(starts[k,]);z<-local_walk(f,starts[k,],lp,ns[k],5)
  scouts[[k]]<-z$h;ends[k,]<-z$x;lps[k]<-z$lp
 }
 scout_pool<-do.call(rbind,scouts);ord<-order(lps,decreasing=TRUE)
 alloc<-rep(bm%/%6L,6L);if(bm%%6L>0) alloc[seq_len(bm%%6L)]<-alloc[seq_len(bm%%6L)]+1L
 followers<-list();pioneers<-list();previous<-matrix(numeric(0),0,D)
 for(k in 1:3) {
  j<-ord[(k-1L)%%2L+1L];x<-ends[j,];lp<-lps[j];n<-alloc[k]
  if(k>2){x<-x+rnorm(D,sd=.7);lp<-f$eval(x);n<-n-1L}
  followers[[k]]<-local_walk(f,x,lp,n,1.15)$h
 }
 for(k in 1:3) {
  n<-alloc[k+3L]
  if(policy=='followers') {
   j<-ord[(k-1L)%%2L+1L];z<-local_walk(f,ends[j,],lps[j],n,1.15)
  } else {
   nc<-min(12L,max(2L,n%/%20L))
   cand<-matrix(rnorm(nc*D,sd=8),nc,D);lp<-apply(cand,1,f$eval)
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
 starts<-matrix(rnorm(nc*D,sd=8),nc,D);starts[1,]<-x0
 n<-rep((explore-nc)%/%nc,nc)
 if((explore-nc)%%nc>0)n[seq_len((explore-nc)%%nc)]<-n[seq_len((explore-nc)%%nc)]+1L
 h<-vector('list',nc);ends<-matrix(NA_real_,nc,D);lp_end<-numeric(nc)
 for(k in seq_len(nc)) {
  lp<-f$eval(starts[k,]);n1<-floor(n[k]/2)
  z1<-local_walk(f,starts[k,],lp,n1,3)
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
# For retrospective diagnostics only: classify each state by dominant true
# mixture responsibility (>0.5) and within Mahalanobis radius 4.
mode_label <- function(x) {
 rr<-responsibility(x);j<-which.max(rr)
 d<-sqrt(sum(((x-means[j,])/scales[j])^2))
 if(rr[j]>.5 && d<=4) j else 0L
}
measure <- function(z,r,method) {
 rr<-vapply(seq_len(nrow(z$draws)),function(i) responsibility(z$draws[i,]),numeric(4))
 w<-rowMeans(rr);mu<-colMeans(z$draws)
 data.frame(rep=r,method=method,evals=z$evals,accept=z$accept,
  mean_rmse=sqrt(mean((mu-truth_mean)^2)),weight_l1=sum(abs(w-weights)),
  rare_est=w[4],rare_abs_error=abs(w[4]-.05),rare_detected=as.integer(w[4]>.005))
}
stage_diagnostics <- function(z,r,method) {
 do.call(rbind,lapply(names(z$stages),function(stage) {
  h<-do.call(rbind,z$stages[[stage]])
  lab<-apply(h,1,mode_label)
  data.frame(rep=r,method=method,stage=stage,mode=1:4,
             visited=vapply(1:4,function(j) any(lab==j),logical(1)),
             count=vapply(1:4,function(j) sum(lab==j),integer(1)),states=nrow(h))
 }))
}
B<-as.integer(Sys.getenv('BUDGET','20000'))
R<-as.integer(Sys.getenv('REPS','100'))
stopifnot(B>=3000,R>=2)
methods<-c('FLOCK_Unexplored','FLOCK_Random','FLOCK_Density',
           'FLOCK_Diversity','FLOCK_Novelty','FLOCK_Followers','AdaptiveMixtureMH')
rows<-vector('list',R*length(methods));stage_rows<-vector('list',length(rows));idx<-0L
for(r in seq_len(R)) {
 seed<-900000L+r;set.seed(seed);x0<-rnorm(D,sd=7)
 for(k in seq_along(methods)) {
  set.seed(seed+100000L);method<-methods[k]
  z<-switch(method,
   FLOCK_Unexplored=flock(B,x0,'unexplored'),
   FLOCK_Random=flock(B,x0,'random'),
   FLOCK_Density=flock(B,x0,'density'),
   FLOCK_Diversity=flock(B,x0,'diversity'),
   FLOCK_Novelty=flock(B,x0,'novelty'),
   FLOCK_Followers=flock(B,x0,'followers'),
   AdaptiveMixtureMH=adaptive_mixture_mh(B,x0))
  idx<-idx+1L;rows[[idx]]<-measure(z,r,method)
  stage_rows[[idx]]<-stage_diagnostics(z,r,method)
 }
 if(r%%5L==0L || r==R) cat('Finished',r,'of',R,'replicates\n')
}
res<-do.call(rbind,rows);stage_res<-do.call(rbind,stage_rows)
stopifnot(all(res$evals==B))
write.csv(res,'flock_experiment_3_results.csv',row.names=FALSE)
write.csv(stage_res,'flock_experiment_3_mode_discovery.csv',row.names=FALSE)
summ<-do.call(rbind,lapply(split(res,res$method),function(z) {
 data.frame(method=z$method[1],n=nrow(z),mean_rmse=mean(z$mean_rmse),
  median_rmse=median(z$mean_rmse),sd_rmse=sd(z$mean_rmse),
  mean_weight_l1=mean(z$weight_l1),mean_rare_est=mean(z$rare_est),
  rare_detection_rate=mean(z$rare_detected),mean_accept=mean(z$accept))
}))
write.csv(summ,'flock_experiment_3_summary.csv',row.names=FALSE)
# Paired bootstrap: differences are FLOCK_Unexplored minus comparator;
# negative values favor FLOCK_Unexplored. Whole replicate pairs resampled.
paired<-do.call(rbind,lapply(setdiff(methods,'FLOCK_Unexplored'),function(other){
 do.call(rbind,lapply(c('mean_rmse','weight_l1','rare_abs_error'),function(metric){
  a<-res[res$method=='FLOCK_Unexplored',c('rep',metric)]
  b<-res[res$method==other,c('rep',metric)]
  zz<-merge(a,b,by='rep',suffixes=c('_flock','_other'))
  d<-zz[[paste0(metric,'_flock')]]-zz[[paste0(metric,'_other')]]
  set.seed(12345+length(d)+nchar(other)+nchar(metric))
  bs<-replicate(1000L,mean(sample(d,length(d),replace=TRUE)))
  data.frame(competitor=other,metric=metric,n=length(d),mean_diff=mean(d),
   ci_low=unname(quantile(bs,.025)),ci_high=unname(quantile(bs,.975)),
   fraction_flock_better=mean(d<0))
 }))
}))
write.csv(paired,'flock_experiment_3_paired.csv',row.names=FALSE)
# Stage discovery frequency: fraction of replicates with >=1 state in mode.
disc<-aggregate(visited~method+stage+mode,data=stage_res,FUN=mean)
names(disc)[4]<-'discovery_rate'
write.csv(disc,'flock_experiment_3_discovery_summary.csv',row.names=FALSE)
pdf('flock_experiment_3_diagnostics.pdf',width=12,height=9)
par(mfrow=c(2,2),mar=c(10,4,3,1))
boxplot(mean_rmse~method,res,las=2,main='Posterior mean RMSE',ylab='RMSE')
boxplot(weight_l1~method,res,las=2,main='Mixture weight L1 error',ylab='L1 error')
boxplot(rare_est~method,res,las=2,main='Rare component weight',ylab='Estimated weight');abline(h=.05,col='red',lty=2)
plot(1,1,type='n',xlim=c(1,4),ylim=c(0,1),xaxt='n',xlab='True mixture component',ylab='Discovery rate',main='Final-stage mode discovery')
axis(1,1:4)
cols<-seq_along(methods)
for(k in seq_along(methods)) {
 v<-disc[disc$method==methods[k]&disc$stage=='final',]
 lines(v$mode,v$discovery_rate,type='b',col=cols[k],pch=k)
}
legend('bottomleft',legend=methods,col=cols,lty=1,pch=seq_along(methods),cex=.65)
dev.off()
print(summ,row.names=FALSE,digits=4)
cat('\nSaved results, summary, paired comparisons, mode discovery, and diagnostics PDF.\n')