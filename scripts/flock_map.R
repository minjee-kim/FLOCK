# FLOCK map v0.1 — standalone base R implementation
# Source this file, then see /mnt/data/flock_map_demo.R for an example.
# Exploratory map only: NOT a posterior mass estimator or a sampler.

flock_map <- function(d, min_eigen = 1e-6) {
  stopifnot(length(d)==1L, d>=1, d==as.integer(d), min_eigen>0)
  structure(list(d=as.integer(d), min_eigen=min_eigen,
                 regions=list(), merges=list()), class="flock_map")
}

flock_region <- function(x, id, min_eigen=1e-6) {
  x <- as.matrix(x); storage.mode(x) <- "double"
  stopifnot(nrow(x)>=1L, ncol(x)>=1L, all(is.finite(x)),
            length(id)==1L, nzchar(as.character(id)))
  d <- ncol(x)
  mu <- colMeans(x)
  S <- if (nrow(x)>1L) stats::cov(x) else matrix(0,d,d)
  S <- matrix(S,d,d)
  list(id=as.character(id), n=nrow(x), mu=as.numeric(mu),
       covariance=S, observations=x, provenance=as.character(id))
}

flock_add <- function(map, x, id) {
  stopifnot(inherits(map,"flock_map"))
  x <- as.matrix(x)
  stopifnot(ncol(x)==map$d, !as.character(id) %in% names(map$regions))
  map$regions[[as.character(id)]] <- flock_region(x,id,map$min_eigen)
  map
}

flock_update <- function(map, id, x) {
  stopifnot(inherits(map,"flock_map"), as.character(id) %in% names(map$regions))
  old <- map$regions[[as.character(id)]]
  x <- as.matrix(x)
  stopifnot(ncol(x)==map$d)
  updated <- flock_region(rbind(old$observations,x),id,map$min_eigen)
  updated$provenance <- old$provenance
  map$regions[[as.character(id)]] <- updated
  map
}

flock_spd <- function(S, floor=1e-6) {
  S <- (S+t(S))/2
  eig <- eigen(S,symmetric=TRUE)
  eig$vectors %*% diag(pmax(eig$values,floor),nrow(S)) %*% t(eig$vectors)
}

flock_symmetric_kl <- function(a,b,floor=1e-6) {
  d <- length(a$mu)
  A <- flock_spd(a$covariance,floor)
  B <- flock_spd(b$covariance,floor)
  delta <- a$mu-b$mu
  # Symmetrized KL between N(mu_a,A) and N(mu_b,B).
  as.numeric(0.25*(sum(diag(solve(B,A)))+sum(diag(solve(A,B))) -
                    2*d + crossprod(delta,(solve(A,delta)+solve(B,delta)))))
}

flock_pair_metrics <- function(map,id1,id2) {
  stopifnot(inherits(map,"flock_map"), id1!=id2,
            id1 %in% names(map$regions),id2 %in% names(map$regions))
  a <- map$regions[[id1]]; b <- map$regions[[id2]]
  P <- flock_spd((a$covariance+b$covariance)/2,map$min_eigen)
  delta <- a$mu-b$mu
  list(mahalanobis=sqrt(as.numeric(crossprod(delta,solve(P,delta)))),
       symmetric_kl=flock_symmetric_kl(a,b,map$min_eigen))
}

# Optional density-valley test: evaluate a log-target along the center-to-center
# segment. A curved connection may fail this test; treat it as conservative.
# If supplied, log_target evaluations are NOT counted automatically.
flock_bridge_check <- function(map,id1,id2,log_target,n_points=9L,
                               max_drop=2) {
  stopifnot(is.function(log_target),n_points>=3L,max_drop>=0)
  a <- map$regions[[id1]]$mu; b <- map$regions[[id2]]$mu
  tt <- seq(0,1,length.out=n_points)
  vals <- vapply(tt,function(t) log_target((1-t)*a+t*b),numeric(1))
  if (!all(is.finite(vals))) return(list(pass=FALSE,drop=Inf,evals=length(vals)))
  # Compare interior log density with the LOWER endpoint, not higher one.
  drop <- min(vals[c(1L,length(vals))])-min(vals)
  list(pass=(drop<=max_drop),drop=drop,evals=length(vals))
}

flock_merge <- function(map,id1,id2,new_id=paste(id1,id2,sep="+"),
                        max_mahalanobis=1.5,max_symmetric_kl=1,
                        log_target=NULL,max_bridge_drop=2) {
  stopifnot(inherits(map,"flock_map"),id1!=id2,
            !new_id %in% setdiff(names(map$regions),c(id1,id2)))
  metrics <- flock_pair_metrics(map,id1,id2)
  bridge <- NULL
  pass <- metrics$mahalanobis<=max_mahalanobis &&
          metrics$symmetric_kl<=max_symmetric_kl
  if (pass && !is.null(log_target)) {
    bridge <- flock_bridge_check(map,id1,id2,log_target,
                                 max_drop=max_bridge_drop)
    pass <- bridge$pass
  }
  if (!pass) return(list(map=map,merged=FALSE,metrics=metrics,bridge=bridge))
  a <- map$regions[[id1]]; b <- map$regions[[id2]]
  new <- flock_region(rbind(a$observations,b$observations),new_id,map$min_eigen)
  new$provenance <- c(a$provenance,b$provenance)
  map$regions[[id1]] <- NULL
  map$regions[[id2]] <- NULL
  map$regions[[new_id]] <- new
  map$merges[[length(map$merges)+1L]] <- list(
    from=c(id1,id2),to=new_id,metrics=metrics,bridge=bridge)
  list(map=map,merged=TRUE,metrics=metrics,bridge=bridge)
}

flock_map_table <- function(map) {
  stopifnot(inherits(map,"flock_map"))
  if (!length(map$regions)) return(data.frame())
  do.call(rbind,lapply(map$regions,function(r)
    data.frame(region=r$id,n=r$n,
               center=paste(signif(r$mu,4),collapse=", "),
               provenance=paste(r$provenance,collapse="+"),
               stringsAsFactors=FALSE)))
}

flock_plot <- function(map,coords=c(1L,2L),main="FLOCK map (2D projection)") {
  stopifnot(inherits(map,"flock_map"),length(coords)==2L,
            all(coords %in% seq_len(map$d)),length(map$regions)>0L)
  xx <- do.call(rbind,lapply(map$regions,function(r) r$observations[,coords,drop=FALSE]))
  plot(xx[,1],xx[,2],pch=16,cex=.45,col="grey65",
       xlab=paste0("theta[",coords[1],"]"),ylab=paste0("theta[",coords[2],"]"),
       main=main)
  for (r in map$regions) {
    points(r$mu[coords[1]],r$mu[coords[2]],pch=19,cex=1.3)
    text(r$mu[coords[1]],r$mu[coords[2]],labels=r$id,pos=3)
    S <- flock_spd(r$covariance[coords,coords,drop=FALSE],map$min_eigen)
    eig <- eigen(S,symmetric=TRUE)
    ang <- seq(0,2*pi,length.out=100L)
    ring <- eig$vectors %*% diag(sqrt(eig$values),2) %*% rbind(cos(ang),sin(ang))
    lines(r$mu[coords[1]]+ring[1,],r$mu[coords[2]]+ring[2,])
  }
  invisible(map)
}