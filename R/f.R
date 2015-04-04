logit <- function(x) return(log(x/(1-x)))
rlogit <- function(x) return(ifelse(x>100,1,exp(x)/(1+exp(x))))

f.eta <- function(x,mu,delta,zm1,z0,zp1,dat,cvg){
  out <- outer(delta,mu,"+")
  out <- x*sum((zm1-zp1)*dat)+sum(cvg*(zm1*log(1+exp(out-x))+zp1*log(1+exp(out+x))))
  return(out)
}

f.mu.1 <- function(x,delta,eta,dat,cvg,zm1,z0,zp1){
  xdelta <- outer(delta,x,"+")
  out <- apply(cvg*(zm1*rlogit(xdelta-eta)+z0*rlogit(xdelta)+zp1*rlogit(xdelta+eta)),2,sum)-apply(dat,2,sum)
  return(out)
}

f.mu.2 <- function(x,delta,eta,dat,cvg,zm1,z0,zp1){
  xdelta <- outer(delta,x,"+")
  out <- apply(cvg*(zm1*rlogit(xdelta-eta)/(1+exp(xdelta-eta))+
                      z0*rlogit(xdelta)/(1+exp(xdelta))+
                      zp1*rlogit(xdelta+eta)/(1+exp(xdelta+eta))),2,sum)
  return(out)
}
####first and second partial derivatives of delta
f.delta.1 <- function(x,mu,eta,dat,cvg,zm1,z0,zp1){
  mux <- outer(x,mu,"+")
  out <- apply(cvg*(zm1*rlogit(mux-eta)+z0*rlogit(mux)
                    +zp1*rlogit(mux+eta)),1,sum)-apply(dat,1,sum)
  return(out)
}
f.delta.2 <- function(x,mu,eta,dat,cvg,zm1,z0,zp1){
  mux <- outer(x,mu,"+")
  out <- apply(cvg*(zm1*rlogit(mux-eta)/(1+exp(mux-eta))+
                      z0*rlogit(mux)/(1+exp(mux))+
                      zp1*rlogit(mux+eta)/(1+exp(mux+eta))),1,sum)
  return(out)
}
####first and second partial derivatives of eta
f.eta.1 <- function(x,mu,delta,dat,cvg,zm1,z0,zp1){
  out <- outer(delta,mu,"+")
  out <- sum((zm1-zp1)*dat)+sum(cvg*(-zm1*rlogit(out-x)+zp1*rlogit(out+x)))
  return(out)
}

f.eta.2 <- function(x,mu,delta,dat,cvg,zm1,z0,zp1){
  out <- outer(delta,mu,"+")
  out <- sum(cvg*(zm1*rlogit(out-x)/(1+exp(out-x))+zp1*rlogit(out+x)/(1+exp(out+x))))
  return(out)
}

z.adj.1 <- function(zm1,z0,zp1,pm1,p0,pp1,th1=0.5,th2=NULL,adj=FALSE,adj.idx=NULL){
  m <- dim(zm1)
  n <- m[2]
  m <- m[1]
  if (is.null(th2))
    th2 <- 1-th1
  p <- cbind(as.vector(pm1),as.vector(p0),as.vector(pp1))
  z <- cbind(as.vector(zm1),as.vector(z0),as.vector(zp1))
  N <- m*n
  ind <- matrix(0,N,3)
  ind[which(p<=th1)] <- 1  
  zm1.adj <- apply(z*ind,1,sum)
  ind <- matrix(0,N,3)
  ind[which((p>th1)&(p<=th2))] <- 1
  z0.adj <- apply(z*ind,1,sum)
  ind <- matrix(0,N,3)
  ind[which(p>th2)] <- 1
  zp1.adj <- apply(z*ind,1,sum)
  adj <- apply(abs(cbind(zm1.adj,z0.adj,zp1.adj)-z),1,sum)>0
  if (any(adj)){
    adj.idx <- sort(unique(ceiling(which(adj)/m)))
    adj <- TRUE
  }
  else{
    adj.idx <- NULL
    adj <- FALSE
  }
  out <- list(z.RR=matrix(zm1.adj,m,n),z.RV=matrix(z0.adj,m,n),z.VV=matrix(zp1.adj,m,n),adj=adj,adj.idx=adj.idx)
  return(out)
}

mu.recal <- function(dat,cvg,zm1,z0,zp1,delta,eta){
  f <- function(x){
    sum(dat)-sum(cvg*(zm1*exp(x+delta-eta)/(1+exp(x+delta-eta))+z0*exp(x+delta)/(1+exp(x+delta))+zp1*exp(x+delta+eta)/(1+exp(x+delta+eta))))
  }
  if (f(-10)*f(10)>0){
    out <- 0
  }
  else{
    out <- uniroot(f,c(-10,10))$root
  }
  return(out)
}

delta.recal <- function(dat,cvg,zm1,z0,zp1,mu,eta){
  f <- function(x){
    sum(dat)-sum(cvg*(zm1*exp(mu+x-eta)/(1+exp(mu+x-eta))+z0*exp(mu+x)/(1+exp(mu+x))+zp1*exp(mu+x+eta)/(1+exp(mu+x+eta))))
  }
  out <- uniroot(f,c(-10,10))$root
  return(out)
}

ecm <- function(mu0,delta0,eta0,pRR0,pRV0,dat,cvg,eps=1e-3,max.steps=500){
  # start from E step given initials of M step
  if ((any(is.na(cvg)))|(any(is.na(cvg)))){
    stop("There is NA in the data, please check")
  }
  m <- dim(dat)
  n <- m[1]
  m <- m[2]
  theta0 <- c(mu0,delta0,eta0)
  pVV0 <- 1-pRR0-pRV0
  a <- outer(delta0,mu0,"+")
  a <- ifelse(a>100,100,a)
  zm1 <- exp((-eta0)*dat-cvg*log((1+exp(a-eta0))/(1+exp(a))))*pRR0
  zp1 <- exp(eta0*dat-cvg*log((1+exp(a+eta0))/(1+exp(a))))*pVV0
  z0 <- pRV0
  inf.index <- which(is.infinite(zm1))
  if (length(inf.index)>0){
    zm1[inf.index] <- 1e+100*pRR0
  }
  inf.index <- which(is.infinite(zp1))
  if (length(inf.index)>0){
    zp1[inf.index] <- 1e+100*pVV0
  }  
  zsum <- zm1+z0+zp1
  zm1 <- zm1/zsum
  z0 <- z0/zsum
  zp1 <- zp1/zsum
  # M step
  mu1 <- mu0-f.mu.1(mu0,delta0,eta0,dat,cvg,zm1,z0,zp1)/f.mu.2(mu0,delta0,eta0,dat,cvg,zm1,z0,zp1)
  delta1 <- delta0-f.delta.1(delta0,mu1,eta0,dat,cvg,zm1,z0,zp1)/f.delta.2(delta0,mu1,eta0,dat,cvg,zm1,z0,zp1)
  eta1 <- eta0-f.eta.1(eta0,mu1,delta1,dat,cvg,zm1,z0,zp1)/f.eta.2(eta0,mu1,delta1,dat,cvg,zm1,z0,zp1)
  pRR1 <- sum(zm1)/(m*n)
  pVV1 <- sum(zp1)/(m*n)
  pRV1 <- 1-pRR1-pVV1
  theta1 <- c(mu1,delta1,eta1)
  s <- 1
  while ((sum((theta1-theta0)^2)>eps)&(s<=max.steps)){
    mu0 <- mu1
    delta0 <- delta1
    eta0 <- eta1
    pRR0 <- pRR1
    pRV0 <- pRV1
    pVV0 <- pVV1
    theta0 <- theta1
    a <- outer(delta0,mu0,"+")
    a <- ifelse(a>100,100,a)
    zm1 <- exp((-eta0)*dat-cvg*log((1+exp(a-eta0))/(1+exp(a))))*pRR0
    zp1 <- exp(eta0*dat-cvg*log((1+exp(a+eta0))/(1+exp(a))))*pVV0
    z0 <- pRV0
    inf.index <- which(is.infinite(zm1))
    if (length(inf.index)>0){
      zm1[inf.index] <- 1e+100*pRR0
    }
    inf.index <- which(is.infinite(zp1))
    if (length(inf.index)>0){
      zp1[inf.index] <- 1e+100*pVV0
    }  
    zsum <- zm1+z0+zp1
    zm1 <- zm1/zsum
    z0 <- z0/zsum
    zp1 <- zp1/zsum
    a <- f.mu.1(mu0,delta0,eta0,dat,cvg,zm1,z0,zp1)
    mu1 <- ifelse(a<1e-3,mu0,mu0-f.mu.1(mu0,delta0,eta0,dat,cvg,zm1,z0,zp1)/f.mu.2(mu0,delta0,eta0,dat,cvg,zm1,z0,zp1))
    a <- f.delta.1(delta0,mu1,eta0,dat,cvg,zm1,z0,zp1)
    delta1 <- ifelse(a<1e-3,delta0,delta0-f.delta.1(delta0,mu1,eta0,dat,cvg,zm1,z0,zp1)/f.delta.2(delta0,mu1,eta0,dat,cvg,zm1,z0,zp1))
    a <- f.eta.1(eta0,mu1,delta1,dat,cvg,zm1,z0,zp1)
    eta1 <- ifelse(eta0,eta0-f.eta.1(eta0,mu1,delta1,dat,cvg,zm1,z0,zp1)/f.eta.2(eta0,mu1,delta1,dat,cvg,zm1,z0,zp1))
    pRR1 <- sum(zm1)/(m*n)
    pVV1 <- sum(zp1)/(m*n)
    pRV1 <- 1-pRR1-pVV1
    theta1 <- c(mu1,delta1,eta1)
    s <- s+1
    print(s)
  }
  a <- (zm1<z0)
  if (all(!a))
    z.est <- (zm1<zp1)*2
  else{
    z.est <- a
    z.est[which(a)] <- (z0<zp1)[which(a)]+1
    z.est[-which(a)] <- (zm1<zp1)[-which(a)]*2
  }
  para.est <- list(mu=mu1,delta=delta1,eta=eta1,probs=c(pRR1,pRV1,pVV1))
  post.probs <- list(z.RR=zm1,z.RV=z0,z.VV=zp1)
  return(list(para.est=para.est,
              post.probs=post.probs,
              geno.est=z.est,
              em.steps=s))
}



loglike <- function(mu,delta,eta,pRR,pRV,dat,cvg){
  pVV <- 1-pRR-pRV
  out <- outer(delta,mu,"+")
  z0 <- out*dat-cvg*log(1+exp(out))# 
  zm1 <- exp((out-eta)*dat-cvg*log(1+exp(out-eta))-z0)*pRR
  zp1 <- exp((out+eta)*dat-cvg*log(1+exp(out+eta))-z0)*pVV
  inf.index <- which(is.infinite(zm1))
  if (length(inf.index)>0){
    zm1[inf.index] <- 1e+100*pRR
  }
  inf.index <- which(is.infinite(zp1))
  if (length(inf.index)>0){
    zp1[inf.index] <- 1e+100*pVV
  } 
  out <- log(zm1+pRV+zp1)+z0
  ##
  mu.like <- apply(out,2,sum)
  delta.like <- apply(out,1,sum)
  return(list(like=sum(out),mu.like=mu.like,delta.like=delta.like))
}


ECM <- function(mu0,delta0,eta0,pRR0,pRV0,dat,cvg,eps=1e-3,max.steps=500,para.by=0.1,check.range=10){
  m <- dim(dat)
  n <- m[1]
  m <- m[2]
  res <- ecm(mu0,delta0,eta0,pRR0,pRV0,dat,cvg,eps=eps,max.steps=max.steps)
  muhat0 <- res$para.est$mu
  muhat0 <- sapply(muhat0,function(x) seq(x-check.range,x+check.range,by=para.by))
  k <- apply(apply(muhat0,1,function(x) loglike(x,res$para.est$delta,res$para.est$eta,res$para.est$probs[1],res$para.est$probs[2],dat,cvg)$mu.like),1,which.max)
  muhat1 <- muhat0[k+(seq(k)-1)*dim(muhat0)[1]]
  if (any(abs(muhat1-res$para.est$mu)>=para.by)){
    res <- ecm(muhat1,res$para.est$delta,res$para.est$eta,res$para.est$probs[1],res$para.est$probs[2],dat,cvg,eps=eps,max.steps=max.steps)
  }
  delta0 <- res$para.est$delta
  delta0 <- sapply(delta0,function(x) seq(x-check.range,x+check.range,by=para.by))
  k <- apply(apply(delta0,1,function(x) loglike(res$para.est$mu,x,res$para.est$eta,res$para.est$probs[1],res$para.est$probs[2],dat,cvg)$delta.like),1,which.max)
  delta1 <- delta0[k+(seq(k)-1)*dim(delta0)[1]]
  if (any(abs(delta1-res$para.est$delta)>=para.by)){
    res <- ecm(res$para.est$mu,delta1,res$para.est$eta,res$para.est$probs[1],res$para.est$probs[2],dat,cvg,eps=eps,max.steps=max.steps)
  }  
  pest <- outer(res$para.est$delta,res$para.est$mu,"+")
  pm1 <- rlogit(pest-res$para.est$eta)
  p0 <- rlogit(pest)
  pp1 <- rlogit(pest+res$para.est$eta)
  dens <- density(c(pm1,p0,pp1))
  x <- dens$x
  y <- dens$y
  a <- which((x>0)&(x<1))
  x <- x[a]
  y <- y[a]
  th1 <- x[which.min(y[which(x<0.5)])]
  th2 <- x[sum(x<0.5)+which.min(y[which(x>=0.5)])]
  zadj <- z.adj.1(res$post.probs$z.RR,res$post.probs$z.RV,res$post.probs$z.VV,pm1,p0,pp1,th1,th2)
  if (zadj$adj){
    adj.idx <- zadj$adj.idx
    zm1 <- zadj[[1]]
    z0 <- zadj[[2]]
    zp1 <- zadj[[3]]
    res$para.est$mu[adj.idx] <- sapply(adj.idx,function(i) mu.recal(dat[,i],cvg[,i],zm1[,i],z0[,i],zp1[,i],res$para.est$delta,res$para.est$eta))
    res$para.est$delta <- sapply(1:dim(dat)[1],function(j) delta.recal(dat[j,],cvg[j,],zm1[j,],z0[j,],zp1[j,],res$para.est$mu,res$para.est$eta))
    res$para.est$eta <- optimize(f.eta,c(2,20),res$para.est$mu,res$para.est$delta,zm1=zm1,z0=z0,zp1=zp1,dat,cvg)$min
    pRR.est <- sum(zm1)/(m*n)
    pVV.est <- sum(zp1)/(m*n)
    pRV.est <- 1-pRR.est-pVV.est
    res$para.est$probs <- c(pRR.est,pRV.est,pVV.est)
    a <- outer(res$para.est$delta,res$para.est$mu,"+")
    zm1 <- exp((-res$para.est$eta)*dat-cvg*log((1+exp(a-res$para.est$eta))/(1+exp(a))))*pRR.est
    zp1 <- exp(res$para.est$eta*dat-cvg*log((1+exp(a+res$para.est$eta))/(1+exp(a))))*pVV.est
    z0 <- pRV.est
    inf.index <- which(is.infinite(zm1))
    if (length(inf.index)>0){
      zm1[inf.index] <- 1e+100*pRR.est
    }
    inf.index <- which(is.infinite(zp1))
    if (length(inf.index)>0){
      zp1[inf.index] <- 1e+100*pVV.est
    }  
    zsum <- zm1+z0+zp1
    zm1 <- zm1/zsum
    z0 <- z0/zsum
    zp1 <- zp1/zsum
    a <- (zm1<z0)
    if (all(!a))
      z.est <- (zm1<zp1)*2
    else{
      z.est <- a
      z.est[which(a)] <- (z0<zp1)[which(a)]+1
      z.est[-which(a)] <- (zm1<zp1)[-which(a)]*2
    }
    res$post.probs <- list(z.RR=zm1,z.RV=z0,z.VV=zp1)
    res$geno.est <- z.est
  }
  return(res)
}