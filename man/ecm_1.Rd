\name{ECM}
\alias{ECM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Genotyping using Next Generation sequencing Data
}
\description{
ECM is the most important function in this package. ECM implements the parameter estimation for the model described in 'Genotyping for Rare Variant Detection Using Next-generation Sequencing Data' for genotyping using Next-generation Sequencing Data.
ECM runs ecm for initial estimation first and make some adjustment to deal with local optimum problem and other real data issues.
}
\usage{
ECM(mu0,delta0,eta0,pRR0,pRV0,dat,cvg,eps,max.steps,para.by,check.range)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mu0}{
a vetor of the same length as number of positions: the initial value of mu.
}
  \item{delta0}{
a vetor of the same length as number of samples: the initial value of delta.
}
  \item{eta0}{
a single value: the initial value of eta. 
}
  \item{pRR0}{
a single value: the initial estimate of the prior probability of RR. 
}
  \item{pRV0}{
a single value: the initial estimate of the prior probability of RV. 
}
  \item{dat}{
a n*m matrix: the ith row, jth column of the matrix represents the non-reference counts of ith sample at jth position. 
}
  \item{cvg}{
a n*m matrix: the ith row, jth column of the matrix represents the depth of ith sample at jth position. 
}
  \item{eps}{
a single value: a threshold to control the convergence criterion. The default is 0.001.
}
  \item{max.steps}{
a single value: the maximum steps to run iterative algorithm to estimate parameters. The default is 500.(Adjustment is needed according to the number of parameters to estimate and the initial value of them) 
}
  \item{para.by}{
a single value: the step size to try to run out of local optimum. The default is 0.1
}
  \item{check.range}{
a single value: the range to search for global optimum. The default is 10.
}
}
\value{
 \item{para.est}{the estimate of position effect, mu, sample effect, delta, and the probability of 3 genotypes, RR, RV and VV}
 \item{post.probs}{3 matrix: the estimate of the posterior probabilities of 3 genotypes for n samples at m positions}
 \item{geno.est}{the estimated genotypes(RR:0; RV:1; VV:2) of n samples at m positions}
 \item{em.steps}{the total steps used to run iterative algorithm}
}
\references{
Na You and Gongyi Huang.(2015) Genotyping for Rare Variant Detection Using Next-generation Sequencing Data.
}
\author{
Na You <youn@mail.sysu.edu.cn> and Gongyi Huang<huanggyi@mail2.sysu.edu.cn>
}
\note{
Compared with ecm, which only implements an ECM algorithm to estimate parameters, this function takes other practical problems into account in the meanwhile, such as local optimum
and the real data problem.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
# generate simulation data
set.seed(2015)
m <- 100
m0 <- m*0.95
m1 <- m-m0
n.ctrl <- 10
n.case <- 20
n <- n.ctrl+n.case

##
M <- 4
lp <- matrix(-M,n,m)
Q <- 0.7
z <- rbind(matrix(0,n.ctrl,m),
           cbind(matrix(0,n.case,m0),matrix(rbinom(n.case*m1,1,Q),n.case,m1)))
b <- which(z==1)
R <- 0.8 # proportion of homozygous SNP######################
w <- rbinom(length(which(z==1)),1,R)
################ z are genotypes
z[b[which(w==0)]] <- 1
z[b[which(w==1)]] <- 2
lp[which(z==1)] <- 0
lp[which(z==2)] <- M
###########33
mu <- rnorm(m,0,0.1)
delta <- rnorm(n,0,0.1)
lp <- lp+outer(delta,mu,"+")
p <- rlogit(lp)
cvg <- matrix(rbinom(m*n,100,0.5),n,m)
dat <- matrix(sapply(1:(m*n),function(i) rbinom(1,cvg[i],p[i])),n,m)
######################parameters estimation
###(initial values(muhat,deltahat,etahat,pRR,pRV) for ecm input)
phat <- dat/cvg
phat[which(cvg==0)] <- NA
lphat <- logit(phat)
lphat[which(phat==0)] <- NA
#call initial genotype
eps <- 1e-6

geno <- matrix(rep(NA,m*n),n,m)
geno[which(phat>1-eps)] <- 2
geno[which(phat<eps)] <- 0
geno[which(lphat<(-0.27))] <- 0 ### 
geno[which(lphat>0.27)] <- 2
geno[which((lphat>=(-0.27))&(lphat<=0.27))] <- 1


##we use genotype information to adjust phat:adj.phat and adj2.phat
adj.phat <- phat
geno.1 <- which(geno==1)
geno.2 <- which(geno==2)
adj.phat[geno.1] <- abs(phat[geno.1]-0.5)
adj.phat[geno.2] <- 1-phat[geno.2]
adj.phat[which(adj.phat==0)] <- NA
adj.lphat<-logit(adj.phat)
#####parameters estimation
etahat <- -median(adj.lphat,na.rm=TRUE)
#etahat <- -min(adj.lphat)
na.delta <- which(apply(adj.lphat,1, function(x) all(is.na(x))))
na.mu <- which(apply(adj.lphat,2, function(x) all(is.na(x))))
deltahat <- apply(adj.lphat+etahat,1,function(x) median(x,na.rm=TRUE))
print(length(which(is.na(deltahat))))
deltahat[which(is.na(deltahat))]=0#no nas
nosample <- adj.lphat+etahat-matrix(rep(deltahat,m),n,m)
muhat <- apply(nosample,2,function(x) median(x,na.rm=TRUE))
print(length(which(is.na(muhat))))
muhat[which(is.na(muhat))]=0#no nas
print(length(na.mu))
print(length(na.delta))
pRR <- mean(geno==0,na.rm=TRUE)
pRV <- mean(geno==1,na.rm=TRUE)
pVV <- mean(geno==2,na.rm=TRUE)
################ECM 
T1 <- proc.time()
res1 <- ECM(muhat,deltahat,etahat,pRR,pRV,dat,cvg,eps=1e-3,max.steps=500,check.range=10)
proc.time()-T1
}

