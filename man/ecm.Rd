\name{ecm}
\alias{ecm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Genotyping and SNP Detection using Next Generation Sequencing Data
}
\description{
This function implements the method described in 'An Empirical Bayes Method for Genotyping and SNP detection Using Multi-sample Next-generation Sequencing Data'.

}
\usage{
ecm(dat,cvg,eps=1e-6,max.steps=500,eps.bisec=1e-6,ini.m=-7,ini.d=-7)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{
a n*m matrix: the ith row, jth column of the matrix represents the non-reference counts of ith sample at jth position. 
}
  \item{cvg}{
a n*m matrix: the ith row, jth column of the matrix represents the depth of ith sample at jth position. 
}
  \item{eps}{
a single value: a threshold to control the convergence criterion. The default is 1e-06.
}
  \item{max.steps}{
a single value: the maximum steps to run iterative algorithm to estimate parameters. The default is 500.(Adjustment is needed according to the number of parameters to estimate and the initial value of them.) 
}
  \item{eps.bisec}{
a single value: a threshold to control the convergence criterion of bisection criterion. The default is 1e-06.
}
  \item{ini.m}{
the initial value of each element of mu. We suggest users to use default -7.
}
  \item{ini.d}{
the initial value of each element of delta. We suggest users to use default -7.
}
}
\details{
This function implements the method described in 'An Empirical Bayes Method for Genotyping and SNP detection Using Multi-sample Next-generation Sequencing Data'. 
According to the paper, users can do genotyping with the estimated genotypes("geno.est"), 
and do SNP detection with the posterior probabilities of RR("post.probs$zRR"), based on the non-reference counts("dat")and depth("cvg").
}
\value{
\item{par.est}{a list including the estimate of position effect(mu), sample effect(delta), and the probability of RR and RV.}
\item{post.probs}{3 matrix: the estimate of the posterior probabilities of 3 genotypes for n samples at m positions.}
\item{steps}{the total steps to run iterative algorithm.}
\item{geno.est}{a n*m matrix: the estimated genotypes(0 for RR, 1 for RV and 2 for VV) of n samples at m positions.}
}
\references{
Na You and Gongyi Huang.(2016) An Empirical Bayes Method for Genotyping and SNP detection Using Multi-sample Next-generation Sequencing Data.
}
\author{
Na You <youn@mail.sysu.edu.cn> and Gongyi Huang<53hgy@163.com>
}
\examples{
#---------generate simulation data-----------
#start:generate simulation data#
set.seed(2016)
m <- 100
m0 <- m*0.95
m1 <- m-m0
n <- 30
Q <- 0.8
z <- cbind(matrix(0,n,m0),matrix(rbinom(n*m1,1,Q),n,m1))
b <- which(z==1)
R <- 0.8 # proportion of homozygous SNP
w <- rbinom(length(which(z==1)),1,R)
# z are genotypes
z[b[which(w==0)]] <- 1
z[b[which(w==1)]] <- 2
mu <- rep(-3,m)# stands for no effect
delta <- rep(-3,n)# stands for no effect
er.p <- -abs(outer(delta,mu,"+"))
p <- rlogit(er.p)
p[which(z==1)] <- 1/2
p[which(z==2)] <- 1-p[which(z==2)]
cvg <- matrix(rbinom(m*n,50,0.5),n,m)
dat <- matrix(sapply(1:(m*n),function(i) rbinom(1,cvg[i],p[i])),n,m)
#end:generate simulation data-#
#-----genotyping and SNP detection----------------
res <- ecm(dat=dat,cvg=cvg)
mean(z!=res$geno.est)#genotyping error
#----------call SNP---------
#start:call snp#
# define a function to calculate power, typeI error and FDR.
cutsnp <- function(fdr,alpha,true){
wh <- (true!=0)
tp <- sum((wh)&(fdr<alpha));
tn <- sum((!wh)&(fdr>=alpha));
fp <- sum((!wh)&(fdr<alpha));
fn <- sum((wh)&(fdr>=alpha));
pw  <- tp/(tp+fn);
t1 <- fp/(fp+tn);
fdr <- fp/(fp+tp);
return(c(TP=tp,TN=tn,FP=fp,FN=fn,power=pw,typeI=t1,FDR=fdr))
}
cutsnp(fdr=res$post.probs$zRR,alpha=0.001,true=z)
cutsnp(fdr=res$post.probs$zRR,alpha=0.01,true=z)
cutsnp(fdr=res$post.probs$zRR,alpha=0.05,true=z)
#end:call snp#
}

