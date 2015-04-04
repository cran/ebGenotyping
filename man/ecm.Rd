\name{ecm}
\alias{ecm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
An initial model of genotyping using Next Generation sequencing Data
}
\description{
ecm implements the initial parameter estimation for the model described in 'Genotyping for Rare Variant Detection Using Next-generation Sequencing Data' for genotyping using Next-generation Sequencing Data. 
}
\usage{
ecm(mu0, delta0, eta0, pRR0, pRV0, dat, cvg, eps = 0.001, max.steps = 500)
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
a single value: the initial estimate of the probability of RR. 
}
  \item{pRV0}{
a single value: the initial estimate of the probability of RV. 
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
This function doesn't take optimum problems into account, so it is likely that the estimators are trapped in a local optimum. 
Another function named ECM in this package is suggested to cope with this problem and other real data issues.
}

