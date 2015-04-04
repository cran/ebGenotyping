\name{loglike}
\alias{loglike}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
The observed data log-likelihood
}
\description{
loglike is the observed data log-likelihood
}
\usage{
loglike(mu, delta, eta, pRR, pRV, dat, cvg)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mu}{
a vetor of the same length as number of positions: the current estimate of mu.
}
  \item{delta}{
a vetor of the same length as number of samples: the current estimate of delta.
}
  \item{eta}{
a single value: the current estimate of eta 
}
  \item{pRR}{
the prior probability of RR
}
  \item{pRV}{
the prior probability of RV
}
  \item{dat}{
a n*m matrix: the ith row, jth column of the matrix represents the non-reference counts of ith sample at jth position. 
}
  \item{cvg}{
a n*m matrix: the ith row, jth column of the matrix represents the depth of ith sample at jth position. 
}
}
\value{
 \item{like}{a single value: the observed data log-likelihood}
 \item{mu.like}{a vetor of the same length as number of positions: the observed data log-likelihood for each positions}
 \item{delta.like}{a vetor of the same length as number of samples: the observed data log-likelihood for each samples}
}
\references{
Na You and Gongyi Huang.(2015) Genotyping for Rare Variant Detection Using Next-generation Sequencing Data.
}
\author{
Na You <youn@mail.sysu.edu.cn> and Gongyi Huang<huanggyi@mail2.sysu.edu.cn>
}

