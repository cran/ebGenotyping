\name{f.mu.1}
\alias{f.mu.1}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
The first-order derivative with respect to mu
}
\description{
f.mu.1 calculates the first-order derivative of Q function with respect to mu.
}
\usage{
f.mu.1(x, delta, eta, dat, cvg, zm1, z0, zp1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
a vetor of the same length as number of positions: the current estimate of mu.
}
  \item{delta}{
a vetor of the same length as number of samples: the current estimate of delta.
}
  \item{eta}{
a single value: the current estimate of eta. 
}
  \item{dat}{
a n*m matrix: the ith row, jth column of the matrix represents the non-reference counts of ith sample at jth position. 
}
  \item{cvg}{
a n*m matrix: the ith row, jth column of the matrix represents the depth of ith sample at jth position. 
}
  \item{zm1}{
a n*m matrix: the ith row, jth column of the matrix represents the posterior probability of RR of the ith sample at jth position.
}
  \item{z0}{
a n*m matrix: the ith row, jth column of the matrix represents the posterior probability of RV of the ith sample at jth position.
}
  \item{zp1}{
a n*m matrix: the ith row, jth column of the matrix represents the posterior probability of VV of the jth sample at jth position.
}
}
\value{
a vetor of the same length as number of positions: the first-order derivative with respect to mu.
}
\references{
Na You and Gongyi Huang.(2015) Genotyping for Rare Variant Detection Using Next-generation Sequencing Data.
}
\author{
Na You <youn@mail.sysu.edu.cn> and Gongyi Huang<huanggyi@mail2.sysu.edu.cn>
}