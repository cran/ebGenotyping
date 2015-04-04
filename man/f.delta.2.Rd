\name{f.delta.2}
\alias{f.delta.2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
The second-order derivative with respect to delta
}
\description{
f.delta.2 calculates the second-order derivative of Q function with respect to delta.
}
\usage{
f.delta.2(x, mu, eta, dat, cvg, zm1, z0, zp1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
a vetor of the same length as number of samples: the current estimate of delta.
}
  \item{mu}{
a vetor of the same length as number of positions: the current estimate of mu.
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
a vetor of the same length as number of samples: the second-order derivative with respect to delta.
}
\references{
Na You and Gongyi Huang.(2015) Genotyping for Rare Variant Detection Using Next-generation Sequencing Data.
}
\author{
Na You <youn@mail.sysu.edu.cn> and Gongyi Huang<huanggyi@mail2.sysu.edu.cn>
}