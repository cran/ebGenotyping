\name{z.adj.1}
\alias{z.adj.1}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Adjustment based on data
}
\description{
z.adj.1 is a only part of main function of this package. It calculates the posterior probability of RR, RV and VV of n samples at m positions according to the real data. 
}
\usage{
z.adj.1(zm1, z0, zp1, pm1, p0, pp1, th1 = 0.5, th2 = NULL, adj = FALSE, adj.idx = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{zm1}{
a n*m matrix: the ith row, jth column of the matrix represents the posterior probability of RR of the ith sample at jth position.
}
  \item{z0}{
a n*m matrix: the ith row, jth column of the matrix represents the posterior probability of RV of the ith sample at jth position.
}
  \item{zp1}{
a n*m matrix: the ith row, jth column of the matrix represents the posterior probability of VV of the jth sample at jth position.
}
  \item{pm1}{
a n*m matrix: an estimate of the probability of observing a non-reference count when the genotype is RR.
}
  \item{p0}{
a n*m matrix: an estimate of the probability of observing a non-reference count when the genotype is RV
}
  \item{pp1}{
a n*m matrix: an estimate of the probability of observing a non-reference count when the genotype is VV.
}
  \item{th1}{
a single value: the cutoff of the probability of observing a non-reference count when the genotype is RR.
}
  \item{th2}{
a single value: the cutoff of the probability of observing a non-reference count when the genotype is VV.
}
  \item{adj}{
the default is FALSE. 
}
  \item{adj.idx}{
the default is NULL.
}
}
\value{
  \item{zm1}{
a n*m matrix: the ith row, jth column of the matrix represents the posterior probability of RR of the ith sample at jth position.
}
  \item{z0}{
a n*m matrix: the ith row, jth column of the matrix represents the posterior probability of RV of the ith sample at jth position.
}
  \item{zp1}{
a n*m matrix: the ith row, jth column of the matrix represents the posterior probability of VV of the jth sample at jth position.
}
 \item{adj}{
 final adjustment status which denotes whether the adjustment is done. 
} 
 \item{adj.idx}{
 the sample ID, which needs to be adjusted.
}
}
\references{
Na You and Gongyi Huang.(2015) Genotyping for Rare Variant Detection Using Next-generation Sequencing Data.
}
\author{
Na You <youn@mail.sysu.edu.cn> and Gongyi Huang<huanggyi@mail2.sysu.edu.cn>
}

