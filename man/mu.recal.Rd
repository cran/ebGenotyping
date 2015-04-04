\name{mu.recal}
\alias{mu.recal}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
      One step maximization of mu when given other parameters
}
\description{
mu.recal calculates the optimal estimator of the ith position effect when given other parameters.
}
\usage{
mu.recal(dat, cvg, zm1, z0, zp1, delta, eta)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{
a vetor of the same length as number of samples: the non-reference counts of n samples at the ith position. 
}
  \item{cvg}{
a vetor of the same length as number of samples: the coverage of n samples at the ith position.
}
  \item{zm1}{
a vetor of the same length as number of samples: posterior probability of RR of n samples at the ith position.
}
  \item{z0}{
a vetor of the same length as number of samples: posterior probability of RV of n samples at the ith position.
}
  \item{zp1}{
a vetor of the same length as number of samples: posterior probability of VV of n samples at the ith position.
}
  \item{delta}{
a vetor of the same length as number of samples: the estimator of delta(sample effects).
}
  \item{eta}{
a single value representing the estimator of eta.
}
}

\value{
a single value representing one step maximization of mu of the ith position.
}
\references{
Na You and Gongyi Huang.(2015) Genotyping for Rare Variant Detection Using Next-generation Sequencing Data.
}
\author{
Na You <youn@mail.sysu.edu.cn> and Gongyi Huang<huanggyi@mail2.sysu.edu.cn>
}

