\name{estep}
\alias{estep}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
E step 
}
\description{
This function calculates the E step of ECM algorithm for the model described in 'An Empirical Bayes Method for Genotyping and SNP detection Using Multi-sample Next-generation Sequencing Data'.

}
\usage{
estep(mu, delta, pm1, p0, dat, cvg)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mu}{
a vetor of the same length as number of positions: the position effect.
}
  \item{delta}{
a vetor of the same length as number of samples: the sample effect.
}
  \item{pm1}{
a single value,which is larger than 0 and less than 1: the probability of RR. 
}
  \item{p0}{
a single value,which is larger than 0 and less than 1: the probability of RV. 
}
  \item{dat}{
a n*m matrix: the ith row, jth column of the matrix represents the non-reference counts of ith sample at jth position. 
}
  \item{cvg}{
a n*m matrix: the ith row, jth column of the matrix represents the depth of ith sample at jth position. 
}
}
\details{
The value of mu and delta must satisfy that each element of outer(delta,mu,"+") must less than zero. This is the requirement of the model 
described in paper "Genotyping for Rare Variant Detection Using Next-generation Sequencing Data."
}
\value{
\item{zRR}{a n*m matrix: the posterior probabilities of genotype RR for n samples at m positions}
\item{zRV}{a n*m matrix: the posterior probabilities of genotype RV for n samples at m positions}
\item{zVV}{a n*m matrix: the posterior probabilities of genotype VV for n samples at m positions}
}
\references{
Na You and Gongyi Huang.(2016) An Empirical Bayes Method for Genotyping and SNP detection Using Multi-sample Next-generation Sequencing Data.
}
\author{
Na You <youn@mail.sysu.edu.cn> and Gongyi Huang<53hgy@163.com>
}
\note{
The most important function in this package is "ecm". "estep" is a function called by "ecm" to realize one E step in the whole process of iteration in "ecm".
}
