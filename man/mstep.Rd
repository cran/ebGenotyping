\name{mstep}
\alias{mstep}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
CM steps
}
\description{
This function calculates the CM steps of ECM algorithm for the model described in 'An Empirical Bayes Method for Genotyping and SNP detection Using Multi-sample Next-generation Sequencing Data'.

}
\usage{
mstep(mu0, delta0, zm1, z0, zp1, dat, cvg, eps = 1e-06)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mu0}{
a vetor of the same length as number of positions: the initial value of position effect mu.
}
  \item{delta0}{
a vetor of the same length as number of samples: the initial value of sample effect delta.
}
  \item{zm1}{
the output of estep: the posterior probabilities of genotype RR for n samples at m positions
}
  \item{z0}{
the output of estep: the posterior probabilities of genotype RV for n samples at m positions
}
  \item{zp1}{
the output of estep: the posterior probabilities of genotype VV for n samples at m positions
}
  \item{dat}{
a n*m matrix: the ith row, jth column of the matrix represents the non-reference counts of ith sample at jth position. 
}
  \item{cvg}{
a n*m matrix: the ith row, jth column of the matrix represents the depth of ith sample at jth position. 
}
  \item{eps}{
a single value: a threshold to control the convergence criterion. The default is 1e-06.
}
}
\details{
The value of mu and delta must satisfy that each element of outer(delta,mu,"+") must less than zero. This is the requirement of the model 
described in 'An Empirical Bayes Method for Genotyping and SNP detection Using Multi-sample Next-generation Sequencing Data'.
}
\value{
\item{mu}{the optimal value of mu in current CM steps}
\item{delta}{the optimal value of delta in current CM steps}
\item{pRR}{the optimal value of the probability of RR in current CM steps}
\item{pRV}{the optimal value of the probability of RV in current CM steps}
}
\references{
Na You and Gongyi Huang.(2016) An Empirical Bayes Method for Genotyping and SNP detection Using Multi-sample Next-generation Sequencing Data.
}
\author{
Na You <youn@mail.sysu.edu.cn> and Gongyi Huang<53hgy@163.com>
}
\note{
The most important function in this package is "ecm". "mstep" is a function called by "ecm" to realize one M step(several CM steps) in the whole process of iteration in "ecm".
}

