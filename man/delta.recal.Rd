\name{delta.recal}
\alias{delta.recal}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
      One step maximization of delta when given other parameters
}
\description{
delta.recal calculates the optimal estimator of certain sample effect when given other parameters.
}
\usage{
delta.recal(dat, cvg, zm1, z0, zp1, mu, eta)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{
a vetor of the same length as number of positions: the non-reference counts of the jth sample at m positions. 
}
  \item{cvg}{
a vetor of the same length as number of positions: the coverage of jth sample at m positions. 
}
  \item{zm1}{
a vetor of the same length as number of positions: posterior probability of RR of the jth sample.
}
  \item{z0}{
a vetor of the same length as number of positions: posterior probability of RV of the jth sample.
}
  \item{zp1}{
a vetor of the same length as number of positions: posterior probability of VV of the jth sample.
}
  \item{mu}{
a vetor of the same length as number of positions: the estimator of mu(position effects).
}
  \item{eta}{
a single value representing the estimator of eta.
}
}
\value{
a single value representing one step maximization of delta of this sample.
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Na You and Gongyi Huang.(2015) Genotyping for Rare Variant Detection Using Next-generation Sequencing Data.
}
\author{
Na You <youn@mail.sysu.edu.cn> and Gongyi Huang<huanggyi@mail2.sysu.edu.cn>
}

