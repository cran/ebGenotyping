\name{logit}
\alias{logit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Logit Transformation
}
\description{
This function is for calculating the logit transformation:ln(x/(1-x))
}
\usage{
logit(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A numeric vector, whose elements are all greater than 0 and less than 1.
}
}
\value{
ln(x/(1-x))
}
\references{
Na You and Gongyi Huang.(2016) An Empirical Bayes Method for Genotyping and SNP detection Using Multi-sample Next-generation Sequencing Data.
}
\author{
Na You <youn@mail.sysu.edu.cn> and Gongyi Huang<53hgy@163.com>
}
\examples{
logit(0.5)
}
