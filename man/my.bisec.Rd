\name{my.bisec}
\alias{my.bisec}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Bisection method to find the root
}
\description{
This function is to apply bisection method to find the root of a function f.
}
\usage{
my.bisec(f, int.l, int.u, eps = 1e-06)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{f}{
the function for which the root is sought.
}
  \item{int.l}{
a vector containing the lower bound of the interval to be searched for the root. The length of the vector is the same as that of the input of function f.
}
  \item{int.u}{
a vector containing the upper bound of the interval to be searched for the root. The length of the vector is the same as that of the input of function f.
}
  \item{eps}{
a single value: a threshold to control the convergence criterion. The default is 1e-06.
}
}
\details{
Both int.l and int.u must be specified: the upper bound must be strictly larger than the lower bound. 

The function f must be well defined without invalid output(NA, nan, Inf, etc). 

The length of the input of function f, the output of function f, int.l and int.u must be the same.

}
\value{
a vector containing the root of the function. If there is no root in the interval (int.l, int.u), lower bound of the interval will be returned.
}
\references{
Na You and Gongyi Huang.(2016) An Empirical Bayes Method for Genotyping and SNP detection Using Multi-sample Next-generation Sequencing Data.
}
\author{
Na You <youn@mail.sysu.edu.cn> and Gongyi Huang<53hgy@163.com>
}

\examples{
f <- function(x){
a <- 1:10
return(x-a)
}
my.bisec(f=f, int.l=rep(-1,10), int.u=rep(11,10), eps = 1e-08)
}

