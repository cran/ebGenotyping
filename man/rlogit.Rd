\name{rlogit}
\alias{rlogit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Inverse Logit Transformation
}
\description{
This function is for calculating the inverse logit transformation:exp(x)/(1+exp(x))
}
\usage{
rlogit(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A numeric vector
}
}
\details{
In order to avoid overflow, we define the function like this: rlogit <- function(x) return(ifelse(x>100,1,exp(x)/(1+exp(x))))
}
\value{
exp(x)/(1+exp(x))
}
\references{
Na You and Gongyi Huang.(2016) An Empirical Bayes Method for Genotyping and SNP detection Using Multi-sample Next-generation Sequencing Data.
}
\author{
Na You <youn@mail.sysu.edu.cn> and Gongyi Huang<53hgy@163.com>
}
\examples{
rlogit(-3)
}