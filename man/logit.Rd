\name{logit}
\alias{logit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
       Logit Transformation
}
\description{
Logit Transformation:ln(x/(1-x))
}
\usage{
logit(x)
}
\arguments{
  \item{x}{
A probability, which is greater than 0 and less than 1.
}
}
\value{
ln(x/(1-x))
}
\references{
Na You and Gongyi Huang.(2015) Genotyping for Rare Variant Detection Using Next-generation Sequencing Data.
}
\author{
Na You <youn@mail.sysu.edu.cn> and Gongyi Huang<huanggyi@mail2.sysu.edu.cn>
}
\examples{
logit(0.5)
}

