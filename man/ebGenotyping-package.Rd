\name{ebGenotyping-package}
\alias{ebGenotyping-package}
\alias{ebGenotyping}
\docType{package}
\title{
Genotyping and SNP Detection using Next Generation Sequencing Data
}
\description{
Genotyping the population using next generation sequencing data is essentially important for the rare variant detection. In order to distinguish the genomic structural variation from sequencing error, we propose a statistical model which involves the genotype effect through a latent variable to depict the distribution of non-reference allele frequency data among different samples and different genome loci, while decomposing the sequencing error into sample effect and positional effect. An ECM algorithm is implemented to estimate the model parameters, and then the genotypes and SNPs are inferred based on the empirical Bayes method.
}
\details{
\tabular{ll}{
Package: \tab ebGenotyping\cr
Type: \tab Package\cr
Version: \tab 2.0.1\cr
Date: \tab 2016-04-07\cr
License: \tab GPL-2 \cr
}
The most important function is ecm, which is used to establish the model described in 'An Empirical Bayes Method for Genotyping and SNP detection Using Multi-sample Next-generation Sequencing Data' to do genotyping and SNP detection using NGS data.

}
\author{
Na You <youn@mail.sysu.edu.cn> and Gongyi Huang<53hgy@163.com>
}
\references{
Na You and Gongyi Huang.(2016) An Empirical Bayes Method for Genotyping and SNP detection Using Multi-sample Next-generation Sequencing Data.
}
