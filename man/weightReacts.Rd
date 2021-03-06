% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/weightReacts.R
\name{weightReacts}
\alias{weightReacts}
\title{Compute weights of rescue reactions}
\usage{
weightReacts(model, mc.cores = 1, gene.num = 1, draw.num = 1000)
}
\arguments{
\item{model}{An object of class \code{modelorg} indicating the weighted \code{rescue} model obtained from the rescue process.}

\item{mc.cores}{The number of cores to use (at least 1), i.e. at most how many child processes will be run simultaneously. Default: 1.}

\item{gene.num}{The number of genes to remove. If 1, \code{oneGeneDel} will be performed and draw.num will be ignored. Default: 1.}

\item{draw.num}{The number of random draws. Default: 1000.}
}
\value{
A vector of weights for rescue reactions and an object of class \code{phylo} for colored plot of fitness weighting schema.
}
\description{
The function weightReacts computes the weights of rescue reactions.
}
\examples{
data(Ec_core)
mod <- rescue(Ec_core, target=0.1)
weightReacts(changeObjFunc(mod$rescue, react=rownames(mod$coef),
                          obj_coef=mod$coef))
}
