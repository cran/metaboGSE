% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rescue.R
\name{rescue}
\alias{rescue}
\title{Rescue a model}
\usage{
rescue(
  model,
  target,
  react = NULL,
  weight.type = "r",
  timeout = 12,
  prefix.rescue = NA,
  prefix.rescued = NA,
  rescue.threshold = 1e-05
)
}
\arguments{
\item{model}{An object of class \code{modelorg}.}

\item{target}{A numeric vector for growth target.}

\item{react}{A numeric vector or a character vector containing reaction id's. Default: reactions in objective function.}

\item{weight.type}{A character indicating which type of weighting to use in model objective modification, i: 1, r: 1/coefficient, s: 1/sqrt(coefficient). Default: r.}

\item{timeout}{The maximum time in seconds to allow for LP call to return. Default: 12.}

\item{prefix.rescue}{A string indicating the prefix of output rescue model. Default: no output.}

\item{prefix.rescued}{A string indicating the prefix of output rescued model. Default: no output.}

\item{rescue.threshold}{A numeric value indicating the threshold to consider a rescue. Default: 1e-5.}
}
\value{
The rescue and rescued models, as well as the coefficient set to rescue reactions.
\code{SYBIL_SETTINGS("OPT_DIRECTION")} is set to "min".
}
\description{
The function rescues a given model.
}
\examples{
data(Ec_core)
rescue(Ec_core, target=0.1)
}
