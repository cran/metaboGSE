% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/metaboGSE.R
\name{area}
\alias{area}
\title{Compute area between two curves}
\usage{
area(y1, y2, x = NULL)
}
\arguments{
\item{y1}{A vector of y-coordinates for the first curve.}

\item{y2}{A vector of y-coordinates for the second curve.}

\item{x}{A vector of increasing x-coordinates. Default, an equidistant vector of length length(y1) from 0 to 1.}
}
\value{
The area
}
\description{
This function computes the area formed by two curves.
}
\keyword{internal}
