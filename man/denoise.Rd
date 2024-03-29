% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{denoise}
\alias{denoise}
\title{denoise}
\usage{
denoise(points, k = 16, prob = 0.3)
}
\arguments{
\item{points}{a data frame with columns x,y,z.}

\item{k}{integer, number of mixture components for the GMM}

\item{prob}{probability level in the range [0,1] to identify high density regions}
}
\value{
a point set
}
\description{
Point density is estimated using a Gaussian mixture model and points in low
density regions are considered as noise and removed.
}
