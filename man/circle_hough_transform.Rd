% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{circle_hough_transform}
\alias{circle_hough_transform}
\title{Circle Hough transform}
\usage{
circle_hough_transform(
  pixels,
  rmin,
  rmax,
  threshold,
  resolution = 360,
  min.separation = rmin/4,
  ncpu = 1
)
}
\arguments{
\item{pixels}{input data, either a matrix representing a 2D image or a data frame of signal coordinates with columns x, y.
For images, background is expected to be 0 and signal to have positive values.}

\item{rmin}{minimum search radius.}

\item{rmax}{maximum search radius.}

\item{threshold}{score threshold between 0 and 1.}

\item{resolution}{number of steps in the circle transform (default: 360). This represents the maximum number of votes a point can get.}

\item{min.separation}{distance between circle centres below which overlapping circles are considered the same and merged (default to 0.25*rmin)}

\item{ncpu}{number of threads to use to speed up computation (default: 1)}
}
\value{
a data frame with columns x, y, r and score
}
\description{
Extract coordinates of the centres of circles from a 2D image using the Hough transform
}
\examples{
point.set <- data.frame(x = c(-9.8,-5.2,12.5,2.5,4.5,1.3,-0.2,0.4,9.3,-1.4,0.5,-1.1,-7.7),
                        y = c(-4.2,1.5,-0.5,12,-3,-7.2,10.9,6.7,-1.3,10,6.7,-6.2,2.9))
circles <- circle_hough_transform(pixels = point.set, rmin = 3, rmax = 6, resolution = 100,
                                  threshold = 0.1, ncpu = 1)
}
