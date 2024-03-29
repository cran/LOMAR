% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data_input.R
\name{point_sets_from_locs}
\alias{point_sets_from_locs}
\title{point_sets_from_locs}
\usage{
point_sets_from_locs(
  locs = NULL,
  channels = NULL,
  min.cardinality = NULL,
  max.cardinality = NULL,
  crop.size = NULL,
  keep.locprec = TRUE,
  sample.size = NULL,
  ignore.site = FALSE,
  cluster.points = FALSE,
  eps = NULL,
  minPts = NULL
)
}
\arguments{
\item{locs, }{a data frame with columns x[nm], y[nm], z[nm] and optionally site[numbers], locprec[nm] and locprecz[nm], other columns are ignored.}

\item{channels}{vector of integers indicating which channel(s) of a multicolour experiment to extract point sets from.}

\item{min.cardinality}{filter out point sets with less than the specified number of points.}

\item{max.cardinality}{filter out point sets with more than the specified number of points.}

\item{crop.size}{remove points from a set if they are further away than the specified distance from the center of the set.}

\item{keep.locprec}{logical (default:TRUE). Whether to keep locprec information for each point.}

\item{sample.size}{returns this number of randomly selected point sets. Selects the point sets after applying eventual filtering.}

\item{ignore.site}{logical (default: FALSE), set to TRUE if point set membership is not present or needed.}

\item{cluster.points}{logical (default: FALSE), whether to cluster the points using DBSCAN (only if ignore.site is also TRUE).}

\item{eps}{DBSCAN parameter, size of the epsilon neighbourhood}

\item{minPts}{DBSCAN parameter, number of minimum points in the eps region}
}
\value{
a list of matrices with columns x,y,z, optionally locprec and name set to the value of the site column (if applicable).
}
\description{
Extracts list of point sets from a data frame of single molecule localization coordinates.
 By default, uses point set membership indicated in the site column.
}
\examples{
data.file <- system.file("test_data", "simulated_NUP107_data.csv", package = "LOMAR",
 mustWork = TRUE)
locs <- locs_from_csv(file = data.file, locprec.filter = 20)
point.sets <- point_sets_from_locs(locs, keep.locprec = TRUE, min.cardinality = 15)
}
