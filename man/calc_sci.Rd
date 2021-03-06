% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sci.R
\name{calc_sci}
\alias{calc_sci}
\title{calc_sci}
\usage{
calc_sci(imgdf, mnidf, normdf = NULL, weights = 1, rng = 1)
}
\arguments{
\item{imgdf}{K x 2 data frame of K seeds and their respective brain
images to calculate the SCI on.}

\item{mnidf}{N x 4 data frame of K seeds and N MNI coordinates.}

\item{normdf}{N x 2 data frame to use for normalization.}

\item{weights}{N x 1 vector to use for weighting.}

\item{rng}{1 x 1 vector for how many voxels to additionally include
in the ROI extraction. The default value of 1 means: target voxel
+/- 1 voxel, i.e., a box of 3 x 3 x 3 = 27 voxels in total}
}
\description{
This function calculates the striatal connectivity index (sci). In
the original paper (Sarpal 2016, Am J Psychiatry), there were K=12
seeds and N=91 connections. The function, however, can calculate the
sci for various combinations of K seeds and N connections.
}
\examples{
calc_sci()
}
\keyword{connectivity,}
\keyword{correlation}
\keyword{functional}
\keyword{striatum,}
