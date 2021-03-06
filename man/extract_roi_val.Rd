% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sci.R
\name{extract_roi_val}
\alias{extract_roi_val}
\title{extract_roi_val}
\usage{
extract_roi_val(vox, img, rng = 1)
}
\arguments{
\item{vox}{N x 3 matrix of voxel coordinates to use for regions of
interest.}

\item{img}{a brain image filename.}

\item{rng}{1 x 1 vector for the additional range of voxels to
extract.}
}
\description{
This function extracts the roi value of a given coordinate (range)
using AFNI's 3dmaskave. It is the workhorse of the sci-package.
}
\examples{
extract_roi_val()
}
\keyword{AFNI}
\keyword{brain}
\keyword{extraction,}
\keyword{imaging,}
\keyword{voxel,}
