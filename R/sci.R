#' calc_sci 
#'
#' This function calculates the striatal connectivity index (sci).
#'
#' @param img brain image to calculate the sci on.
#' @param vox N x 3 matrix of voxel coordinates to use for regions of
#'     interest.
#' @param normdf N x 2 data frame to use for normalization.
#' @param weights N x 1 vector to use for weighting.
#' @keywords striatum, functional connectivity, correlation 
#' @export
#' @examples
#' calc_sci()
calc_sci <- function(img, vox, normdf=NULL,
                     weights=1) {
  #
  # calculate the striatal connectivity index
  
  # extract all 91 values from the rois
  vals <- as.numeric(sapply(as.list(data.frame(t(vox))),
                            function(x) extract_roi_val(x, img=img)))

  # normalize and weight the raw sci values
  if (!is.null(normdf)) {
    valsz <- normalize_sci(vals, normdf)
  } else {
    valsz <- vals
  }
  
  # weight by whatever is given as weights (defaults to 1)
  valszw <- weight_sci(valsz, weights) 

  # compute the two flavors of sci
  sci <- sum(valszw)
  mu_sci <- sci/nrow(vox)
  return(list("sci"=sci, "sci_r"=mu_sci))
}

#' is_high
#'
#' This function tests if a given SCI is higher than the cutoff value
#' (which defaults to 3.8 as in Sarpal 2016, Am J Psychiatry).
#
#' @param sci striatal connectivity index 
#' @param cutoff striatal connectivity index cutoff (default is 3.8) 
#' @keywords striatum, functional connectivity, correlation 
#' @export
#' @examples
#' is_high(sci)
is_high <- function(sci, cutoff=3.8) {
  #
  # tests if the sci is higher than the cutoff 
  return(sci > cutoff)
}

#' is_low
#'
#' This function tests if a given SCI is equal or lower than the
#' cutoff value (which defaults to 3.8 as in Sarpal 2016,
#' Am J Psychiatry).
#
#' @param sci striatal connectivity index 
#' @param cutoff striatal connectivity index cutoff (default is 3.8) 
#' @keywords striatum, functional connectivity, correlation 
#' @export
#' @examples
#' is_low(sci)
is_low <- function(sci, cutoff=3.8) {
  #
  # tests if the sci is equal or lower than the cutoff
  return(!is_high(sci=sci, cutoff=cutoff))
}

#' load_params
#'
#' This function loads the default parameters needed for the
#' SCI calculation, including the default 91 MNI coordinates, the
#' normalization data frame, and the weights vector.
#' @keywords striatum, functional connectivity, correlation 
#' @export
#' @examples
#' load_params()
load_params <- function() {
  #
  # loads the default sci parameters
  mni <- read.csv(system.file("data", "sci_mni.csv", package="sci"))
  weights <- read.csv(system.file("data", "sci_weights.csv",
                                  package="sci"))
  normdf <- read.csv(system.file("data", "sci_norm.csv", package="sci"))
  return(list("mni"=mni, "weights"=weights, "normdf"=normdf))
}


#' normalize_sci 
#'
#' This function normalizes the N x 1 individual SCI connectivity values
#' using a N x 2 normalization data frame (which should comprise columns
#' 'mean' and 'sd').
#' @keywords normalization, striatal connectivity
#' @param vals individual raw sci connectivity values
#' @param normdf normalization data frame
#' @export
#' @examples
#' normalize_sci(vals, normdf)
normalize_sci <- function(vals, normdf) {
  #
  # normalize the sci values
  valsz <- (vals-normdf$mean)/normdf$sd
  return(valsz)
}

#' weight_sci 
#'
#' This function weights the N x 1 individual SCI connectivity values
#' using a N x 1 vector of weights.
#' @keywords striatal connectivity, weighting
#' @param vals individual sci connectivity values
#' @param weights vector of weights
#' @export
#' @examples
#' weight_sci(vals, weights)
weight_sci <- function(vals, weights) {
  #
  # weight the sci values
  valsw <- vals * weights
  return(valsw)
}


#' mni2vox
#'
#' This function transforms MNI coordinates into voxel coordinates.
#' @keywords neuroimaging, MNI space, voxel space
#' @param mni MNI coordinates
#' @export
#' @examples
#' mni2vox(mni)
mni2vox <- function(mni) {
  #
  # transform mni 2 voxel coordinates

  vox <- matrix(ncol=3, nrow=nrow(mni))
  # recipe to transform mni coordinates to voxels
  vox[, 1] <- (45-(mni[, 1]/2))
  vox[, 2] <- (63+(mni[, 2]/2))
  vox[, 3] <- (36+(mni[, 3]/2))
  return(vox)
}

#' has_afni
#'
#' This function checks if afni is installed (at the standard location).
#' @keywords AFNI, neuroimaging, UNIX 
#' @export
#' @examples
#' has_afni()
has_afni <- function() {
  #
  # checks if afni is installed

  # standard path to afni
  pth <- "/usr/local/opt/afni"

  cmd <- paste0("export PATH=$PATH:", pth, " && which 3dmaskave")
  out <- system(cmd, intern=FALSE)
  #print(out)
  if (out==1) {
    return(FALSE)
  }
  else {
    return(TRUE)
  }
}

#' which_afni
#'
#' This function returns the path to AFNI
#'
#' @export
#' @examples
#' which_afni()
which_afni <- function() {
  #
  # returns the path to afni
  pth <- system("which afni", intern=TRUE)
  return(pth)
}

#' parse_afni_cmd
#'
#' This function returns the full AFNI command 3dmaskave
#'
#' @param params string of options passed to 3dmaskave (defaults to
#'     -ibox)
#' @export
#' @examples
#' parse_afni_cmd(params)
parse_afni_cmd <- function(params=" -ibox") {
  #
  # returns the path to afni
  return(paste0(which_afni(), params))
}


#' extract_roi_val 
#'
#' This function extracts the roi value of a given coordinate (range)
#' using AFNI's 3dmaskave. It is the workhorse of the sci-package.
#' @param vox N x 3 matrix of voxel coordinates to use for regions of
#'     interest.
#' @param img a brain image filename.
#' @keywords voxel, extraction, brain imaging, AFNI
#' @export
#' @examples
#' extract_roi_val()
extract_roi_val <- function(vox, img) {
  #
  # extract voxel roi value using afni's 3dmaskave

  # check if afni is installed
  if (!has_afni()) {
    cat("afni not installed!\n")
    return(-1)
  }

  # defaut afni location
  afni_bin <- "/usr/local/opt/afni/afni"

  afni_3dmaskave <- paste0(ifelse(file.exists(afni_bin),
                           "/usr/local/opt/afni/3dmaskave",
                           gsub("[\r\n]", "", which_afni())), " -ibox")

  rng <- paste0(vox[1]-1, ":", vox[1]+1, " ", 
                vox[2]-1, ":", vox[2]+1, " ",
                vox[3]-1, ":", vox[3]+1, " ")
  cmd <- paste0(afni_3dmaskave, " ", rng, img)
  #print(cmd)
  # run the command
  out <- system(cmd, intern=TRUE)

  # split and transform the output string
  val <- as.numeric(strsplit(out, split=" ")[[1]][1])

  return(val)
}
