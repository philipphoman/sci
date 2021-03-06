#' calc_sci 
#'
#' This function calculates the striatal connectivity index (sci). In
#' the original paper (Sarpal 2016, Am J Psychiatry), there were K=12
#' seeds and N=91 connections. The function, however, can calculate the
#' sci for various combinations of K seeds and N connections.
#'
#' @param imgdf K x 2 data frame of K seeds and their respective brain
#'     images to calculate the SCI on.
#' @param mnidf N x 4 data frame of K seeds and N MNI coordinates.
#' @param normdf N x 2 data frame to use for normalization.
#' @param weights N x 1 vector to use for weighting.
#' @param rng 1 x 1 vector for how many voxels to additionally include
#'     in the ROI extraction. The default value of 1 means: target voxel
#'     +/- 1 voxel, i.e., a box of 3 x 3 x 3 = 27 voxels in total
#' @keywords striatum, functional connectivity, correlation
#' @export
#' @examples
#' calc_sci()
calc_sci <- function(imgdf, mnidf, normdf=NULL,
                     weights=1, rng=1) {
  #
  # calculate the striatal connectivity index

  # get all seeds
  seeds <- unique(mnidf$seed)

  # extract all values from the rois
  vals <- unlist(lapply(seeds, function(x) {
    sapply(as.list(data.frame(t(get_seed_connections(mnidf, x)))),
           function(y) extract_roi_val(mni2vox(y),
                                       imgdf$img[imgdf$seed==x],
                                       rng))
  }))

    
  #vals <- as.numeric(sapply(as.list(data.frame(t(vox))),
  #                          function(x) extract_roi_val(x, img=img)))

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
  mu_sci <- sci/nrow(mnidf)
  return(list("sci"=sci, "sci_r"=mu_sci, "valszw"=valszw,
              "seed"=mnidf$seed))
}

#' get_seed_connections
#'
#' This function returns a K x 3 matrix of roi coordinates for a seed
#' region
#'
#' @param mnidf N x 4 data frame of regions of seed regions and their
#'     respective connectivity coordinates
#' @param seed a string indicating the seed for which connectivity
#'     coordinates should be returned
#' @keywords functional connectivity, seed regions
#' @export
get_seed_connections <- function(mnidf, seed) {
  return(as.matrix(mnidf[mnidf$seed==seed, 2:4]))
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
#' normalization data frame, the weights vector, and an illustrative
#' set of 12 whole brain images for each seed.
#' @keywords striatum, functional connectivity, correlation 
#' @export
#' @examples
#' load_params()
load_params <- function() {
  #
  # loads the default sci parameters
  mnidf <- readr::read_csv(system.file("data", "sci_mni.csv",
                                     package="sci"))

  weights <- readr::read_csv(system.file("data", "sci_weights.csv",
                                         package="sci"))

  normdf <- readr::read_csv(system.file("data", "sci_norm.csv",
                                        package="sci"))

  #imgdf <- readr::read_csv(system.file("data", "sci_img_example.csv",
  #                                     package="sci"))
  seeds <- unique(mnidf$seed)

  imgdf <- data.frame(seed=seeds,
                      img=system.file("data/nii/2mm",
                                      paste0(seeds, ".nii.gz"),
                                      package="sci"))
                                                  
  #sapply(imgdf$img, function(x) system(paste0("gunzip ", x))) 

  return(list("mnidf"=mnidf, "weights"=weights, "normdf"=normdf,
              "imgdf"=imgdf))
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

  #vox <- matrix(ncol=3, nrow=nrow(mni))
  vox <- vector("numeric", 3)
  # recipe to transform mni coordinates to voxels
  vox[1] <- (45-(mni[1]/2))
  vox[2] <- (63+(mni[2]/2))
  vox[3] <- (36+(mni[3]/2))
  return(vox)
}

#' add_afni_pth_cmd
#'
#' This function adds the standard AFNI path to PATH
#' @export
add_afni_pth_cmd <- function(pth="/usr/local/opt/afni") {
  #
  # adds the standard afni path to PATH
  return(paste0("PATH=$PATH:", pth))
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
  pth <- system(paste0(add_afni_pth_cmd(), " && which 3dmaskave"),
                       intern=TRUE)
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
parse_afni_cmd <- function(params=" -quiet -ibox") {
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
#' @param rng 1 x 1 vector for the additional range of voxels to
#'     extract.
#' @keywords voxel, extraction, brain imaging, AFNI
#' @export
#' @examples
#' extract_roi_val()
extract_roi_val <- function(vox, img, rng=1) {
  #
  # extract voxel roi value using afni's 3dmaskave

  # check if afni is installed
  if (!has_afni()) {
    cat("afni not installed!\n")
    return(-1)
  }

  # defaut afni location
  afni_bin <- "/usr/local/opt/afni/afni"

  #afni_3dmaskave <- paste0(ifelse(file.exists(afni_bin),
  #                         "/usr/local/opt/afni/3dmaskave",
  #                         gsub("[\r\n]", "", which_afni())), " -ibox")


  afni_3dmaskave <- parse_afni_cmd()
    
  box <- paste0(vox[1]-rng, ":", vox[1]+rng, " ", 
                vox[2]-rng, ":", vox[2]+rng, " ",
                vox[3]-rng, ":", vox[3]+rng, " ")
  cmd <- paste0(afni_3dmaskave, " ", box, img)
  #print(cmd)
  # run the command
  out <- system(cmd, intern=TRUE)

  # split and transform the output string
  val <- as.numeric(strsplit(out, split=" ")[[1]][1])

  return(val)
}

#' run_example
#'
#' This function runs the SCI calculation using the example parameters
#' provided with the package
#' @export
run_example <- function() {
  params <- load_params()
  out <- calc_sci(params$imgdf, params$mnidf, params$normdf,
                  params$weights)
  return(out)
}
