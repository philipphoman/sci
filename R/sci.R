#' calc_sci 
#'
#' This function calculates the striatal connectivity index.
#' @param img brain image to calculate the sci on.
#' @param vox N x 3 - matrix of voxel coordinates to use for regions of
#'     interest.
#' @param hcdf N x 2 - data frame to use for normalization.
#' @param loadings N x 1 - vector to use for weighting.
#' @keywords striatum, functional connectivity, correlation 
#' @export
#' @examples
#' calc_sci()
calc_sci <- function(img=img, vox=vox, hcdf=hcdf,
                     loadings=sci_loadings[, 1]) {
  #
  # calculate the striatal connectivity index
  
  # extract all 91 values from the rois
  vals <- as.numeric(sapply(as.list(data.frame(t(vox))),
                            function(x) extract_roi_val(x, img=img)))

  # normalize and weight the raw sci values
  valsz <- normalize_sci(vals, hcdf) * sci_loadings[, 1]

  sci <- sum(valsz)
  mu_sci <- sci/nrow(vox)
  return(list(sci, mu_sci))
}

normalize_sci <- function(vals, hcdf) {
  #
  # normalize the sci values
  valsz <- (vals-hcdf$mean)/hcdf$sd
  return(valsz)
}


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
  out <- system(cmd, intern=TRUE)
  if (!out=="") {
    return(TRUE)
  }
  else {
    return(FALSE)
  }
}


extract_roi_val <- function(vox, img) {
  #
  # extract voxel roi value using afni's 3dmaskave

  afni_3dmaskave <- "/usr/local/opt/afni/3dmaskave -ibox"
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
