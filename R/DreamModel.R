#' Data used for the DREAM3 challenge
#'
#' This data object contains the model used in the package vignette, already 
#' loaded and formatted as a Model object. This is to be used with the data 
#' in "CNOListDREAM"
#'
#' J. Saez-Rodriguez, L. G. Alexopoulos, J. Epperlein, R. Samaga, D. A.
#' Lauffenburger, S. Klamt and P. K. Sorger. Discrete logic modeling as a means to
#' link protein signaling networks with functional analysis of mammalian signal
#' transduction, Molecular Systems Biology, 5:331, 2009.
#' 
#' Prill RJ, Marbach D, Saez-Rodriguez J, Sorger PK, Alexopoulos LG, Xue X,
#' Clarke ND, Altan-Bonnet G, and Stolovitzky G. Towards a rigorous assessment of
#' systems biology models: the DREAM3 challenges. PLoS One, 5(2):e9202, 2010.
#' @format DreamModel is a list with fields "reacID" (character vector), "namesSpecies" (character vector), "interMat" (numerical matrix), "notMat"(numerical matrix).
#' @usage DreamModel
#' @aliases LiverDREAM, DreamModel
"DreamModel"
