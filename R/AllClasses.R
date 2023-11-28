
#-------------------------------------------------------------------------------
#' @title Bridge: An S4 class for Bridge analysis
#'
#' @slot cogids Object of class \code{"data.frame"}, a data frame with COG data.
#' @slot tree Object of class \code{"phylo"}, a given species tree.
#' @slot spbranches Object of class \code{"data.frame"}, a data frame 
#' listing branches of a given species tree.
#' @slot orthoroot Object of class \code{"data.frame"}, a data.frame with 
#' results from the 'runBridge' function.
#' @slot orthoct Object of class \code{"data.frame"}, a data.frame with 
#' results from the 'Bridge' constructor. 
#' @slot status Object of class \code{"character"}, a string specifying the 
#' status of the Bridge object based on the available methods.
#'
#' @method runBridge \code{\link{runBridge}}
#' @method getBridge \code{\link{getBridge}}
#' @method plotBridge \code{\link{plotBridge}}
#' @aliases Bridge-class
#' @return An S4 class object.
#' @section Constructor:
#' see \code{\link{newBridge}} constructor.
#' @exportClass Bridge
#'
setClass(
    "Bridge",
    representation(
        cogids = "data.frame",
        tree = "ANY",
        spbranches = "data.frame",
        orthoroot = "data.frame",
        orthoct = "data.frame",
        status = "character"
    ),
    prototype = list(
        cogids = data.frame(),
        tree = "ANY",
        spbranches = data.frame(),
        orthoct = data.frame(),
        orthoroot = data.frame(),
        status = character()
    )
)
