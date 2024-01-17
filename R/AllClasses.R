
#-------------------------------------------------------------------------------
#' @title Bridge: An S4 class for Bridge analysis
#'
#' @slot refsp Object of class \code{"character"}, a single string 
#' specifying a reference species used to sort the \code{"phylo"} object.
#' @slot ogids Object of class \code{"character"}, a vector of strings 
#' listing OGs evaluated by the Bridge algorithm.
#' @slot tree Object of class \code{"phylo"}, a given species tree.
#' @slot spbranches Object of class \code{"data.frame"}, a data frame 
#' listing branches of the \code{"phylo"} object.
#' @slot orthocount Object of class \code{"data.frame"}, a data.frame with 
#' results from the 'Bridge' constructor. 
#' @slot orthostats Object of class \code{"data.frame"}, a data.frame with 
#' results from the 'runBridge' function.
#' @slot branchprobs Object of class \code{"matrix"}, a numeric matrix with 
#' results from the 'runBridge' function.
#' @slot bridgeprobs Object of class \code{"matrix"}, a numeric matrix with 
#' results from the 'runBridge' function.
#' @slot rootscore Object of class \code{"matrix"}, a numeric matrix with 
#' results from the 'runBridge' function.
#' @slot branches Object of class \code{"numeric"}, a numeric vector with 
#' results from the 'runBridge' function.
#' @slot roots Object of class \code{"numeric"}, a numeric vector with 
#' results from the 'runBridge' function.
#' @slot misc A list with intermediate objects for downstream methods.
#' @slot status Object of class \code{"character"}, a string specifying the 
#' status of the Bridge object based on the available methods.
#'
#' @method runBridge \code{\link{runBridge}}
#' @method getBridge \code{\link{getBridge}}
#' @method plotBridgeTree \code{\link{plotBridgeTree}}
#' @method plotBridgeStats \code{\link{plotBridgeStats}}
#' @aliases Bridge-class
#' @return An S4 class object.
#' @section Constructor:
#' see \code{\link{newBridge}} constructor.
#' @exportClass Bridge
#'
setClass(
    "Bridge",
    representation(
        refsp = "character",
        ogids = "character",
        tree = "ANY",
        spbranches = "data.frame",
        orthocount = "data.frame",
        orthostats = "data.frame",
        branchprobs = "matrix",
        bridgeprobs = "matrix",
        rootscore = "matrix",
        branches = "numeric",
        roots = "numeric",
        misc = "list",
        status = "character"
    ),
    prototype = list(
        refsp = character(),
        ogids = character(),
        tree = "ANY",
        spbranches = data.frame(),
        orthocount = data.frame(),
        orthostats = data.frame(),
        branchprobs = matrix(0),
        bridgeprobs = matrix(0),
        rootscore = matrix(0),
        branches = numeric(),
        roots = numeric(),
        misc = list(),
        status = character()
    )
)
setValidity("Bridge", function(object) {
    if (!is(object@tree, "phylo")){
        return("slot 'tree' should contain a 'phylo' class object.")
    }
    slot_lengths <- c(
        Ntip(object@tree), 
        nrow(object@spbranches), 
        nrow(object@orthocount))
    if (length(unique(slot_lengths)) != 1){
        msg <- paste0("n. rows of slots 'spbranches' and 'orthocount' ", 
            "differ from n. tips of slot 'tree'.")
        return(msg)
    }
    TRUE
})
