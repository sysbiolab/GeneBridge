
#-------------------------------------------------------------------------------
setGeneric("runBridge", signature = "gbr",
    function(gbr, ...) {
        standardGeneric("runBridge")
    }, package = "GeneBridge"
)

#-------------------------------------------------------------------------------
setGeneric("getBridge", signature = "gbr",
    function(gbr, what = "status") {
        standardGeneric("getBridge")
    }, package = "GeneBridge"
)

#-------------------------------------------------------------------------------
setGeneric("plotBridge", signature = "gbr",
    function(gbr, ...) {
        standardGeneric("plotBridge")
    }, package = "GeneBridge"
)
