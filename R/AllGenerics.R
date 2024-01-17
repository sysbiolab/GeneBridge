
#-------------------------------------------------------------------------------
setGeneric("runBridge", signature = "gbr",
    function(gbr, ...) {
        standardGeneric("runBridge")
    }, package = "GeneBridge"
)

#-------------------------------------------------------------------------------
setGeneric("runPermutation", signature = "gbr",
    function(gbr, ...) {
        standardGeneric("runPermutation")
    }, package = "GeneBridge"
)

#-------------------------------------------------------------------------------
setGeneric("getBridge", signature = "gbr",
    function(gbr, what = "status") {
        standardGeneric("getBridge")
    }, package = "GeneBridge"
)

#-------------------------------------------------------------------------------
setGeneric("plotBridgeTree", signature = "gbr",
    function(gbr, ...) {
        standardGeneric("plotBridgeTree")
    }, package = "GeneBridge"
)

#-------------------------------------------------------------------------------
setGeneric("plotBridgeStats", signature = "gbr",
    function(gbr, ...) {
        standardGeneric("plotBridgeStats")
    }, package = "GeneBridge"
)

#-------------------------------------------------------------------------------
setGeneric("plotBridgeSimulation", signature = "gbr",
    function(gbr, ...) {
        standardGeneric("plotBridgeSimulation")
    }, package = "GeneBridge"
)
