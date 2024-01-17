# Unit tests fot OGR-class methods
test_gbr <- function(){
  data("ogdata")
  gbr <- newBridge(ogdata=ogdata, phyloTree=phyloTree, 
    refsp="9606", ogids=ogids)
  gbr <- runBridge(gbr, verbose=FALSE)
  gbr <- runPermutation(gbr, nPermutations=50, verbose=FALSE)
  res <- getBridge(gbr, what="results")
  checkTrue(is.data.frame(res) && ncol(res)==5)
}
