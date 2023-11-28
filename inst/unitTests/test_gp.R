# Unit tests fot OGR-class methods
test_ogr <- function(){
  data("cogData")
  gbr <- newBridge(cogdata=cogdata, phyloTree=phyloTree, spid="9606", cogids=cogids)
  gbr <- runBridge(gbr, nPermutations=50, verbose=FALSE)
  res <- getBridge(gbr,what="results")
  checkTrue(is.data.frame(res) && ncol(res)==4)
}
