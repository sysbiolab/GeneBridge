# Unit tests fot OGR-class methods
test_gbr <- function(){
  data("gene_bridge_data")
  gbr <- newBridge(ogdata=gene_bridge_data$ogdata, 
    phyloTree=gene_bridge_data$phyloTree, 
    refsp="9606", ogids=gene_bridge_data$ogids)
  gbr <- runBridge(gbr, verbose=FALSE)
  gbr <- runPermutation(gbr, nPermutations=50, verbose=FALSE)
  res <- getBridge(gbr, what="results")
  checkTrue(is.data.frame(res) && ncol(res)==5)
}
