# Tests for GeneBridge-class methods
test_that("Check Bridge-class methods", {
    data("gene_bridge_data")
    gbr <- newBridge(ogdata=gene_bridge_data$ogdata, 
        phyloTree=gene_bridge_data$phyloTree, 
        refsp="9606", ogids=gene_bridge_data$ogids)
    expect_true(is(gbr, "Bridge"))
})
