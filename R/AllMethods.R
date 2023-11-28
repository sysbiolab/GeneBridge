
#-------------------------------------------------------------------------------
#' @title Constructor of Bridge-class objects
#'
#' @description \code{Bridge} is a constructor of Bridge-class objects.
#'
#' @param cogdata A data frame with COG data.
#' @param phyloTree Object of class \code{phylo}.
#' @param spid A single string or integer specifying the reference species 
#' to be used in the Bridge algorithm. This species should be listed in the 
#' \code{phyloTree}.
#' @param cogids An optional data frame with COG annotation. Alternatively, 
#' this can be a character vector with COG IDs.
#' @param verbose A single logical value specifying to display detailed 
#' messages (when \code{verbose=TRUE}) or not (when \code{verbose=FALSE}).
#' @return A preprocessed \linkS4class{Bridge} class object.
#' @author Sysbiolab.
#' @seealso \code{\link{runBridge}}
#' @examples
#' # Load datasets used for demonstration
#' data(cogData)
#'
#' # Create an object of class 'Bridge' for H. sapiens
#' gbr <- newBridge(cogdata, phyloTree, spid="9606", cogids=cogids)
#' 
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom stats pnorm aggregate ecdf p.adjust sd uniroot reorder
#' @importFrom data.table as.data.table
#' @aliases newBridge
#' @export
newBridge <- function(cogdata, phyloTree, spid, cogids = NULL,
    verbose = TRUE) {
    
    .validate.args("singleLogical", "verbose", verbose)
    if (verbose) message("-Preprocessing input data...\n")
    #--- validate args
    cogdata <- .validate.bridge.args(name = "cogdata", para = cogdata)
    cogids <- .validate.bridge.args(name = "cogids", para = cogids)
    .validate.bridge.args(name = "phyloTree", para = phyloTree)
    .validate.args("singleStringOrInteger", "spid", spid)
    spid <- as.character(spid)
    #--- check phyloTree ids in cogdata
    if (any(!phyloTree$tip.label %in% cogdata$ssp_id)) {
        stop("all id(s) in 'phyloTree' should be listed in 'cogdata'.")
    }
    #---rotate phyloTree to set spid at the top
    tip <- which(phyloTree$tip.label == spid)
    if (phyloTree$edge[Nedge(phyloTree), 2] != tip) {
        phyloTree <- rotatePhyloTree(phyloTree, spid)
        if (phyloTree$edge[Nedge(phyloTree), 2] != tip) {
            warning("spid seems not placed at the top of the phyloTree.")
        }
    }
    #---compute spbranches from phyloTree
    spbranches <- spBranches(phyloTree = phyloTree, spid = spid)
    spbranches <- .validate.bridge.args(name = "spbranches", para = spbranches)
    #---get cogids
    if (is.null(cogids)) {
        cogids <- cogdata$cog_id[cogdata$ssp_id == spid]
        cogids <- unique(as.character(cogids))
        cogids <- cogids[!is.na(cogids)]
        cogids <- cogids[cogids != ""]
        cogids <- data.frame(cog_id = cogids, stringsAsFactors = FALSE,
        row.names = cogids)
    } else {
        if (any(!cogids$cog_id %in% cogdata$cog_id)) {
            stop("one or more 'cogids' not listed in 'cogdata'.")
        }
    }
    #---remove non-usefull data
    cogdata <- cogdata[cogdata$cog_id %in% cogids$cog_id, ]
    checkcogid <- unique(cogdata$cog_id[which(cogdata$ssp_id == spid)])
    idx <- cogids$cog_id %in% checkcogid
    if (any(!idx)) {
        checkcogid <- cogids$cog_id[!idx]
        cogids <- cogids[idx, , drop = FALSE]
        tp <- paste(checkcogid, collapse = ",")
        warning("'spid' not listed in one or more 'cogids':\n", tp, 
        call. = FALSE)
    }
    #---compute orthoct
    orthoct <- .orthoCount(cogdata = cogdata, cogvec = cogids$cog_id, 
        sspvec = spbranches$ssp_id, verbose = verbose)
    orthoct <- data.frame(ssp_id = row.names(orthoct), orthoct, 
        stringsAsFactors = FALSE)
    orthoct <- merge(spbranches, orthoct, by = "ssp_id", all = TRUE)
    rownames(orthoct) <- orthoct[, 1]
    orthoct <- orthoct[, -c(1, 2)]
    orthoct <- orthoct[sort.list(orthoct[, spid]), ]
    #---return new 'Bridge'
    status <- c(Preprocess = "[x]", Rooting = "[ ]")
    object <- new("Bridge", cogids = cogids, tree = phyloTree, 
        spbranches = spbranches, orthoct = orthoct, status=status)
    
    return(object)
}

#-------------------------------------------------------------------------------
#' @title Run gene root inference with the Bridge algorithm
#'
#' @description The \code{runBridge} function assesses the 
#' evolutionary root of a gene based on the distribution of its orthologs.
#'
#' @param gbr A preprocessed \linkS4class{Bridge} class object.
#' @param penalty A single number specifying the penalty used in the 
#' rooting algorithm (see details).
#' @param cutoff A single number in [0,1] specifying the cutoff used in the 
#' BR statistics (see details).
#' @param nPermutations A single integer specifying the number of permutations
#' used to compute a null distribution for the inferred roots 
#' in the species tree.
#' @param pAdjustMethod A single string value specifying the p-value 
#' adjustment method to be used (see 'p.adjust' for details).
#' @param verbose A single logical value specifying to display detailed 
#' messages (when \code{verbose=TRUE}) or not (when \code{verbose=FALSE}).
#' @return A preprocessed \linkS4class{Bridge} class object, including 
#' results from the rooting algorithm.
#' @author Sysbiolab.
#' @details
#' The 'runBridge' function addresses the problem of finding the evolutionary 
#' root of a feature in an phylogenetic tree. The method infers the probability
#' that such feature was present in the Last Common Ancestor (LCA) of a given 
#' lineage. Events like horizontal gene transfer, gene deletion, and de novo 
#' gene formation add noise to vertical heritage patterns. The 'runBridge' 
#' function assesses the presence and absence of the orthologs in the extant
#' species of the phylogenetic tree in order to build a probability 
#' distribution, which is used to identify vertical heritage patterns.
#' 
#' The 'penalty' argument weighs gene gain and loss; penalty=1 indicates equal 
#' probability; penalty > 1 indicates higher probability of gene loss while 
#' penalty < 1 indicates higher probability of gene gain (penalty value should
#' be greater than zero; default penalty=2). 
#' 
#' After the probability distribution is built for a given lineage, then a 
#' rooting algorithm is used to search the LCA that provides the best vertical 
#' heritage pattern. The rooting algorithms finds the 
#' optimum point that splits the probability distribution into two 
#' components: one enriched with the queried feature (supporting vertical 
#' heritage in the lineage) and another with low evidence in favour of the 
#' feature's presence. The cutoff sets the tolerance for the discrimination 
#' between the two components (default cutoff=0.3). Based on the optimization 
#' settings, then the 'runBridge' function computes the inconsistency score 
#' called 'Dscore', which assesses the significance of the inferred root 
#' against a null distribution derived by permutation analysis.
#' 
#' @seealso \code{\link{newBridge}}
#' @examples
#' # Load datasets used for demonstration
#' data(cogData)
#'
#' # Create an object of class 'Bridge' for H. sapiens
#' gbr <- newBridge(cogdata, phyloTree, spid="9606", cogids=cogids)
#' 
#' # Run the Bridge algorithm
#' gbr <- runBridge(gbr, nPermutations=100)
#' res <- getBridge(gbr, what="results")
#'
#' \donttest{
#' # Option: parallel version with SNOW package!
#' library(snow)
#' options(cluster=makeCluster(2, "SOCK"))
#' gbr <- runBridge(gbr, nPermutations=100)
#' stopCluster(getOption("cluster")) 
#' }
#' 
#' @import methods
#' @importFrom snow parSapply parLapply makeCluster stopCluster clusterExport
#' @importFrom ape drop.tip rotateConstr node.height node.depth Nedge Ntip
#' @importFrom ape mrca dist.nodes node_height_clado
#' @docType methods
#' @rdname runBridge-methods
#' @aliases runBridge
#' @export
#'
setMethod(
    "runBridge", "Bridge",
    function(gbr, penalty = 2, cutoff = 0.3, nPermutations = 1000, 
        pAdjustMethod = "bonferroni", verbose = TRUE) {
        
        if (gbr@status["Preprocess"] != "[x]") {
            stop("NOTE: input 'gbr' needs preprocessing!")
        }

        #--- main checks
        .validate.bridge.args(name = "penalty", para = penalty)
        .validate.bridge.args(name = "cutoff", para = cutoff)
        .validate.bridge.args(name = "nPermutations", para = nPermutations)
        .validate.bridge.args(name = "pAdjustMethod", para = pAdjustMethod)
        .validate.args("singleLogical", "verbose", verbose)
        if (.isParallel()) {
            if (verbose) {
                message("-Performing rooting analysis (parallel version)...\n")
            }
        } else {
            if (verbose) message("-Performing rooting analysis...\n")
        }
        
        if (verbose) 
            message("--For ", nrow(gbr@cogids), " orthologous groups...\n")
        
        # run
        rootprobs <- aggregate(gbr@orthoct[, -1],
        by = list(gbr@orthoct[, 1]),
        FUN = .getprobs, penalty = penalty
        )
        rownames(rootprobs) <- rootprobs[, 1]
        rootprobs <- rootprobs[, -1]

        #--- for testing
        # rootprobs[,]<-sample(c(0,1),prod(dim(rootprobs)),replace=TRUE)
        # #nr<-nrow(rootprobs)
        # for(i in seq_len(ncol(rootprobs))){
        #   rootprobs[,i]<-as.numeric( rnorm(nr,sd=(1/(seq_len(nr)))) > (1/nr))
        # }
        # plot(seq_len(nr),rowSums(rootprobs))
        # rootprobs[1,]<-1

        #---
        orthoroot <- .runBridge(rootprobs, cutoff)
        sspCount <- colSums(gbr@orthoct[, -1])

        #--- get stats
        stats <- .runPermutation(orthoroot, rootprobs, nPermutations, verbose)
        stats$adjpvals <- p.adjust(stats$pvals, method = pAdjustMethod)

        #--- update orthoroot
        orthoroot <- cbind(orthoroot,
            Zscore = stats$zscore, Pvalue = stats$pvals,
            AdjPvalue = stats$adjpvals, sspCount = sspCount
        )
        orthoroot <- orthoroot[sort.list(orthoroot$Root, decreasing = FALSE), ]
        orthorootsort <- NULL
        #---
        for (i in unique(orthoroot$Root)) {
            tp <- orthoroot[orthoroot$Root == i, ]
            tp <- tp[sort.list(tp$Pvalue), ]
            orthorootsort <- rbind(if (!is.null(orthorootsort)) 
                orthorootsort, tp)
        }
        gbr@orthoroot <- orthorootsort
        gbr@status["Rooting"] <- "[x]"
        return(gbr)
    }
)

#-------------------------------------------------------------------------------
#' @title Accessors for fetching slots from a Bridge object
#'
#' @description \code{getBridge} retrives information from
#' individual slots available in a Bridge object.
#'
#' @param gbr A preprocessed \linkS4class{Bridge} class object.
#' @param what A single character value specifying which information should 
#' be retrieved from the slots.
#' Options: 'cogids', 'spbranches', 'orthoroot', 'status', 'results', 'tree'.
#' @return Content from slots in the \linkS4class{Bridge} object.
#' @examples
#' # Load datasets used for demonstration
#' data(cogData)
#'
#' # Create an object of class 'Bridge' for H. sapiens
#' gbr <- newBridge(cogdata, phyloTree, spid="9606", cogids=cogids)
#'
#' # Get the 'status' slot in gbr
#' getBridge(gbr, what = 'status')
#'
#' @import methods
#' @docType methods
#' @rdname getBridge-methods
#' @aliases getBridge
#' @export
setMethod(
    "getBridge", "Bridge", function(gbr, what = "status") {
    
        ## -----check input arguments
        .validate.args("singleString", "what", what)
        opts <- c("cogids", "spbranches", "orthoroot", 
            "status", "results", "tree")
        if (!what %in% opts) {
            opts <- paste0(opts, collapse = ", ")
            stop("'what' must be one of:\n", opts, call. = FALSE)
        }
        ## -----get query
        query <- NULL
        if (what == "cogids") {
            query <- gbr@cogids
        } else if (what == "spbranches") {
            query <- gbr@spbranches
        } else if (what == "tree") {
            query <- gbr@tree
        } else if (what == "orthoroot") {
            query <- gbr@orthoroot
        } else if (what == "status") {
            query <- gbr@status
        } else if (what == "results") {
            if (gbr@status["Rooting"] != "[x]") {
                warning("input 'gbr' needs 'Rooting' evaluation!")
                query <- data.frame()
            } else {
                query <- data.frame(
                    Root = gbr@orthoroot$Root,
                    Dscore = round(gbr@orthoroot$Dscore, 2),
                    Pvalue = signif(gbr@orthoroot$Pvalue, 3),
                    AdjPvalue = signif(gbr@orthoroot$AdjPvalue, 3),
                    stringsAsFactors = FALSE
                )
                rownames(query) <- rownames(gbr@orthoroot)
            }
        }
        return(query)
    }
)

#-------------------------------------------------------------------------------
setMethod(
    "show", "Bridge", function(object) {
        cat("A Bridge object:\n")
        message("--status:")
        print(getBridge(object, what = "status"), quote = FALSE)
    }
)
