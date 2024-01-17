
#-------------------------------------------------------------------------------
#' @title Constructor of Bridge-class objects
#'
#' @description \code{Bridge} is a constructor of Bridge-class objects.
#'
#' @param ogdata A data frame with COG data.
#' @param phyloTree Object of class \code{phylo}.
#' @param refsp A single string or integer specifying the reference species 
#' to be used in the Bridge algorithm. The \code{refsp} must be listed in 
#' the \code{tip.label} of the \code{phylo} object.
#' @param ogids An optional character vector listing which OGs should be  
#' evaluated by the Bridge algorithm.
#' @param verbose A single logical value specifying to display detailed 
#' messages (when \code{verbose=TRUE}) or not (when \code{verbose=FALSE}).
#' @return A preprocessed \linkS4class{Bridge} class object.
#' @author Sysbiolab.
#' @seealso \code{\link{runBridge}}
#' @examples
#' # Load datasets used for demonstration
#' data(ogdata)
#'
#' # Create an object of class 'Bridge' for H. sapiens
#' gbr <- newBridge(ogdata, phyloTree, refsp="9606", ogids=ogids)
#' 
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom stats pnorm aggregate ecdf p.adjust sd uniroot reorder
#' @importFrom data.table as.data.table
#' @importFrom ape rotateConstr Nedge Ntip
#' @importFrom ape mrca dist.nodes
#' @aliases newBridge
#' @export
newBridge <- function(ogdata, phyloTree, refsp, ogids = NULL,
    verbose = TRUE) {
    
    .validate.args("singleLogical", "verbose", verbose)
    if (verbose) message("-Preprocessing input data...")
    #--- validate args
    ogdata <- .validate.bridge.args(name = "ogdata", para = ogdata)
    phyloTree <- .validate.bridge.args(name = "phyloTree", para = phyloTree)
    if(!is.null(ogids)) .validate.args("allCharacter", "ogids", ogids)
    .validate.args("singleStringOrInteger", "refsp", refsp)
    refsp <- as.character(refsp)
    #--- check refsp in phyloTree
    if (!refsp %in% phyloTree$tip.label) {
        msg <- "'refsp' not listed in 'phyloTree' tip labels."
        stop(msg, call. = FALSE)
    }
    #--- check phyloTree ids in ogdata
    if (any(!phyloTree$tip.label %in% ogdata$ssp_id)) {
        msg <- "some 'phyloTree' tip labels are no listed in the 'ogdata'."
        warning(msg, call. = FALSE)
    }
    #---rotate phyloTree to set refsp at the top
    phyloTree <- .rotatePhyloTree(phyloTree, refsp)
    #---compute spbranches from phyloTree
    spbranches <- .spBranches(phyloTree = phyloTree, refsp = refsp)
    branches <- sort(unique(spbranches$branch))
    #---get ogids
    if (is.null(ogids)) {
        ogids <- ogdata$og_id[ogdata$ssp_id == refsp]
        ogids <- unique(as.character(ogids))
        ogids <- ogids[!is.na(ogids)]
        ogids <- ogids[ogids != ""]
    } else {
        if (any(!ogids %in% ogdata$og_id)) {
            msg <- "one or more 'ogids' not listed in 'ogdata'."
            stop(msg, call. = FALSE)
        }
    }
    #---remove non-usefull data
    ogdata <- ogdata[ogdata$og_id %in% ogids, ]
    checkogid <- unique(ogdata$og_id[which(ogdata$ssp_id == refsp)])
    idx <- ogids %in% checkogid
    if (any(!idx)) {
        checkogid <- ogids[!idx]
        ogids <- ogids[idx]
        if(length(ogids)==0){
            msg <- "'refsp' not listed in the input OGs."
            stop(msg, call. = FALSE)
        }
        msg <- paste0("'refsp' is not listed in ", 
            length(checkogid), " OGs.")
        warning(msg, call. = FALSE)
    }
    if (verbose) message("--For ", length(ogids), " orthologous groups...")
    #---compute orthocount
    orthocount <- .orthoCount(ogdata = ogdata, ogids = ogids, 
        spids = spbranches$ssp_id, verbose = verbose)
    orthocount <- data.frame(ssp_id = row.names(orthocount), orthocount, 
        stringsAsFactors = FALSE)
    orthocount <- merge(spbranches, orthocount, by = "ssp_id", all = TRUE)
    rownames(orthocount) <- orthocount[, 1]
    orthocount <- orthocount[, -c(1, 2)]
    orthocount <- orthocount[sort.list(orthocount$branch), ]
    #---return new 'Bridge'
    status <- c(Preprocess = "[x]", Bridge = "[ ]", Permutation = "[ ]")
    object <- new("Bridge", 
        refsp = refsp,
        ogids = ogids, 
        tree = phyloTree, 
        spbranches = spbranches, 
        orthocount = orthocount, 
        branches = branches,
        status = status)
    return(object)
}

#-------------------------------------------------------------------------------
#' @title Run gene root inference with the Bridge algorithm
#'
#' @description The \code{runBridge} function assesses the 
#' evolutionary root of a gene based on the distribution of its orthologs.
#'
#' @param gbr A preprocessed \linkS4class{Bridge} class object.
#' @param penalty A single number (>0) specifying the penalty used in the 
#' rooting algorithm (see details).
#' @param threshold A single number in `[0,1)`, specifying a local threshold 
#' used in the Bridge algorithm (see details).
#' @param verbose A single logical value specifying to display detailed 
#' messages (when \code{verbose=TRUE}) or not (when \code{verbose=FALSE}).
#' @return A preprocessed \linkS4class{Bridge} class object, including 
#' results from the Bridge algorithm.
#' @author Sysbiolab.
#' @details
#' The 'Bridge' algorithm addresses the problem of finding the evolutionary 
#' root of a feature in an phylogenetic tree. The method infers the probability
#' that such feature was present in the Last Common Ancestor (LCA) of a given 
#' lineage. Events like horizontal gene transfer, gene deletion, and de novo 
#' gene formation add noise to vertical heritage patterns. The 'runBridge' 
#' function assesses the presence and absence of orthologs in the extant
#' species of a phylogenetic tree in order to build a probability 
#' distribution, which is used to identify vertical heritage patterns.
#' 
#' The 'penalty' argument weighs gene gain and loss; penalty=1 indicates equal 
#' probability; penalty > 1 indicates higher probability of gene loss while 
#' penalty < 1 indicates higher probability of gene gain (penalty value should
#' be greater than zero; default \code{penalty = 2}).
#' 
#' After the probability distribution is calculated for a given lineage, 
#' then a rooting algorithm is used to search the LCA that provides the best  
#' vertical heritage pattern. The rooting algorithm finds the optimal point 
#' that splits the probability distribution into two components: one enriched
#' with the queried feature (supporting vertical heritage in the lineage)
#' and another with low evidence in favor of the feature's presence.
#' The 'threshold' sets the tolerance for the discrimination between the 
#' two components (default \code{threshold = 0.5}).
#' 
#' @seealso \code{\link{newBridge}}
#' @examples
#' # Load datasets used for demonstration
#' data(ogdata)
#'
#' # Create an object of class 'Bridge' for H. sapiens
#' gbr <- newBridge(ogdata, phyloTree, refsp="9606", ogids=ogids)
#' 
#' # Run the Bridge algorithm
#' gbr <- runBridge(gbr)
#' 
#' @import methods
#' @importFrom stats weighted.mean
#' @docType methods
#' @rdname runBridge-methods
#' @aliases runBridge
#' @export
#'
setMethod(
    "runBridge", "Bridge",
    function(gbr, penalty = 2, threshold = 0.5, verbose = TRUE) {
        
        if (!.checkStatus(gbr, "Preprocess")) {
            stop("NOTE: input 'gbr' needs preprocessing.")
        }
        
        #--- main checks
        .validate.args("singleNumber","penalty", penalty)
        .validate.args("singleNumber","threshold", threshold)
        .validate.args("singleLogical", "verbose", verbose)
        if(penalty<0) stop("'penalty' should be >0")
        if(threshold<0 | threshold>1) stop("'threshold' should be in [0,1)")
        if(penalty==0) penalty <- penalty + 1e-10
        if(threshold==0) threshold <- threshold + 1e-10
        if (verbose){
            message("-Performing rooting analysis...")
            message("--For ", length(gbr@ogids), " orthologous groups...")
        }
        
        phyloTree <- getBridge(gbr, what = "tree")
        spbranches <- getBridge(gbr, what = "spbranches")
        orthocount <- getBridge(gbr, what = "orthocount")

        #--- get branch probs
        branchprobs <- .branchProbs(orthocount, penalty)
        
        #--- run Bridge
        bridgeprobs <- .runBridge(spbranches, branchprobs, threshold)
        
        #--- update status
        gbr@branchprobs <- branchprobs
        gbr@bridgeprobs <- bridgeprobs$probs
        gbr@rootscore <- bridgeprobs$score
        gbr@roots <- bridgeprobs$roots
        gbr@misc$pars$penalty <- penalty
        gbr@misc$pars$threshold <- threshold
        gbr <- .updateStatus(gbr, "Bridge")
        return(gbr)
    }
)

#-------------------------------------------------------------------------------
#' @title Check consistency of root placement
#'
#' @description The \code{runPermutation} function assesses the consistency
#' of the roots inferred by the Bridge algorithm.
#'
#' @param gbr A preprocessed \linkS4class{Bridge} class object.
#' @param nPermutations A single integer (>1) specifying the number of
#' permutations used to compute a null distribution for the inferred roots 
#' in the species tree.
#' @param pAdjustMethod A single string specifying the p-value 
#' adjustment method to be used (see 'p.adjust' for details).
#' @param verbose A single logical value specifying to display detailed 
#' messages (when \code{verbose=TRUE}) or not (when \code{verbose=FALSE}).
#' @return A preprocessed \linkS4class{Bridge} class object, including 
#' results from the \code{runPermutation} function.
#' @author Sysbiolab.
#' @details
#' The 'runPermutation' function computes an inconsistency score, 
#' 'Dscore', which is used to assess the significance of the inferred 
#' root against a null distribution obtained by permutation analysis.
#' 
#' @seealso \code{\link{newBridge}}
#' @examples
#' # Load datasets used for demonstration
#' data(ogdata)
#'
#' # Create an object of class 'Bridge' for H. sapiens
#' gbr <- newBridge(ogdata, phyloTree, refsp="9606", ogids=ogids)
#' 
#' # Run the Bridge algorithm
#' gbr <- runBridge(gbr)
#' 
#' # Assess root consistency
#' gbr <- runPermutation(gbr, nPermutations=100)
#' res <- getBridge(gbr, what="results")
#'
#' \donttest{
#' # Option: parallel version with the SNOW package
#' library(snow)
#' options(cluster=makeCluster(2, "SOCK"))
#' gbr <- runPermutation(gbr, nPermutations=100)
#' stopCluster(getOption("cluster")) 
#' }
#' 
#' @import methods
#' @importFrom snow parSapply parLapply makeCluster stopCluster clusterExport
#' @importFrom stats pt var pbinom integrate
#' @docType methods
#' @rdname runPermutation-methods
#' @aliases runPermutation
#' @export
#'
setMethod(
    "runPermutation", "Bridge",
    function(gbr, nPermutations = 1000, pAdjustMethod = "bonferroni", 
        verbose = TRUE) {
        
        if (!.checkStatus(gbr, "Bridge")) {
            stop("'gbr' object should be evaluated by 'runBridge'.")
        }
        
        #--- main checks
        .validate.args("singleInteger","nPermutations", nPermutations)
        .validate.bridge.args(name = "pAdjustMethod", para = pAdjustMethod)
        .validate.args("singleLogical", "verbose", verbose)
        if(nPermutations < 1) stop("'nPermutations' should be >1")
        
        #--- get bridge data
        ogids <- getBridge(gbr, what = "ogids")
        branchprobs <- getBridge(gbr, what = "branchprobs")
        branches <- getBridge(gbr, what = "branches")
        roots <- getBridge(gbr, what = "roots")
        
        #--- run permutation
        if (verbose){
            if (.isParallel()){
                message("-Assessing root uncertainty (parallel version)...")
            } else {
                message("-Assessing root uncertainty...")
            }
            message("--For ", length(ogids), " orthologous groups...")
        }
        orthostats <- .permAsymptotic(branchprobs, roots, nPermutations, 
            pAdjustMethod, verbose)
        
        #--- sort orthostats by Root and Pvalue
        idx <- order(orthostats$Root, orthostats$Pvalue, decreasing = FALSE)
        orthostats <- orthostats[idx, ]
        
        #--- update status
        gbr@orthostats <- orthostats
        gbr@misc$pars$nPermutations <- nPermutations
        gbr@misc$pars$pAdjustMethod <- pAdjustMethod
        gbr <- .updateStatus(gbr, "Permutation")
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
#' Options: 'status', 'results', 'tree', 'ogids', 'orthoroots',
#' 'branches', 'roots'.
#' @return Content from slots in the \linkS4class{Bridge} object.
#' @examples
#' # Load datasets used for demonstration
#' data(ogdata)
#'
#' # Create an object of class 'Bridge' for H. sapiens
#' gbr <- newBridge(ogdata, phyloTree, refsp="9606", ogids=ogids)
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
        opts <- c("status", "results", "tree", "ogids",  
            "branches", "roots", "orthostats",
            "spbranches", "orthocount", "branchprobs", 
            "bridgeprobs", "rootscore", "misc", "simulation")
        if (!what %in% opts) {
            opts <- paste0(opts, collapse = ", ")
            stop("'what' must be one of:\n", opts, call. = FALSE)
        }
        ## -----get query
        query <- NULL
        if (what == "ogids") {
            query <- gbr@ogids
        } else if (what == "spbranches") {
            query <- gbr@spbranches
        } else if (what == "tree") {
            query <- gbr@tree
        } else if (what == "roots") {
            query <- gbr@roots
        } else if (what == "branches") {
            query <- gbr@branches
        } else if (what == "orthostats") {
            query <- gbr@orthostats
        } else if (what == "orthocount") {
            query <- gbr@orthocount
        } else if (what == "branchprobs") {
            query <- gbr@branchprobs
        } else if (what == "bridgeprobs") {
            query <- gbr@bridgeprobs
        } else if (what == "misc") {
            query <- gbr@misc
        } else if (what == "simulation") {
            query <- .get.simulation(gbr)
        } else if (what == "rootscore") {
            query <- gbr@rootscore
        } else if (what == "status") {
            query <- gbr@status
        } else if (what == "results") {
            if (!.checkStatus(gbr, "Permutation")) {
                msg <- paste("For full report, input 'gbr' needs",
                    "'runPermutation' evaluation!")
                warning(msg, call. = FALSE)
                query <- data.frame(Root = gbr@roots,
                    stringsAsFactors = FALSE)
            } else {
                query <- data.frame(
                    Root = gbr@orthostats$Root,
                    Dscore = round(gbr@orthostats$Dscore, 4),
                    Statistic = round(gbr@orthostats$Statistic, 3),
                    Pvalue = signif(gbr@orthostats$Pvalue, 3),
                    AdjPvalue = signif(gbr@orthostats$AdjPvalue, 3),
                    stringsAsFactors = FALSE
                )
                rownames(query) <- rownames(gbr@orthostats)
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

################################################################################
### Internal accessors
################################################################################
.updateMisc <- function(gbr, obj, name) {
    gbr@misc[[name]] <- obj
    return(gbr)
}
.updateStatus <- function(gbr, name, check = TRUE) {
    gbr@status[name] <- ifelse(check, "[x]", "[ ]")
    return(gbr)
}
.checkStatus <- function(gbr, name) {
    status <- getBridge(gbr, what = "status")
    if(name %in% names(status)){
        sts <- status[name] == "[x]"
    } else {
        sts <- FALSE
    }
    sts
}
