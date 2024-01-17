
#-------------------------------------------------------------------------------
#' @title Simulating Random Orthologous Groups (ROGs)
#'
#' @description The \code{simulateRogs} function generates random orthologous 
#' groups by modeling gain and loss events within a given 'phylo' tree object.
#' The simulation starts by assigning a random root to each ROG. Following 
#' that, a vertical heritage is propagated through the tree, extending up to a
#' designated reference species. The resulting structure is then subjected to 
#' random gain and loss events. The primary goal of the simulation is to assess
#' variations from ideal vertical heritage patterns that include the reference 
#' species.
#'
#' @param phyloTree Object of class \code{phylo}.
#' @param n.rogs A single integer (>0) specifying the number of ROGs to be 
#' simulated.
#' @param gain A single numerical value in [0, 1]. Probability of gene gain
#' events, modeled by the \code{\link{rbinom}} funtion.
#' @param loss A single numerical value in [0, 1]. Probability of gene loss 
#' events, modeled by the \code{\link{rbinom}} funtion.
#' @param event.trials A single integer (>0) specifying a fixed number of 
#' trials associated with each possible outcome (success or failure)  
#' in modeling gain and loss events. Each trial is independent, so the total
#' number of random events increases with more trials.
#' @param interaction.trials A single integer (>0) specifying a fixed number 
#' of trials associated with the event order. Each interaction trial 
#' randomizes the event order, so the event types interact along the 
#' simulation.
#' @param ref.sp A single string specifying which \code{phyloTree} tip will 
#' be used as the reference species. When ref.sp = 'top_tip', the reference 
#' species will be the one at the top of the provided \code{phylo} object; and
#' when ref.sp = 'max_space', the reference species will be the one that 
#' maximizes the number the possible roots. The reference species will 
#' establish a "root space", encompassing all ancestors of the reference
#' species in the tree.
#' @param return.list A single logical value specifying to return results
#' into a \code{list} object. In this case, the list will include all
#' data objects required to run the \code{\link{newBridge}} constructor.
#' @param verbose A single logical value specifying to display detailed 
#' messages (when \code{verbose=TRUE}) or not (when \code{verbose=FALSE}).
#' @return Either a \code{Bridge} object or a list with data objects required 
#' to build a \code{Bridge} object.
#' @author Sysbiolab.
#' @seealso \code{\link{newBridge}}
#' @examples
#' 
#' # Make a random tree with 'ape'
#' library(ape)
#' phyloTree <- rtree(n=25, br=NULL)
#' 
#' # Simulate 10 ROGs
#' rogs <- simulateRogs(phyloTree, n.rogs = 10, return.list = TRUE)
#' 
#' @name simulateRogs
#' @importFrom stats rbinom
#' @aliases simulateRogs
#' @export
#' 
simulateRogs <- function(phyloTree, n.rogs = 100, gain = 0.2, loss = 0.4, 
    event.trials = 5, interaction.trials = event.trials * 2, 
    ref.sp = c("max_space", "top_tip"), 
    return.list = FALSE, verbose = TRUE){
    
    .validate.args("singleLogical", "verbose", verbose)
    if (verbose) message("-Simulating Random Orthologous Groups...")
    
    #--- validate args
    tree <- .validate.bridge.args(name="phyloTree", para=phyloTree)
    .validate.args("singleInteger", "n.rogs", n.rogs)
    .validate.args("singleNumber", "gain", gain)
    .validate.args("singleNumber", "loss", loss)
    .validate.args("singleInteger", "event.trials", event.trials)
    .validate.args("singleInteger", "interaction.trials", interaction.trials)
    .validate.args("singleLogical", "return.list", return.list)
    ref.sp <- match.arg(ref.sp)
    #--- validate arg values
    if(n.rogs < 1) stop("'n.rogs' should be >=1")
    if(gain < 0 | gain >1) stop("'gain' should be in [0,1]")
    if(loss < 0 | loss >1) stop("'loss' should be in [0,1]")
    if(event.trials < 1) stop("'event.trials' should be >=1")
    if(interaction.trials < 1) stop("'interaction.trials' should be >=1")
    #--- get refsp
    if(ref.sp=="top_tip" || is(tree,"symmtree")){
        top.tip <- tree$edge[Nedge(tree), 2]
        refsp <- tree$tip.label[top.tip]
    } else {
        lcas <- mrca(tree)
        nlcas <- unlist(lapply(apply(lcas, 1, table), length))
        refsp <- names(nlcas)[which.max(nlcas)]
        #---rotate tree to place refsp at the top
        tree <- .rotatePhyloTree(tree, refsp)
    }
    
    #--- get branches map
    spbranches <- .spBranches(phyloTree = tree, refsp = refsp)
    
    #--- assign a random root to each rog and propagate heritage to the tree
    x.rogs <- .expand.nrogs(n.rogs, loss, event.trials, interaction.trials)
    simdata <- .simulate.heritage(spbranches, x.rogs)
    refroots <- .get.rogroots(spbranches, simdata)
    
    #--- simulate HGTs, gains, and losses
    simdata <- .simulate.gain.loss(simdata, gain, loss, event.trials, 
        interaction.trials)
    sim <- .filter.rogs(simdata, refroots, n.rogs)
    simdata <- sim$data; refroots <- sim$refroots
    rogdata <- .get.rogs(simdata, verbose)
    
    #--- return results
    pars <- list(n.rogs = n.rogs, gain = gain, 
        loss = loss, event.trials = event.trials, 
        interaction.trials = interaction.trials,
        ref.sp = ref.sp)
    if(return.list){
        rogs <- list(data = rogdata, tree = tree, refsp = refsp, 
            ref_roots = refroots, pars = pars)
        return(rogs)
    } else {
        gbr <- newBridge(ogdata=rogdata, phyloTree=tree, 
            refsp=refsp, verbose=verbose)
        rogs <- list(ref_roots = refroots, pars = pars)
        gbr <- .updateMisc(gbr, rogs, "rogs")
        gbr <- .updateStatus(gbr, "Simulation")
        return(gbr)
    }
}
#-------------------------------------------------------------------------------
.get.rogroots <- function(spbranches, simdata){
    brs <- spbranches$branch
    rog_roots <- vapply(colnames(simdata), function(rog){
        lcas <- brs[simdata[,rog]>0]
        max(lcas)
    }, numeric(1))
    return(rog_roots)
}
#-------------------------------------------------------------------------------
.simulate.heritage <- function(spbranches, n.rogs){
    roots <- spbranches$branch
    lcas <- unique(roots)
    sim <- vapply(seq_len(n.rogs), function(i){
        rt <- sample(lcas, size = 1)
        as.numeric(roots<=rt)
    },numeric(length(roots)))
    colnames(sim) <- paste0("ROG",seq_len(n.rogs))
    rownames(sim) <- rownames(spbranches)
    return(sim)
}
#-------------------------------------------------------------------------------
.simulate.gain.loss <- function(simdata, gain, loss, 
    event.trials, interaction.trials){
    simfuns <- list(f1=.simulate.gain, f2=.simulate.loss)
    simargs <- list(a1=gain/interaction.trials, a2=loss/interaction.trials)
    gl_game <- simdata
    gl_game[] <- 0
    total.trials <- event.trials + interaction.trials - 1
    for(i in seq_len(total.trials)){
        j <- sample(c(1, 2))
        gl_game <- do.call(simfuns[[j[1]]], list(gl_game, simargs[[j[1]]]) )
        gl_game <- do.call(simfuns[[j[2]]], list(gl_game, simargs[[j[2]]]) )
    }
    simdata <- simdata + gl_game
    simdata[simdata<0] <- 0
    return(simdata)
}
.simulate.gain <- function(gl_game, gain){
    nr <- nrow(gl_game)
    nc <- ncol(gl_game)
    gains <- vapply(seq_len(nc), function(i){
        rbinom(n = nr, size = 1, prob = gain)
    }, numeric(nr))
    gl_game <- gl_game + gains
    return(gl_game)
}
.simulate.loss <- function(gl_game, loss){
    nr <- nrow(gl_game)
    nc <- ncol(gl_game)
    losses <- vapply(seq_len(nc), function(i){
        rbinom(n = nr, size = 1, prob = loss)
    }, numeric(nr))
    gl_game <- gl_game - losses
    return(gl_game)
}
#-------------------------------------------------------------------------------
.get.rogs <- function(simdata, verbose){
    rog_ids <- colnames(simdata)
    ssp_ids <- rownames(simdata)
    rogs <- NULL
    if (verbose) pb <- txtProgressBar(style = 3)
    for(i in seq_along(rog_ids)){
        if (verbose) setTxtProgressBar(pb, i / length(rog_ids))
        id <- rog_ids[i]
        its <- simdata[,id]
        rog <- lapply(seq_along(its), function(j){
            n <- its[j]
            if(n>0){
                sp <- paste0("sp", ssp_ids[j])
                pts <- paste0("prot", seq_len(its[j]))
                pts <- paste(rep(sp, its[j]), pts, id, sep=".")
                pts <- cbind(protein_id=pts, ssp_id=ssp_ids[j], og_id=id)
            } else {
                pts <- character()
                pts <- cbind(protein_id=pts, ssp_id=pts, og_id=pts)
            }
            return(pts)
        })
        rog <- do.call(rbind, rog)
        rog <- data.frame(rog)
        rownames(rog) <- NULL
        rogs <- rbind(rogs, rog)
    }
    if (verbose) close(pb)
    rownames(rogs) <- NULL
    return(rogs)
}
#-------------------------------------------------------------------------------
# Note: ".expand.nrogs" sets more rogs to be created, so the 'refsp'
# will be represented in at least 'n.rogs' after loss events. However,
# this may not be achieved if 'loss' and 'trials' combination
# simulates a large number of losses; in that case, the 'refsp' may be 
# represented in '< n.rogs'.
.expand.nrogs <- function(n.rogs, loss, event.trials, interaction.trials){
    trials <- event.trials + interaction.trials
    p <- 1 - (1 - loss)^trials
    x.rogs <- n.rogs/2 + ( n.rogs * (1 + 2 * p) )
    x.rogs <- ceiling(x.rogs)
    return(x.rogs)
}
.filter.rogs <- function(simdata, refroots, n.rogs){
    sel <- sort(simdata[1,] > 0, decreasing = TRUE)
    sel <- names(sel)[seq_len(n.rogs)]
    simdata <- simdata[, sel, drop=FALSE]
    refroots <- refroots[sel]
    newnames <- paste0("ROG", seq_len(n.rogs))
    names(refroots) <- colnames(simdata) <- newnames
    sim <- list(data=simdata, refroots=refroots)
    return(sim)
}

#-------------------------------------------------------------------------------
#' @title Symmetric Branching Tree Generator
#'
#' @description The \code{symmetricBranching} function generates random binary 
#' trees with branches exhibiting various symmetry patterns, all mapping to 
#' the ancestral nodes of a reference tip. This generator is intended for 
#' assessing potential biases in null distributions.
#' 
#' @param n.branches A single integer (>=2) specifying the number of branches
#' to be simulated.
#' @param n.tips A single integer (>=1) or vector of integers specifying the 
#' number of tips in each branch. When 'n.tips' is a vector, its length 
#' should be equal to 'n.branches'.
#' @return A \code{phylo} object.
#' @author Sysbiolab.
#' @seealso \code{\link{simulateRogs}}
#' @examples
#' # Make a tree with 3 branches, 2 tips each, all branches
#' # anchored to the ancestral nodes of a single reference tip.
#' phy <- symmetricBranching(n.branches = 3 , n.tips = 2)
#' # plot(phy, type = "cladogram")
#' 
#' @name symmetricBranching
#' @importFrom ape rtree bind.tree stree
#' @aliases symmetricBranching
#' @export
#' 
symmetricBranching <- function(n.branches = 10, n.tips = seq(1, n.branches)){
    .validate.args("singleInteger", "n.branches", n.branches)
    .validate.args("allInteger", "n.tips", n.tips)
    if(any(n.branches < 1)) stop("'n.branches' should be >=2")
    if(any(n.tips < 1)) stop("'n.tips' should be >=1")
    if(length(n.tips) == 1){
        n.tips <- rep(n.tips, n.branches)
    } else {
        msg <- paste0("When 'n.tips' is a vector, its length should be\n",
            "equal to 'n.branches'; otherwise it will be recycled.")
        if(length(n.tips) > n.branches){
            n.tips <- n.tips[seq_len(n.branches)]
            warning(msg, call. = FALSE)
        } else if(length(n.tips) < n.branches){
            n.tips <- rep(n.tips, length.out = n.branches)
            warning(msg, call. = FALSE)
        }
    }
    tree <- .symmBranching(n.branches, n.tips)
    class(tree) <- c(class(tree), "symmtree")
    return(tree)
}
.symmBranching <- function(n.branches, n.tips){
    tree <- .make.branches(n.branches)
    n.tips <- rev(n.tips); flag <- "rm"
    tips <- .get.tips(c(n.tips, 1), flag)
    for(i in seq_len(length(tips))){
        tree <- bind.tree(tree, tips[[i]], 1)
    }
    tree <- drop.tip(tree, which(tree$tip.label==flag))
    labs <- rev(seq_len(Ntip(tree)-1))
    tree$tip.label <- c(paste0("t", labs),"ref")
    return(tree)
}
.make.branches <- function(n.branches) {
    type <- ifelse(n.branches > 1, "left", "balanced")
    tree <- stree(n.branches + 1, type = type)
    return(tree)
}
.get.tips <- function(n.tips, flag = "rm"){
    tips <- lapply(n.tips, function(n){
        if(n>=2){
            tr <- stree(n + 1,"left")
        } else {
            tr <- rtree(n + 1, br = NULL)
        }
        tr$tip.label[1] <- flag
        tr
    })
    return(tips)
}

################################################################################
### Package documentation
################################################################################
#' @keywords internal
#' @title GeneBridge: large-scale evolutionary analysis of orthologous groups
#'
#' @description
#' GeneBridge is designed to assess the distribution of orthologous groups in 
#' a given species tree. It implements the Bridge algorithm to determine the 
#' evolutionary root of genes in large-scale evolutionary analysis.
#'
#' @details
#'
#' \tabular{ll}{
#' Package: \tab GeneBridge\cr
#' Type: \tab Software\cr
#' License: \tab GPL-3\cr
#' Maintainer: \tab Mauro Castro \email{mauro.a.castro@@gmail.com}\cr
#' }
#'
#' @section Index:
#' \tabular{ll}{
#' \link{newBridge}:
#' \tab Constructor of Bridge-class objects.\cr
#' \link{runBridge}:
#' \tab Run gene root inference with the Bridge algorithm.\cr
#' \link{getBridge}:
#' \tab Accessors for fetching slots from a Bridge object.\cr
#' }
#' Further information is available in the vignettes by typing
#' \code{vignette('GeneBridge')}. Documented topics are also available in
#' HTML by typing \code{help.start()} and selecting the GeneBridge package
#' from the menu.
#'
#'
"_PACKAGE"
#> [1] '_PACKAGE'

################################################################################
### Documentation for 'ogdata'
################################################################################
#' @title A pre-processed dataset for the GeneBridge package
#'
#' @description Dataset used to demonstrate GeneBridge main functions.
#'
#' @format A data frame containing orthologous groups (OGs).
#'
#' @usage data(ogdata)
#'
#' @source This package.
#' 
#' @details
#' 
#' The dataset consists of 4 R objects to be used for demonstration purposes 
#' in the geneplast vignette.
#' 
#' \describe{
#' 
#' \item{ogdata}{
#' A data frame with three columns listing orthology 
#' annotation retrieved from the STRING database (http://string-db.org/), 
#' release 9.1. Column 1 = Ensembl protein ID; column 2 = NCBI species ID; 
#' column 3 = OG ID. Note: This dataset is to be used for demonstration 
#' purposes only as it represents a subset of the STRING database; in order 
#' to reduce the dataset size, orthology annotation was mapped to genome 
#' stability genes (Castro et al.).
#' }
#' \item{sspids}{
#' A data frame containing the species listed in STRING database 
#' (http://string-db.org/), release 9.1. Column 1 = NCBI species ID; 
#' column 2 = NCBI species name;column 3 = species domain (eukaryotes).
#' }
#' \item{ogids}{
#' A vector listing IDs of the OGs available in the 'ogdata' object.
#' }
#' \item{phyloTree}{
#' An phylo-class object for the eukaryotes listed in the STRING database, 
#' release 9.1.
#' }
#' }
#' 
#' @references
#' Franceschini et al. STRING v9.1: protein-protein interaction networks, 
#' with increased coverage and integration. Nucleic Acids Research, 
#' 41(Database issue):D808-15, 2013. doi:10.1093/nar/gks1094.
#' 
#' @docType data
#' @keywords ogdata
#' @name ogdata
#' @aliases ogdata
#' @return A pre-processed dataset.
#' @examples
#' data(ogdata)
NULL

