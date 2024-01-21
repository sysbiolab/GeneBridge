
#-------------------------------------------------------------------------------
.orthoCount <- function(ogdata, ogids, spids, verbose) {
    if (verbose) pb <- txtProgressBar(style = 3)
    len1 <- length(ogids)
    len2 <- length(spids)
    orthotable <- vapply(seq_len(len1), function(i) {
        if (verbose) setTxtProgressBar(pb, i / len1)
        dt <- ogdata[which(ogdata$og_id == ogids[i]), ]
        vapply(spids, function(sp) {
            sum(dt$ssp_id == sp)
        }, numeric(1))
    }, numeric(len2))
    if (verbose) close(pb)
    rownames(orthotable) <- spids
    colnames(orthotable) <- ogids
    orthotable[orthotable > 0] <- 1
    return(orthotable)
}

#-------------------------------------------------------------------------------
.runBridge <- function(spbranches, branchprobs, threshold) {
    #-- adjust threshold for sampling errors from unbalanced branch sizes;
    #-- each branch error will add to the relative error of the series (rers)
    thr <- threshold
    ntips <- as.numeric(table(spbranches$branch))
    stder <- sqrt( (thr * (1 - thr)) / ntips )
    ## independ. errors assumption doesn't apply
    ## thr <- thr - (thr * stder)/2
    stder[1] <- NA # 1st branch is given
    rers <- vapply(seq_along(stder), function(i){
        mean( stder[i] - stder[-i] , na.rm = TRUE)
    }, numeric(1))
    rers[1] <- 0 # 1st branch is given
    thr <- thr - (thr * rers)/2
    #-- start Bridge
    nroot <- ncol(branchprobs)
    probs <- score <- branchprobs
    roots <- NULL
    for (i in seq_len(nrow(branchprobs))) {
        pbranch <- branchprobs[i, ]
        pbridge <- vapply(seq_len(nroot), function(m){
            pr <- vapply(seq(m, nroot), function(n) {
                # mean(pbranch[seq.int(m, n)])
                # faster than mean:
                sum(pbranch[seq.int(m, n)])/(n - m + 1)
            }, numeric(1))
            max(pr)
        }, numeric(1))
        #-- find k-branches below threshold
        k <- pbridge < thr ## decision boundary
        #-- assign root
        rt <- .assignRoot(k)
        is_bridge <- pbridge > pbranch
        rt <- .unrootBridge(is_bridge, rt)
        #-- assign 'pbridge' a cumulative score (not used for rooting)
        rscore <- .assignRootScore(pbridge)
        #--get results
        probs[i, ] <- pbridge
        score[i, ] <- rscore
        roots <- c(roots, rt)
    }
    names(roots) <- rownames(branchprobs)
    bridgeprobs <- list(probs=probs, score=score, roots=roots)
    return(bridgeprobs)
}

#-------------------------------------------------------------------------------
.assignRoot <- function(k){
    if(any(k)){
        rt <- which(k)[1] - 1
    } else {
        rt <- length(k)
    }
    return(rt)
}

#-------------------------------------------------------------------------------
# m = 5; n = 20
# pb <- c(rep(1, m), rep(0.1, n-m))
# pr <- .assignRootScore(pb)
# plot(seq_along(pr), pr)
# plot(seq_along(pbridge), pbridge)
.assignRootScore <- function(pbridge){
    pb <- c(1, pbridge, 0)
    m <- cumsum(pb)
    m <- (1 + m) / (1 + sum(m) )
    n <- rev( cumsum( rev(1 - pb) ) )
    n <- (1 + n) / ( 1 + sum(n) )
    pr <- m + n
    pr <- pr[ -c(1, length(pr) ) ]
    pr <- pr / max(pr)
    return(pr)
}

#-------------------------------------------------------------------------------
# check if a bridge was not assigned to the root: the root is an 
# "anchor point", not the other way around.
.unrootBridge <- function(is_bridge, rt){
    while(is_bridge[rt] && rt>1){
        rt <- rt - 1
    }
    return(rt)
}

#-------------------------------------------------------------------------------
.permAsymptotic <- function(branchprobs, roots, nPermutations,
    pAdjustMethod, verbose) {
    # get root space
    roots <- roots[rownames(branchprobs)]
    nogs <- nrow(branchprobs)
    nroots <- ncol(branchprobs)
    rext <- (nroots*2)-1
    rootspace <- t(vapply(seq_len(nogs), function(i){
        .getRootSpace(proot = branchprobs[i, ], rt = roots[i])
    }, numeric(rext)))
    # get Dscore
    dscore <- vapply(seq_len(nogs), function(i) {
        .getAdjDscore(rspace = rootspace[i, ])
    }, numeric(1))
    orthostats <- data.frame(Root=roots, Dscore=dscore)
    # function to compute null Dscore from random rootspace
    .nullDscore <- function(rootspace) {
        vapply(seq_len(nogs), function(i) {
            .getAdjDscore(rspace = sample(rootspace[i, ]))
        }, numeric(1))
    }
    # compute null
    if (.isParallel()) {
        cl <- getOption("cluster")
        clusterExport(cl=cl, envir = environment(),
            list=list(".nullDscore", "rootspace", ".getAdjDscore"))
        nulldist <- parSapply(cl, seq_len(nPermutations), function(i) {
            .nullDscore(rootspace)
        }, simplify = FALSE)
    } else {
        if (verbose) pb <- txtProgressBar(style = 3)
        nulldist <- vapply(seq_len(nPermutations), function(i) {
            if (verbose) setTxtProgressBar(pb, i / nPermutations)
            .nullDscore(rootspace)
        }, numeric(nogs))
        if (verbose) close(pb)
    }
    # get null stats
    if(nogs==1){
        nu_m <- mean(nulldist)
        nu_s <- sd(nulldist)
    } else {
        nu_m <- apply(nulldist, 1, mean)
        nu_s <- apply(nulldist, 1, sd)
    }
    # comput z-statistics
    stat <- (orthostats$Dscore - nu_m) / nu_s
    pvals <- pnorm(stat, lower.tail = FALSE)
    # add stats to orthostats
    orthostats$Statistic <- stat
    orthostats$Pvalue <- pvals
    orthostats$AdjPvalue <- p.adjust(pvals, method = pAdjustMethod)
    return(orthostats)
}

#-------------------------------------------------------------------------------
#-- Dscore is adjusted to stabilize the root space as the root 
#-- approaches the boundaries.
.getAdjDscore <- function(rspace) {
    rt <- (length(rspace) + 1)/2
    # get ensemble M
    M <- mean(rspace[seq_len(rt)]) 
    # get ensemble N
    N <- mean(rspace[ seq((rt + 1), length(rspace)) ])
    # get d-score
    D <- M - N
    return(D)
}

#-------------------------------------------------------------------------------
.getRootSpace <- function(proot, rt) {
    # adjust the root space by adding an 'outgroup', 
    # pseudo counts that extends boundary conditions and
    # place the root at the center of 'rspace'
    names(proot) <- NULL
    len <- length(proot)
    m <- len - rt
    n <- (len - m) - 1
    outm <- 1 - mean(proot[ seq((rt + 1), len) ])
    outn <- 1 - mean(proot[seq_len(rt)])
    rspace <- c(rep(outm, m), proot, rep(outn, n))
    return(rspace)
}

#-------------------------------------------------------------------------------
.branchProbs <- function(orthocount, penalty){
    .getprobs <- function(x, penalty) {
        loss <- length(x) - sum(x)
        gain <- sum(x) * penalty
        gain / (gain + loss)
    }
    branchprobs <- aggregate(orthocount[, -1, drop=FALSE],
        by = list(orthocount[, 1]),
        FUN = .getprobs, penalty = penalty
    )
    rownames(branchprobs) <- branchprobs[, 1]
    branchprobs <- branchprobs[, -1, drop=FALSE]
    branchprobs <- as.matrix(branchprobs)
    branchprobs <- t(branchprobs)
    return(branchprobs)
}

#-------------------------------------------------------------------------------
#--- Return LCAs and branches from a 'phylo' obj
.spBranches <- function(phyloTree, refsp) {
    if (all(phyloTree$tip.label != refsp)) {
        stop("NOTE: 'refsp' should be listed in 'phyloTree'!")
    }
    lcas <- .getLCAs(phyloTree)
    spbranches <- .getBranches(phyloTree, lcas)
    #---set refsp as the 1st branch
    spbranches$branch <- spbranches$branch + 1
    spbranches[refsp, "branch"] <- 1
    #---set obj to 'sspbranches' format
    spbranches <- spbranches[, c("ssp_id", "ssp_name", "branch")]
    idx <- match(spbranches$ssp_id, phyloTree$tip.label)
    spbranches <- spbranches[order(spbranches$branch, -idx), ]
    rownames(spbranches) <- spbranches$ssp_id
    return(spbranches)
}
.getLCAs <- function(phyloTree) {
    ntips <- length(phyloTree$tip.label)
    edgetree <- phyloTree$edge
    tip <- edgetree[nrow(edgetree), 2]
    LCAs <- edgetree[which(edgetree[, 2] == tip), 1]
    while (LCAs[length(LCAs)] > (ntips + 1)) {
        idx <- which(edgetree[, 2] == LCAs[length(LCAs)])
        res <- edgetree[idx, 1]
        LCAs <- c(LCAs, res)
    }
    return(LCAs)
}
.getBranches <- function(phyloTree, LCAs) {
    ntips <- length(phyloTree$tip.label)
    edgetree <- phyloTree$edge
    edges <- seq_len(nrow(edgetree))
    edgetree <- cbind(edgetree, NA)
    for (loc in LCAs) {
        idx1 <- which(edgetree[, 2] == loc)
        idx2 <- which(edges > idx1)
        edgetree[edges[idx2], 3] <- loc
        edges <- edges[-idx2]
    }
    edgetree[is.na(edgetree[, 3]), 3] <- edgetree[1, 1]
    branches <- edgetree[edgetree[, 2] <= ntips, 2:3]
    ids <- phyloTree$tip.label[branches[, 1]]
    if (!is.null(phyloTree$tip.alias)) {
        alias <- phyloTree$tip.alias[branches[, 1]]
    } else {
        alias <- ids
    }
    branches <- data.frame(ids, alias, branches, stringsAsFactors = FALSE)
    colnames(branches) <- c("ssp_id", "ssp_name", "tip", "lca")
    gps <- as.integer(as.factor(rank(-branches$lca)))
    branches$branch <- gps
    rownames(branches) <- branches$ssp_id
    return(branches)
}
.tipOrder <- function(phyloTree) {
    tporder <- phyloTree$edge[, 2]
    tporder <- tporder[tporder <= Ntip(phyloTree)]
    tporder <- as.character(phyloTree$tip.label[tporder])
    return(tporder)
}
.rotatePhyloTree <- function(phyloTree, refsp) {
    top.tip <- which(phyloTree$tip.label == refsp)
    if (phyloTree$edge[Nedge(phyloTree), 2] == top.tip) {
        return(phyloTree)
    }
    mrcas <- mrca(phyloTree)[, refsp]
    phyloTree$edge.length <- rep(1, Nedge(phyloTree))
    tgroup <- dist.nodes(phyloTree)[, top.tip]
    phyloTree$edge.length <- NULL
    tgroup <- tgroup[mrcas]
    names(tgroup) <- names(mrcas)
    ct <- 1; tp <- tgroup
    for (i in sort(unique(tgroup))) {
        tgroup[tp == i] <- ct
        ct <- ct + 1
    }
    tord <- names(sort(tgroup, decreasing = TRUE))
    phyloTree <- rotateConstr(phyloTree, constraint = tord)
    if (phyloTree$edge[Nedge(phyloTree), 2] != top.tip) {
        msg <- "'refsp' could not be placed at the top of the phyloTree."
        stop(msg, call. = FALSE)
    }
    return(phyloTree)
}

#-------------------------------------------------------------------------------
.isParallel <- function() {
    b1 <- "package:snow" %in% search()
    b2 <- tryCatch(
        {
            cl <- getOption("cluster")
            cl.check <- FALSE
            if (is(cl, "cluster")) {
                idx <- seq_len(length(cl))
                cl.check <- vapply(idx, function(i) {
                    isOpen(cl[[i]]$con)
                }, logical(1))
                cl.check <- all(cl.check)
            }
            cl.check
        },
        error = function(e) { FALSE }
    )
    all(c(b1, b2))
}

#-------------------------------------------------------------------------------
.get.simulation <- function(gbr){
    if (!.checkStatus(gbr, "Simulation")) {
        stop("'gbr' object should have 'simulation' data.")
    }
    branches <- getBridge(gbr, what="branches")
    prediction <- getBridge(gbr, what="roots")
    misc <- getBridge(gbr, what="misc")
    reference <- misc$rogs$ref_roots
    reference <- reference[names(prediction)]
    prediction <- factor(prediction, levels = branches)
    reference <- factor(reference, levels = branches)
    simdata <- list(classes=branches, reference=reference, prediction=prediction)
    return(simdata)
}
