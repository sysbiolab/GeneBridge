
#-------------------------------------------------------------------------------
.orthoCount <- function(cogdata, cogvec, sspvec, verbose) {
    if (verbose) pb <- txtProgressBar(style = 3)
    len <- length(cogvec)
    orthotable <- sapply(seq_len(len), function(i) {
        if (verbose) setTxtProgressBar(pb, i / len)
        dt <- cogdata[which(cogdata$cog_id == cogvec[i]), ]
        sapply(sspvec, function(sp) {
            sum(dt$ssp_id == sp)
        })
    })
    if (verbose) close(pb)
    rownames(orthotable) <- sspvec
    colnames(orthotable) <- cogvec
    orthotable[orthotable > 0] <- 1
    return(orthotable)
}

#-------------------------------------------------------------------------------
.runBridge <- function(rootprobs, cutoff) {
    # compute roots
    Root <- NULL
    Dscore <- NULL
    D <- NULL
    nroot <- nrow(rootprobs) + 1
    for (i in seq_len(ncol(rootprobs))) {
        proot <- c(rootprobs[, i], 0)
        rt <- NA
        dr <- 1
        while (is.na(rt)) {
            k <- which(proot < cutoff)[1]
            res <- sapply((k + 1):nroot, function(j) {
                mean(proot[k:j])
            })
            bridge <- which(res > cutoff)
            if (length(bridge) > 0) {
                dr <- res[bridge[1]]
                proot[seq_len(k + bridge[1])] <- 1
            } else {
                rt <- k - 1
            }
        }
        proot <- c(rootprobs[, i], 0)
        dscore <- mean(proot[seq_len(rt)]) - mean(proot[(rt + 1):(nroot)])
        Root <- c(Root, rt)
        D <- c(D, dr)
        Dscore <- c(Dscore, dscore)
    }
    orthoroot <- data.frame(Root = Root, D = D, Dscore = Dscore,
        stringsAsFactors = FALSE)
    rownames(orthoroot) <- colnames(rootprobs)
    return(orthoroot)
}

#-------------------------------------------------------------------------------
.runPermutation <- function(orthoroot, rootprobs, nPermutations, verbose) {
    # function to compute Dscore from random rootprobs
    .nullRootProbs <- function(rootprobs, orthoroot) {
        sapply(seq_len(ncol(rootprobs)), function(i) {
            .getAdjDscore(proot = rootprobs[, i], 
            rt = orthoroot[i, "Root"], TRUE)
        })
    }
    # compute null
    if (.isParallel()) {
        cl <- getOption("cluster")
        clusterExport(cl, list(
            ".nullRootProbs", "rootprobs", "orthoroot", ".getAdjDscore"
        ), envir = environment())
        nulldist <- parSapply(cl, seq_len(nPermutations), function(i) {
            .nullRootProbs(rootprobs, orthoroot)
        })
    } else {
        if (verbose) pb <- txtProgressBar(style = 3)
        nulldist <- sapply(seq_len(nPermutations), function(i) {
            if (verbose) setTxtProgressBar(pb, i / nPermutations)
            .nullRootProbs(rootprobs, orthoroot)
        })
        if (verbose) close(pb)
    }
    adjDscore <- sapply(seq_len(ncol(rootprobs)), function(i) {
        .getAdjDscore(rootprobs[, i], rt = orthoroot[i, "Root"])
    })
    # z-transformation
    xmd <- apply(nulldist, 1, mean)
    xsd <- apply(nulldist, 1, sd)
    zscore <- (adjDscore - xmd) / xsd
    pvals <- pnorm(zscore, lower.tail = FALSE)
    return(list(zscore = zscore, pvals = pvals))
}

#-------------------------------------------------------------------------------
.getAdjDscore <- function(proot, rt, perm = FALSE) {
    nroot <- length(proot)
    # pspace <- c(rep(1, nroot),rep(0, nroot))
    pspace <- c(seq_len(nroot * 2)) %% 2
    ct <- 1 + nroot - rt
    mask <- ct:(ct + nroot - 1)
    if (perm) {
        pspace[mask] <- sample(pspace)[mask]
    } else {
        pspace[mask] <- proot
    }
    rt <- nroot
    nroot <- length(pspace)
    mean(pspace[seq_len(rt)]) - mean(pspace[(rt + 1):(nroot)])
}
###---code for testing
# proot <- c(1,1,1,1,1,1,1,1,0,0,0,0,0,0,0);rt=8
# proot <- c(1,0,1,0,1,0,1,0,1,0,1,0,1,0,1);rt=8
# .getAdjDscore(proot,rt)
# null <- sapply(seq_len(1000),function(i){
#   .getAdjDscore(proot=proot,rt=rt,TRUE)
# })
# max(null)
# plot(density(null))
#---
# testdnoise <- function(nroot=25,nPermutations=1000,noise=0){
#   sapply(seq_len(nroot),function(rt){
#     proot<-rep(0,nroot);proot[seq_len(rt)]=1
#     if(noise>0) proot[sample(nroot,noise)]<-sample(c(0,1),noise,replace=TRUE)
#     .getAdjDscore(proot,rt)-.getAdjDscore(proot,rt,TRUE)
#   })
# }
# res <- testdnoise(nroot=25,nPermutations=1000,noise=25); plot(density(res))
## check nulldist across rootprobs for the best distributions
# testddist <- function(nroot=25,nPermutations=1000){
#   testsd <- function(proot,rt,nPermutations){
#     xd <- sapply(seq_len(nPermutations),function(i){
#       .getAdjDscore(proot,rt,TRUE)}
#       )
#     xd <- quantile(xd,probs=1-1/nPermutations,names=FALSE)
#     .getAdjDscore(proot,rt)/xd
#   }
#   sapply(seq_len(nroot),function(rt){
#     proot<-rep(0,nroot);proot[seq_len(rt)]=1
#     testsd(proot,rt,nPermutations)
#   })
# }
# testddist(nroot=25,nPermutations=1000)

#-------------------------------------------------------------------------------
.getprobs <- function(x, penalty) {
    loss <- length(x) - sum(x)
    gain <- sum(x) * penalty
    gain / (gain + loss)
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
#--- Return LCAs and branches from a 'phylo' obj
spBranches <- function(phyloTree, spid) {
    if (all(phyloTree$tip.label != spid)) {
        stop("NOTE: 'spid' should be listed in 'phyloTree'!")
    }
    lcas <- getLCAs(phyloTree)
    spbranches <- getBranches(phyloTree, lcas)
    #---set spid as the 1st branch
    spbranches$branch <- spbranches$branch + 1
    spbranches[spid, "branch"] <- 1
    #---set obj to 'sspbranches' format
    spbranches <- spbranches[, c("ssp_id", "ssp_name", "branch")]
    spbranches <- spbranches[sort.list(spbranches$branch), ]
    colnames(spbranches)[3] <- spid
    rownames(spbranches) <- spbranches$ssp_id
    return(spbranches)
}
getLCAs <- function(phyloTree) {
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
getBranches <- function(phyloTree, LCAs) {
    ntips <- length(phyloTree$tip.label)
    edgetree <- phyloTree$edge
    edges <- seq_len(nrow(edgetree))
    lcas <- rev(which(edgetree[, 2] %in% LCAs))
    edgetree <- cbind(edgetree, NA)
    for (loc in LCAs) {
        idx1 <- which(edgetree[, 2] == loc)
        idx2 <- which(edges > idx1)
        edgetree[edges[idx2], 3] <- loc
        edges <- edges[-idx2]
    }
    edgetree[is.na(edgetree[, 3]), 3] <- edgetree[1, 1]
    branches <- edgetree[edgetree[, 2] <= ntips, 2:3]
    #---
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
tipOrder <- function(phyloTree) {
    tporder <- phyloTree$edge[, 2]
    tporder <- tporder[tporder <= Ntip(phyloTree)]
    tporder <- as.character(phyloTree$tip.label[tporder])
    return(tporder)
}
rotatePhyloTree <- function(phyloTree, spid) {
    tip <- which(phyloTree$tip.label == spid)
    lcas <- mrca(phyloTree)[, spid]
    phyloTree$edge.length <- rep(1, nrow(phyloTree$edge))
    tgroup <- dist.nodes(phyloTree)[, tip]
    tgroup <- tgroup[lcas]
    names(tgroup) <- names(lcas)
    #---
    ct <- 1
    tp <- tgroup
    for (i in sort(unique(tgroup))) {
        tgroup[tp == i] <- ct
        ct <- ct + 1
    }
    #---
    tord <- rev(rank(rev(tgroup), ties.method = "first"))
    #---
    phyloTree <- rotateConstr(phyloTree, names(sort(tord, decreasing = TRUE)))
    tord <- tord[tipOrder(phyloTree)]
    #---
    tgroup <- tgroup[names(tord)]
    phyloTree$tip.group <- tgroup
    #---
    lcas <- lcas[names(tord)]
    # atualiza lca do spid, troca pelo nodo mais proximo
    tp <- phyloTree$edge[, 2]
    lcas[spid] <- phyloTree$edge[which(tp == tip), 1]
    phyloTree$tip.lcas <- lcas
    #---
    return(phyloTree)
}
# getLCAs <- function(phyloTree){
#     edgetree<-phyloTree$edge
#     LCAs<-edgetree[1,1]
#     while( !is.na(LCAs[length(LCAs)]) ){
#         idx<-which(edgetree[,1] == LCAs[length(LCAs)] )[2]
#         res<-edgetree[idx,2]
#         LCAs<-c(LCAs,res)
#     }
#     LCAs<-LCAs[seq_len(length(LCAs)-1)]
#     return(LCAs)
# }
