#-------------------------------------------------------------------------------
#' @title plotBridgeSimulation
#'
#' @description The \code{plotBridgeSimulation} generated a ggplot2 object for 
#' statistics inferred by the Bridge algorithm.
#'
#' @param gbr A preprocessed \linkS4class{Bridge} class object.
#' @return A ggplot object.
#' @examples
#' # Make a 'phylo' tree
#' tree <- symmetricBranching(n.branches = 10)
#' 
#' # Simulated 100 ROGs
#' gbr <- simulateRogs(phyloTree = tree, n.rogs = 100)
#'
#' # Run the Bridge algorithm
#' gbr <- runBridge(gbr)
#' 
#' # Check predictions
#' plotBridgeSimulation(gbr)
#' 
#' @import methods
#' @importFrom ggplot2 ggplot aes theme labs geom_point scale_fill_manual
#' @importFrom ggplot2 geom_hline annotate theme_minimal theme_bw xlab ylab
#' @importFrom ggplot2 element_blank element_text geom_tile scale_x_discrete
#' @importFrom ggplot2 scale_fill_gradient2 geom_text scale_colour_gradient2
#' @importFrom ggrepel geom_text_repel 
#' @docType methods
#' @rdname plotBridgeSimulation-methods
#' @aliases plotBridgeSimulation
#' @export
setMethod("plotBridgeSimulation", "Bridge", 
    function(gbr){
        
        if (!.checkStatus(gbr, "Simulation")) {
            stop("'gbr' object should have 'simulation' data.")
        }
        if (!.checkStatus(gbr, "Bridge")) {
            stop("'gbr' object should be evaluated by 'runBridge'.")
        }
        
        gg <- .confusion.mtx.simulation(gbr)
        
        plot(gg)
        
        return(invisible(gg))
        
    }
)

.confusion.mtx.simulation <- function(gbr){
    sim <- getBridge(gbr, what="simulation")
    cm_d <- table(Reference=sim$reference, 
        Prediction=sim$prediction)
    cm_f <- cm_d/rowSums(cm_d)
    cm_d <- as.data.frame(cm_d[nrow(cm_d):1,])
    cm_f <- as.data.frame(cm_f[nrow(cm_f):1,])
    cm_d$Counts <- cm_d$Freq
    cm_d$Freq <- cm_f$Freq*100
    cm_d$Counts[cm_d$Counts == 0] <- ""
    cm_d$Freq[cm_d$Freq==0] <- NA
    Counts <- Freq <- Prediction <- Reference <- NULL
    ggp <-  ggplot(data = cm_d, aes(x = Prediction , 
        y =  Reference, fill = Freq)) +
        geom_tile() + scale_x_discrete(position = "top") +
        scale_fill_gradient2(guide = FALSE , 
            low="grey97", mid = "#56B1F7", high = "#132B43",
            midpoint = 50, na.value = 'grey97') +
        geom_text(aes(label = Counts, colour=Freq), size = 2.8) + 
        scale_colour_gradient2(
            low="grey50", mid="grey85", high="white",
            midpoint=50) +
        xlab("Root Predictions") + ylab("Reference") + 
        theme_bw() +
        theme(
            aspect.ratio = 1, 
            axis.text.x = element_text(size = 8),
            axis.text.y = element_text(size = 8),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            legend.position = "none",
            panel.border = element_blank(),
            plot.background = element_blank(),
            axis.line = element_blank(),
        )
    return(ggp)
}

#-------------------------------------------------------------------------------
#' @title plotBridgeStats
#'
#' @description The \code{plotBridgeStats} generated a ggplot2 object for 
#' statistics inferred by the Bridge algorithm.
#'
#' @param gbr A preprocessed \linkS4class{Bridge} class object.
#' @param thr A pvalue threshold.
#' @return A ggplot object.
#' @examples
#' # Load datasets used for demonstration
#' data(ogdata)
#'
#' # Create an object of class 'Bridge' for H. sapiens
#' gbr <- newBridge(ogdata, phyloTree, refsp="9606", ogids=ogids)
#'
#' # Run the Bridge algorithm
#' gbr <- runBridge(gbr)
#' gbr <- runPermutation(gbr, nPermutations=100)
#' plotBridgeStats(gbr)
#' 
#' @import methods
#' @importFrom ggplot2 ggplot aes theme labs geom_point scale_fill_manual
#' @importFrom ggplot2 geom_hline annotate theme_minimal theme_bw xlab ylab
#' @importFrom ggplot2 element_blank element_text geom_tile scale_x_discrete
#' @importFrom ggplot2 scale_fill_gradient2 geom_text scale_colour_gradient2
#' @importFrom ggrepel geom_text_repel 
#' @importFrom stats quantile
#' @docType methods
#' @rdname plotBridgeStats-methods
#' @aliases plotBridgeStats
#' @export
setMethod("plotBridgeStats", "Bridge", 
    function(gbr, thr = 0.01){
        
        if (!.checkStatus(gbr, "Permutation")) {
            stop("'gbr' object should be evaluated by 'runPermutation'.")
        }
        
        .validate.args("singleNumberOrNA", "thr", thr)
        
        if(!is.na(thr)){
            if(thr < 0 | thr >1) stop("'thr' should be in [0,1]")
            gg <- .dscoreStats1(gbr, thr)
        } else {
            gg <- .dscoreStats2(gbr)
        }
        
        plot(gg)
        
        return(invisible(gg))
        
    }
)
.dscoreStats1 <- function(gbr, thr){
    # Get Bridge stats
    res_df <- getBridge(gbr, what="results")
    pars <- getBridge(gbr, what="misc")$pars
    adj.thr <- p.threshold(res_df$Pvalue, thr, pars$pAdjustMethod)
    # Set labels to top10 outliers
    ntop <- 10
    if(nrow(res_df) < 3*ntop) ntop <- nrow(res_df)/3
    qth <- ntop/nrow(res_df)
    res_df$Outliers <- rownames(res_df)
    sort.list(res_df$Statistic)
    outl <- quantile(res_df$Pvalue, c(qth, 1-qth))
    outl <- res_df$Pvalue < outl[1] | res_df$Pvalue > outl[2]
    res_df$Outliers[!outl] <- ""
    res_df$Signif <- as.numeric(res_df$AdjPvalue<thr) + 1
    res_df$Signif <- factor(c("no", "yes")[res_df$Signif], 
        levels = c("yes", "no"))
    # Set colors and a title
    cls1 <- c(yes="green4", no="grey40")
    cls2 <- c("grey40","magenta", "red")
    title <- paste("Phyletic pattern consistency of root\nplacements for n =",
        nrow(res_df),"OGs")
    lab <- paste("Adj.P <", thr)
    adj.thr <- -log10(adj.thr)
    rg <- range(res_df$Dscore)
    rp <- pretty(-log10(res_df$Pvalue))
    if(adj.thr >= mean(rp)){
        vjust <- 1.5
        hjust <- 0.15
        labx <- rg[1]
    } else {
        vjust <- -1
        hjust <- 0.85
        labx <- rg[2]
    }
    # Make a ggplot
    Signif <- Dscore <- Pvalue <- Outliers <- NULL
    gg <- ggplot(res_df, aes(x=Dscore, y = -log10(Pvalue),
        label=Outliers, fill=Signif)) + labs(title= title) + 
        geom_hline(yintercept = adj.thr, color=cls2[3], 
            linewidth = 0.4, linetype = 2) + 
        annotate("text", x=labx, y = adj.thr, label= lab, 
            color=cls2[3], hjust = hjust, vjust = vjust, size = 3) +
        geom_point(size = 1.2, shape = 21, color=cls1[2], alpha=0.5) + 
        scale_fill_manual(values = cls1) +
        theme_minimal() + theme(aspect.ratio = 1, legend.position = "none") +
        geom_text_repel(colour=cls2[2], box.padding = 0.5, 
            max.overlaps = nrow(res_df), size = 2, 
            segment.colour=cls2[1])
    return(gg)
}
.dscoreStats2 <- function(gbr){
    # Get Bridge stats
    res_df <- getBridge(gbr, what="results")
    # Set labels to outliers
    res_df$Outliers <- rownames(res_df)
    outl <- quantile(res_df$AdjPvalue, c(0.05, 0.95))
    outl <- res_df$AdjPvalue < outl[1] | res_df$AdjPvalue > outl[2]
    res_df$Outliers[!outl] <- ""
    # Set colors and a title
    cls <- c("green4", "grey40", "magenta", "red")
    title <- paste("Uncertainty of root placement\nfor n =", 
        nrow(res_df),"OGs")
    # Make a ggplot
    Signif <- Dscore <- AdjPvalue <- Outliers <- NULL
    gg <- ggplot(res_df, aes(x=Dscore, y = -log10(AdjPvalue),
        label=Outliers)) + labs(title= title) + 
        geom_point(size = 2, shape = 21, fill = cls[1], 
            color = cls[2], alpha = 0.5) +
        theme_minimal() + theme(aspect.ratio = 1, legend.position = "none") +
        geom_text_repel(colour=cls[3], box.padding = 0.5, 
            max.overlaps = nrow(res_df), size = 2, 
            segment.colour=cls[2])
    return(gg)
}
p.threshold <- function (pvals, alpha=0.05, method="BH") {
    pvals <- sort(pvals)
    padj <- p.adjust(pvals, method = method)
    thset <- which(padj <= alpha)
    if(length(thset)>0){
        mx1 <- mx2 <- which.max(thset)
        if( mx2 < length(padj) ) mx2 <- mx2 + 1
        th <- (pvals[mx1] + min(pvals[mx2], alpha) ) / 2
    } else {
        th <- min(c(alpha, pvals))
    }
    return(th)
}

#-------------------------------------------------------------------------------
#' @title plotBridgeTree
#'
#' @description The \code{plotBridgeTree} generated a tree diagram for the inferred 
#' evolutionary root of a given OG.
#'
#' @param gbr A preprocessed \linkS4class{Bridge} class object.
#' @param whichOG A single string indicating the OG to be plotted.
#' @param fname A name for the file in which the plot will be generated.
#' @param width A single number specifying plot 'width' in inches.
#' @param height A single number specifying plot 'height' in inches.
#' @param cex.lab A single number specifying expansion factor for tip labels.
#' @param cex.nodes A single number specifying expansion factor for 
#' node symbols.
#' @param adj.tips A vector with two numbers.
#' @param lab.offset A single number to offset labels.
#' @param col.tips A vector (length=2) specifying colors of the tips.
#' @param col.edges A vector (length=2) specifying colors of the edges.
#' @param col.root A single color specifying the color of the inferred root.
#' @param plot.sspnames A single logical value specifying whether 'ssp' names 
#' should be plotted.
#' @param plot.lcas A single logical value specifying whether a species tree 
#' should be generated mapping the positions of all possible roots.
#' @param plot.subtree A single logical value specifying whether a 
#' sub-tree (that includes the root) should be plotted.
#' @param plot.pdf A single logical value specifying to generate a PDF file.
#' @return Either a graphics or a pdf file.
#' @examples
#' # Load datasets used for demonstration
#' data(ogdata)
#'
#' # Create an object of class 'Bridge' for H. sapiens
#' gbr <- newBridge(ogdata, phyloTree, refsp="9606", ogids=ogids)
#'
#' # Run the Bridge algorithm
#' gbr <- runBridge(gbr)
#' gbr <- runPermutation(gbr, nPermutations=100)
#' 
#' # This example plots 'NOG40170' in the phyloTree
#' plotBridgeTree(gbr, whichOG="NOG40170")
#' 
#' @import methods
#' @importFrom grDevices col2rgb colorRampPalette dev.off pdf graphics.off
#' @importFrom graphics legend par plot plot.default 
#' @importFrom graphics points segments strwidth text
#' @importFrom ape node.height node.depth drop.tip
#' @docType methods
#' @rdname plotBridgeTree-methods
#' @aliases plotBridgeTree
#' @export
setMethod("plotBridgeTree", "Bridge", 
    function(gbr, whichOG, fname = "gproot", width = 4.5, height = 6.5,
    cex.lab = 0.3, cex.nodes = 0.6, adj.tips = c(1.2, 0.5), lab.offset = 1.5,
    col.tips = c("green", "blue"), col.edges = c("black", "grey"),
    col.root = "red2", plot.sspnames = TRUE, plot.lcas = FALSE, 
    plot.subtree = FALSE, plot.pdf = FALSE) {
    
    .validate.args("singleString", "fname", fname)
    .validate.args("singleNumber", "width", width)
    .validate.args("singleNumber", "height", height)
    .validate.args("singleNumber", "cex.lab", cex.lab)
    .validate.args("singleNumber", "cex.nodes", cex.nodes)
    .validate.bridge.args(name = "adj.tips", para = adj.tips)
    .validate.args("singleNumber", "lab.offset", lab.offset)
    .validate.colors("allColors", "col.tips", col.tips)
    .validate.colors("allColors", "col.edges", col.edges)
    .validate.colors("singleColor", "col.root", col.root)
    .validate.args("singleLogical", "plot.sspnames", plot.sspnames)
    .validate.args("singleLogical", "plot.subtree", plot.subtree)
    .validate.args("singleLogical", "plot.lcas", plot.lcas)
    .validate.args("singleLogical", "plot.pdf", plot.pdf)
    if (!plot.lcas) {
        if(missing(whichOG)){
            plot.lcas <- TRUE
        } else {
            if (!.checkStatus(gbr, "Permutation")) {
                warning("'gbr' object not evaluated by 'runPermutation'.", 
                    call. = FALSE)
            }
            .validate.args("singleString", "whichOG", whichOG)
        }
    }
    pargs <- list(
        cex.lab = cex.lab, cex.nodes = cex.nodes, adj.tips = adj.tips,
        lab.offset = lab.offset, col.tips = col.tips, 
        col.edges = col.edges, col.root = col.root, 
        plot.sspnames = plot.sspnames, plot.subtree = plot.subtree
    )
    phyloTree <- getBridge(gbr, what = "tree")
    spbranches <- getBridge(gbr, what = "spbranches")
    if(!plot.pdf){
        old.par <- par(no.readonly = TRUE)
        on.exit(par(old.par))
    }
    if (plot.lcas) { 
        #---plot LCAS
        refsp <- spbranches$ssp_id[1]
        refspname <- spbranches$ssp_name[1]
        if (plot.sspnames) {
            phyloTree$tip.label <- phyloTree$tip.alias
        }
        #---get OG root position
        lcas <- .getLCAs(phyloTree)
        #---plot
        if(plot.pdf){
            fname <- paste(fname, "_", refsp, "LCAs.pdf", sep = "")
            pdf(file = fname, width = width, height = height)
        }
        .plot.lcas(phyloTree, pargs = pargs, lcas = lcas, 
            refsp = refsp, refspname = refspname)
        if(plot.pdf){ 
            invisible(dev.off())
            message("PDF file ", fname, " has been generated!", sep = "'")
        }
    } else {
        orthocount <- getBridge(gbr, what = "orthocount")
        results <- getBridge(gbr, what = "results")
        #---get, tree, root of whichOG and other relevant results
        if (!whichOG %in% rownames(results)) {
            msg <- paste(whichOG, "should be listed in the 'gbr' object.")
            stop(msg)
        }
        root <- results[whichOG, "Root"]
        orthocount <- orthocount[ , whichOG, drop = FALSE]
        if(ncol(results)>1){
            lgres <- as.numeric(results[whichOG, 
                c("Dscore", "Pvalue", "AdjPvalue")])
            lgres <- c(format(lgres[1], digits = 3, nsmall=1), 
                format(lgres[2:3], digits = 3, scientific = TRUE))
            lgres <- paste(c("Dscore", "Pvalue", "AdjPvalue"), lgres, sep = " = ") 
        } else {
            lgres <- NULL
        }
        refsp <- spbranches$ssp_id[1]
        #---aline orthocount and spbranches with phyloTree
        spbranches <- spbranches[as.character(phyloTree$tip.label), , 
            drop = FALSE]
        orthocount <- orthocount[as.character(phyloTree$tip.label), , 
            drop = FALSE]
        if (plot.sspnames) {
            rownames(spbranches) <- spbranches$ssp_name
            rownames(orthocount) <- spbranches$ssp_name
            phyloTree$tip.label <- phyloTree$tip.alias
        }
        fname <- paste(fname, "_", whichOG, "_", refsp, "LCAs.pdf", sep = "")
        if (plot.subtree) {
            #---cut tree in the OG root
            msg1 <- "Note: the OG root is placed in a subtree, which "
            msg2 <- "lists only part of the OGs."
            warning(msg1, msg2, call. = FALSE)
            subbranch <- spbranches[spbranches[, refsp] <= root, ]
            idx <- !phyloTree$tip.label %in% rownames(subbranch)
            subPhyloTree <- drop.tip(phyloTree, 
                tip = phyloTree$tip.label[idx])
            subct <- orthocount[subPhyloTree$tip.label, , drop = FALSE]
            #---get OG root position
            rtnode <- subPhyloTree$edge[1, 1]
            #---plot
            if(plot.pdf) pdf(file = fname, width = width, height = height)
            .plot.phylo(subPhyloTree, whichOG = whichOG, 
                orthocount = subct, rtnode = rtnode, lgres = lgres, 
                pargs = pargs, root = root, is.last.root = TRUE)
            if(plot.pdf){
                dev.off()
                message("PDF file ", fname, " has been generated!", sep = "'")
            }
        } else {
            #---plot all tree
            #---get OG root position
            lcas <- .getLCAs(phyloTree)
            if(root>1){
                rtnode <- lcas[root - 1]
            } else {
                rtnode <- phyloTree$edge[nrow(phyloTree$edge), 2]
            }
            is.last.root <- ifelse(root > length(lcas), TRUE, FALSE)
            #---plot
            if(plot.pdf) pdf(file = fname, width = width, height = height)
            .plot.phylo(phyloTree, whichOG = whichOG, orthocount = orthocount, 
                rtnode = rtnode, lgres = lgres, pargs = pargs, 
                root = root, is.last.root = is.last.root)
            if(plot.pdf){
                invisible(dev.off())
                message("PDF file ", fname, " has been generated!", sep = "'")
            }
        }
    }

    }
)

#-------------------------------------------------------------------------------
.plot.phylo <- function(x, whichOG, orthocount, rtnode, lgres, pargs,
    root, is.last.root) {
    #---number of segments below the inferred root
    nsegs <- which(x$edge[, 2] == rtnode)
    if (length(nsegs) == 0) nsegs <- 0
    #---set x,y positions
    yy <- node.height(x, clado.style = TRUE)
    xx <- node.depth(x, method = 1) - 1
    xx <- max(xx) - xx
    Ntip <- length(x$tip.label)
    Nedge <- dim(x$edge)[1]
    #---set plot limits
    y.lim <- c(1, Ntip)
    x.lim <- c(0, NA)
    pin1 <- par("pin")[1]
    strWi <- strwidth(x$tip.label, "inches", cex = pargs$cex.lab)
    xx.tips <- xx[seq_len(Ntip)] * 1.04
    alp <- try( uniroot(function(a) max(a * xx.tips + strWi) - 
                pin1,c(0, 1e+06))$root, silent = TRUE )
    if (is.character(alp)) {
        tmp <- max(xx.tips) * 1.5
    } else {
        tmp <- max(xx.tips + strWi / alp)
    }
    tmp <- tmp + pargs$lab.offset
    x.lim[2] <- tmp
    #---get cols
    col.edges <- pargs$col.edges
    col.root <- pargs$col.root
    col.tips <- rev(pargs$col.tips)
    col.root[2] <- colorRampPalette(c(col.root, "white"))(10)[8]
    col.tips[3] <- colorRampPalette(c(col.tips[1], "white"))(10)[8]
    col.tips[4] <- colorRampPalette(c(col.tips[2], "white"))(10)[8]
    #---plot segments
    par(mai = c(0, 0.1, 0, 0))
    plot.default(0, type = "n", xlim = x.lim, ylim = y.lim,
        xlab = "", ylab = "", axes = FALSE, asp = NA)
    segments(xx[x$edge[, 1]], yy[x$edge[, 1]], xx[x$edge[, 2]],
        yy[x$edge[, 2]],
        col = c(rep(col.edges[2], nsegs), rep(col.edges[1], Nedge - nsegs)), 
        lwd = 1, lty = 1)
    #---plot root
    points(xx[rtnode], yy[rtnode], pch = 23, col = col.root[1],
        bg = col.root[2], cex = pargs$cex.nodes * 1.5)
    if(is.last.root){
        text(xx[rtnode], yy[rtnode], labels = paste("root", root),
            cex = 0.5, pos = 3, offset = 0.5)
    } else {
        text(xx[rtnode], yy[rtnode], labels = paste("root", root),
            cex = 0.5, pos = 2, offset = 0.2)
    }
    #---plot tips
    idx <- orthocount[, 1] + 1
    points(xx[seq_len(Ntip)] + pargs$adj.tips[1] - 0.5, 
        yy[seq_len(Ntip)] + pargs$adj.tips[2] - 0.5,
        pch = 21, cex = pargs$cex.nodes, lwd = 0.5,
        col = col.tips[seq_len(2)][idx], bg = col.tips[3:4][idx]
    )
    #---plot labels
    text(xx[seq_len(Ntip)] + pargs$lab.offset, yy[seq_len(Ntip)], x$tip.label,
        adj = 0, font = 3, srt = 0, cex = pargs$cex.lab, col = "black"
    )
    #---plot legend
    legend("topleft", legend = c(whichOG), cex = 1.0, bty = "n")
    if(!is.null(lgres)){
        legend("bottomleft",
            legend = c("OG presence", "OG absence", "Inferred root", NA, lgres),
            col = c(col.tips[2], col.tips[1], col.root[1], NA, NA, NA, NA),
            pt.bg = c(col.tips[4], col.tips[3], col.root[2], NA, NA, NA, NA),
            pch = c(21, 21, 23, NA, NA, NA, NA), cex = 0.6, pt.cex = 0.8,
            inset = 0.05, bty = "n"
        )
    }
}

#-------------------------------------------------------------------------------
.plot.lcas <- function(x, pargs, lcas, refsp, refspname) {
    #---set x,y positions
    yy <- node.height(x, clado.style = TRUE)
    # yy <- .node.height(x)
    xx <- node.depth(x, method = 1) - 1
    xx <- max(xx) - xx
    Ntip <- length(x$tip.label)
    Nedge <- dim(x$edge)[1]

    #---set plot limits
    y.lim <- c(1, Ntip)
    x.lim <- c(0, NA)
    pin1 <- par("pin")[1]
    strWi <- strwidth(x$tip.label, "inches", cex = pargs$cex.lab)
    xx.tips <- xx[seq_len(Ntip)] * 1.04
    alp <- try(
        uniroot(function(a) {
            max(a * xx.tips + strWi) -
                pin1
        }, c(0, 1e+06))$root,
        silent = TRUE
    )
    if (is.character(alp)) {
        tmp <- max(xx.tips) * 1.5
    } else {
        tmp <- max(xx.tips + strWi / alp)
    }
    tmp <- tmp + pargs$lab.offset
    x.lim[2] <- tmp

    #---get cols
    col.edges <- pargs$col.edges
    col.root <- pargs$col.root
    col.tips <- rev(pargs$col.tips)
    col.root[2] <- "white"
    col.tips[3] <- colorRampPalette(c(col.tips[1], "white"))(10)[8]
    col.tips[4] <- colorRampPalette(c(col.tips[2], "white"))(10)[8]

    #---plot segments
    par(mai = rep(0, 4))
    plot.default(0,
        type = "n", xlim = x.lim, ylim = y.lim,
        xlab = "", ylab = "", axes = FALSE, asp = NA
    )
    segments(xx[x$edge[, 1]], yy[x$edge[, 1]], xx[x$edge[, 2]],
        yy[x$edge[, 2]],
        col = col.edges[2], lwd = 1, lty = 1
    )

    #---plot root
    text(xx[lcas], yy[lcas],
        labels = paste(2:(length(lcas) + 1)), cex = 0.4,
        pos = 3, offset = 0.1, srt = 0
    )

    #---plot tips
    points(xx[seq_len(Ntip)] + pargs$adj.tips[1] - 0.5, yy[seq_len(Ntip)] +
        pargs$adj.tips[2] - 0.5,
    pch = 21, cex = pargs$cex.nodes, lwd = 0.5,
    col = col.tips[1], bg = col.tips[3]
    )

    #---plot labels
    text(xx[seq_len(Ntip)] + pargs$lab.offset, yy[seq_len(Ntip)], x$tip.label,
        adj = 0, font = 3, srt = 0, cex = pargs$cex.lab, col = "black"
    )

    #---plot legend
    if(refspname==refsp){
        lab <- paste("REF: ", refspname, sep = "")
    } else {
        lab <- paste("REF: ", refspname, " (", refsp, ")", sep = "")
    }
    legend("topleft", legend = lab, cex = 0.8, bty = "n")
}
