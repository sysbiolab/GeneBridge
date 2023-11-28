
#-------------------------------------------------------------------------------
#' @title plotBridge
#'
#' @description The \code{plotBridge} generated a tree diagram for the inferred 
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
#' @param plot.subtree A single logical value specifying whether a 
#' sub-tree (that includes the root) should be plotted.
#' @param plot.lcas A single logical value specifying whether a species tree 
#' should be generated mapping the positions of all possible roots.
#' @return A pdf file.
#' @examples
#' # Load datasets used for demonstration
#' data(cogData)
#'
#' # Create an object of class 'Bridge' for H. sapiens
#' gbr <- newBridge(cogdata, phyloTree, spid="9606", cogids=cogids)
#'
#' # Run the Bridge algorithm
#' gbr <- runBridge(gbr, nPermutations=100)
#' 
#' # This example plots 'NOG40170' in the phyloTree
#' plotBridge(gbr, whichOG="NOG40170")
#' 
#' @import methods
#' @importFrom grDevices col2rgb colorRampPalette dev.off pdf
#' @importFrom graphics legend par plot plot.default 
#' @importFrom graphics points segments strwidth text
#' @docType methods
#' @rdname plotBridge-methods
#' @aliases plotBridge
#' @export
setMethod("plotBridge", "Bridge", 
    function(gbr, whichOG, fname = "gproot", width = 4.5, height = 6.5,
    cex.lab = 0.3, cex.nodes = 0.6, adj.tips = c(1, 0.5), lab.offset = 1.5,
    col.tips = c("green2", "grey"), col.edges = c("black", "grey"),
    col.root = "red", plot.sspnames = TRUE, plot.subtree = FALSE,
    plot.lcas = FALSE) {
        
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
    if (!plot.lcas) {
        if (gbr@status["Rooting"] != "[x]") {
            stop("'gbr' object should be evaluated by 'runBridge'.")
        }
        .validate.args("singleString", "whichOG", whichOG)
    }
    if (!is(phyloTree,"phylo")) {
        stop("'phyloTree' should be an object of class 'phylo'.")
    }
    pargs <- list(
        cex.lab = cex.lab, cex.nodes = cex.nodes, adj.tips = adj.tips,
        lab.offset = lab.offset, col.tips = col.tips, 
        col.edges = col.edges, col.root = col.root, 
        plot.sspnames = plot.sspnames, plot.subtree = plot.subtree
    )
    if (plot.lcas) { 
        #---plot LCAS
        phyloTree <- gbr@tree
        spbranches <- gbr@spbranches
        refsp <- spbranches$ssp_id[1]
        refspname <- spbranches$ssp_name[1]
        if (plot.sspnames) {
            phyloTree$tip.label <- phyloTree$tip.alias
        }
        #---get OG root position
        lcas <- getLCAs(phyloTree)
        #---plot
        fname <- paste(fname, "_", refsp, "LCAs.pdf", sep = "")
        pdf(file = fname, width = width, height = height)
        .plot.lcas(phyloTree,
            pargs = pargs, lcas = lcas, refsp = refsp,
            refspname = refspname
        )
        invisible(dev.off())
        cat("PDF file ", fname, " has been generated!", sep = "'")
    } else {
        #---get, tree, root of whichOG and other relevant results
        phyloTree <- gbr@tree
        orthoroot <- gbr@orthoroot
        if (!whichOG %in% rownames(orthoroot)) {
            msg <- paste(whichOG, "should be listed in the 'gbr' object.")
            stop(msg)
        }
        root <- orthoroot[whichOG, "Root"]
        orthoct <- gbr@orthoct[, whichOG, drop = FALSE]
        lgres <- as.numeric(
            orthoroot[whichOG, c("Dscore", "Pvalue", "AdjPvalue")])
        lgres <- c(format(lgres[1], digits = 3), format(lgres[2:3],
            digits = 3, scientific = TRUE
        ))
        lgres <- paste(c("Dscore", "Pvalue", "AdjPvalue"), 
            lgres, sep = " = ")
        spbranches <- gbr@spbranches
        refsp <- spbranches$ssp_id[1]
        #---aline orthoct and spbranches with phyloTree
        spbranches <- spbranches[as.character(phyloTree$tip.label), , 
            drop = FALSE]
        orthoct <- orthoct[as.character(phyloTree$tip.label), , 
            drop = FALSE]
        if (plot.sspnames) {
            rownames(spbranches) <- spbranches$ssp_name
            rownames(orthoct) <- spbranches$ssp_name
            phyloTree$tip.label <- phyloTree$tip.alias
        }
        fname <- paste(fname, "_", whichOG, "_", refsp, "LCAs.pdf", 
            sep = "")
        if (plot.subtree) {
            #---cut tree in the OG root
            msg1 <- "Note: the OG root is placed in a subtree.\n"
            msg2 <- "Only part of the OG distribution is represented."
            warning(msg1, msg2)
            #---filtra tudo acima da raiz!
            subbranch <- spbranches[spbranches[, refsp] <= root, ]
            idx <- !phyloTree$tip.label %in% rownames(subbranch)
            subPhyloTree <- drop.tip(phyloTree, 
                tip = phyloTree$tip.label[idx])
            subct <- orthoct[subPhyloTree$tip.label, , drop = FALSE]
            #---get OG root position
            rtnode <- subPhyloTree$edge[1, 1]
            #---plot
            pdf(file = fname, width = width, height = height)
            .plot.phylo(subPhyloTree,
                whichOG = whichOG, orthoct = subct, rtnode = rtnode,
                lgres = lgres, pargs = pargs
            )
            dev.off()
            cat("PDF file ", fname, " has been generated!", sep = "'")
        } else {
            #---plot all tree
            #---get OG root position
            lcas <- getLCAs(phyloTree)
            rtnode <- lcas[root - 1]
            #---plot
            pdf(file = fname, width = width, height = height)
            .plot.phylo(phyloTree,
                whichOG = whichOG, orthoct = orthoct, rtnode = rtnode,
                lgres = lgres, pargs = pargs, root = root, lcas = lcas
            )
            invisible(dev.off())
            cat("PDF file ", fname, " has been generated!", sep = "'")
        }
    }
}
)

#-------------------------------------------------------------------------------
.plot.phylo <- function(x, whichOG, orthoct, rtnode, lgres, 
    pargs, root, lcas) {
    #---number of segments below the inferred root
    nsegs <- which(x$edge[, 2] == rtnode)
    if (length(nsegs) == 0) nsegs <- 0
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
        uniroot(
            function(a) max(a * xx.tips + strWi) - pin1,
            c(0, 1e+06)
        )$root,
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
        col = c(
            rep(col.edges[2], nsegs),
            rep(col.edges[1], Nedge - nsegs)
        ), lwd = 1, lty = 1
    )
    #---plot root
    points(xx[rtnode], yy[rtnode],
        pch = 23, col = col.root[1],
        bg = col.root[2], cex = pargs$cex.nodes * 1.3
    )
    text(xx[rtnode], yy[rtnode],
        labels = paste("root", root),
        cex = 0.5, pos = 2, offset = 0.5
    )
    #---plot tips
    idx <- orthoct[, 1] + 1
    points(
        xx[seq_len(Ntip)] + pargs$adj.tips[1] -
            0.5, yy[seq_len(Ntip)] + pargs$adj.tips[2] - 0.5,
        pch = 21, cex = pargs$cex.nodes, lwd = 0.5,
        col = col.tips[seq_len(2)][idx], bg = col.tips[3:4][idx]
    )
    #---plot labels
    text(xx[seq_len(Ntip)] + pargs$lab.offset, yy[seq_len(Ntip)], x$tip.label,
        adj = 0, font = 3, srt = 0, cex = pargs$cex.lab, col = "black"
    )
    #---plot legend
    legend("topleft", legend = c(whichOG), cex = 1.0, bty = "n")
    legend("bottomleft",
        legend = c("OG presence", "OG absence", "Inferred root", NA, lgres),
        col = c(col.tips[2], col.tips[1], col.root[1], NA, NA, NA, NA),
        pt.bg = c(col.tips[4], col.tips[3], col.root[2], NA, NA, NA, NA),
        pch = c(21, 21, 23, NA, NA, NA, NA), cex = 0.6, pt.cex = 0.8,
        inset = 0.05, bty = "n"
    )
}

#-------------------------------------------------------------------------------
# fix a bug in "node.height" from ape
# .node.height <- function(phy){
#   n <- length(phy$tip.label)
#   m <- phy$Nnode
#   N <- dim(phy$edge)[1]
#   phy <- reorder(phy)
#   yy <- numeric(n + m)
#   e2 <- phy$edge[, 2]
#   yy[e2[e2 <= n]] <- seq_len(n)
#   phy <- reorder(phy, order = "postorder")
#   e1 <- phy$edge[, 1]
#   e2 <- phy$edge[, 2]
#   .C(node_height_clado, as.integer(n), as.integer(e1),
#      as.integer(e2), as.integer(N), double(n + m), as.double(yy))[[6]]
# }

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
        labels = paste("root", 2:(length(lcas) + 1)), cex = 0.25,
        pos = 2, offset = 0.1, srt = -45
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
    legend("topleft",
        legend = paste("REF: ", refspname, " (", refsp, ")", sep = ""),
        cex = 0.8, bty = "n"
    )
}
