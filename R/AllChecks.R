
#-------------------------------------------------------------------------------
.validate.args <- function(check, name, para) {
    if (check == "numeric_vec") {
        msg <- paste0("'", name, "' should be a numeric vector.")
        if (!.is_numericVector(para)) stop(msg, call. = FALSE)
    } else if (check == "integer_vec") {
        msg <- paste0("'", name, "' should be an integer vector.")
        if (!.is_integerVector(para)) stop(msg, call. = FALSE)
    } else if (check == "character_vec") {
        msg <- paste0("'", name, "' should be a character vector.")
        if (!.is_characterVector(para)) stop(msg, call. = FALSE)
    } else if (check == "numeric_mtx") {
        msg <- paste0("'", name, "' should be a numeric matrix")
        if (!is.numeric(para) || !is.matrix(para)) stop(msg, call. = FALSE)
    } else if (check == "allCharacter") {
        msg <- paste0("'", name, "' should be a vector of strings.")
        if (!.all_characterValues(para)) stop(msg, call. = FALSE)
    } else if (check == "allCharacterOrInteger") {
        msg <- paste("'", name, "'should be a vector of strings or integers.")
        if (! (.all_characterValues(para) | .all_integerValues(para) ) ) 
            stop(msg, call. = FALSE)
    } else if (check == "allCharacterOrNa") {
        msg <- paste0("'", name, "' should be a vector of strings.")
        if (!.all_characterValues(para, notNA=FALSE)) stop(msg, call. = FALSE)
    } else if (check == "allBinary") {
        msg <- paste0("'", name, "' should be a vector of binary values.")
        if (!.all_binaryValues(para)) stop(msg, call. = FALSE)
    } else if (check == "allInteger") {
        msg <- paste0("'", name, "' should be a vector of integer values.")
        if (!.all_integerValues(para)) stop(msg, call. = FALSE)
    } else if (check == "singleStringOrInteger") {
        msg <- paste0("'", name, "' should be a single string or integer.")
        if (!.is_singleString(para) && !.is_singleInteger(para))
            stop(msg, call. = FALSE) 
    } else if (check == "singleString") {
        msg <- paste0("'", name, "' should be a single string.")
        if (!.is_singleString(para)) stop(msg, call. = FALSE)
    } else if (check == "singleStringOrNa") {
        msg <- paste0("'", name, "' should be a single string.")
        if (!.is_singleString(para, notNA = FALSE)) stop(msg, call. = FALSE)
    } else if (check == "singleInteger") {
        msg <- paste0("'", name, "' should be a single integer value.")
        if (!.is_singleInteger(para)) stop(msg, call. = FALSE)
    } else if (check == "singleIntegerOrNa") {
        msg <- paste0("'", name, "' should be a single integer value.")
        if (!.is_singleInteger(para, notNA = FALSE)) stop(msg, call. = FALSE)
    } else if (check == "singleNumber") {
        msg <- paste0("'", name, "' should be a single numeric value.")
        if (!.is_singleNumber(para)) stop(msg, call. = FALSE)
    } else if (check == "singleNumberOrNA") {
        msg <- paste0("'", name, "' should be a single numeric value.")
        if (!.is_singleNumber(para, notNA = FALSE)) stop(msg, call. = FALSE)
    } else if (check == "function") {
        msg <- paste0("'", name, "' should be a function.")
        if (!is.function(para)) stop(msg, call. = FALSE)
    } else if (check == "singleLogical") {
        msg <- paste0("'", name, "' should be a single logical value.")
        if (!.is_singleLogical(para)) stop(msg, call. = FALSE)
    }
}

#-------------------------------------------------------------------------------
.validate.colors <- function(check, name, para) {
    if (name == "singleColor") {
        if (!.is_singleColor(para)) {
            msg <- paste0("'", name, "' should be a single color.")
            stop(msg, call. = FALSE)
        }
    } else if (name == "allColors") {
        if (!.is_color(para)) {
            msg <- paste0("'", name, "' should be a vector with colors.")
            stop(msg, call. = FALSE)
        }
    }
}

#-------------------------------------------------------------------------------
.validate.bridge.args <- function(name, para) {
    if (name == "ogdata") {
        if ((!is.matrix(para) && !is.data.frame(para)) || ncol(para) < 3) {
            stop("'ogdata' object should be a data frame with length >=3",
                call. = FALSE
            )
        }
        clpars <- c("protein_id", "ssp_id", "og_id")
        clname <- tolower(colnames(para))
        colnames(para) <- clname
        if (!all(clpars %in% clname)) {
            msg <- paste0("'ogdata' colnames should include: ",
                paste(clpars, collapse = ", "))
            idx <- c(grep(c("protein"), clname), grep(c("ssp"), clname), 
                grep(c("og"), clname))
            if(any(is.na(idx)) || length(idx)!=3){
                stop(msg, call. = FALSE)
            } else {
                warning(msg, call. = FALSE)
            }
            colnames(para)[idx] <- clpars
            clname <- colnames(para)
        }
        para <- para[, c(clpars, clname[which(!clname %in% clpars)])]
        para[, 1] <- as.character(para[, 1])
        para[, 2] <- as.character(para[, 2])
        para[, 3] <- as.character(para[, 3])
        for (i in seq_len(ncol(para))) {
            if (!is.numeric(para[, i]) || !is.integer(para[, i])) {
                para[, i] <- as.character(para[, i])
            }
        }
        if (is.matrix(para)) {
            para <- as.data.frame(para,
                stringsAsFactors = FALSE,
                check.names = FALSE
            )
        }
        for (i in seq_len(3)) {
            idx <- !is.na(para[, i]) & !para[, i] == ""
            para <- para[which(idx), ]
        }
        if (nrow(para) <= 3) {
            stop("'ogdata' has no useful data!\n", call. = FALSE)
        }
        return(para)
    } else if (name == "phyloTree") {
        if (!is(para, "phylo")) {
            stop("'phyloTree' should be an object of class 'phylo'")
        }
        if(is.null(para$tip.label)){
            stop("'phyloTree' should have 'tip.label' names.")
        }
        para$tip.label <- as.character(para$tip.label)
        if(is.null(para$tip.alias)){
            para$tip.alias <- para$tip.label
        } else {
            para$tip.alias <- as.character(para$tip.label)
        }
        para$edge.length <- NULL
        return(para)
    } else if (name == "pAdjustMethod") {
        .validate.args("singleString","pAdjustMethod", para)
        opt <- c("holm", "hochberg", "hommel", "bonferroni", "BH", 
            "BY", "fdr", "none")
        if (!(para %in% opt)) {
            msg <- paste0("'pAdjustMethod' should be any one of\n",
                paste0(opt, collapse = ", "))
            stop(msg, call. = FALSE)
        }
    } else if (name == "adj.tips") {
        .validate.args("numeric_vec","adj.tips", para)
        if (length(para) != 2) {
            stop("'adj.tips' should be a numeric vector of length 2.", 
                call. = FALSE)
        }
    }
}

#-------------------------------------------------------------------
.is_singleNumber <- function(para, notNA = TRUE) {
    lg <- length(para) == 1L && 
        (is.integer(para) || is.numeric(para) || is.na(para))
    if(lg && notNA) lg <- !is.na(para)
    return(lg)
}
.is_singleInteger <- function(para, notNA = TRUE) {
    lg <- length(para) == 1L && 
        (is.integer(para) || is.numeric(para) || is.na(para))
    if (lg && !is.na(para)) {
        para <- abs(para)
        lg <- (para - round(para)) == 0
    }
    if(lg && notNA) lg <- !is.na(para)
    return(lg)
}
.is_singleString <- function(para, notNA = TRUE) {
    lg <- length(para) == 1L && (is.character(para) || is.na(para))
    if(lg && notNA) lg <- !is.na(para)
    return(lg)
}
.is_singleLogical <- function(para, notNA = TRUE) {
    lg <- length(para) == 1L && (is.logical(para) || is.na(para))
    if(lg && notNA) lg <- !is.na(para)
    return(lg)
}
.all_binaryValues <- function(para, notNA = TRUE) {
    lg <- all(para %in% c(0, 1, NA))
    if(lg && notNA) lg <- !any(is.na(para))
    return(lg)
}
.all_integerValues <- function(para, notNA = TRUE) {
    lg <- is.integer(para) || is.numeric(para) || all(is.na(para))
    if (lg) {
        para <- abs(para)
        lg <- all( (para - round(para)) == 0, na.rm=TRUE)
    }
    if(lg && notNA) lg <- !any(is.na(para))
    return(lg)
}
.all_numericValues <- function(para, notNA = TRUE) {
    lg <- is.numeric(para) || all(is.na(para))
    if(lg && notNA) lg <- !any(is.na(para))
    return(lg)
}
.all_characterValues <- function(para, notNA = TRUE) {
    lg <- is.character(para) || all(is.na(para))
    if(lg && notNA) lg <- !any(is.na(para))
    return(lg)
}
.is_numericVector <- function(para){
    is.vector(para) && .all_numericValues(para)
}
.is_integerVector <- function(para){
    is.vector(para) && .all_integerValues(para)
}
.is_characterVector <- function(para){
    is.vector(para) && .all_characterValues(para)
}
.is_color <- function(x) {
    res <- try(col2rgb(x), silent = TRUE)
    return(!"try-error" %in% class(res))
}
.is_singleColor <- function(para) {
    .is_color(para) && length(para) == 1L && !is.na(para)
}

