
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
        msg <- paste0("'", name, "' should be a vector with strings.")
        if (!.all_characterValues(para)) stop(msg, call. = FALSE)
    } else if (check == "allCharacterOrInteger") {
        msg <- paste("'", name, "'should be a vector of strings of integers.")
        if (! (.all_characterValues(para) | .all_integerValues(para) ) ) 
            stop(msg, call. = FALSE)
    } else if (check == "allCharacterOrNa") {
        msg <- paste0("'", name, "' should be a vector with strings.")
        if (!.all_characterValues(para, notNA=FALSE)) stop(msg, call. = FALSE)
    } else if (check == "allBinary") {
        msg <- paste0("'", name, "' should be a vector with binary values.")
        if (!.all_binaryValues(para)) stop(msg, call. = FALSE)
    } else if (check == "allInteger") {
        msg <- paste0("'", name, "' should be a vector with integer values.")
        if (!.all_integerValues(para)) stop(msg, call. = FALSE)
    } else if (check == "singleStringOrInteger") {
        msg <- paste0("'", name, "' should be a single string or integer.")
        if (!.is_singleString(para) && !.is_singleInteger(para))
            stop(msg, call. = FALSE) 
    } else if (check == "singleString") {
        msg <- paste0("'", name, "' should be a single string.")
        if (!.is_singleString(para)) stop(msg, call. = FALSE)
    } else if (check == "singleInteger") {
        msg <- paste0("'", name, "' should be a single integer value.")
        if (!.is_singleInteger(para)) stop(msg, call. = FALSE)
    } else if (check == "singleNumber") {
        msg <- paste0("'", name, "' should be a single numeric value.")
        if (!.is_singleNumber(para)) stop(msg, call. = FALSE)
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
    if (name == "cogdata") {
        if ((!is.matrix(para) && !is.data.frame(para)) || ncol(para) < 3) {
            stop("'cogdata' object should be a data frame with length >=3",
                call. = FALSE
            )
        }
        clpars <- c("protein_id", "ssp_id", "cog_id")
        clname <- tolower(colnames(para))
        if (!all(clpars %in% clname)) {
            stop("'cogdata' colnames should include: ", 
                paste(clpars, collapse = ", "), call. = FALSE
            )
        }
        colnames(para) <- clname
        para <- para[, c(clpars, clname[which(!clname %in% clpars)])]
        para[, 1] <- as.character(para[, 1])
        para[, 2] <- as.character(para[, 2])
        para[, 3] <- as.character(para[, 3])
        for (i in 1:ncol(para)) {
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
        for (i in 1:3) {
            idx <- !is.na(para[, i]) & !para[, i] == ""
            para <- para[which(idx), ]
        }
        if (nrow(para) <= 3) {
            stop("'cogdata' has no useful data!\n", call. = FALSE)
        }
        return(para)
    } else if (name == "cogids") {
        if (is.null(para)) {
            return(para)
        } else {
            b1 <- is.character(para)
            b2 <- is.matrix(para) || is.data.frame(para)
            if (!b1 && !b2) {
                stop("'cogids' should be a vector of data frame.",
                    call. = FALSE
                )
            }
            if (b1) {
                para <- unique(as.character(para))
                para <- para[!is.na(para)]
                para <- para[para != ""]
                para <- sort(para)
                para <- data.frame(cog_id = para, stringsAsFactors = FALSE, 
                    check.names = FALSE)
            } else {
                clpars <- c("cog_id")
                clname <- tolower(colnames(para))
                if (!all(clpars %in% clname)) {
                    stop("'cogids' colnames should include: '", clpars, "'!", 
                        call. = FALSE)
                }
                colnames(para) <- clname
                para <- para[, c(clpars, clname[which(!clname %in% clpars)]), 
                    drop = FALSE]
                para[, 1] <- as.character(para[, 1])
                for (i in 1:ncol(para)) {
                    if (!is.numeric(para[, i]) || !is.integer(para[, i])) {
                        para[, i] <- as.character(para[, i])
                    }
                }
                para <- data.frame(para, stringsAsFactors = FALSE, 
                    check.names = FALSE)
                uni <- unique(para[, 1])
                uni <- uni[!is.na(uni)]
                uni <- uni[uni != ""]
                uni <- sort(uni)
                para <- para[match(uni, para[, 1]), , drop = FALSE]
            }
            rownames(para) <- para[, 1]
            return(para)
        }
    } else if (name == "sspids") {
        if (is.null(para)) {
            return(para)
        } else {
            b1 <- is.character(para) || is.integer(para) || is.numeric(para)
            b2 <- is.matrix(para) || is.data.frame(para)
            if (!b1 && !b2) {
                stop("'sspids' should be a vector of characters or dataframe.",
                    call. = FALSE
                )
            }
            if (b1) {
                para <- unique(as.character(para))
                para <- para[!is.na(para)]
                para <- para[para != ""]
                para <- sort(para)
                para <- data.frame(ssp_id = para, stringsAsFactors = FALSE)
            } else {
                clpars <- c("ssp_id")
                clname <- tolower(colnames(para))
                if (!all(clpars %in% clname)) {
                    stop("'cogdata' colnames should include: ", 
                        paste(clpars, collapse = ", "),
                        call. = FALSE
                    )
                }
                colnames(para) <- clname
                para <- para[, c(clpars, clname[which(!clname %in% clpars)]), 
                    drop = FALSE]
                para[, 1] <- as.character(para[, 1])
                para[, 2] <- as.character(para[, 2])
                for (i in 1:ncol(para)) {
                    if (!is.numeric(para[, i]) || !is.integer(para[, i])) {
                        para[, i] <- as.character(para[, i])
                    }
                }
                para <- data.frame(para, stringsAsFactors = FALSE)
                uni <- unique(para[, 1])
                uni <- uni[!is.na(uni)]
                uni <- uni[uni != ""]
                uni <- sort(uni)
                para <- para[match(uni, para[, 1]), , drop = FALSE]
            }
            rownames(para) <- para[, 1]
            return(para)
        }
    } else if (name == "spbranches") {
        b1 <- is.matrix(para) || is.data.frame(para)
        if (!b1) {
            stop("'spbranches' should be a data.frame with spp branches.",
                call. = FALSE
            )
        }
        clpars <- c("ssp_id", "ssp_name")
        clname <- tolower(colnames(para))
        clNOTpars <- clname[which(!clname %in% clpars)]
        if (!all(clpars %in% clname)) {
            stop("'spbranches' colnames should include: ", 
                paste(clpars, collapse = " AND "),
                call. = FALSE
            )
        }
        colnames(para) <- clname
        para <- para[, c(clpars, clNOTpars), drop = FALSE]
        para[, 1] <- as.character(para[, 1])
        para[, 2] <- as.character(para[, 2])
        for (i in 1:ncol(para)) {
            if (!is.numeric(para[, i]) && !is.integer(para[, i])) {
                para[, i] <- as.character(para[, i])
            }
        }
        para <- data.frame(para, stringsAsFactors = FALSE, check.names = FALSE)
        if (any(duplicated(para[, 1]))) {
            stop("NOTE: 'spbranches' should have unique spp ids!")
        }
        if (any(duplicated(para[, 2]))) {
            stop("NOTE: 'spbranches' should have unique spp names!")
        }
        rownames(para) <- para[, 1]
        return(para)
    } else if (name == "phyloTree") {
        if (!is(para, "phylo")) {
            stop("'phyloTree' should be an object of class 'phylo'")
        }
    } else if (name == "penalty") {
        .validate.args("singleNumber","penalty", para)
        if (para <= 0) {
            stop("'penalty' should be numeric value >0", 
                call. = FALSE)
        }   
    } else if (name == "cutoff") {
        .validate.args("singleNumber","cutoff", para)
        if (para > 1 || para < 0) {
            stop("'cutoff' should be numeric value >=0 and <=1", 
                call. = FALSE)
        }
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
.is_singleNumber <- function(para) {
    (is.integer(para) || is.numeric(para)) &&
        length(para) == 1L && !is.na(para)
}
.is_singleInteger <- function(para) {
    lg <- (is.integer(para) || is.numeric(para)) &&
        length(para) == 1L && !is.na(para)
    if (lg) {
        para <- abs(para)
        lg <- abs(para - round(para)) <= para
    }
    return(lg)
}
.is_singleString <- function(para) {
    is.character(para) && length(para) == 1L && !is.na(para)
}
.is_singleLogical <- function(para) {
    is.logical(para) && length(para) == 1L && !is.na(para)
}
.all_binaryValues <- function(para) {
    all(para %in% c(0, 1, NA))
}
.all_integerValues <- function(para, notNA = TRUE) {
    lg <- is.integer(para) || is.numeric(para) || all(is.na(para))
    if (lg) {
        para <- abs(para)
        lg <- all(abs(para - round(para)) <= para, na.rm=TRUE)
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

