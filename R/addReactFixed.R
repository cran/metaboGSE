#' A fixed version of sybil::addReact
#' 
#' The function fixes sybil::addReact to avoid deprecated functions rBind and cBind.
#' @param model See \code{\link{addReact}}.
#' @param id See \code{\link{addReact}}.
#' @param met See \code{\link{addReact}}.
#' @param Scoef See \code{\link{addReact}}.
#' @param reversible See \code{\link{addReact}}.
#' @param lb See \code{\link{addReact}}.
#' @param ub See \code{\link{addReact}}.
#' @param obj See \code{\link{addReact}}.
#' @param subSystem See \code{\link{addReact}}.
#' @param gprAssoc See \code{\link{addReact}}.
#' @param reactName See \code{\link{addReact}}.
#' @param metName See \code{\link{addReact}}.
#' @param metComp See \code{\link{addReact}}.
#' @return An object of class \code{modelorg}.
#' @import sybil 
#' @importFrom Matrix Matrix
#' @keywords internal
addReactFixed <- function (model, id, met, Scoef, reversible = FALSE, lb = 0, 
    ub = SYBIL_SETTINGS("MAXIMUM"), obj = 0, subSystem = NA, 
    gprAssoc = NA, reactName = NA, metName = NA, metComp = NA) 
{
    if (!is(model, "modelorg")) {
        stop("needs an object of class modelorg!")
    }
    stopifnot(checkVersion(model))
    if (length(met) != length(Scoef)) {
        stop("arguments 'met' and 'Scoef' must have the same length")
    }
    if (length(id) > 1) {
        stop("add/change one reaction")
    }
    if (((ub > 0) && (lb < 0)) && (!isTRUE(reversible))) {
        Crev <- TRUE
        warning(paste("'lb' and 'ub' are signed different,", 
            "set reversible to 'TRUE'"))
    }
    else {
        Crev <- reversible
    }
    colInd <- match(id, react_id(model))
    addCol <- FALSE
    nCols <- react_num(model)
    if (is.na(colInd)) {
        colInd <- react_num(model) + 1
        addCol <- TRUE
        nCols <- nCols + 1
    }
    rowInd <- match(met, met_id(model))
    newM <- which(is.na(rowInd))
    nRows <- met_num(model)
    nNewRows <- length(newM)
    addRow <- FALSE
    for (i in seq(along = newM)) {
        addRow <- TRUE
        nRows <- nRows + 1
        rowInd[newM[i]] <- nRows
    }
    if ((isTRUE(addCol)) || (isTRUE(addRow))) {
        newmet_num <- met_num(model)
        newmet_id <- met_id(model)
        newmet_name <- met_name(model)
        newmet_comp <- met_comp(model)
        newmet_single <- met_single(model)
        newmet_de <- met_de(model)
        newreact_num <- react_num(model)
        newreact_rev <- react_rev(model)
        newreact_id <- react_id(model)
        newreact_name <- react_name(model)
        newreact_single <- react_single(model)
        newreact_de <- react_de(model)
        newlowbnd <- lowbnd(model)
        newuppbnd <- uppbnd(model)
        newobj_coef <- obj_coef(model)
        newgprRules <- gprRules(model)
        newgenes <- genes(model)
        newgpr <- gpr(model)
        newallGenes <- allGenes(model)
        newrxnGeneMat <- rxnGeneMat(model)
        newsubSys <- subSys(model)
        newS <- S(model)
        newMetAttr <- met_attr(model)
        newReactAttr <- react_attr(model)
        newCompAttr <- comp_attr(model)
        newModAttr <- mod_attr(model)
        if (isTRUE(addRow)) {
            newmet_num <- nRows
            newmet_id <- append(met_id(model), met[newM])
            if (any(is.na(metName))) {
                newmet_name <- append(met_name(model), met[newM])
            }
            else {
                newmet_name <- append(met_name(model), metName[newM])
            }
            if (any(is.na(metComp))) {
                newmet_comp <- append(met_comp(model), rep(NA, 
                  nNewRows))
            }
            else {
                if (is(metComp, "numeric")) {
                  newmet_comp <- append(met_comp(model), metComp[newM])
                }
                else {
                  newmet_comp <- append(met_comp(model), match(metComp[newM], 
                    mod_compart(model)))
                }
            }
            newmet_single <- append(met_single(model), rep(NA, 
                nNewRows))
            newmet_de <- append(met_de(model), rep(NA, nNewRows))
            newRows <- Matrix::Matrix(0, nrow = nNewRows, ncol = react_num(model))
            newS <- rbind(newS, newRows)
            if (ncol(newMetAttr) > 0) {
                newMetAttr[nrow(newMetAttr) + 1:nNewRows, ] <- NA
            }
        }
        if (isTRUE(addCol)) {
            newreact_num <- nCols
            newreact_id <- append(react_id(model), id)
            if (is.na(reactName)) {
                newreact_name <- append(react_name(model), id)
            }
            else {
                newreact_name <- append(react_name(model), reactName)
            }
            newreact_single <- append(react_single(model), NA)
            newreact_de <- append(react_de(model), NA)
            newreact_rev <- append(react_rev(model), Crev)
            newlowbnd <- append(lowbnd(model), lb)
            newuppbnd <- append(uppbnd(model), ub)
            newobj_coef <- append(obj_coef(model), obj)
            newS <- cbind(newS, rep(0, nrow(newS)))
            if (ncol(newReactAttr) > 0) {
                newReactAttr[nrow(newReactAttr) + 1, ] <- NA
            }
            if (any(is.na(subSystem))) {
                ss <- subSys(model)
                if (ncol(ss) == 0) {
                  dim(ss) <- c(nrow(ss) + 1, ncol(ss))
                  newsubSys <- ss
                }
                else {
                  newsubSys <- rbind(ss, rep(FALSE, ncol(subSys(model))))
                }
            }
            else {
                if (is(subSystem, "logical")) {
                  newsubSys <- rbind(subSys(model), subSystem)
                }
                else {
                  nSubsRow <- colnames(subSys(model)) %in% subSystem
                  newsubSys <- rbind(subSys(model), nSubsRow)
                }
            }
            if (ncol(rxnGeneMat(model)) > 0) {
                newrxnGeneMat <- rbind(rxnGeneMat(model), rep(FALSE, 
                  ncol(rxnGeneMat(model))))
            }
            else {
                newrxnGeneMat <- rxnGeneMat(model)
                dim(newrxnGeneMat) <- c(nrow(newrxnGeneMat) + 
                  1, ncol(newrxnGeneMat))
            }
            if ((is.na(gprAssoc)) || (gprAssoc == "")) {
                if ((length(gprRules(model)) > 0)) {
                  newgprRules <- append(gprRules(model), "")
                  newgenes <- append(genes(model), "")
                  newgpr <- append(gpr(model), "")
                }
            }
            else {
                gene_rule <- parseBooleanCopy(gprAssoc)
                geneInd <- match(gene_rule$gene, allGenes(model))
                new_gene <- which(is.na(geneInd))
                if (length(new_gene) > 0) {
                  newallGenes <- append(allGenes(model), gene_rule[["gene"]][new_gene])
                  geneInd <- match(gene_rule[["gene"]], newallGenes)
                  if (ncol(newrxnGeneMat) == 0) {
                    newrxnGeneMat <- Matrix::Matrix(FALSE, nCols, 
                      max(geneInd))
                  }
                  else {
                    for (i in seq(along = gene_rule[["gene"]][new_gene])) {
                      newrxnGeneMat <- cbind(newrxnGeneMat, rep(FALSE, 
                        nrow(newrxnGeneMat)))
                    }
                  }
                }
                newrxnGeneMat[nCols, geneInd] <- TRUE
                newgpr <- append(gpr(model), gprAssoc)
                newgenes <- append(genes(model), list(gene_rule$gene))
                newrule <- gene_rule$rule
                newgprRules <- append(gprRules(model), newrule)
            }
        }
        newS[, colInd] <- 0
        newS[rowInd, colInd] <- Scoef
        if (is(model, "modelorg_irrev")) {
            mod_out <- modelorg_irrev(mod_id(model), mod_name(model))
            irrev(mod_out) <- TRUE
            matchrev(mod_out) <- append(matchrev(model), 0L)
            revReactId <- as.integer(max(irrev2rev(model)) + 
                1)
            irrev2rev(mod_out) <- append(irrev2rev(model), revReactId)
            rev2irrev(mod_out) <- rbind(rev2irrev(model), c(nCols, 
                nCols))
        }
        else {
            mod_out <- modelorg(mod_id(model), mod_name(model))
        }
        mod_desc(mod_out) <- mod_desc(model)
        mod_compart(mod_out) <- mod_compart(model)
        met_num(mod_out) <- as.integer(newmet_num)
        met_id(mod_out) <- newmet_id
        met_name(mod_out) <- newmet_name
        met_comp(mod_out) <- as.integer(newmet_comp)
        met_single(mod_out) <- newmet_single
        met_de(mod_out) <- newmet_de
        react_num(mod_out) <- as.integer(newreact_num)
        react_rev(mod_out) <- newreact_rev
        react_id(mod_out) <- newreact_id
        react_name(mod_out) <- newreact_name
        react_single(mod_out) <- newreact_single
        react_de(mod_out) <- newreact_de
        lowbnd(mod_out) <- newlowbnd
        uppbnd(mod_out) <- newuppbnd
        obj_coef(mod_out) <- newobj_coef
        gprRules(mod_out) <- newgprRules
        genes(mod_out) <- newgenes
        gpr(mod_out) <- newgpr
        allGenes(mod_out) <- newallGenes
        rxnGeneMat(mod_out) <- newrxnGeneMat
        subSys(mod_out) <- newsubSys
        S(mod_out) <- newS
        react_attr(mod_out) <- newReactAttr
        met_attr(mod_out) <- newMetAttr
        comp_attr(mod_out) <- newCompAttr
        mod_attr(mod_out) <- newModAttr
    }
    else {
        mod_out <- model
        react_rev(mod_out)[colInd] <- Crev
        lowbnd(mod_out)[colInd] <- lb
        uppbnd(mod_out)[colInd] <- ub
        obj_coef(mod_out)[colInd] <- obj
        S(mod_out)[, colInd] <- 0
        S(mod_out)[rowInd, colInd] <- Scoef
    }
    check <- validObject(mod_out, test = TRUE)
    if (check != TRUE) {
        msg <- paste("Validity check failed:", check, sep = "\n    ")
        warning(msg)
    }
    return(mod_out)
}


#' Copy of sybil:::.parseBoolean
#'
#' This is the copy of sybil:::.parseBoolean, which is used to fix sybil::addReact
#' @param gprRule See \code{sybil:::.parseBoolean}.
#' @param token See \code{sybil:::.parseBoolean}.
#' @return See \code{sybil:::.parseBoolean}.
#' @import sybil
#' @keywords internal
parseBooleanCopy <- function (gprRule, tokens = "()&|~") 
{
    if (is.na(gprRule) || (gprRule == "")) {
        return(list(gene = "", rule = ""))
    }
    if (grepl("\\s*\\(\\s*\\)\\s*", gprRule)) {
        warning("found empty expression rule: '( )'. check if this is intended.")
        return(list(gene = "", rule = ""))
    }
    gpr <- gsub("and ", "& ", gprRule, ignore.case = TRUE)
    gpr <- gsub("or ", "| ", gpr, ignore.case = TRUE)
    gpr <- gsub("not ", "~ ", gpr, ignore.case = TRUE)
    gpr <- gsub("[", "", gpr, fixed = TRUE)
    gpr <- gsub("]", "", gpr, fixed = TRUE)
    rule <- gpr
    genes_tmp <- strsplit(gpr, paste("[", tokens, "]", sep = ""))
    genes_tmp <- gsub("(^\\s+)|(\\s+$)", "", genes_tmp[[1]], 
                      perl = TRUE)
    not_empty <- which(nchar(genes_tmp) > 0)
    genes <- genes_tmp[not_empty]
    num_genes <- length(genes)
    gene_uniq <- unique(genes)
    newTok <- match(genes, gene_uniq)
    newTok <- sapply(newTok, function(x) paste("x[", x, "]", 
                                               sep = ""))
    for (i in 1:num_genes) {
        rule <- sub(genes[i], newTok[i], rule, fixed = TRUE)
    }
    return(list(gene = gene_uniq, rule = rule))
}
