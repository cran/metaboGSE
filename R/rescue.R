#' Relax model external source
#'
#' This function returns the \emph{generic} form of a given model, where external reactions are rendered bidirectional and unrestricted.
#' @param model An object of class \code{modelorg}.
#' @return The \emph{generic} model.
#' @import sybil methods
#' @keywords internal
generic <- function(model) {
    if (!is(model, "modelorg")) {
        stop("needs an object of class modelorg!")
    }
    if (!checkVersion(model)) {
        stop("model is of wrong version!")
    }
    
    extremum <- max(abs(lowbnd(model)), abs(uppbnd(model)))
    exc <- findExchReact(model)
    if (!is.null(exc)) {
        model <- changeBounds(model, (react_id(exc))[which(uptake(exc) == F)], lb=-extremum, ub=extremum)
        rev  <- react_rev(model)
        names(rev) <- react_id(model)
        rev[(react_id(exc))[which(uptake(exc) == F)]] <- TRUE
        react_rev(model) <- rev
    }
    
    return (model)
}


#' Compute appropriate tolerance for the given model
#'
#' The function computes a tolerance value which is appropriate to a model. This is an empirical procedure whose results must be considered with caution.
#' @param model An object of class \code{modelorg}.
#' @param mc.cores The number of cores to use (at least 2), i.e. at most how many child processes will be run simultaneously. Default: 2.
#' @param solver A character indicating the solver to be used in \code{\link{sybil}}. Default: \code{SYBIL_SETTINGS("SOLVER")}.
#' @param method A character indicating the method to be used in \code{\link{sybil}}. Default: \code{SYBIL_SETTINGS("METHOD")}.
#' @return Tolerance
#' @import sybil parallel
#' @keywords internal
diagnosticTolerance <- function(model, mc.cores = 2,
                                solver = SYBIL_SETTINGS("SOLVER"), 
                                method = SYBIL_SETTINGS("METHOD")) {
    if (!is(model, "modelorg")) {
        stop("needs an object of class modelorg!")
    }
    if (!checkVersion(model)) {
        stop("model is of wrong version!")
    }
    if (mc.cores < 2) {
        stop("mc.cores should be at least 2.")
    }
    if (SYBIL_SETTINGS("SOLVER") != solver) {
        SYBIL_SETTINGS("SOLVER", solver)
        cat("SYBIL_SETTINGS(SOLVER) has been set to", SYBIL_SETTINGS("SOLVER"), "\n")
    }
    if (SYBIL_SETTINGS("METHOD") != method) {
        SYBIL_SETTINGS("METHOD", method)
        cat("SYBIL_SETTINGS(METHOD) has been set to", SYBIL_SETTINGS("METHOD"), "\n")
    }

    fva     <- suppressMessages(multiDel(model,
                                         nProc=mc.cores,
                                         todo='fluxVar',
                                         del1=react_id(model)))
    fva.lim <- data.frame(min=unlist(mclapply(fva, mc.cores=mc.cores, minSol, lp_obj)),
                          max=unlist(mclapply(fva, mc.cores=mc.cores, maxSol, lp_obj)),
                          row.names=react_id(model),
                          stringsAsFactors=F)

    deadend <- deadEndMetabolites(model)
    if (length(deadend$der) == 0) {
        tol <- 0
    } else {
        tol <- max(abs(fva.lim[deadend$der, ])) * (1 + .Machine$double.eps)
    }
    smat    <- as.matrix(S(model))
    tol     <- max(tol,
                   10^max(-10,
                          min(-6,
                              floor(log10(min(abs(smat[smat!=0]))/max(abs(smat[smat!=0])))) - 1)))

    return (tol)
}


#' Rescue a model
#'
#' The function rescues a given model.
#' @param model An object of class \code{modelorg}.
#' @param target A numeric vector for growth target.
#' @param react A numeric vector or a character vector containing reaction id's. Default: reactions in objective function.
#' @param weight.type A character indicating which type of weighting to use in model objective modification, i: 1, r: 1/coefficient, s: 1/sqrt(coefficient). Default: r.
#' @param timeout The maximum time in seconds to allow for LP call to return. Default: 12.
#' @param prefix.rescue A string indicating the prefix of output rescue model. Default: no output.
#' @param prefix.rescued A string indicating the prefix of output rescued model. Default: no output.
#' @param rescue.threshold A numeric value indicating the threshold to consider a rescue. Default: 1e-5.
#' @return The rescue and rescued models, as well as the coefficient set to rescue reactions. SYBIL_SETTINGS("OPT_DIRECTION") is set as "min".
#' @import sybil
#' @importFrom sys eval_safe
#' @export
#' @examples 
#' data(Ec_core)
#' rescue(Ec_core, target=0.1)
rescue <- function(model, target, react = NULL, weight.type = 'r', timeout = 12,
                   prefix.rescue = NA, prefix.rescued = NA, rescue.threshold = 1e-5) {
    if (!is(model, "modelorg")) {
        stop("needs an object of class modelorg!")
    }
    if (!checkVersion(model)) {
        stop("model is of wrong version!")
    }
    if (is.null(react)) {
        obj.ind <- which(obj_coef(model) == 1)
    } else if (is.character(react)) {
        obj.ind <- match(react, react_id(model))
    } else if (is.numeric(react)) {
        if (all(react == round(react))) {
            obj.ind <- react
        } else {
            stop("argument react must be integer!")
        }
    } else {
        stop("argument react must be numeric or character!")
    }
    if (length(target) != length(obj.ind)) {
        stop("target must have the same length as react!")
    }

    options(stringsAsFactors=F)
    if (SYBIL_SETTINGS("TOLERANCE") != 1e-12) {
        SYBIL_SETTINGS("TOLERANCE", 1e-12)
        cat("SYBIL_SETTINGS(TOLERANCE) has been set to", SYBIL_SETTINGS("TOLERANCE"), "\n")
    }
    if (SYBIL_SETTINGS("OPT_DIRECTION") != "min") {
        SYBIL_SETTINGS("OPT_DIRECTION", "min")
        cat("SYBIL_SETTINGS(OPT_DIRECTION) has been set to", SYBIL_SETTINGS("OPT_DIRECTION"), "\n")
    }

    ##- initialize temporary model
    model.rescue   <- model
    smat           <- as.matrix(S(model))
    colnames(smat) <- react_id(model)
    rownames(smat) <- met_id(model)
    met.obj        <- smat[, obj.ind, drop=F]
    ## both substrates and products
    #met.obj.ind    <- which(rowSums(abs(met.obj)) != 0)
    ## not rescue BIOMASS, NEFA, DROPLET
    #biomass.ind    <- grep(names(met.obj.ind), pattern="BIOMASS|biomass|NEFA|DROPLET|MNXM2\\[", value=F, perl=T)
    #if (length(biomass.ind) > 0)
        #met.obj.ind    <- met.obj.ind[-biomass.ind]
    ## only substrates
    met.obj.ind    <- which(rowSums(met.obj) < 0)
    met.comp       <- met_comp(model)
    ref.ind        <- which.max(apply(met.obj[met.obj.ind, , drop=F], 2, function(x) {
        length(which(x != 0))
    }))  # the reaction with most metabolites among obj.ind: BIOMASS reaction

    ##- add HELP and RECO reactions to metabolites of targeted reactions
    coef <- data.frame(character(), rep(double(), length(obj.ind)), stringsAsFactors=F)
    for (ri in 1:length(obj.ind)) {
        for (mi in 1:length(met.obj.ind)) {
            met <- names(met.obj.ind)[mi]
            if (met.obj[met.obj.ind[mi], ri] != 0) {
                if (sign(met.obj[met.obj.ind[mi], ri]) == sign(met.obj[met.obj.ind[mi], ref.ind])) {
                    reco.type <- "_PUSH_"
                } else {
                    reco.type <- "_PULL_"
                }
                comp.reco <- paste("RECO",
                                   reco.type,
                                   substr(met,
                                          unlist(regexpr(pattern='[', met, fixed=T)) + 1,
                                          unlist(regexpr(pattern=']', met, fixed=T)) - 1),
                                   sep='')
                met.reco  <- sub(met, pattern="[", replacement=paste("[RECO", reco.type, sep=''), fixed=T)

                mod_compart(model.rescue) <- union(mod_compart(model.rescue), comp.reco)
                rea.suffix <- paste(sub(sub(met, pattern='[', replacement='_', fixed=T),
                                        pattern=']', replacement='', fixed=T),
                                    sep='')

                ##- HELP reactions
                if (met.obj[met.obj.ind[mi], ri] < 0) {
                    mets     <- c(met, met.reco)
                    metcomps <- c(met.comp[met.obj.ind[mi]], which(mod_compart(model.rescue) == comp.reco))
                } else {
                    mets     <- c(met.reco, met)
                    metcomps <- c(which(mod_compart(model.rescue) == comp.reco), met.comp[met.obj.ind[mi]])
                }
                if (! paste("HELP", reco.type, rea.suffix, sep='') %in% react_id(model.rescue)) {
                    model.rescue <- addReactFixed(model   = model.rescue,
                                                  id         = paste("HELP", reco.type, rea.suffix, sep=''),
                                                  met        = mets,
                                                  Scoef      = c(-1, 1),
                                                  reversible = F,
                                                  lb         = 0,
                                                  ub         = SYBIL_SETTINGS("MAXIMUM"),
                                                  obj        = 0,
                                                  gprAssoc   = "",
                                                  metName    = mets,
                                                  metComp    = metcomps
                                                  )
                }

                ##- RECO reactions
                if (! paste("RECO", reco.type, rea.suffix, sep='') %in% react_id(model.rescue)) {
                    model.rescue <- addReactFixed(model      = model.rescue,
                                                  id         = paste("RECO", reco.type, rea.suffix, sep=''),
                                                  met        = c(met.reco),
                                                  Scoef      = ifelse(met.obj[met.obj.ind[mi], ri] < 0, c(1), c(-1)),
                                                  reversible = F,
                                                  lb         = 0,
                                                  ub         = SYBIL_SETTINGS("MAXIMUM"),
                                                  obj        = 0,
                                                  gprAssoc   = "",
                                                  metName    = c(met.reco),
                                                  metComp    = c(which(mod_compart(model.rescue) == comp.reco))
                                                  )
                }

                if (! paste("RECO", reco.type, rea.suffix, sep='') %in% coef[, 1]) {
                    coef <- rbind(coef,
                                  c(paste("RECO", reco.type, rea.suffix, sep=''),
                                    abs(met.obj[met.obj.ind[mi], ])))
                }

                ##- modify targeted reactions in rescue model
                S(model.rescue)[met.obj.ind[mi], obj.ind[ri]] <- 0
                S(model.rescue)[match(met.reco, met_id(model.rescue)),
                                obj.ind[ri]] <- met.obj[met.obj.ind[mi], ri]
            }
        }
        model.rescue <- changeBounds(model.rescue, react=obj.ind[ri], lb=target[ri], ub=target[ri])
    }

    ##- write rescue model
    if (!is.na(prefix.rescue)) {
        modelorg2tsv(model.rescue, prefix=prefix.rescue, quote=F)
    }

    ##- rescue coefficients of metabolites of targeted reactions
    coefs <- matrix(apply(as.matrix(coef[, -1]), 2, as.numeric),
                    ncol=length(obj.ind),
                    dimnames=list(coef[, 1], react_id(model)[obj.ind]))
    if (weight.type == 'o') {
        coef.weight <- rep(1, nrow(coefs))
    } else {
        if (weight.type == 'r') {
            coef.weight <- 1/coefs
        } else if (weight.type == 's') {
            coef.weight <- 1/sqrt(coefs)
        }
        coef.weight <- coef.weight[, react_id(model)[which(obj_coef(model) == 1)], drop=F]
    }

    ##- change objective function on RECO reaction fluxes
    model.weight <- changeObjFunc(model.rescue,
                                  react=rownames(coefs),
                                  obj_coef=coef.weight)
    if (SYBIL_SETTINGS("SOLVER") == "glpkAPI") {
        fba <- optimizeProb(model.weight, algorithm='fba', retOptSol=T, lpdir='min',
                            solverParm=list(TM_LIM=1000*timeout))
    } else if (SYBIL_SETTINGS("SOLVER") == "lpSolveAPI") {
        fba <- optimizeProb(model.weight, algorithm='fba', retOptSol=T, lpdir='min',
                            solverParm=list(timeout=timeout))
    } else if (SYBIL_SETTINGS("SOLVER") == "cplexAPI") {
        fba <- optimizeProb(model.weight, algorithm='fba', retOptSol=T, lpdir='min',
                            solverParm=list(CPX_PARAM_TILIM=timeout))
    } else {
        fba <- tryCatch(
            sys::eval_safe(optimizeProb(model.weight, algorithm='fba',
                                        retOptSol=T, lpdir='min'),
                           timeout=timeout),
            error=function(e) {
                fba <- NULL
            })
    }
    if (!is.null(fba))
        cos <- checkOptSol(fba)
    if (is.null(fba) ||
        exit_code(cos)!=0 ||
        (SYBIL_SETTINGS("SOLVER") == "glpkAPI" && status_code(cos)!=5) ||
        (SYBIL_SETTINGS("SOLVER") %in% c("clpAPI", "lpSolveAPI") && status_code(cos)!=0)
        )
                                        #if (is.null(fba) || !checkOptSol(fba, onlywarn=T))
        stop("FBA failed. Please try with another SOLVER-METHOD!")
    
    wtflux <- setNames(getFluxDist(fba), react_id(model.weight))
    
    ##- necessary reactions for rescue
    flux.reco <- wtflux[rownames(coefs)]
    threshold <- (rescue.threshold * target %*% t(coefs))[1, ]
    reco.keep <- names(flux.reco)[which(flux.reco > threshold)]
    help.keep <- sub(reco.keep, pattern="RECO", replacement="HELP")
    
    ##- initialize rescued model
    model.rescued         <- model
    smat.rescue           <- as.matrix(S(model.rescue))
    colnames(smat.rescue) <- react_id(model.rescue)
    rownames(smat.rescue) <- met_id(model.rescue)

    ##- add rescued reactions to the given model
    for (rea in c(reco.keep, help.keep)) {
        mets <- smat.rescue[, rea, drop=T]
        mind <- which(mets != 0)
        rind <- which(react_id(model.rescue) == rea)
        metcomps <- mod_compart(model.rescue)[met_comp(model.rescue)[mind]]
        mod_compart(model.rescued) <- union(mod_compart(model.rescued), metcomps)
        ub   <- ifelse(regexpr(rea, pattern="^RECO_", perl=T) > 0,
                       abs(flux.reco[rea]),
                       uppbnd(model.rescue)[rind])
        model.rescued <- addReactFixed(model    = model.rescued,
                                  id         = rea,
                                  met        = names(mets)[mind],
                                  Scoef      = mets[mind],
                                  reversible = react_rev(model.rescue)[rind],
                                  lb         = lowbnd(model.rescue)[rind],
                                  ub         = ub,
                                  obj        = 0,
                                  gprAssoc   = "",
                                  metName    = names(mets)[mind],
                                  metComp    = match(metcomps, mod_compart(model.rescued))
                                  )

        ##- modify targeted reactions in rescued model
        if (regexpr(rea, pattern="^RECO_", perl=T) > 0) {
            if (regexpr(rea, pattern="_PUSH_", perl=T) > 0) {
                mind.reco   <- match(names(mets)[mind], met_id(model.rescued))
                mind.toreco <- match(sub(names(mets)[mind], pattern='RECO_PUSH_', replacement=''),
                                     met_id(model.rescued))
                S(model.rescued)[mind.reco, obj.ind]    <- S(model.rescued)[mind.toreco, obj.ind]
                S(model.rescued)[mind.toreco, obj.ind]  <- 0
            } else {
                mind.reco   <- match(names(mets)[mind], met_id(model.rescued))
                mind.toreco <- match(sub(names(mind), pattern='RECO_PULL_', replacement=''),
                                     met_id(model.rescued))
                coef.toreco <- met.obj[met_id(model.rescued)[mind.toreco], ]
                pull.ind    <- which(sign(coef.toreco) != sign(coef.toreco[ref.ind]))
                S(model.rescued)[mind.reco, pull.ind]   <- S(model.rescued)[mind.toreco, pull.ind]
                S(model.rescued)[mind.toreco, pull.ind] <- 0
            }
        }
    }
    for (ri in 1:length(obj.ind)) {
        model.rescued <- changeBounds(model.rescued, react=obj.ind[ri],
                                      lb=lowbnd(model)[obj.ind[ri]],
                                      ub=uppbnd(model)[obj.ind[ri]])
    }

    ##- write rescued model
    if (!is.na(prefix.rescued)) {
        modelorg2tsv(model.rescued, prefix=prefix.rescued, quote=F)
    }
    
    return (list(rescued=model.rescued,
                 rescue=model.rescue,
                 coef=coef.weight))
}
