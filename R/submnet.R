#' Structure of Class "scoreGeneDel"
#'
#' Structure of the class \code{scoreGeneDel}. Objects of this class are returned by the function submnet.
#' @param model An object of class \code{modelorg} indicating the weighted \code{rescue} model obtained from the rescue process.
#' @param condition The experimental condition ID.
#' @param fitness.random Random-based fitness with weighting scheme.
#' @param fitness.ranks Ranks-based fitness with weighting scheme.
#' @param fitness.id.random Random-based fitness without weighting scheme.
#' @param fitness.id.ranks Ranks-based fitness without weighting scheme.
#' @param ess.gene Percentages of essential genes. The computation of essentiality is deprecated in this version.
#' @param ess.reaction Percentages of essential reactions. The computation of essentiality is deprecated in this version.
#' @param gene.del Number of deleted genes.
#' @param gene.sets Gene sets.
#' @param ratio.GS Percentages of remaining genes in each gene set.
#' @param sub.genes Remaining genes in submodels after propagation.
#' @param sub.reacs Remaining reactions in submodels after propagation.
#' @param rescue.met Fraction of every rescued metabolite among random draws.
#' @examples 
#' data(yarliSubmnets)
#' attributes(yarliSubmnets[[1]])
#' @export
scoreGeneDel <- function(model = NULL, condition = NA,
                         fitness.random = NULL, fitness.ranks = NULL,
                         fitness.id.random = NULL,  fitness.id.ranks = NULL,
                         ess.gene = NULL, ess.reaction = NULL, gene.del = NULL, gene.sets = NULL, 
                         ratio.GS = NULL, sub.genes = NULL, sub.reacs = NULL, rescue.met = NULL) {
    res <- list(
        model             = model,
        condition         = condition,
        fitness.random    = fitness.random,
        fitness.ranks     = fitness.ranks,
        fitness.id.random = fitness.id.random,
        fitness.id.ranks  = fitness.id.ranks,
        ess.gene          = ess.gene,
        ess.reaction      = ess.reaction,
        gene.del          = gene.del,
        gene.sets         = gene.sets,
        ratio.GS          = ratio.GS,
        sub.genes         = sub.genes,
        sub.reacs         = sub.reacs,
        rescue.met        = rescue.met
    )
    class(res) <- "scoreGeneDel"

    return(res)
}


#' Fitness of gene removal-based submodels with different gene rankings
#'
#' This function computes the fitness of submodels by removing genes in different gene rankings.
#' @param model An object of class \code{modelorg} indicating the weighted \code{rescue} model obtained from the rescue process.
#' @param ranks A list of data frames of scores for ranking genes, with gene per row, e.g. data.frame(pkm=pkm expression, rel=relative expression).
#' @param rescue.weight A vector of rescue reaction weights. Default: NULL, the weights are computed from the given model with gene.num=1.
#' @param step An integer indicating the step in numbers of genes to remove. Default: 1, gene-by-gene removal. 
#' When there are many genes in the model, the step is multiplied by an exponent of 2 for later removals. 
#' This is to reduce the computing time for non-informative sub-models at the end of the series.
#' @param draw.num Number of random draws. Default: 0.
#' @param obj.react A string indicating objective reaction ID. Default: reaction producing BIOMASS.
#' @param mc.cores The number of cores to use (at least 1), i.e. at most how many child processes will be run simultaneously. Default: 1.
# #' @param MILP Logical value indicating whether MILP is used. Default: FALSE.
# #' @param timeout The maximum time in seconds to allow for MILP call to return. Default: 12.
# #' @param verboseMode An integer value indicating the amount of output to stdout: 0: nothing, 1: MILP status messages. Default: 0.
#' @param tol The maximum value to be considered null. Default: \code{SYBIL_SETTINGS("TOLERANCE")}.
#' @param solver \code{\link{sybil}} solver. Default: \code{SYBIL_SETTINGS("SOLVER")}.
#' @param method \code{\link{sybil}} method. Default: \code{SYBIL_SETTINGS("METHOD")}.
#' @return An object of class \code{scoreGeneDel} for the submodel construction simulation.
#' @import sybil stats
# #' @importFrom glpkAPI GLP_ON
# #' @importFrom sys eval_safe
#' @examples 
#' data(Ec_core)
#' mod <- rescue(Ec_core, target=0.1)
#' mod.weight <- changeObjFunc(mod$rescue, react=rownames(mod$coef), obj_coef=mod$coef)
#' ranks <- list(rep.1=data.frame(expr=setNames(rnorm(length(sybil::allGenes(mod.weight)),
#'                                              mean=5, sd=4), sybil::allGenes(mod.weight))),
#'               rep.2=data.frame(expr=setNames(rnorm(length(sybil::allGenes(mod.weight)),
#'                                              mean=5, sd=4.1), sybil::allGenes(mod.weight))))
#' fn <- fitness(model=mod.weight, ranks=ranks, step=200, draw.num=1)
#' @export
fitness <- function(model, ranks, rescue.weight = NULL, step = 1, draw.num = 0, obj.react = NA,
                    mc.cores = 1, #MILP = FALSE, timeout = 12, verboseMode = 0,
                    tol = SYBIL_SETTINGS("TOLERANCE"),
                    solver = SYBIL_SETTINGS("SOLVER"),
                    method = SYBIL_SETTINGS("METHOD")) {
    ##- settings ----
    if (!is(model, "modelorg")) {
        stop("needs an object of class modelorg!")
    }
    if (!checkVersion(model)) {
        stop("model is of wrong version!")
    }
    if (is.null(ranks)) {
        stop("no rankings are provided!")
    }
    if (is.null(names(ranks))) {
        stop("no names of ranks are found!")
    }
    ## MILP <- FALSE
    ## if (MILP) {
    ##     if (SYBIL_SETTINGS("SOLVER") != "glpkAPI") {
    ##         SYBIL_SETTINGS("SOLVER", "glpkAPI")
    ##         cat("SYBIL_SETTINGS(SOLVER) has been set to", SYBIL_SETTINGS("SOLVER"), "\n")
    ##     }
    ## }
    if (SYBIL_SETTINGS("SOLVER") != solver) {
        SYBIL_SETTINGS("SOLVER", solver)
        cat("SYBIL_SETTINGS(SOLVER) has been set to", SYBIL_SETTINGS("SOLVER"), "\n")
    }
    if (SYBIL_SETTINGS("METHOD") != method) {
        SYBIL_SETTINGS("METHOD", method)
        cat("SYBIL_SETTINGS(METHOD) has been set to", SYBIL_SETTINGS("METHOD"), "\n")
    }
    if (SYBIL_SETTINGS("OPT_DIRECTION") != "min") {
        SYBIL_SETTINGS("OPT_DIRECTION", "min")
        cat("SYBIL_SETTINGS(OPT_DIRECTION) has been set to", SYBIL_SETTINGS("OPT_DIRECTION"), "\n")
    }
    if (SYBIL_SETTINGS("TOLERANCE") != tol) {
        SYBIL_SETTINGS("TOLERANCE", tol)
        cat("SYBIL_SETTINGS(TOLERANCE) has been set to", SYBIL_SETTINGS("TOLERANCE"), "\n")
    }
    options(stringsAsFactors=F)
    RNGkind("L'Ecuyer-CMRG")
    set.seed(1000)
    
    reps      <- names(ranks)
    rep.num   <- length(reps)
    reps.char <- do.call(cbind, strsplit(reps, ''))
    reps.commonend <- min(which(apply(reps.char, 1, function(rc) {length(unique(rc))}) > 1L)) - 1L
    while (!grepl("[a-zA-Z]", reps.char[reps.commonend, 1L])) {
        reps.commonend <- reps.commonend - 1L
    }
    condition <- paste(reps.char[1L:reps.commonend, 1L], collapse='')
    genes     <- rownames(ranks[[1]])
    if (sum(genes != sybil::allGenes(model)) > 0L || 
        (!is.null(ranks) && sum(rownames(ranks[[1]]) != genes) > 0L)) {
        stop("score: genes in model, expr and ranks do not match")
    }
    gene.num  <- length(genes)
    if (gene.num > 1000) {
        gene.num.draw <- c(seq(0L, floor(gene.num/4), step),
                           seq(floor(gene.num/4)+1, floor(gene.num/2), step*2),
                           seq(floor(gene.num/2)+1, floor(3*gene.num/4), step*4),
                           seq(floor(3*gene.num/4)+1, gene.num, step*8))
    } else if (gene.num > 500) {
        gene.num.draw <- c(seq(0L, floor(gene.num/2), step),
                           seq(floor(gene.num/2)+1, gene.num, step*2))
    } else {
        gene.num.draw <- seq(0L, gene.num, step)
    }
    if (is.null(rescue.weight)) {
        rescue.weight <- (weightReacts(model, mc.cores=mc.cores, gene.num=1))$weight
    }
    if (length(setdiff(names(rescue.weight), react_id(model))) > 0) {
        stop("names of rescue.weight not found in react_id(model)!")
    }
    mc.cores2 <- max(1L, as.integer(floor(mc.cores/length(gene.num.draw))))
    mc.cores1 <- max(1L, as.integer(floor(mc.cores/mc.cores2)))
    recos     <- grep('RECO', react_id(model), perl=T, value=T)
    draw.lim  <- gene.num
    draw.num  <- as.integer(draw.num)
    ##-----
    
    ##- initial fba
    fba.weight <- optimizeProb(model, algorithm="fba", retOptSol=T, lpdir='min')
    if (!checkOptSol(fba.weight, onlywarn=T)) stop("cannot perform FBA for input model!")
    wtflux     <- getFluxDist(fba.weight)
    if (is.na(obj.react)) {
	if (length(grep("biomass", tolower(met_id(model)), value=F)) > 0)
            obj.react <- react_id(model)[which(S(model)[grep("biomass",
                                                         tolower(met_id(model)),
                                                         value=F), ]
                                           > 0)]
	else
            obj.react <- react_id(model)[grep("biomass", tolower(react_id(model)))]
    }
    if (length(obj.react) < 1) {
        stop("cannot determine obj.react producing BIOMASS!")	
    } else if (length(obj.react) > 1) {
	stop("too many obj.react producing BIOMASS!")
    } else if (uppbnd(model)[which(react_id(model)==obj.react)]-lowbnd(model)[which(react_id(model)==obj.react)] > tol) {
        stop("obj.react does not have a fixed flux!")
    } else {
        obj.fixed <- uppbnd(model)[which(react_id(model)==obj.react)]
    }
    ##- compute scores for different types of gene removal ----
    recoscores <- mclapply(1L:length(gene.num.draw), mc.cores=mc.cores1, mc.set.seed=T, function(i) {
        draw.num1 <- 1L
      
        ##- draw based on ranks, only once ----
        draw.ranks <- NULL
        if (!is.null(ranks)) {
            draw.ranks <- mclapply(1L:rep.num, mc.cores=1L, function(l) {
                mclapply(1L:ncol(ranks[[l]]), mc.cores=mc.cores2, function(k) {
                    rank        <- ranks[[l]][, k, drop=T]
                    names(rank) <- rownames(ranks[[l]])
                    rank.sort <- sort(rank, decreasing=F)
                    draw.rank <- matrix(0L, nrow=gene.num, ncol=draw.num1)
                    if (draw.lim < gene.num.draw[i]) {
                        ##- draw the first 'draw.lim' genes of lowest expression and
                        ##- 'gene.num.draw[i] - draw.lim' random genes
                        if (gene.num.draw[i] > 0L) {
                            draw1rowindex <- matrix(
                                replicate(draw.num1,
                                          which(genes %in%
                                                    names(rank.sort)[
                                                        c(1L:draw.lim,
                                                          sample((draw.lim+1L):gene.num,
                                                                 gene.num.draw[i]-draw.lim,
                                                                 replace=F))])),
                                ncol=draw.num1)
                        } else {
                            draw1rowindex <- matrix(0L, nrow=0L, ncol=draw.num1)
                        }
                    } else {
                        ##- draw the first 'gene.num.draw[i]' genes of lowest expression
                        if (gene.num.draw[i] > 0L) {
                            draw1rowindex <- matrix(
                                replicate(draw.num1,
                                          which(genes %in%
                                                    names(rank.sort)[1L:gene.num.draw[i]])),
                                ncol=draw.num1)
                        } else {
                            draw1rowindex <- matrix(0L, nrow=0L, ncol=draw.num1)
                        }
                    }
                    draw.rank <- sapply(1L:draw.num1, function(j) {
                        draw.rank[draw1rowindex[, j], j] <- 1L
                        return (draw.rank[, j])
                    })
                    
                    return (draw.rank)
                })
            })
        }
        ##-----
                
        ##- draw randomly, 'draw.num' times ----
        draw.random <- matrix(0L, nrow=gene.num, ncol=draw.num)
        if (draw.num > 0L) {
            #- draw 'gene.num.draw[i]' genes randomly
            if (gene.num.draw[i] > 0L) {
                draw1rowindex <- matrix(
                    replicate(draw.num,
                              sample(1L:gene.num,
                                     gene.num.draw[i],
                                     replace=F)),
                    ncol=draw.num)
            } else {
                draw1rowindex <- matrix(0L, nrow=0L, ncol=draw.num)
            }
            draw.random <- sapply(1L:draw.num, function(j) {
                draw.random[draw1rowindex[, j], j] <- 1L
                return (draw.random[, j])
            })
        }
        ##-----

        ##- combine draws into list ----
        draw.list <- c(lapply(apply(draw.random, 2, list), unlist),
                       unlist(lapply(draw.ranks, function(dr) {
                           sapply(dr, function(x) {
                               lapply(apply(x,2,list), unlist)
                           })
                       }), recursive=F))
        ##-----
        
        ##- genes deleted in each draw ----
        genes.del.param <- mclapply(draw.list, mc.cores=mc.cores2, function(x) {
            genes[which(as.numeric(x) == 1L)]
        })
        names(genes.del.param) <- c(if (draw.num > 0) paste0("random.", 1L:draw.num) else NULL,
                                    sapply(1L:length(ranks), function(r) {
                                        paste0(colnames(ranks[[r]]), '.', names(ranks)[r])}))
        gene.del <- gene.num.draw[i]
        ##-----

        ## if (!MILP) {
            ##< FBA LP 0
            fluxlist.param <- mclapply(genes.del.param, mc.cores=1, function(x) {
                if (length(x)==0) x <- NULL
                obj.coef <- setNames(obj_coef(model), react_id(model))
                obj.coef.reco <- obj.coef[names(rescue.weight)]/obj.fixed*rescue.weight
                fba.x <- optimizeProb(changeObjFunc(model,
                                                    react=names(obj.coef.reco),
                                                    obj_coef=obj.coef.reco),
                                      gene=x, lb=0, ub=0,
                                      algorithm="fba",
                                      lpdir="min")
                if (checkOptSol(fba.x, onlywarn=T)) {                      
                    return (list(flux=setNames(getFluxDist(fba.x), react_id(model)),
                                 obj=lp_obj(fba.x)))
                }
                return (list(flux=setNames(rep(NaN, react_num(model)), react_id(model)),
                             obj=NaN))
            })
            fluxmat.param <- do.call(rbind, mclapply(fluxlist.param, mc.cores=mc.cores2, function(flp) {
                flp$flux
            }))
            ##>
            
            ## ##< FBA LP 1
            ## fluxmat.param <- do.call(rbind, mclapply(genes.del.param, mc.cores=1, function(x) {
            ##     if (length(x)==0) x <- NULL
            ##     fba.x <- optimizeProb(model,
            ##                           gene=x, lb=0, ub=0,
            ##                           algorithm="fba",
            ##                           lpdir="min")
            ##     if (checkOptSol(fba.x, onlywarn=T)) {                      
            ##     	return (setNames(getFluxDist(fba.x), react_id(model)))
            ##     }
            ##     return (setNames(rep(NaN, react_num(model)), react_id(model)))
            ## }))
            ## ##>
            
            ## ##< FBA LP 2
            ## ##- delete genes ----
            ## gd.param <- mclapply(genes.del.param, mc.cores=mc.cores2, function(x) {
            ##     if (length(x) == 0) {
            ##         return (fba.weight)
            ##     }
            ##     gd <- suppressMessages(geneDeletion(model, x, combinations=length(x),
            ##                                         lb=NULL, ub=NULL,
            ##                                         poCmd=list("getFluxDist"), verboseMode=1))
            ##     if (checkOptSol(gd, onlywarn=T)) {
            ##         return (gd)
            ##     }
            ##     return (NaN)
            ## })
            ## ##-----
            
            ## ##- then obtain the corresponding fluxes ----
            ## fluxmat.param <- matrix(
            ##     unlist(
            ##         lapply(
            ##             mclapply(gd.param, mc.cores=mc.cores2, function(x) {
            ##                 if (is(x, "optsol_genedel") || is(x, "optsol_optimizeProb")) {
            ##                     return (postProc(x))
            ##                 }
            ##                 return (NaN)
            ##             }),
            ##             function(y) {
            ##                 if (is(y, "ppProc")) {
            ##                     return (sybil::pa(y))
            ##                 }
            ##                 return (list(list(rep(NaN, reac.num))))
            ##             }),
            ##         use.names=F),
            ##     ncol=length(reactions), byrow=T)
            ## colnames(fluxmat.param) <- reactions
            ## ##-----
            ## ##>
        ## } else {      
        ## ##< MILP
        ##     fluxlist.param <- mclapply(1L:length(genes.del.param), mc.cores=mc.cores2, function(j) {
        ##         x <- genes.del.param[[j]]
        ##                                 #print(paste("gene", gene.del, j))
        ##         if (length(x) == 0) {
        ##             return (list(flux=setNames(wtflux, react_id(model)),
        ##                          obj=lp_obj(fba.weight)))
        ##         } else {
        ##             fba.x <- optimizeProb(model,
        ##                                   gene=x, lb=0, ub=0,
        ##                                   algorithm="fba",
        ##                                   lpdir="min", retOptSol=TRUE)
        ##             if (!checkOptSol(fba.x, onlywarn=T)) {
        ##                 return (list(flux=setNames(rep(NaN, react_num(model)), react_id(model)),
        ##                              obj=NA))
        ##             }
        ##             mip <- tryCatch(
        ##                 sys::eval_safe(
        ##                     suppressWarnings(suppressMessages(
        ##                         optimizeProb(
        ##                             model, algorithm='room',
        ##                             wtflux=wtflux,
        ##                             gene=x, lb=0, ub=0,
        ##                             minRescue=lp_obj(fba.x) + 0.001,
        ##                             rescueWeight=rescue.weight,
        ##                             solverParm=list(PRESOLVE=glpkAPI::GLP_ON),
        ##                                 #delta=0, epsilon=0,
        ##                             LPvariant=F))),
        ##                     timeout=timeout),
                        ## error=function(e1) {
                        ##     print(paste0("first-attempt error catched: ", condition, ", gene.del: ", length(x), ", draw: ", j, ", trying delta=0, epsilon=0.01"))
                        ##     mip <- tryCatch(
                        ##         sys::eval_safe(suppressWarnings(suppressMessages(
                        ##             optimizeProb(
                        ##                 model, algorithm='room',
                        ##                 wtflux=wtflux,
                        ##                 gene=x, lb=0, ub=0,
                        ##                 minRescue=lp_obj(fba.x) + 0.01,
                        ##                 solverParm=list(PRESOLVE=glpkAPI::GLP_ON),
                        ##                 delta=0, epsilon=0, LPvariant=F))),
                        ##                   timeout=timeout),
                        ##         error=function(e2) {
                        ##             print(paste0("second-attempt error catched: ", condition, ", gene.del: ", length(x), ", draw: ", j, ", trying delta=0.01, epsilon=0.01"))
                        ##             mip <- tryCatch(
                        ##                 sys::eval_safe(suppressWarnings(suppressMessages(
                        ##                     optimizeProb(
                        ##                         model, algorithm='room',
                        ##                         wtflux=wtflux,
                        ##                         gene=x, lb=0, ub=0,
                        ##                         minRescue=lp_obj(fba.x)*1.01 + 0.01,
                        ##                         solverParm=list(PRESOLVE=glpkAPI::GLP_ON),
                        ##                         delta=0, epsilon=0, LPvariant=F))),
                        ##                           timeout=timeout),
        ##                 error=function(e3) {
        ##                     if (verboseMode > 0) 
        ##                         print(paste0("MILP-attempt error catched: ", condition, ", gene.del: ", length(x), ", draw: ", j, ", trying with LPvariant"))
        ##                     mip <- tryCatch(
        ##                         sys::eval_safe(
        ##                             suppressWarnings(suppressMessages(
        ##                                 optimizeProb(
        ##                                     model, algorithm='room',
        ##                                     wtflux=wtflux,
        ##                                     gene=x, lb=0, ub=0,
        ##                                     minRescue=lp_obj(fba.x) + 0.001,
        ##                                     rescueWeight=rescue.weight,
        ##                                     solverParm=list(PRESOLVE=glpkAPI::GLP_ON),
        ##                                     delta=0, epsilon=0,
        ##                                     LPvariant=T))),
        ##                             timeout=timeout),
        ##                         error=function(e4) {
        ##                             if (verboseMode > 0) 
        ##                                 print("MILP LPvariant failed. NULL is returned.")
        ##                             mip <- NULL
        ##                         })
        ##                 })
                    
        ##             if (!is.null(mip) && suppressWarnings(checkOptSol(mip, onlywarn=T))) {
        ##                 return (list(flux=setNames(getFluxDist(mip), react_id(model)),
        ##                              obj=lp_obj(mip)))
        ##             }
        ##             return (list(flux=setNames(rep(NaN, react_num(model)), react_id(model)),
        ##                          obj=NA))
        ##         }
        ##     })
        ##     fluxmat.param <- do.call(rbind, mclapply(fluxlist.param, mc.cores=mc.cores2, function(flp) {
        ##         flp$flux
        ##     }))
        ## }
        
        ##- fluxes of RECO reactions ----
        recoflux.param <- fluxmat.param[, recos, drop=F]
        recoflux.param[is.nan(recoflux.param)] <- SYBIL_SETTINGS("MAXIMUM")
        ##-----
        
        ##- existing status of RECO reactions ----
        recoexist.param <- abs(recoflux.param) > tol
        recoexist.param <- apply(recoexist.param, c(1,2), as.numeric)
        if (sum(colnames(recoexist.param) != names(rescue.weight)) > 0) {
            stop("Reco weight: Different reco orders!")
        }
        ##-----
        
        ##- computing the model fitness ----
        #- with reaction weights
        ## fitness    <- setNames(c(recoexist.param %*% rescue.weight),
        ##                        c(if (draw.num > 0) paste0("F.random.",   1L:draw.num) else NULL,
        ##                          sapply(reps, function(reps.i) {
        ##                              paste("F",    colnames(ranks[[1]]), reps.i, sep='.')
        ##                          })))       
        fitness <- setNames(unlist(mclapply(fluxlist.param, mc.cores=mc.cores2, function(flp) {
            flp$obj
        })), c(if (draw.num > 0) paste0("F.random.",   1L:draw.num) else NULL,
               sapply(reps, function(reps.i) {
                   paste("F",    colnames(ranks[[1]]), reps.i, sep='.')
               })))
        
        #- without reaction weights
        fitness.id <- setNames(c(recoexist.param %*% rep(1/length(rescue.weight), length(rescue.weight))),
                               c(if (draw.num > 0) paste0("F.id.random.", 1L:draw.num) else NULL,
                                 sapply(reps, function(reps.i) {
                                     paste("F.id", colnames(ranks[[1]]), reps.i, sep='.')
                                 })))
        ##-----

        return (list(score=c(
                         ifelse(1 - fitness > 0,    1 - fitness,    0),
                         ifelse(1 - fitness.id > 0, 1 - fitness.id, 0),
                         if (draw.num > 0) apply(recoexist.param[1L:draw.num, , drop=F], 2, mean) else
                         apply(recoexist.param[1, , drop=F], 2, mean)),
                     gene.del=gene.del,
                     genes.del=genes.del.param,
                     ranks.name=colnames(ranks[[1]]),
                     draw.num=draw.num
                     )
                )
    })
    return (recoscores)
}


#' Identify the best ranking
#'
#' This function computes the performance indices of different rankings compared to the random ranking for gene removal and identify the best ranking
#' @param fns List of fitness objects.
#' @return The performance indices for all rankings and the best ranking.
#' @examples
#' data(Ec_core)
#' mod <- rescue(Ec_core, target=0.1)
#' mod.weight <- changeObjFunc(mod$rescue, react=rownames(mod$coef), obj_coef=mod$coef)
#' ranks <- list(rep.1=data.frame(expr=setNames(rnorm(length(sybil::allGenes(mod.weight)),
#'                                              mean=5, sd=4), sybil::allGenes(mod.weight))),
#'               rep.2=data.frame(expr=setNames(rnorm(length(sybil::allGenes(mod.weight)),
#'                                              mean=5, sd=4.1), sybil::allGenes(mod.weight))))
#' fn <- fitness(model=mod.weight, ranks=ranks, step=200, draw.num=1)
#' bestRanking(list(fn))
#' @export
bestRanking <- function(fns) {
    ranks.name <- fns[[1]][[1]]$ranks.name
    perfIndex <- do.call(cbind, mclapply(fns, mc.cores=length(fns), function(fn) {
        param.num <- length(fn)
        reps <- gsub(grep(paste0("F.",ranks.name[1], "."), names(fn[[1]]$score), value=T, fixed=T),
                     pattern=paste0("F\\.",ranks.name[1], "\\."), replacement='', perl=T) 
        draw.num <- fn[[1]]$draw.num
        if (draw.num == 0) {
            weights <- rep(1/param.num, param.num)
        } else {
            F.random.means <- sapply(1L:length(fn), function(k) {
                mean(fn[[k]]$score[paste0("F.random.", 1L:draw.num)])
            })
            weights <- F.random.means/sum(F.random.means)
        }
        ##- percentages of random draws better (i.e. higher fitness) than ranked draws ----
        r2random <- do.call(cbind, lapply(reps, function(rep) {
            ranks2random <- sapply(ranks.name, function(rn) {
                100 * sum(sapply(1L:length(fn), function(k) {
                    sum(fn[[k]]$score[paste0("F.", rn, '.', rep)] < fn[[k]]$score[paste0("F.random.", 1L:fn[[k]]$draw.num)]) / draw.num * weights[k] 
                }))
            })
            return (setNames(ranks2random, ranks.name))
        }))
	colnames(r2random) <- reps
        return (r2random)
    }))
    return (list(perfIndex=perfIndex,
                 rank.best=ranks.name[as.numeric(names(which.max(table(apply(perfIndex, 2, which.min)))))]))
}


#' Simulation of gene removal-based submodel series with a given ranking
#'
#' This function simulates the construction of a series of submodels by removing genes in a given ranking.
#' @param model An object of class \code{modelorg} indicating the weighted \code{rescue} model obtained from the rescue process.
#' @param fn An object returned by the fitness function.
#' @param rank.best Name of a ranking among simulated ones. Default: "expr".
#' @param gene.sets Named list of gene sets for gene set enrichment analysis. Default: NULL,
#' depletion fraction of gene sets should be further computed for gene set enrichment analysis.
# #' @param essential A logical value indicating whether essentiality analysis will be fulfilled. Default: FALSE.
#' @param mc.cores The number of cores to use (at least 1), i.e. at most how many child processes will be run simultaneously. Default: 1.
#' @param obj.react A string indicating objective reaction ID. Default: reaction producing BIOMASS.
#' @param tol The maximum value to be considered null. Default: \code{SYBIL_SETTINGS("TOLERANCE")}.
#' @return An object of class \code{scoreGeneDel} for the submodel construction simulation.
#' @import sybil stats
#' @examples 
#' data(Ec_core)
#' mod <- rescue(Ec_core, target=0.1)
#' mod.weight <- changeObjFunc(mod$rescue, react=rownames(mod$coef), obj_coef=mod$coef)
#' ranks <- list(rep.1=data.frame(expr=setNames(rnorm(length(sybil::allGenes(mod.weight)),
#'                                              mean=5, sd=4), sybil::allGenes(mod.weight))),
#'               rep.2=data.frame(expr=setNames(rnorm(length(sybil::allGenes(mod.weight)),
#'                                              mean=5, sd=4.1), sybil::allGenes(mod.weight))))
#' fn <- fitness(model=mod.weight, ranks=ranks, step=200, draw.num=1)
#' gene.sets <- list(X1=head(sybil::allGenes(mod.weight)), X2=tail(sybil::allGenes(mod.weight)))
#' sgd <- submnet(model=mod.weight, fn=fn, rank.best="expr",
#'                obj.react="Biomass_Ecoli_core_w_GAM", gene.sets=gene.sets)
#' @export
submnet <- function(model, fn, rank.best = "expr",
                    gene.sets = NULL, mc.cores = 1, obj.react = NA,
                    tol = SYBIL_SETTINGS("TOLERANCE")) {
    ##- settings ----
    if (!is(model, "modelorg")) {
        stop("needs an object of class modelorg!")
    }
    if (!checkVersion(model)) {
        stop("model is of wrong version!")
    }

    if (SYBIL_SETTINGS("OPT_DIRECTION") != "min") {
        SYBIL_SETTINGS("OPT_DIRECTION", "min")
        cat("SYBIL_SETTINGS(OPT_DIRECTION) has been set to", SYBIL_SETTINGS("OPT_DIRECTION"), "\n")
    }
    if (SYBIL_SETTINGS("TOLERANCE") != tol) {
        SYBIL_SETTINGS("TOLERANCE", tol)
        cat("SYBIL_SETTINGS(TOLERANCE) has been set to", SYBIL_SETTINGS("TOLERANCE"), "\n")
    }
    options(stringsAsFactors=F)

    rank.best.id <- grep(paste0("^", rank.best, "\\."), names(fn[[1]]$genes.del), value=T, perl=T)
    genes.del <- mclapply(fn, mc.cores=mc.cores, function(fni) {
        fni$genes.del[rank.best.id]
    })
    reps <- gsub(rank.best.id, pattern=paste0(rank.best, "\\."), replacement='', perl=T)
    rep.num <- length(reps)
    reps.char <- do.call(cbind, strsplit(reps, ''))
    reps.commonend <- min(which(apply(reps.char, 1, function(rc) {length(unique(rc))}) > 1L)) - 1L
    while (!grepl("[a-zA-Z]", reps.char[reps.commonend, 1L])) {
        reps.commonend <- reps.commonend - 1L
    }
    condition <- paste(reps.char[1L:reps.commonend, 1L], collapse='')
    mc.cores2   <- max(1L, as.integer(floor(mc.cores/length(fn))))
    mc.cores1   <- max(1L, as.integer(floor(mc.cores/mc.cores2)))
    recos       <- grep('^RECO', react_id(model), perl=T, value=T)
    iter.fba    <- 40
    essential   <- F
    if (is.null(gene.sets)) {
        warning("gene.sets is NULL, depletion fraction of gene sets should be further computed for gene set enrichment analysis.")
    }
    ##-----
        
    ##- get initial obj.react ----
    if (is.na(obj.react)) {
	if (length(grep("biomass", tolower(met_id(model)), value=F)) > 0)
            obj.react <- react_id(model)[which(S(model)[grep("biomass",
                                                         tolower(met_id(model)),
                                                         value=F), ]
                                           > 0)]
	else 
	    obj.react <- react_id(model)[grep("biomass", tolower(react_id(model)))]
    } else if (!is.character(obj.react)) {
        stop("argument obj.react must be character!")
    }
    if (length(obj.react) < 1) {
        stop("cannot determine obj.react producing BIOMASS!")
    } else if (length(obj.react) > 1) {
        stop("too many obj.react producing BIOMASS!")
    }

    ##-----

    recoscores <- mclapply(1L:length(fn), mc.cores=mc.cores1, function(i) {
        ##- submodels for each replicate ----
        ratios.GS <- mclapply(1L:rep.num, mc.cores=1, function(j) {
            ratio.GS <- NULL
            esg <- character(0)
            esr <- character(0)
            
            ##- build a submodel while removing genes.del[[i]][[j]] genes
            submod <- rmGenes(model=model,
                              genes=genes.del[[i]][[j]])

            ##- build the submodel by propagation ----
            if (suppressWarnings(checkOptSol(optimizeProb(submod, algorithm="fba", retOptSol=T),
                                             onlywarn=T))) {
                submod.fd.1 <- obj.react
                submod.fd.0 <- setdiff(react_id(submod), submod.fd.1)
                CONTINUE <- T
                iter <- 0
                ##- the sequential FBAs are used to identify reactions with a non-null flux, then
                ##- unnecessary to perform FVA on them
                while (CONTINUE && iter < iter.fba) {
                    submod.obj  <- changeObjFunc(submod, 
                                                 react=submod.fd.1,
                                                 obj_coef=rep(1, length(submod.fd.1)))
                    submod.fd   <- setNames(getFluxDist(optimizeProb(submod.obj, algorithm="fba", 
                                                                     lpdir="min", retOptSol=T)), 
                                            react_id(submod.obj))
                    submod.fd.1 <- union(submod.fd.1, names(submod.fd[abs(submod.fd) >= 1e-4]))
                    CONTINUE    <- (length(intersect(submod.fd.0, submod.fd.1)) > 0)
                    submod.fd.0 <- setdiff(submod.fd.0, submod.fd.1)
                    iter        <- iter + 1
                }
                ##- FVA should not be performed with objective on RECOs since RECOs would be blocked
                ##- reset the objective function to the initial (BIOMASS by default) before FVA
                submod.obj <- changeObjFunc(submod, react=obj.react, obj_coef=1)
                if (mc.cores2 > 2) {
                    subfva.obj <- suppressMessages(multiDel(submod.obj,
                                                            nProc=mc.cores2,
                                                            todo="fluxVar",
                                                            del1=submod.fd.0))
                    reacs.blocked <- setNames(unlist(mclapply(subfva.obj, mc.cores=mc.cores2, blReact)),
                                              submod.fd.0)
                } else {
                    subfva.obj <- suppressMessages(fluxVar(submod.obj, react=submod.fd.0, verboseMode=0))
                    reacs.blocked <- setNames(blReact(subfva.obj), submod.fd.0)
                }
                reacs.blocked.id <- names(reacs.blocked[reacs.blocked==T])
                submod <- rmReact(submod, reacs.blocked.id, rm_met=T)
                ##-----
            }

            gene.sub <- sybil::allGenes(submod)
            reac.sub <- sybil::react_id(submod)
            ##- gene set behavior ----
            if (!is.null(gene.sets)) {
                ratio.GS <- mclapply(gene.sets, mc.cores=mc.cores2, function(x) {
                    gs <- length(intersect(x, gene.sub)) / length(x)
                    return (gs)
                })
            }
            ##-----
            ratio.GS$gene.sub <- gene.sub
            ratio.GS$reac.sub <- reac.sub
                
            ##- essential genes and reactions in submodels ----
            if (essential && suppressWarnings(checkOptSol(optimizeProb(submod, algorithm="fba", retOptSol=T),
                                                          onlywarn=T))) {
                ## NB: dual off to infinity error comes from the two operations below
                ##     this is due to futile cycles (observed in CANAL),
                ##     which make the FVA dependent on SYBIL MAXIMUM!
                subfba <- optimizeProb(submod)
                
                ##- find essential genes ----
                if (length(sybil::allGenes(submod)) > 0) {
                    ogd <- suppressMessages(multiDel(submod,
                                                     nProc=mc.cores2,
                                                     todo="oneGeneDel",
                                                     checkOptSolObj=T))
                    if (sum(sapply(ogd, is.null)) == 0) {
                        lpokg <- unlist(lapply(ogd, lp_ok)) #check for crash
                        esg <- sybil::allGenes(submod)[lpokg == 0 &
                                                       unlist(mclapply(ogd, mc.cores=mc.cores2, lp_obj))
                                                       > lp_obj(subfba) + SYBIL_SETTINGS("TOLERANCE")]
                    }
                                        # lpokg <- unlist(lapply(1:length(ogd), function(ii) {
                        #     if (!is.null(ogd[[ii]])) {
                        #         return (lp_ok(ogd[[ii]]))
                        #     }
                        #     if (names(ogd)[ii] != "") {
                        #         r1 <- strsplit(gsub(names(ogd)[ii], pattern="[\\(\\]]", replacement="", perl=T), ",")[[1]][1]
                        #         r2 <- strsplit(gsub(names(ogd)[ii-1], pattern="[\\(\\]]", replacement="", perl=T), ",")[[1]][2]
                        #         return (rep(1, as.numeric(r1)-as.numeric(r2)))
                        #     } else if (names(ogd)[ii-1] != "") {
                        #         r2 <- strsplit(gsub(names(ogd)[ii-1], pattern="[\\(\\]]", replacement="", perl=T), ",")[[1]][2]
                        #         return (rep(1, length(sybil::allGenes(submod))-r2+1))
                        #     }
                        #     ## TODO: more complex situations may happen with NULL results from dual off to infinity!
                        #     return (character(0))
                        # }))
                        # esg <- sybil::allGenes(submod)[lpokg == 0 &
                        #                                unlist(mclapply(1:length(ogd), mc.cores=mc.cores2, function(ii) {
                        #                                    if (!is.null(ogd[[ii]])) {
                        #                                        return (lp_obj(ogd[[ii]]))
                        #                                    }
                        #                                    if (names(ogd)[ii] != "") {
                        #                                        r1 <- strsplit(gsub(names(ogd)[ii],
                        #                                                            pattern="[\\(\\]]",
                        #                                                            replacement="", perl=T), ",")[[1]][1]
                        #                                        r2 <- strsplit(gsub(names(ogd)[ii-1],
                        #                                                            pattern="[\\(\\]]",
                        #                                                            replacement="", perl=T), ",")[[1]][2]
                        #                                        return (rep(1+lp_obj(subfba), as.numeric(r1)-as.numeric(r2)))
                        #                                    } else if (names(ogd)[ii-1] != "") {
                        #                                        r2 <- strsplit(gsub(names(ogd)[ii-1],
                        #                                                            pattern="[\\(\\]]",
                        #                                                            replacement="", perl=T), ",")[[1]][2]
                        #                                        return (rep(1+lp_obj(subfba), length(sybil::allGenes(submod))-r2+1))
                        #                                    }
                        #                                    ## TODO: more complex situations may happen with NULL results from dual off to infinity!
                        #                                    return (character(0))
                        #                                }))
                        #                                > lp_obj(subfba) + SYBIL_SETTINGS("TOLERANCE")]
                }
                ##-----
                
                ##- find essential reactions ----
                if (react_num(submod) > 0) {
                    ofd <- suppressMessages(multiDel(submod,
                                                     nProc=mc.cores2,
                                                     todo="oneFluxDel",
                                                     del1=react_id(submod),
                                                     checkOptSolObj=T))
                    if (sum(sapply(ofd, is.null)) == 0) {
                        lpokr <- unlist(lapply(ofd, lp_ok)) #check for crash
                        esr <- react_id(submod)[lpokr == 0 &
                                                unlist(mclapply(ofd, mc.cores=mc.cores2, lp_obj))
                                                > lp_obj(subfba) + SYBIL_SETTINGS("TOLERANCE")]
                    }
                                        # lpokr <- unlist(lapply(1:length(ofd), function(ii) {
                        #     if (!is.null(ofd[[ii]])) {
                        #         return (lp_ok(ofd[[ii]]))
                        #     }
                        #     if (names(ofd)[ii] != "") {
                        #         r1 <- strsplit(gsub(names(ofd)[ii], pattern="[\\(\\]]", replacement="", perl=T), ",")[[1]][1]
                        #         r2 <- strsplit(gsub(names(ofd)[ii-1], pattern="[\\(\\]]", replacement="", perl=T), ",")[[1]][2]
                        #         return (rep(1, as.numeric(r1)-as.numeric(r2)))
                        #     } else if (names(ofd)[ii-1] != "") {
                        #         r2 <- strsplit(gsub(names(ofd)[ii-1], pattern="[\\(\\]]", replacement="", perl=T), ",")[[1]][2]
                        #         return (rep(1, react_num(submod)-r2+1))
                        #     }
                        #     ## TODO: more complex situations may happen with NULL results from dual off to infinity!
                        #     return (character(0))                   
                        # }))
                        # esr <- react_id(submod)[lpokr == 0 &
                        #                         unlist(mclapply(1:length(ofd), mc.cores=mc.cores2, function(ii) {
                        #                             if (!is.null(ofd[[ii]])) {
                        #                                 return (lp_obj(ofd[[ii]]))
                        #                             }
                        #                             if (names(ofd)[ii] != "") {
                        #                                 r1 <- strsplit(gsub(names(ofd)[ii],
                        #                                                     pattern="[\\(\\]]",
                        #                                                     replacement="", perl=T), ",")[[1]][1]
                        #                                 r2 <- strsplit(gsub(names(ofd)[ii-1],
                        #                                                     pattern="[\\(\\]]",
                        #                                                     replacement="", perl=T), ",")[[1]][2]
                        #                                 return (rep(1+lp_obj(subfba), as.numeric(r1)-as.numeric(r2)))
                        #                             } else if (names(ofd)[ii-1] != "") {
                        #                                 r2 <- strsplit(gsub(names(ofd)[ii-1],
                        #                                                     pattern="[\\(\\]]",
                        #                                                     replacement="", perl=T), ",")[[1]][2]
                        #                                 return (rep(1+lp_obj(subfba), react_num(submod)-r2+1))
                        #                             }
                        #                             ## TODO: more complex situations may happen with NULL results from dual off to infinity!
                        #                             return (character(0))
                        #                         }))
                        #                         > lp_obj(subfba) + SYBIL_SETTINGS("TOLERANCE")]
                }
                ##-----
            }
            ratio.GS$esg <- esg
            ratio.GS$esr <- esr
            
            return (ratio.GS)
        })
        genes.sub <- setNames(mclapply(ratios.GS, mc.cores=mc.cores2, function(rgs) {
            rgs$gene.sub
        }), reps)
        reacs.sub <- setNames(mclapply(ratios.GS, mc.cores=mc.cores2, function(rgs) {
            rgs$reac.sub
        }), reps)
        esgs <- setNames(mclapply(ratios.GS, mc.cores=mc.cores2, function(rgs) {
            rgs$esg
        }), paste0("esg.", reps))
        esrs <- setNames(mclapply(ratios.GS, mc.cores=mc.cores2, function(rgs) {
            rgs$esr
        }), paste0("esr.", reps))
        ratios.GS <- setNames(mclapply(ratios.GS, mc.cores=mc.cores2, function(rgs) {
            unlist(rgs[-c(length(rgs)-3L:0L)])
        }), paste0("gs.", reps))
        ##-----

        return (list(score=c(
                         fn[[i]]$score,
                         sapply(esgs, length)/length(sybil::allGenes(model)),
                         sapply(esrs, length)/react_num(model),
                         gene.del=fn[[i]]$gene.del,
                         unlist(ratios.GS)),
                     genes.sub=genes.sub,
                     reacs.sub=reacs.sub))
    })
    sub.genes <- mclapply(recoscores, mc.cores=mc.cores1, function(rcs) {
        rcs$genes.sub
    })
    sub.reacs <- mclapply(recoscores, mc.cores=mc.cores1, function(rcs) {
        rcs$reacs.sub
    })
    recoscore <- do.call(cbind, mclapply(recoscores, mc.cores=mc.cores1, function(rcs) {
        rcs$score
    }))
    colnames(recoscore) <- recoscore["gene.del", , drop=T]
    ##-----

    ranks.name <- union(rank.best, fn[[1]]$ranks.name) # move the best ranking to the first one
    ranks.num  <- length(ranks.name)
    draw.num   <- fn[[1]]$draw.num
    GS.num     <- length(gene.sets)
    res <- scoreGeneDel(
        model             = model,
        condition         = condition,
        fitness.random    = recoscore[grep("^F.random", rownames(recoscore)), , drop=F],
        fitness.ranks     = setNames(lapply(1L:rep.num, function(j) {
            rsj <- recoscore[paste("F", ranks.name, reps[j], sep='.'), , drop=F]
            rownames(rsj) <- ranks.name
            return (rsj)
        }), reps),
        fitness.id.random = recoscore[grep("^F.id.random", rownames(recoscore)), , drop=F],
        fitness.id.ranks  = setNames(lapply(1L:rep.num, function(j) {
            rsj <- recoscore[paste("F.id", ranks.name, reps[j], sep='.'), , drop=F]
            rownames(rsj) <- ranks.name
            return (rsj)
        }), reps),
        ess.gene          = recoscore[grep("^esg", rownames(recoscore)), , drop=F],
        ess.reaction      = recoscore[grep("^esr", rownames(recoscore)), , drop=F],
        gene.del          = recoscore["gene.del", , drop=F],
        gene.sets         = gene.sets,
        ratio.GS          = setNames(lapply(1L:rep.num, function(j) {
            rgsj <- recoscore[paste("gs", reps[j], names(gene.sets), sep='.'), , drop=F]
            rownames(rgsj) <- names(gene.sets)
            return (rgsj)
        }), reps),
        sub.genes         = sub.genes,
        sub.reacs         = sub.reacs,
        rescue.met        = recoscore[recos, , drop=F]
        )
    return (res)
}


#' Plot fitness of submodels built by gene removal in a condition
#'
#' This function plots the fitness of submodels built by gene removal in a condition with different rankings.
#' @param sgd An object of class \code{scoreGeneDel}.
#' @param mc.cores The number of cores to use (at least 1), i.e. at most how many child processes will be run simultaneously. Default: 1.
#' @param ranks.name Names of gene expression ranking. Default: NULL.
#' @param njt An object of class \code{phylo} for colored plot of fitness weighting schema resulting from \code{weightReacts}. Default: NULL.
#' @param cols Colors for conditions. Default: rainbow colors.
#' @param ltys Line types for conditions. Default: incrementing line types in R.
# #' @param cutoff A numeric value for the cutoff of gene removal. 
# #' If \code{cutoff} <= 1, apply on removal fitness scores; if \code{cutoff} > 1, apply on number of removed genes. 
# #' Default: NA, \code{cutoff} at fitness^opt.degree/remaining_genes peak.
# #' @param opt.percentage The minimum percentage of genes to remove, used when \code{cutoff = NA}. Default: 10\%.
# #' @param opt.degree A numeric value for the power of optimization score, used when \code{cutoff = NA}. Default: 2.
#' @import sybil grDevices graphics
#' @examples
#' data(yarliSubmnets)
#' simulateSubmnet(yarliSubmnets$DN)
#' @export
simulateSubmnet <- function(sgd, mc.cores = 1, ranks.name = NULL, njt = NULL, cols = NULL, ltys = NULL) {
    ##- settings ----
    draw.num   <- nrow(sgd$fitness.random)
    if (draw.num < 1L) return (0)
    gene.num   <- length(sybil::allGenes(sgd$model))
    gene.del   <- as.vector(sgd$gene.del)
    param.num  <- length(gene.del)
    rep.num    <- length(sgd$fitness.ranks)
    legend.pos <- "topright"
    if (is.null(ranks.name)) {
        ranks.name  <- rownames(sgd$fitness.ranks[[1]])
    }
    ranks.name <- c("random", ranks.name)
    condition  <- sgd$condition
    misc.plot  <- F
    cutoff     <- gene.num*0.1 
    opt.percentage <- 0.1
    opt.degree <- 2
    if (length(cols) == 0) {
        cols         <- c("black", rainbow(nrow(sgd$fitness.ranks[[1]])))
        if (length(cols) == 8L) {
            cols[3]  <- "darkorange1"
            cols[4]  <- "gold3"
        }
    }
    if (length(cols) != nrow(sgd$fitness.ranks[[1]])+1) {
        stop("Number of colors does not match number of conditions.")
    }
    if (length(ltys) == 0) {
        ltys         <- c(4, 1, rep(5, nrow(sgd$fitness.ranks[[1]])-1))
    }
    if (length(ltys) != nrow(sgd$fitness.ranks[[1]])+1) {
        stop("Number of line types does not match number of conditions.")
    }

    ##-----
    
    ##- average fitness in random and ranks draws ----
    fitness.means <- mclapply(1L:rep.num, mc.cores=mc.cores, function(i) {
        fitness.mean    <- cbind(apply(sgd$fitness.random, 2, function(x) {mean(na.omit(x))}),
                                 t(sgd$fitness.ranks[[i]]))
        colnames(fitness.mean) <- ranks.name
        rownames(fitness.mean) <- gene.del
        return (fitness.mean)
    })
    fitness.id.means <- mclapply(1L:rep.num, mc.cores=mc.cores, function(i) {
        fitness.id.mean <- cbind(apply(sgd$fitness.id.random, 2, function(x) {mean(na.omit(x))}),
                                 t(sgd$fitness.id.ranks[[i]]))
        colnames(fitness.id.mean) <- ranks.name
        rownames(fitness.id.mean) <- gene.del
        return (fitness.id.mean)
    })
    ##-----
    
    # ##- standard deviation in random ----
    # fitness.sd    <- apply(sgd$fitness.random, 2, function(x) {sd(na.omit(x))})
    # fitness.id.sd <- apply(sgd$fitness.id.random, 2, function(x) {sd(na.omit(x))})
    # #print(paste("SD mean:", min(fitness.sd), max(fitness.sd), sep=' '))
    # #print(paste("SD id mean:", min(fitness.id.sd), max(fitness.id.sd), sep=' '))
    # fitness.cv    <- fitness.sd / fitness.mean[, 1]
    # fitness.id.cv <- fitness.id.sd / fitness.id.mean[, 1]
    # fitness.cv[is.na(fitness.cv)]       <- 0
    # fitness.id.cv[is.na(fitness.id.cv)] <- 0
    # #print(paste("CV mean:", min(fitness.cv), max(fitness.cv), sep=' '))
    # #print(paste("CV id mean:", min(fitness.id.cv), max(fitness.id.cv), sep=' '))
    # ##-----
    
    # ##- fitness over number of remaining genes ----
    # ratio.mat <- cbind(apply(sgd$fitness.random^opt.degree/rep(gene.num-sgd$gene.del,
    #                                                           each=nrow(sgd$fitness.random)),
    #                          2,
    #                          function(x) {mean(na.omit(x))}),
    #                    # apply(sgd$fitness.ranked^opt.degree/(gene.num-sgd$gene.del),
    #                    #       2,
    #                    #       function(x) {mean(na.omit(x))}),
    #                    matrix(unlist(lapply(sgd$fitness.ranks,
    #                                         function(x) {
    #                                             apply(x^opt.degree/(gene.num-sgd$gene.del),
    #                                                   2,
    #                                                   function(x) {mean(na.omit(x))})
    #                                         })), ncol=length(sgd$fitness.ranks), byrow=F))
    # colnames(ratio.mat) <- expr.names
    # ##-----
    
    ##- percentages of random draws better (i.e. higher fitness) than ranked draws ----
    if (sum(sgd$fitness.random) == 0) {
        weights <- rep(1/param.num, param.num)
    } else {
        weights <- apply(sgd$fitness.random, 2, mean) / sum(apply(sgd$fitness.random, 2, mean))
    }
    r2random <- lapply(1L:rep.num, function(i) {
        ranks2random  <- apply(sgd$fitness.ranks[[i]], 1, function(x) {
            100 * sum(sapply(1L:param.num, function(k) {
                    sum(x[k] < sgd$fitness.random[, k, drop=T]) / draw.num * weights[k]
                }))
        })
        return (setNames(c(100, ranks2random), ranks.name))
    })
    ##-----

    ##- make figures ----
    colors.light <- rep(0.4, 3)
    ##- coordinates for two-percentile plots ----
    pbs <- c(20, 80)
    pbs.coords <- sapply(1L:length(pbs), function(k) {
        xcoords <- gene.del
        ycoords <- apply(sgd$fitness.random, 2, function(x) {
            quantile(na.omit(x), probs=pbs/100)[k]
        })
        return (c(xcoords, ycoords))
    })
    pbs.id.coords <- sapply(1L:length(pbs), function(k) {
        xcoords <- gene.del
        ycoords <- apply(1-sgd$fitness.id.random, 2, function(x) {
            quantile(na.omit(x), probs=pbs/100)[k]
        })
        return (c(xcoords, ycoords))
    })
    
    ##- rows ~ #removed genes
    ##- 2 first columns ~ (xcoords, ycoords) of pbs[1]
    ##- 2 last columns  ~ (xcoords, ycoords) of pbs[2]
    pbs.coords <- matrix(pbs.coords, ncol=2L*length(pbs))
    rownames(pbs.coords) <- gene.del
    pbs.id.coords <- matrix(pbs.id.coords, ncol=2L*length(pbs))
    rownames(pbs.id.coords) <- gene.del
    ##-----
    
    ##- fitness of ranks versus random ----
    pdf(file=paste(condition, "smooth.pdf", sep='_'), width=9, height=6)
    cex.lab    <- 2
    cex.axis   <- 1.6
    cex.leg    <- 1.8
    line.lab   <- 3
    cex.panel  <- 2.8
    at.panel   <- -55
    line.panel <- -2
    invisible(sapply(1L:rep.num, function(i) {
        par(mar=c(0.1,0.1,0.5,0.5), oma=c(4.5,4.5,0,0))
        smoothScatter(rep(gene.del, each=nrow(sgd$fitness.random)),
                      as.vector(sgd$fitness.random),
                      #main    = paste("Quality of ", condition, " submodel", sep=''),
                      #colramp = colorRampPalette(c("white", "gray80", "gray40")),
                      xlab     = "",
                      ylab     = "",
                      cex.axis = cex.axis,
                      col      = cols[1],
                      ylim     = c(0, 1))
        matpoints(gene.del, fitness.means[[i]],
                  type = 'l',
                  lwd  = 2,
                  lty  = ltys,
                  col  = cols)
        lines(polygon(c(pbs.coords[, 1], rev(pbs.coords[, 1+2])),
                      c(pbs.coords[, 2], rev(pbs.coords[, 2+2])),
                      col    = rgb(colors.light[1], colors.light[2], colors.light[3], 0.5),
                      border = NA)
        )
        #lines(gene.del, sgd$ess.gene, lty=3, lwd=2, col="coral")
        #lines(gene.del, sgd$ess.reaction, lty=3, lwd=2, col="cyan2")
        mtext("fitness", side=2, line=line.lab, outer=F, cex=cex.lab)
        mtext("#removed genes", side=1, line=line.lab, cex=cex.lab)
        mtext(ifelse(condition=="WN", "A",
                     ifelse(condition=="WH", "B",
                            ifelse(condition=="SN", "C",
                                   ifelse(condition=="SH", "D",
                                          ifelse(condition=="DN", "E",
                                                 ifelse(condition=="UN", "F", "")))))),
              side=3, line=line.panel, outer=F, at=at.panel, cex=cex.panel, font=2)
        legend(legend.pos,
               c(names(r2random[[i]])[1L],
                 paste(names(r2random[[i]][-1L]), ": ",
                       round(r2random[[i]][-1L], digits=2),
                       "%",
                       sep='')),
               lty=ltys, lwd=2, col=cols, cex=cex.leg)
    }))
    invisible(dev.off())
    ##-----
    
    ##- miscellaneous plots ----
    if (misc.plot) {
        # ##- fitness over numbers of remaining genes vs numbers of deleted genes ----
        # pdf(file=paste0(condition, "opt.pdf"), width=10, height=6)
        # matplot(gene.del, ratio.mat,
        #         xlab = "#removed genes",
        #         ylab = paste("fitness^", opt.degree, " / #remaining genes", sep=''),
        #         type = 'l',
        #         lty  = ltys,
        #         col  = cols)
        # legend(legend.pos, names(ltys), lty=ltys, col=cols, cex=0.8)
        # invisible(dev.off())
        # ##-----
        
        ##- only random fitness ----
        pdf(file=paste(condition, "weightrandom.pdf", sep='_'), width=9, height=6)
        cex.lab  <- 2
        cex.axis <- 1.6
        cex.leg  <- 1.8
        line.lab <- 3
        par(mar=c(0.1,0.1,0.5,0.5), oma=c(4.5,4.5,0,0))
        smoothScatter(rep(gene.del, each=nrow(sgd$fitness.random)),
                      as.vector(sgd$fitness.random),
                      xlab     = "",
                      ylab     = "",
                      cex.axis = cex.axis,
                      col      = cols[1],
                      ylim     = c(0, 1))
        matpoints(gene.del, fitness.means[[1]][,1],
                  type = 'l',
                  lwd  = 2,
                  lty  = ltys,
                  col  = cols)
        lines(polygon(c(pbs.coords[, 1], rev(pbs.coords[, 1+2])),
                      c(pbs.coords[, 2], rev(pbs.coords[, 2+2])),
                      col    = rgb(colors.light[1], colors.light[2], colors.light[3], 0.5),
                      border = NA)
        )
        mtext("fitness", side=2, line=line.lab, outer=F, cex=cex.lab)
        mtext("#removed genes", side=1, line=line.lab, cex=cex.lab)
        invisible(dev.off())
        ##-----
    }
    ##-----
    
    ##- fitness of random gene removal without and with weighting schema ----
    pdf(file=paste(condition, "random.pdf", sep='_'), width=24, height=16)
    cex.main   <- 2.4
    cex.lab    <- 2.5
    cex.axis   <- 2.9
    cex.leg    <- 1.8
    cex.panel  <- 3
    cex.clust  <- 3
    at.panel   <- -57
    line.panel <- 0
    line.lab   <- 5
    line.main  <- 1
    ncol.leg   <- 3
    
    par(mar=c(1,0.1,3,2), oma=c(5.5,8,1.5,0), bg="gray100")
    if (is.null(njt)) {
        layout(matrix(c(rep(1,ncol.leg),
                        rep(2,ncol.leg),
                        c(3:(3+ncol.leg-1))),
                      3, ncol.leg, byrow=T))
    } else {
        layout(rbind(c(1,1,1,1,1,3,3,3,3,3,4), c(2,2,2,2,2,3,3,3,3,3,4)))
    }
    
    #fitness.id.sd.mean <- signif(mean(fitness.id.sd), 3)
    smoothScatter(rep(gene.del, each=nrow(sgd$fitness.id.random)),
                  1-as.vector(sgd$fitness.id.random),
                  #main    = bquote(list("Without weighting scheme", bar(sigma) == .(fitness.id.sd.mean))),
                  #colramp = colorRampPalette(c("white", "blue", "red")),
                  xlab     = "",
                  ylab     = "",
                  xaxt     = "n",
                  cex.axis = cex.axis,
                  col      = cols[1],
                  ylim     = c(0, 1))
    mtext("A", side=3, line=line.panel, outer=F, at=at.panel, cex=cex.panel, font=2)
    mtext("fraction of rescued reactions", side=2, line=line.lab, outer=F, cex=cex.lab)
    # mtext("Without weighting scheme",
    #       side=3, line=line.main, outer=F, cex=cex.main, font=2)
    points(gene.del, 1-fitness.id.means[[2]][, 1],
           type = 'l',
           lwd  = 2,
           lty  = 1,
           col  = 1)
    lines(polygon(c(pbs.id.coords[, 1], rev(pbs.id.coords[, 1+2])),
                  c(pbs.id.coords[, 2], rev(pbs.id.coords[, 2+2])),
                  col    = rgb(colors.light[1], colors.light[2], colors.light[3], 0.5),
                  border = NA)
    )
    
    # fitness.sd.mean <- signif(mean(fitness.sd), 3)
    # print(fitness.sd.mean)
    # smoothScatter(rep(as.vector(sgd$gene.del), each=nrow(sgd$fitness.random)),
    #               as.vector(sgd$fitness.random),
    #               #main    = bquote(list("With weighting scheme", bar(sigma) == .(fitness.sd.mean))),
    #               #colramp = colorRampPalette(c("white", "gray80", "gray40")),
    #               xlab    = "",
    #               ylab    = "",
    #               xaxt    = "n",
    #               cex.axis = cex.axis,
    #               ylim    = c(0, 1))
    # mtext("B", side=3, line=line.panel, outer=F, at=at.panel, cex=cex.panel, font=2)
    # mtext("fitness", side=2, line=line.lab, outer=F, cex=cex.lab)
    # # mtext("With weighting scheme",
    # #       side=3, line=line.main, outer=F, cex=cex.main, font=2)
    # points(genesDel.mean, fitness.mean[, 1],
    #        type = "l",
    #        lwd  = 2,
    #        lty  = 1,
    #        col  = 1)
    # lines(polygon(c(pbs.coords[, 1], rev(pbs.coords[, 1+2])),
    #               c(pbs.coords[, 2], rev(pbs.coords[, 2+2])),
    #               col    = rgb(colors.light[1], colors.light[2], colors.light[3], 0.5),
    #               border = NA)
    # )
    
    metplot <- rownames(sgd$rescue.met)
    nrow.leg <- ceiling(length(metplot)/ncol.leg)
    rescue.met.id <- gsub(metplot, pattern="RECO_PUSH_", replacement="")
    joinids <- sapply(rescue.met.id, function(metid) {
        joinid <- gregexpr("_", metid)[[1]]
        return (joinid[length(joinid)])
    })
    substr(rescue.met.id, joinids, joinids) <- '['
    rescue.met.id <- gsub(rescue.met.id, pattern="$", replacement="]")
    rescue.met.name <- sapply(rescue.met.id, function(metid) {
        gsub(gsub(strsplit(met_name(sgd$model)[which(met_id(sgd$model) == metid)], split=', ')[[1]][1],
                  pattern=" (yeast specific)", replacement="", fixed=T),
             pattern="L-", replacement="")
    })
    colors.met <- rainbow(length(metplot)+5, s=1, v=0.7, alpha=1)
        
    ##- sort curve at 50% rescue to create colors red->violet from left to right
    rescuedhalf <- apply(sgd$rescue.met, 1, function(met) {
        min(c(length(met), which(met > 0.5)))
    })
    metplot       <- metplot[order(rescuedhalf, decreasing=F)]
    rescue.met.name <- rescue.met.name[order(rescuedhalf, decreasing=F)]
    names(colors.met) <- rescue.met.name
    ltys <- setNames(c(rep(1L:6L, times=floor(length(metplot)/6)), 1L:(length(metplot)-6*floor(length(metplot)/6))),
                     rescue.met.name)
    ss <- apply(sgd$rescue.met[metplot, , drop=F], 1, function(met) {
        smooth.spline(gene.del, met)
    })
    matplot(matrix(unlist(lapply(ss, function(s) s$x)), ncol=length(ss)),
            matrix(unlist(lapply(ss, function(s) s$y)), ncol=length(ss)),
            type = "l",
            lwd  = 2,
            lty  = ltys,
            xlab = "",
            ylab = "",
            #xaxt = 'n',
            cex.axis = cex.axis,
            col  = colors.met)
    mtext("B", side=3, line=line.panel, outer=F, at=at.panel, cex=cex.panel, font=2)
    mtext("rescued average", side=2, line=line.lab, cex=cex.lab, outer=F)
    mtext("#removed genes",  side=1, line=line.lab, cex=cex.lab, outer=F)
    
    if (is.null(njt)) {
        sapply(0:(ncol.leg-1), function(li) {
            plot.new()
            legend("bottom",
                   rescue.met.name[(li*nrow.leg+1):min(length(metplot), (li*nrow.leg+nrow.leg))],
                   col=colors.met[(li*nrow.leg+1):min(length(metplot), (li*nrow.leg+nrow.leg))],
                   lty=ltys[(li*nrow.leg+1):min(length(metplot), (li*nrow.leg+nrow.leg))],
                   bty='n', lwd=2, cex=cex.leg)
        })
    }
    if (is.null(njt)) {
        invisible(dev.off())
    } else {
        ##-----
        ##- plot cluster of RECO reactions with weights and colors ----
        #pdf(paste0(condition, "_recoclusterweight.pdf"), width=8, height=12)
        #layout(rbind(c(1,1,2), c(1,1,2), c(1,1,2)))
        tip.color <- colors.met[njt$tip.label]
        tip.label <- njt$tip.label
        njt$tip.label <- paste(njt$tip.label, njt$weight, sep=" # ")
        plot(njt, cex=cex.clust, tip.color=tip.color, align.tip.label=T, no.margin=T)
        tip.order <- rev(order.dendrogram(as.dendrogram(as.hclust(njt))))
        mtext("C", side=3, line=line.panel-3, outer=F, at=0, cex=cex.panel, font=2)
        plot.new()
        sapply(1L:length(tip.label[tip.order]), function(itl) {
            tl <- tip.label[tip.order][itl]
            ytl <- 1.02 - itl*0.0226
            points(c(0.4, 0.8), rep(ytl, 2),
                   col=colors.met[tl],
                   type='l',
                   lty=ltys[tl],
                   lwd=3)
            #text(1.2, ytl, tl, col=colors.met[tl])
        })
        
        #mtext("B", side=3, line=line.panel, outer=F, at=at.panel, cex=cex.panel, font=2)
        invisible(dev.off())
        save(colors.met, file=paste0(condition, "_colors.met.Rdata"))
    }
    ##-----
    
    # ##- standard deviation and coefficient of variation ----
    # pdf(file=paste(condition, "sdcv.pdf", sep='_'), width=10, height=10)
    # par(mfrow=c(2,1))
    # matplot(cbind(genesDel.mean, genesDel.mean), cbind(fitness.sd, fitness.id.sd),
    #         xlab="#removed genes",
    #         ylab="standard deviation",
    #         type='l', col=c(1,2), lty=c(1,2))
    # matplot(cbind(genesDel.mean, genesDel.mean), cbind(fitness.cv, fitness.id.cv),
    #         xlab="#removed genes",
    #         ylab="coefficient of variation",
    #         type='l', col=c(1,2), lty=c(1,2))
    # invisible(dev.off())
    # ##-----
    ##-----
    
    ##- number of genes to cut ----
    if (is.na(cutoff)) {                # optimal cut
        # genedel <- apply(ratio.mat, 2, function(x) {
        #     gene.del[floor(length(x)*opt.percentage) +
        #                       which.max(x[(floor(length(x)*opt.percentage)+1):(floor(length(x)/2))])]
        # })
        genedel <- NULL
    } else {
        if (cutoff > 1) {               # cut on number of genes
            genedel <- lapply(1L:rep.num, function(i) {
                setNames(rep(floor(cutoff), ncol(fitness.means[[i]])), ranks.name)
            })
        } else {                        # cut on fitness
            genedel <- lapply(1L:rep.num, function(i) {
                apply(fitness.means[[i]], 2, function(x) {
                    gene.del[min(which(x < cutoff)) - 1]
                })
            })
        }
    }
    ##-----

    return (NULL)
}


#' Generate a submodel by removing genes
#'
#' This functions creates a submodel by removing genes from a given model. It is similar to \code{deleteModelGenes} from the COBRA Toolbox.
#' @param model An object of class \code{modelorg}.
#' @param genes A vector of genes to remove.
#' @return The submodel.
#' @import sybil
#' @importFrom utils packageVersion
#' @export
rmGenes <- function(model, genes) {
    if (!is(model, "modelorg")) {
        stop("needs an object of class modelorg!")
    }
    if (!checkVersion(model)) {
        stop("model is of wrong version!")
    }
    if (packageVersion("sybil") < '2.0.1') {
        stop("Bug found in sybil, please update to sybil >= 2.0.4!")
    }
    if (length(genes) == 0) {
        return (model)
    }
    
    ## debug <- packageVersion("sybil") < '2.0.1'
    
    #- reactions affected by removing genes
    reacs.ind <- geneDel(model, genes, checkId=T)
    if (!is.null(reacs.ind)) {
        reacs  <- react_id(checkReactId(model, reacs.ind))
        submod <- rmReact(model, react=reacs, rm_met=T)
        # if (debug) {
        #     modelorg2tsv(submod, prefix, makeClosedNetwork=T, quote=F)
        #     submod <- suppressMessages(
        #         readTSVmod(prefix=prefix,
        #                    quoteChar="",
        #                    mergeMet=F, balanceReact=F, remUnusedMetReact=F))
        #     file.remove(list.files(pattern=prefix))
        # }
    } else {
        submod <- model
    }

    genes.del <- intersect(genes, sybil::allGenes(submod))
    while (length(genes.del) > 0) {
        gen <- genes.del[1]
        gen.ind       <- which(sybil::allGenes(submod) == gen)
        reacs.del.ind <- which(rxnGeneMat(submod)[, gen.ind, drop=T] != 0)
        reacs.del     <- react_id(checkReactId(submod, reacs.del.ind))

        #- replace iteratively each reaction from reacs.del in submod
        for (rea in reacs.del) {
            mat       <- S(submod)
            rea.ind   <- which(react_id(submod) == rea)
            rind      <- which(mat[, rea.ind] > 0)
            lind      <- which(mat[, rea.ind] < 0)
            rmet.id   <- met_id(submod)[rind]
            lmet.id   <- met_id(submod)[lind]
            rmet.name <- met_name(submod)[rind]
            lmet.name <- met_name(submod)[lind]
            rmet.comp <- met_comp(submod)[rind]
            lmet.comp <- met_comp(submod)[lind]
            gprs      <- unlist(strsplit(gsub(gsub(gpr(submod)[rea.ind],
                                                   pattern='(',
                                                   replacement='',
                                                   fixed=T),
                                              pattern=")",
                                              replacement="",
                                              fixed=T),
                                         split=" or ",
                                         fixed=T))
            gprs      <- gprs[-c(grep(gprs, pattern=gen))]
            gpr.exp   <- paste('(',
                               paste(gprs, collapse=" or "),
                               ')',
                               sep='')
            submod    <- addReact(model      = submod,
                                  id         = "tmp",
                                  met        = c(lmet.id, rmet.id),
                                  Scoef      = mat[c(lind, rind), rea.ind],
                                  reversible = react_rev(submod)[rea.ind],
                                  lb         = lowbnd(submod)[rea.ind],
                                  ub         = uppbnd(submod)[rea.ind],
                                  obj        = obj_coef(submod)[rea.ind],
                                  gprAssoc   = gpr.exp,
                                  metName    = c(lmet.name, rmet.name),
                                  metComp    = c(lmet.comp, rmet.comp)
            )
            submod    <- rmReact(submod, react=rea, rm_met=T)
            react_id(submod)[react_num(submod)] <- rea
            # if (debug) {
            #     modelorg2tsv(submod, prefix=prefix, makeClosedNetwork=T, quote=F)
            #     submod    <- suppressMessages(
            #         readTSVmod(prefix=prefix,
            #                    quoteChar="",
            #                    mergeMet=F, balanceReact=F, remUnusedMetReact=F))
            #     file.remove(list.files(pattern=prefix))
            # }
        }
        genes.del <- intersect(genes, sybil::allGenes(submod))
    }
    return (submod)
}
