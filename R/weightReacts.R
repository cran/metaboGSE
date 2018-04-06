#' Compute distances of rescue reactions
#'
#' The function rescueDist computes the distances (similarities) between rescue reactions.
#' @param model An object of class \code{modelorg}.
# #' @param type Type of knock-out simulation: \code{g}(ene) or \code{r}(eaction). Default: "\code{g}".
#' @param mc.cores The number of cores to use (at least 1), i.e. at most how many child processes will be run simultaneously. Default: 1.
#' @param tol The maximum value to be considered null. Default: \code{SYBIL_SETTINGS("TOLERANCE")}.
#' @param gene.num The number of genes to remove. If 1, \code{oneGeneDel} will be performed and draw.num will be ignored. Default: 1.
#' @param draw.num The number of random draws. Default: 1000. It is ignored when gene.nume = 1.
#' @return An object of class \code{dist}, containing distances between rescue reactions in the given model.
#' @import sybil
#' @keywords internal
rescueDist <- function(model, mc.cores = 1, gene.num = 1, draw.num = 1000, tol = SYBIL_SETTINGS("TOLERANCE")) {
    if (!is(model, "modelorg")) {
        stop("needs an object of class modelorg!")
    }
    if (!checkVersion(model)) {
        stop("model is of wrong version!")
    }
    type     <- 'g'
    genes    <- sybil::allGenes(model)
    reacs    <- react_id(model)
    reac.num <- react_num(model)
    recos    <- grep('RECO', reacs, perl=T, value=T)
    nonrecos <- setdiff(reacs, recos)
    
    if (type == 'g') {
        set.seed(1000)
        draw.mat <- replicate(draw.num, sample(genes, size=gene.num))
        
        ##- delete each gene from the model
        if (mc.cores == 1) {
            if (gene.num == 1) {
                od <- suppressMessages(oneGeneDel(model, genes,
                                                  poCmd=list("getFluxDist"),
                                                  checkOptSolObj=F, verboseMode=1))
                ##- corresponding fluxes
                flux.list <- unlist(sybil::pa(postProc(od)), recursive=F)
                lpok <- lp_ok(od)
            } else {
                ods <- lapply(1L:ncol(draw.mat), function(i) {
                    od <- suppressMessages(geneDeletion(model, draw.mat[, i], combinations=gene.num,
                                                        poCmd=list("getFluxDist"),
                                                        checkOptSolObj=F, verboseMode=1))
                    fluxes <- unlist(sybil::pa(postProc(od)), recursive=F)
                    return (list(fluxes[[1]], lp_ok(od)))
                })
                ##- corresponding fluxes
                flux.list <- lapply(ods, function(x) {x[[1]]})
                lpok <- unlist(lapply(ods, function(x) {x[[2]]}))
            }
        } else {
            if (gene.num == 1) {
                od <- suppressMessages(multiDel(model,
                                                nProc=mc.cores,
                                                todo="oneGeneDel",
                                                poCmd=list("getFluxDist"),
                                                checkOptSolObj=F))
                ##- corresponding fluxes
                flux.list <- unlist(mclapply(od, mc.cores=mc.cores, function(x) {
                    unlist(sybil::pa(postProc(x)), recursive=F)
                }), recursive=F)
                lpok <- unlist(lapply(od, lp_ok))
            } else {
                ods <- mclapply(1L:ncol(draw.mat), mc.cores=mc.cores, function(i) {
                    od <- suppressMessages(geneDeletion(model, draw.mat[, i], combinations=gene.num,
                                                        poCmd=list("getFluxDist"),
                                                        checkOptSolObj=F, verboseMode=1))
                    fluxes <- unlist(sybil::pa(postProc(od)), recursive=F)
                    return (list(fluxes[[1]], lp_ok(od)))
                })
                ##- corresponding fluxes
                flux.list <- lapply(ods, function(x) {x[[1]]})
                lpok <- unlist(lapply(ods, function(x) {x[[2]]}))
            }
        }
        if (1 %in% lpok) {
            flux.list[which(lpok == 1)] <- rep(list(rep(1, reac.num)), length(which(lpok == 1)))
        }
        if (gene.num == 1) {
            flux.mat <- matrix(unlist((flux.list), use.names=T),
                               ncol=length(reacs), byrow=T, dimnames=list(genes, reacs))
        } else {
            flux.mat <- matrix(unlist((flux.list), use.names=T),
                               ncol=length(reacs), byrow=T, dimnames=list(1L:draw.num, reacs))
        }
    } else if (type == 'r') {
        ##- delete each reaction other than RECO from the model
        if (mc.cores == 1) {
            od <- suppressMessages(oneFluxDel(model, react=nonrecos,
                                              lb=rep(0, length(nonrecos)),
                                              ub=rep(0, length(nonrecos)),
                                              poCmd=list("getFluxDist"),
                                              checkOptSolObj=F, verboseMode=1))
            ##- corresponding fluxes
            flux.list <- unlist(sybil::pa(postProc(od)), recursive=F)
            lpok <- lp_ok(od)
        } else {
            od <- suppressMessages(multiDel(model,
                                            nProc=mc.cores,
                                            todo="oneFluxDel",
                                            del1=nonrecos,
                                            poCmd=list("getFluxDist"),
                                            checkOptSolObj=F))
            ##- corresponding fluxes
            flux.list <- unlist(mclapply(od, mc.cores=mc.cores, function(x) {
                unlist(sybil::pa(postProc(x)), recursive=F)
            }), recursive=F)
            lpok <- unlist(lapply(od, lp_ok))
        }
        if (1 %in% lpok) {
            flux.list[which(lpok == 1)] <- rep(list(rep(1, reac.num)), length(which(lpok == 1)))
        }
        flux.mat <- matrix(unlist((flux.list), use.names=T),
                           ncol=length(reacs), byrow=T, dimnames=list(nonrecos, reacs))
    } else {
        stop("Unknown type: 'g' or 'r' required.")
    }
    
    if (length(which(colnames(flux.mat) != reacs)) > 0) {
        stop("Colnames differ from reactions.")
    }
    # if (type == 'g') {
    #     if (length(which(rownames(flux.mat) != genes)) > 0) {
    #         stop("Rownames differ from genes.")
    #     }
    # } else {
    #     if (length(which(rownames(flux.mat) != nonrecos)) > 0) {
    #         stop("Rownames differ from reactions.")
    #     }
    # }
    
    ##- fluxes of RECO reactions
    flux.recos <- flux.mat[, recos, drop=F]
    
    ##- existing status of RECO fluxes = necessity of RECO reactions
    exist.recos <- apply(abs(flux.recos) > tol, 2, as.numeric)
    # print(sum(exist.recos))
    # pdf(paste0("heatmap_", gene.num, "_", draw.num, ".pdf"), width=10, height=12)
    # heatmap.2(exist.recos[rowSums(exist.recos) > 0, , drop=F], scale="none", margins = c(8, 8),
    #           trace="none", col=colorRampPalette(c("white", "darkblue"))(n = 2))
    # dev.off()
    
    ##- distance matrix of RECO existing status
    return (dist(t(exist.recos)))
}


#' Compute weights of rescue reactions
#'
#' The function weightReacts computes the weights of rescue reactions.
#' @param model An object of class \code{modelorg} indicating the weighted \code{rescue} model obtained from the rescue process.
# #' @param type Type of tree: \code{h}(ierachical) clustering (average), \code{i}(soweight), \code{n}(eibor-joining). Default: \code{h}.
#' @param mc.cores The number of cores to use (at least 1), i.e. at most how many child processes will be run simultaneously. Default: 1.
#' @param gene.num The number of genes to remove. If 1, \code{oneGeneDel} will be performed and draw.num will be ignored. Default: 1.
#' @param draw.num The number of random draws. Default: 1000.
#' @return A vector of weights for rescue reactions.
#' @import ape
#' @examples 
#' data(Ec_core)
#' mod <- rescue(Ec_core, target=0.1)
#' weightReacts(changeObjFunc(mod$rescue, react=rownames(mod$coef), obj_coef=mod$coef))
#' @export
weightReacts <- function(model, mc.cores = 1, gene.num = 1, draw.num = 1000) {
    if (!is(model, "modelorg")) {
        stop("needs an object of class modelorg!")
    }
    if (!checkVersion(model)) {
        stop("model is of wrong version!")
    }
    type <- 'h'
    dist <- rescueDist(model, mc.cores=mc.cores, gene.num=gene.num, draw.num=draw.num)
    stopifnot(!is.null(dist))
    distance <- as.matrix(dist)
    recos <- rownames(distance)
    
    ##- associated genes for reco metabolites
    rescue.met.id <- gsub(gsub(recos, pattern="RECO_PUSH_", replacement=""),
                          pattern="$", replacement="]")
    comp.pos <- regexpr("\\_[^\\_]*$", rescue.met.id)
    substr(rescue.met.id, comp.pos, comp.pos) <- '['
    rescue.met.name <- sapply(rescue.met.id, function(metid) {
        gsub(gsub(strsplit(met_name(model)[which(met_id(model) == metid)], split=', ')[[1]][1],
                  pattern=" (yeast specific)", replacement="", fixed=T),
             pattern="L-", replacement="")
    })
    GENEINFO <- F
    if (GENEINFO) {
        Smat <- as.matrix(S(model))
        rownames(Smat) <- met_id(model)
        colnames(Smat) <- react_id(model)
        rescue.met.reacs <- apply(Smat[rescue.met.id, ], 1, function(metcoef) {
            react_id(model)[which(abs(metcoef) > SYBIL_SETTINGS("TOLERANCE"))]
        })
        rxnGenes <- as.matrix(rxnGeneMat(model))
        rownames(rxnGenes) <- react_id(model)
        colnames(rxnGenes) <- sybil::allGenes(model)
        
        rescue.met.genes <- lapply(rescue.met.reacs, function(metreacs) {
            unlist(lapply(genes(model)[which(react_id(model) %in% metreacs)], function(g) {
                g[g != ""]
            }))
        })
        rescue.met.genes.num <- lapply(rescue.met.genes, length)
    }
    
    ##- set names for plotting clustering
    rownames(distance) <- rescue.met.name
    colnames(distance) <- rownames(distance)
    
    reco.num <- nrow(distance)
    tol      <- 1/max(100, reco.num)^2
    
    if (type == 'i') {
        return (rep(1/reco.num, reco.num))
    }
    if (max(distance) == 0) {
        return (rep(1/reco.num, reco.num))
    }
    
    if (type == 'n') {
        njt <- bionj(distance)                              #- neighbor-joining tree
        ##-- root at the furthest node
        njtroot <- which(
            rownames(distance) == names(
                which.max(apply(distance[rowSums(distance == 0) == 1, , drop=F],
                                1,
                                sum))
            )
        )
        if (length(njtroot) == 1) {
            njt <- root(njt, njtroot, resolve.root=T)
        }
    } else if (type == 'h') {
        njt <- as.phylo(hclust(as.dist(distance), method='average'))
    } else {
        stop("Type unknown: h, i or n required.")
    }
    
    ## pdf(paste0("/tmp/recocluster_", gene.num, "_", draw.num, ".pdf"), width=8, height=12)
    ## plot(njt, cex=0.8)
    ## dev.off()
    
    tmpedge <- njt$edge                                     #- all edges in njt
    tmpedgelength <- njt$edge.length                        #- length of edges
    while (length(unique(tmpedge[, 1])) > 1) {
        idnjt.recos <- which(tmpedge[, 2] <= reco.num)      #- index of leaves
        
        ##-- keep only recos leaves that are not siblings with an internal node
        idnjt.recos.leaves         <- idnjt.recos[which(
            !(tmpedge[idnjt.recos, 1] %in% tmpedge[tmpedge[, 2] > reco.num, 1])
        )]
        weight.recos.leaves        <- tmpedgelength[idnjt.recos.leaves]
        names(weight.recos.leaves) <- njt$tip.label[tmpedge[idnjt.recos.leaves, 2]]
        ids.remove  <- NULL
        edges.add   <- NULL
        lengths.add <- NULL
        for (i in unique(tmpedge[idnjt.recos.leaves, 1])) {
            id.remove    <- which(tmpedge[, 2] == i)
            ids.remove   <- c(ids.remove, id.remove, which(tmpedge[, 1] == i))
            node.parent  <- tmpedge[id.remove, 1]
            weight.added <- tmpedgelength[id.remove]
            weight.sum   <- sum(
                weight.recos.leaves[njt$tip.label[tmpedge[idnjt.recos.leaves,
                                                          ][tmpedge[idnjt.recos.leaves, 1
                                                                    ] == i, 2
                                                            ]
                                                  ]
                                    ]
            )
            
            for (j in tmpedge[idnjt.recos.leaves, ][tmpedge[idnjt.recos.leaves, 1] == i, 2]) {
                if (weight.sum > tol) {
                    weight.recos.leaves[njt$tip.label[j]] <-
                        weight.recos.leaves[njt$tip.label[j]] +
                        weight.added*weight.recos.leaves[njt$tip.label[j]] / weight.sum
                } else {
                    weight.recos.leaves[njt$tip.label[j]] <-
                        weight.recos.leaves[njt$tip.label[j]] +
                        weight.added / length(tmpedge[idnjt.recos.leaves,
                                                      ][tmpedge[idnjt.recos.leaves, 1] == i, 2])
                }
                edges.add   <- rbind(edges.add, c(node.parent, j))
                lengths.add <- rbind(lengths.add, weight.recos.leaves[njt$tip.label[j]])
            }
        }
        tmpedge <- tmpedge[-ids.remove, ]
        tmpedge <- rbind(tmpedge, edges.add)
        tmpedgelength <- tmpedgelength[-ids.remove]
        tmpedgelength <- c(tmpedgelength, lengths.add)
    }
    weight.recos        <- tmpedgelength
    names(weight.recos) <- njt$tip.label[tmpedge[, 2]]
    
    weight.recos.norm   <- weight.recos/sum(weight.recos)   #- normalize weights
    weight.recos.norm   <- weight.recos.norm[rownames(distance)]
    
    pdf(paste0("/tmp/recoclusterweight_", gene.num, "_", draw.num, substr(mod_key(model), 1,5), ".pdf"), width=8, height=12)
    njt$weight <- signif(weight.recos.norm[njt$tip.label],3)
    save(njt, file=paste0("/tmp/njt_", substr(mod_key(model), 1,5), ".RData"))
    njt$tip.label <- paste(njt$tip.label, njt$weight, sep=" # ")
    plot(njt, cex=0.8)
    dev.off()
    
    names(weight.recos.norm) <- recos
    
    return (weight.recos.norm)
}
