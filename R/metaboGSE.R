#' Description of gene set IDs
#'
#' This function returns the description of given gene set IDs, e.g. GO, KEGG.
#' @param x A vector of gene set IDs.
#' @param desc.data A named vector of descriptions of all studied gene set IDs. Default: NULL, 
#' AnnotationDbi is used if x is a GO term ID. KEGGREST will be called with internet connection required if x is a KEGG pathway ID.
#' @return Description
#' @importFrom AnnotationDbi Term
#' @examples
#' pwDesc("GO:0006696")
#' pwDesc("genesetX", desc.data=setNames("processX", "genesetX"))
#' \donttest{
#' pwDesc("hsa04930")
#' }
#' @export
pwDesc <- function(x, desc.data = NULL) {
    if (!is.null(desc.data)) {
        return (sub(desc.data[x], 
                    pattern=" - Cricetulus griseus (Chinese hamster)",
                    replacement="", fixed=T))
    }
    if (length(grep("^GO:", x)) == length(x)) {
        return (Term(x))
    }
    if (length(grep("^GO:", x)) > 0) {
        warning("Wrong format of GO ID found.")
        return (rep("", length(x)))
    }
    if (!requireNamespace("KEGGREST", quietly = TRUE)) {
        stop("Please install KEGGREST: source('https://bioconductor.org/biocLite.R'); biocLite('KEGGREST')")
    }
    return (sapply(x, function(xx) {
        kegg.path <- tryCatch(KEGGREST::keggList(xx), error=function(e) {""})
        names(kegg.path) <- NULL
        return (sub(kegg.path, 
                    pattern=" - Cricetulus griseus (Chinese hamster)",
                    replacement="", fixed=T))
    }, USE.NAMES=T))
}


#' Abbreviate GO term description
#'
#' This function produces abbreviations in GO term description.
#' @param x A string describing a GO term
#' @return The abbreviation of x
#' @keywords internal
abbr <- function(x) {
    return (
        gsub(
            gsub(
                gsub(
                    gsub(
                        gsub(
                            gsub(
                                gsub(
                                    gsub(x,
                                         pattern="metabolic process", replacement="mp"),
                                    pattern="biosynthetic process", replacement="bp"),
                                pattern="catabolic process", replacement="cp"),
                            pattern="transmembrane transport", replacement="tt"),
                        pattern="amino acid", replacement="aa"),
                    pattern="response", replacement="res."),
                pattern="signaling pathway", replacement="sp"),
            pattern="regulation", replacement="reg.")
        )
}


#' Compute area between two curves
#' 
#' This function computes the area formed by two curves.
#' @param y1 A vector of y-coordinates for the first curve.
#' @param y2 A vector of y-coordinates for the second curve.
#' @param x A vector of increasing x-coordinates. Default, an equidistant vector of length length(y1) from 0 to 1.
#' @return The area
#' @keywords internal
area <- function(y1, y2, x = NULL) {
    if (length(y1) != length(y2)) {
        stop("Curves of different lengths.")
    }
    if (is.null(x)) {
        x <- seq(0, 1, length.out = length(y1))
    } else if (length(y1) != length(x)) {
        stop("Length of x does not match to curves.")
    }
    xdiff <- diff(x)
    diff.ij <- y1 - y2
    area.int <- sum(sapply(1L:(length(x)-1L), function(k) {
        if (sign(diff.ij[k]) == sign(diff.ij[k+1]) || sign(diff.ij[k])*sign(diff.ij[k+1])==0) {
            return (abs(diff.ij[k]+diff.ij[k+1])*xdiff[k]/2)
        } else {
            return ((diff.ij[k]^2+diff.ij[k+1]^2)*xdiff[k]/(2*(abs(diff.ij[k])+abs(diff.ij[k+1]))))
        }
    }))
    return (area.int)
}


#' Compute maxium area between curves
#'
#' This function computes the maximum area between multiple curves.
#' @param y A numeric matrix for all curves with each curve per column.
#' @param x A vector of increasing x-coordinates. Default, an equidistant vector of length \code{nrow(y)} from 0 to 1.
#' @return The maximum area
#' @keywords internal
maxArea <- function(y, x = NULL) {
    if (is.null(x)) {
        x <- seq(0, 1, length.out = nrow(y))
    }
    stopifnot(length(x)==nrow(y))
    areas <- sapply(1L:(ncol(y)-1L), function(i) {
        sapply((i+1L):ncol(y), function(j) {
            area(y[, j], y[, i], x)
        })
    })
    return (max(unlist(areas)))
}


#' Gene set enrichment analysis
#'
#' This function performs the gene set enrichment analysis.
#' @param scores A list of \code{scoreGeneDel} objects.
# #' @param expr A matrix of expression data (log2 counts), with gene per row and experiment per column.
#' @param gene.sets A named list of gene sets for gene set enrichment analysis, or a vector of gene set IDs computed in \code{scores}. Default: NULL, all gene sets from \code{scores}.
#' @param method Statistical testing method \code{c("perm", "survival")}. Default: \code{"perm"}. \code{"survial"} may be used for exploration.
#' @param test Type of test c(\code{"likelihood", "logrank", "wald"}), when method = \code{"survival"}. Default: \code{"likelihood"}.
#' @param nperm Number of permutations for testing, when method = \code{"perm"}. Default: 10000.
#' @param nrand Number of draws for random gene set generation, when method = \code{"perm"}. Default: 1000.
#' @param mc.cores The number of cores to use (at least 1), i.e. at most how many child processes will be run simultaneously. Default: 1.
#' @param posthoc A logical value indicating if pairwise tests are performed. Default: TRUE.
#' @param contrast A logical value indicating if the Newick-based contrast will be computed. Default: FALSE.
#' @param prefix A string indicating prefix of output plots. Default: NA, no plot.
#' @param desc.data A vector of descriptions of a priori KEGG pathway IDs. Default: NULL,
#' KEGGREST will be called with internet connection required if gene.sets is KEGG pathway.
#' @param cols Colors for conditions. Default: rainbow colors.
#' @param ltys Line types for conditions. Default: incrementing line types in R.
#' @return Gene set enrichment information
#' @import parallel utils grDevices graphics
#' @examples
#' data(yarliSubmnets)
#' \donttest{
#' metaboGSE(yarliSubmnets[c('SH','SN')], gene.sets = "GO:0006696",
#'          method="perm", nperm=10, nrand=10)
#' }
#' @export
metaboGSE <- function(scores, gene.sets = NULL, method = "perm", test = NA,
                      nperm = 1000, nrand = 1000, mc.cores = 1,
                      posthoc = TRUE, contrast = FALSE, prefix = NA, desc.data = NULL,
                      cols = NULL, ltys = NULL) {
    ##- settings ----
    ties         <- mean
    conds        <- unlist(lapply(scores, function(sgd) sgd$condition), use.names=F)
    cond.num     <- length(conds)
    samples      <- mclapply(scores, mc.cores=mc.cores, function(sgd) {names(sgd$ratio.GS)})
    sample.num   <- sapply(samples, length)
    samples.num  <- sum(sample.num)
    if (method == "perm") {
        nperm    <- as.integer(nperm)
        nrand    <- as.integer(nrand)
        ncombi   <- factorial(samples.num)/(prod(sapply(sample.num, factorial))*prod(sapply(table(sample.num), factorial)))
        ncombi2  <- ncombi*factorial(samples.num)/prod(sapply(sample.num, factorial)) - factorial(cond.num) + 1
    }
    if (length(cols) == 0) {
        cols        <- c(rainbow(cond.num))
        if (length(cols) == 7L) {
            cols[2] <- "darkorange1"
            cols[3] <- "gold3"
        }
    }
    if (length(cols) != cond.num) {
        stop("Number of colors does not match number of conditions.")
    }
    cols.rep     <- rep(cols, times=sapply(samples, length))
    if (length(ltys) == 0) {
        ltys     <- 1L:cond.num
    }
    if (length(ltys) != cond.num) {
        stop("Number of line types does not match number of conditions.")
    }
    ltys.rep     <- rep(ltys, times=sapply(samples, length))
    cex.lab      <- 2.2
    cex.leg      <- 2.5
    cex.axis     <- 2.6
    cex.panel    <- 2.8
    line.panel   <- -1
    line.lab     <- 4
    at.panel     <- -100
    legend.pos   <- "topright"
    func         <- mean
    gene.del     <- as.vector(scores[[1L]]$gene.del)
    pw.posthoc   <- NULL
    gs.newick.ec <- ""
    group1.ec    <- ""
    group2.ec    <- ""
    if (method == "survival") {
        if (!requireNamespace("survival", quietly = TRUE)) {
            stop("Please install survival: install.packages('survival')")
        }
        if (is.na(test)) {
            test <- "likelihood"
        }
    }
#    if (contrast) {
#        if (!requireNamespace("ctc", quietly = TRUE)) {
#            stop("Please install ctc: source('https://bioconductor.org/biocLite.R'); biocLite('ctc')")
#        }
#    }
    RNGkind("L'Ecuyer-CMRG")
    set.seed(1000)
    ##-----
    
    ##- gene set enrichment ----
    if (is.null(gene.sets)) {
        gene.sets <- scores[[1L]]$gene.sets
    } else if (is.character(gene.sets)) {
        gene.sets <- scores[[1L]]$gene.sets[gene.sets]
    }
    if (1 %in% lengths(gene.sets)) {
        warning("gene set of length 1 will be ignored!")
    }
    gene.sets <- gene.sets[lengths(gene.sets) > 1]
    fitness.ranked <- t(do.call(rbind, mclapply(scores, mc.cores=mc.cores, function(sgd) {
        do.call(rbind, lapply(sgd$fitness.ranks, function(fn) {
            fn[1,]
        }))
    })))
    fitness.ranked.pcg <- fitness.ranked * (1-gene.del/max(gene.del))
    GS.ulen <- c(3:50) #sort(unique(sapply(c(scores[[1L]]$gene.sets, gene.sets), length)))
    
    if (method == "perm") {
        ##- AUCs created by random gene set of size GS.ulen
        rGS.auc <- mclapply(GS.ulen, mc.cores=mc.cores, mc.set.seed=T, function(gsul) {
            rgss <- replicate(nrand, sample(gene.del, size=gsul, replace=(gsul>length(gene.del))))
            rgss.auc <- apply(rgss, 2, function(rgs) {
                gs.frac <- 1 - cumsum(gene.del %in% rgs)/gsul
                gs.frac.u <- unique(gs.frac)
                idx.out <- as.vector(sapply(gs.frac.u, function(gfu) {
                    gfu.id <- which(gs.frac==gfu)
                    return (c(min(gfu.id), max(gfu.id)))
                }))
                apply(fitness.ranked.pcg[idx.out,], 2, function(x) {
                    area(y1=rep(0,length(idx.out)), y2=rev(gs.frac[idx.out]), x=rev(x))
                })
            })
            return (rgss.auc)
        })
        names(rGS.auc) <- GS.ulen
        rGS.auc.mean <- mclapply(rGS.auc, mc.cores=mc.cores, function(rgs.auc) {
            rgs.auc.mean <- t(sapply(1L:cond.num, function(k) {
                apply(rgs.auc[samples[[k]], , drop=F], 2, mean)
            }))
            rownames(rgs.auc.mean) <- conds
            return (rgs.auc.mean)
        })
        
        mc.cores1 <- floor(sqrt(mc.cores))
        mc.cores2 <- floor(mc.cores/mc.cores1)
        GS.metric <- mclapply(1L:length(gene.sets), mc.cores=mc.cores1, function(j) {
            gs <- names(gene.sets)[j]
            i  <- match(gs, names(scores[[1L]]$gene.sets))
            ##- degradation of gs w.r.t number of removed genes
            gs.fracs <- do.call(cbind, mclapply(scores, mc.cores=mc.cores2, function(sgd) {
                if (is.na(i)) {
                    gs.frac <- t(sapply(sgd$sub.genes, function(sgdsg) {
                        sapply(sgdsg, function(sg) {
                            length(intersect(gene.sets[[j]], sg)) / length(gene.sets[[j]])
                        }, USE.NAMES=T)
                    }))
                } else {
                    gs.frac <- sapply(sgd$ratio.GS, function(rgs.rep) {
                        rgs.rep[gs, , drop=F]
                    }, USE.NAMES=T)
                }
                return (gs.frac)
            }))
            initial.ratio <- max(gs.fracs[1L,])
            gs.fracs.mean <- do.call(cbind, mclapply(samples, mc.cores=mc.cores2, function(csample) {
                apply(gs.fracs[, csample, drop=F], 1, func)
            }))
            ##- positions where gs.fracs change their value
            idx.out <- mclapply(1L:samples.num, mc.cores=mc.cores2, function(k) {
                gs.frac <- gs.fracs[, k]
                gs.frac.u <- unique(gs.frac)
                as.vector(sapply(gs.frac.u, function(gfu) {
                    gfu.id <- which(gs.frac==gfu)
                    return (c(min(gfu.id), max(gfu.id)))
                }))
            })
            ##- fitness where gs.fracs change their value
            xout <- unique(sort(unlist(mclapply(1L:samples.num, mc.cores=mc.cores2, function(k) {
                fitness.ranked.pcg[idx.out[[k]], k]
            }))))
            xlen <- length(xout)
            ##- interpolation for AUC computation
            interps <- mclapply(1L:samples.num, mc.cores=mc.cores2, function(k) {
                localfit <- approx(fitness.ranked.pcg[, k, drop=T],
                                   gs.fracs[, k, drop=T],
                                   rule=2,
                                   xout=xout,
                                   ties=ties)
                return (localfit)
            })
            names(interps) <- colnames(gs.fracs)
            gs.fracs.itp <- do.call(cbind, mclapply(interps, mc.cores=mc.cores2, function(itp) {itp$y}))
            colnames(gs.fracs.itp) <- colnames(gs.fracs)
            # pdf("curverep.pdf")
            # par(mar=c(5,5,1,1))
            # matplot(xout, gs.fracs.itp, xlim=c(1,0), xlab='fitness*frac_gene_model', ylab="frac_gene_set",
            #         type='l', col=cols.rep, lty=ltys.rep, lwd=2, cex.lab=2, cex.axis=1.5)
            # dev.off()
            gs.fracs.itp.mean <- do.call(cbind, mclapply(samples, mc.cores=mc.cores2, function(csample) {
                apply(gs.fracs.itp[, csample, drop=F], 1, func)
            }))
            # pdf("curvemean.pdf")
            # par(mar=c(5,5,1,1))
            # matplot(xout, gs.fracs.itp.mean, xlim=c(1,0), xlab='fitness*frac_gene_model', ylab="frac_gene_set",
            #         type='l', col=cols, lty=ltys, lwd=2, cex.lab=2, cex.axis=1.5)
            # dev.off()
            ##- AUC for each sample
            gs.fracs.itp.auc <- setNames(unlist(mclapply(1L:samples.num, mc.cores=mc.cores2, function(k) {
                area(y1=rep(0, xlen), y2=gs.fracs.itp[,k], x=xout)
            })), colnames(gs.fracs.itp))
            ##- AUC for each condition
            gs.fracs.itp.auc.mean <- sapply(samples, function(csample) {
                func(gs.fracs.itp.auc[csample])
            })
            ##- compare to random AUC
            rgs.auc.mean <- initial.ratio*do.call(cbind, rGS.auc.mean)
            
            p.Cond <- setNames(unlist(mclapply(conds, mc.cores=mc.cores2, function(cs) {
                (length(which(rgs.auc.mean[cs,] <= gs.fracs.itp.auc.mean[cs]))+1L)/(length(rgs.auc.mean[cs,])+1L)
            })), conds)
            # pdf('auccond.pdf')
            # par(mar=c(5,5,4,1))
            # idx.auc.min <- which.min(p.Cond)
            # idx.auc.max <- which.max(p.Cond)
            # hist(rgs.auc.mean[idx.auc.min,],
            #      main="",#paste0(conds[idx.auc.min], ": AUC(", gs, ") vs random AUCs"),
            #      xlab="AUC", cex.main=2, cex.lab=1.8, cex.axis=1.8)
            # abline(v=gs.fracs.itp.auc.mean[idx.auc.min], col=cols[idx.auc.min], lty=ltys[idx.auc.min], lwd=4)
            # #abline(v=gs.fracs.itp.auc.mean[idx.auc.max], col=cols[idx.auc.max], lty=ltys[idx.auc.max], lwd=4)
            # # sapply(1L:length(p.Cond), function(ip) {
            # #     abline(v=gs.fracs.itp.auc.mean[ip], col=cols[ip], lty=ltys[ip], lwd=4)
            # # })
            # mtext("D", side=3, line=1, outer=F, cex=3, at=-0.07, font=2)
            # dev.off()
            
            if (cond.num == 1L) {
                statistic <- mean(gs.fracs.itp.auc)
                p.Val     <- NA
            } else {
                statistic <- maxArea(y=gs.fracs.itp.mean, x=xout)
                p.Val <- NA
                if (ncombi2 > 1) {
                    rstatistic <- unlist(mclapply(1L:nperm, mc.cores=mc.cores2, mc.set.seed=T, function(k) {
                        ridx <- sample(1L:xlen, size=1L)
                        gs.fracs.itp.perm <- gs.fracs.itp
                        gs.fracs.itp.perm[c(1L:ridx), ] <- gs.fracs.itp[c(1L:ridx), sample(1L:samples.num)]
                        gs.fracs.itp.perm[-c(1L:ridx), ] <- gs.fracs.itp[-c(1L:ridx), sample(1L:samples.num)]
                        gs.fracs.itp.perm.mean <- sapply(samples, function(csample) {
                            apply(gs.fracs.itp.perm[, csample, drop=F], 1, func)
                        })
                        return (maxArea(y=gs.fracs.itp.perm.mean, x=xout))
                    }))
                    p.Val <- (length(which(rstatistic >= statistic))+1L)/(length(rstatistic)+1L)
                }
                
                ##- pairwise tests ----
                if (posthoc) {
                    pw.idx <- combn(cond.num, 2)
                    pw.posthoc <- as.data.frame(t(sapply(1L:ncol(pw.idx), function(k) {
                        pw.statistic <- maxArea(y=gs.fracs.itp.mean[, pw.idx[,k]], x=xout)
                        pw.aucdiff <- area(y1=rep(0, xlen), y2=gs.fracs.itp.mean[, pw.idx[1,k]], x=xout) -
                            area(y1=rep(0, xlen), y2=gs.fracs.itp.mean[, pw.idx[2,k]], x=xout)
                        pw.p.Val <- NA
                        pw.sample.num <- sapply(samples[pw.idx[,k]], length)
                        pw.samples.num <- sum(pw.sample.num)
                        pw.ncombi <- factorial(pw.samples.num)/
                            (prod(sapply(pw.sample.num, factorial))*prod(sapply(table(pw.sample.num), factorial)))
                        pw.ncombi2 <- pw.ncombi*factorial(pw.samples.num)/prod(sapply(pw.sample.num, factorial)) - 1
                        if (pw.ncombi2 > 1) {
                            pw.gfi <- gs.fracs.itp[, unlist(samples[pw.idx[,k]]), drop=F]
                            pw.rstatistic <- unlist(mclapply(1L:nperm, mc.cores=mc.cores2, mc.set.seed=T, function(l) {
                                ridx <- sample(1L:xlen, size=1)
                                pw.gfi.perm <- pw.gfi
                                pw.gfi.perm[c(1L:ridx), ] <- pw.gfi[c(1L:ridx), sample(1L:ncol(pw.gfi))]
                                pw.gfi.perm[-c(1L:ridx), ] <- pw.gfi[-c(1L:ridx), sample(1L:ncol(pw.gfi))]
                                
                                pw.gfi.perm.mean <- sapply(samples[pw.idx[,k]], function(csample) {
                                    apply(pw.gfi.perm[, csample, drop=F], 1, func)
                                })
                                return (maxArea(y=pw.gfi.perm.mean, x=xout))
                            }))
                            pw.p.Val <- (length(which(pw.rstatistic >= pw.statistic))+1L)/(length(pw.rstatistic)+1L)
                        }
                        return (c(cond1     = conds[pw.idx[1,k]],
                                  cond2     = conds[pw.idx[2,k]],
                                  statistic = signif(pw.statistic, 2),
                                  aucdiff   = signif(pw.aucdiff, 2),
                                  p.Val     = pw.p.Val))
                    })), stringsAsFactors=F)
                }
                ##-----
                
                ##- Newick ----
                #if (contrast) {
                #    gs.hc.ec <- hclust(dist(t(gs.fracs.itp.mean), method="manhattan"), method="ward.D")
                #    gs.newick.ec <- gsub(ctc::hc2Newick(gs.hc.ec), pattern=":[.0-9]+|;$", replacement="", perl=T)
                #    groupcut.ec <- sort(cutree(gs.hc.ec, k=2))
                #    group1.ec <- names(which(groupcut.ec==1))
                #    group2.ec <- names(which(groupcut.ec==2))
                #}
                ##-----
            }
            
            
            return (list(res=list(GS.ID       = gs,
                                  Description = unname(pwDesc(gs, desc.data=desc.data)),
                                  statistic   = statistic,
                                  p.Cond      = p.Cond,
                                  p.Val       = p.Val,
                                  auc         = gs.fracs.itp.auc,
                                  pw.posthoc  = pw.posthoc,
                                  contrast    = c(newick = gs.newick.ec,
                                                  group1 = paste(group1.ec, collapse=','),
                                                  group2 = paste(group2.ec, collapse=',')
                                      )),
                         plot=list(
                             gs.fracs=gs.fracs,
                             xout=xout, 
                             gs.fracs.itp.mean=gs.fracs.itp.mean))
            )
        })
        names(GS.metric) <- names(gene.sets)
    }
    
    if (method == "survival") {
        mc.cores1 <- floor(sqrt(mc.cores))
        mc.cores2 <- floor(mc.cores/mc.cores1)
        GS.metric <- mclapply(1L:length(gene.sets), mc.cores=mc.cores1, function(j) {
            gs <- names(gene.sets)[j]
            i  <- match(gs, names(scores[[1L]]$gene.sets))
            gs.len <- length(scores[[1L]]$gene.sets[[i]])
            ##- degradation of gs w.r.t number of removed genes
            gs.fracs <- do.call(cbind, mclapply(scores, mc.cores=mc.cores2, function(sgd) {
                if (is.na(i)) {
                    gs.frac <- t(sapply(sgd$sub.genes, function(sgdsg) {
                        sapply(sgdsg, function(sg) {
                            length(intersect(gene.sets[[j]], sg)) / length(gene.sets[[j]])
                        }, USE.NAMES=T)
                    }))
                } else {
                    gs.frac <- sapply(sgd$ratio.GS, function(rgs.rep) {
                        rgs.rep[gs, , drop=F]
                    }, USE.NAMES=T)
                }
                return (gs.frac)
            }))
            initial.ratio <- max(gs.fracs[1L,])
            gs.fracs.mean <- do.call(cbind, mclapply(samples, mc.cores=mc.cores2, function(csample) {
                apply(gs.fracs[, csample, drop=F], 1, func)
            }))
            ##- positions where gs.fracs change their value
            idx.out <- mclapply(1L:samples.num, mc.cores=mc.cores2, function(k) {
                gs.frac <- gs.fracs[, k]
                gs.frac.u <- unique(gs.frac)
                as.vector(sapply(gs.frac.u, function(gfu) {
                    gfu.id <- which(gs.frac==gfu)
                    return (min(gfu.id))
                }))
            })
            
            surv.data <- as.data.frame(do.call(rbind, mclapply(1L:samples.num, mc.cores=mc.cores2, function(k) {
                cond.idx <- sapply(conds, function(cnd) {grep(cnd, colnames(gs.fracs)[k], value=T)}, USE.NAMES=T)
                cond <- names(cond.idx)[which(lengths(cond.idx)>0)]
                do.call(rbind, lapply(2L:length(idx.out[[k]]), function(l) {
                    idx.l <- idx.out[[k]][l]
                    idv.num <- round(gs.len*(gs.fracs[idx.l-1L,k]-gs.fracs[idx.l,k]), 0L)
                    if (idv.num <= 0) {
                        return (NULL)
                    }
                    return (t(replicate(idv.num,
                                        c(time   = unname((1-fitness.ranked.pcg[idx.l-1L,k])^fitness.ranked.pcg[idx.l-1L,k]),
                                          weight = unname(fitness.ranked.pcg[idx.l-1L,k]),
                                          cond   = cond))))
                }))
            })), stringsAsFactors=F)
            for (ele in c("time", "weight")) {
                surv.data[, ele] <- as.numeric(surv.data[, ele])
            }
            
            if (cond.num == 1) {
                statistic <- NA
                p.Val <- NA
            } else {
                #plot(survfit(Surv(time) ~ cond, data=surv.data), mark.time=T, col=cols, lty=ltys)
                gs.coxph <- survival::coxph(survival::Surv(time) ~ cond,
                                            data=surv.data,
                                            weights=surv.data$weight, 
                                            control=survival::coxph.control(iter.max = 40))
                gs.coxph.sum <- summary(gs.coxph)
                test.idx <-ifelse(test=="likelihood", 9, ifelse(test=="logrank", 10, 12))
                statistic <- unname(gs.coxph.sum[[test.idx]]["test"])
                p.Val <- unname(gs.coxph.sum[[test.idx]]["pvalue"])
                
                ##- pairwise test ----
                if (posthoc) {
                    pw.idx <- combn(cond.num, 2)
                    pw.posthoc <- as.data.frame(t(sapply(1L:ncol(pw.idx), function(k) {
                        pw.surv.data <- surv.data[surv.data$cond %in% conds[pw.idx[,k]], , drop=F]
                        pw.gs.coxph  <- survival::coxph(survival::Surv(time) ~ cond,
                                                        data=pw.surv.data,
                                                        weights=pw.surv.data$weight, 
                                                        control=survival::coxph.control(iter.max = 40))
                        pw.gs.coxph.sum <- summary(pw.gs.coxph)
                        pw.statistic <- unname(pw.gs.coxph.sum[[test.idx]]["test"])
                        pw.p.Val <- unname(pw.gs.coxph.sum[[test.idx]]["pvalue"])
                        
                        return (c(cond1     = conds[pw.idx[1,k]],
                                  cond2     = conds[pw.idx[2,k]],
                                  statistic = signif(pw.statistic, 2),
                                  p.Val     = pw.p.Val))
                    })), stringsAsFactors=F)
                }
                ##-----
                
                ##- Newick ----
                if (contrast) {
                    sfit <- summary(survival::survfit(survival::Surv(time) ~ cond, data=surv.data))
                    sfit.data <- data.frame(time=sfit$time, 
                                            n.risk=sfit$n.risk, 
                                            n.event=sfit$n.event, 
                                            cond=gsub(as.vector(sfit$strata),
                                                      pattern='cond=', replacement=''))
                    conds.dist <- dist(setNames(mclapply(conds, mc.cores=mc.cores2, function(cond) {
                        sfit.cond <- sfit.data[sfit.data$cond==cond, , drop=T]
                        sfit.cond.x <- c(0, sfit.cond$time)
                        sfit.cond.y <-  sfit.cond$n.risk/max(sfit.cond$n.risk)
                        return (sum(sfit.cond.y*diff(sfit.cond.x)))
                    }), conds))
                    
                    gs.hc.ec <- hclust(conds.dist, method='ward.D')
                    #gs.newick.ec <- gsub(ctc::hc2Newick(gs.hc.ec),
                    #                     pattern=":[.0-9]+|e-[0-9]+|;$", replacement="", perl=T)
                    groupcut.ec <- sort(cutree(gs.hc.ec, k=2))
                    group1.ec <- names(which(groupcut.ec==1))
                    group2.ec <- names(which(groupcut.ec==2))
                }
                ##-----
            }
            
            return (list(res=list(GS.ID       = gs,
                                  Description = unname(pwDesc(gs, desc.data=desc.data)),
                                  statistic   = statistic,
                                  p.Val       = p.Val,
                                  pw.posthoc  = pw.posthoc,
                                  contrast    = c(newick = gs.newick.ec,
                                                  group1 = paste(group1.ec, collapse=','),
                                                  group2 = paste(group2.ec, collapse=',')
                                  )),
                         plot=list(
                             gs.fracs=gs.fracs,
                             surv.data=surv.data))
            )
        })
        names(GS.metric) <- names(gene.sets)
    }

    if (!is.na(prefix)) {
        ##- plot ----
        cex.lab    <- 1.2
        cex.main   <- 1.3
        cex.leg    <- 0.9
        cex.axis   <- 0.9
        cex.panel  <- 1.6
        line.panel <- -1
        panels     <- c("A", "B")
        
        pdf(file=paste(prefix, "GSE.pdf", sep='_'), width=10, height=5)
        invisible(sapply(GS.metric, function(gsm) {
            par(mfrow=c(1,2), mar=c(4.2,5.0,0.3,0.5), oma=c(0.1,0.5,5,0.5), xpd=NA)
            matplot(gene.del, gsm$plot$gs.fracs,
                    xlab="#removed_genes",
                    ylab="frac_gene_set",
                    ylim=c(0,1),
                    col=cols.rep, lty=ltys.rep,
                    type='l', lwd=2)
            mtext(paste0(gsm$res$GS.ID, ", ", abbr(pwDesc(gsm$res$GS.ID, desc.data=desc.data))),
                  side=3, line=3, outer=T, cex=cex.main, font=2)
            if (method == 'perm') {
                mtext(paste0("E = ", signif(gsm$res$statistic,2), 
                                        #", min(p) = ", signif(min(gsm$res$p.Cond),2),
                                        #", max(p) = ", signif(max(gsm$res$p.Cond),2),
                             ", p-val = ", signif(gsm$res$p.Val,2)),
                      side=3, line=1, outer=T, cex=cex.main, font=2)
            } else {
                mtext(paste0("E = ", signif(gsm$res$statistic,2), 
                             ", p-val = ", signif(gsm$res$p.Val,2)),
                      side=3, line=1, outer=T, cex=cex.main, font=2)
            }
            mtext(panels[1], side=3, line=line.panel, outer=T, adj=0, cex=cex.panel, font=2)
            legend(ifelse(max(gsm$plot$gs.fracs[min(which(gene.del>0.75*max(gene.del))),]) < 0.5, 
                          'topright', 'bottomleft'),
                   conds, lty=ltys, lwd=2, col=cols, cex=cex.leg)
            if (method == 'perm') {
                matplot(gsm$plot$xout, gsm$plot$gs.fracs.itp.mean,
                        xlab="fitness*frac_gene_model", ylab="frac_gene_set",
                        xlim=c(1,0), ylim=c(0,1),
                        col=cols, lty=ltys,
                        type='l', lwd=2)
                mtext(panels[2], side=3, line=line.panel+0.4, outer=F, at=1.29, cex=cex.panel, font=2)
            } else {
                plot(survival::survfit(survival::Surv(time) ~ cond, data=gsm$plot$surv.data), mark.time=T, 
                     xlim=c(0,1),
                     xlab="Time (1-fitness*frac_gene_model)", ylab="Survival fraction",
                     col=cols, lty=ltys, lwd=2)
                mtext(panels[2], side=3, line=line.panel+0.4, outer=F, at=-0.23, cex=cex.panel, font=2)
            }
        }))
        invisible(dev.off())
        ##-----
    }
    return (GS.metric)
}
