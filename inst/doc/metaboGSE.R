## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(fig.width=12, fig.height=4, fig.path='./', fig.align='center', echo = TRUE)
knitr::opts_knit$set(root.dir = normalizePath('..'))

## ---- eval=TRUE, message=FALSE, warning=FALSE---------------------------------
library(metaboGSE)
SYBIL_SETTINGS("SOLVER", "glpkAPI")
SYBIL_SETTINGS("METHOD", "simplex")
SYBIL_SETTINGS("OPT_DIRECTION", "max")

## ---- eval=TRUE, message=FALSE, warning=FALSE---------------------------------
data(iMK735)
iMK735[1]

## ---- eval=TRUE, message=FALSE, warning=FALSE---------------------------------
data(exprMaguire)
names(exprMaguire)
str(exprMaguire$expr, vec.len=3)

## ---- eval=TRUE, message=FALSE, warning=FALSE---------------------------------
data(yarli2GO)
length(yarli2GO)
str(head(yarli2GO))

## ---- eval=TRUE, message=FALSE, warning=FALSE---------------------------------
data(yarliSubmnets)
names(yarliSubmnets)
data(yarliGSE)
length(yarliGSE)

## ---- eval=TRUE---------------------------------------------------------------
target.ratio <- 0.2
hmodel <- iMK735$hypoxia$model
hobj <- iMK735$hypoxia$obj
hmodel.rescue <- rescue(model = hmodel,
                        target = c(target.ratio*hobj),
                        react = c(which(obj_coef(hmodel) == 1)))

## ---- eval=TRUE---------------------------------------------------------------
hmodel.rescue$rescue
head(hmodel.rescue$coef)

## ---- eval=TRUE---------------------------------------------------------------
SYBIL_SETTINGS("TOLERANCE", 1e-08)

## ---- eval=FALSE, message=FALSE, warning=FALSE--------------------------------
#  ## not run
#  mc.cores <- 10
#  fva <- multiDel(model=hmodel.rescue$rescue,
#                  nProc=mc.cores,
#                  todo="fluxVar",
#                  fixObjVal=F,
#                  del1=react_id(hmodel.rescue$rescue))
#  reacs.blo <- names(which(setNames(unlist(lapply(fva, blReact)),
#                     react_id(hmodel.rescue$rescue))))
#  hmodel.clean <- rmReact(hmodel.rescue$rescue, reacs.blo, rm_met=T)
#  ##

## ---- eval=TRUE---------------------------------------------------------------
hmodel.clean <- iMK735$hypoxia$comp
hmodel.clean

## ---- eval=TRUE---------------------------------------------------------------
SYBIL_SETTINGS("OPT_DIRECTION", "min")
hmodel.weight <- changeObjFunc(hmodel.clean, react=rownames(hmodel.rescue$coef), 
                               obj_coef=hmodel.rescue$coef)
hmodel.weight
optimizeProb(hmodel.weight)

## ---- eval=TRUE, message=FALSE------------------------------------------------
mc.cores <- 1
rescue.weight <- (weightReacts(hmodel.weight, mc.cores=mc.cores, gene.num=1))$weight
str(rescue.weight, vec.len=2)

## ---- eval=TRUE, message=FALSE, warning=FALSE---------------------------------
if (!requireNamespace("topGO", quietly = TRUE)) {
    print("Please install topGO: BiocManager::install('topGO')")
} else {
    require(topGO)
    GO2geneID <- inverseList(yarli2GO)
    length(GO2geneID)
    str(head(GO2geneID), vec.len=3)
    gene.name <- names(yarli2GO)
    gene.list <- factor(as.integer(gene.name %in% sybil::allGenes(hmodel.clean)))
    names(gene.list) <- gene.name
    GOdata <- new("topGOdata",
                  ontology = "BP",
                  nodeSize = 5,
                  allGenes = gene.list,
                  annot    = annFUN.gene2GO,
                  gene2GO  = yarli2GO
                  )
    result <- runTest(GOdata, statistic="fisher", algorithm="weight01")
    table  <- GenTable(GOdata,
                       weight   = result,
                       orderBy  = "weight",
                       numChar  = 10000,
                       topNodes = result@geneData[4]
                       )
    table$weight <- as.numeric(sub("<", "", table$weight))
    table <- table[!is.na(table$weight), ]
    MINSIG <- 3
    MAXSIG <- 50
    WCUTOFF <- 0.1
    GO.interest <- table[table$Significant >= MINSIG & table$Significant <= MAXSIG & 
                         table$weight < WCUTOFF, ]$GO.ID
    GO2geneID.interest.proteome <- genesInTerm(GOdata, GO.interest)
    GO2geneID.interest <- lapply(GO2geneID.interest.proteome, function(git) {
        intersect(sybil::allGenes(hmodel.clean), git)
    })
    length(GO.interest)
    str(head(GO2geneID.interest), vec.len=3)
}

## ---- eval=TRUE---------------------------------------------------------------
cond <- "UH"
step <- 50
draw.num <- 4
reps.i <- grep(cond, colnames(exprMaguire$expr), value=T)
ranks <- lapply(reps.i, function(ri) {
    data.frame(
        # ranks1. voom-normalized expression
        expr = exprMaguire$expr[, ri, drop=T],
        # ranks2. pkm normalized expression
        pkmExpr  = exprMaguire$pkmExpr[, ri, drop=T], 
        # ranks3. relative expression power 1
        relExpr1 = relativeExpr(exprMaguire$expr, power=1)[, ri, drop=T], 
        # ranks4. relative expression power 2
        relExpr2 = relativeExpr(exprMaguire$expr, power=2)[, ri, drop=T],  
        # ranks5. relative expression power 3
        relExpr3 = relativeExpr(exprMaguire$expr, power=3)[, ri, drop=T], 
        # ranks6. reverse expression (the worst)
        revExpr  = 1/(1 + exprMaguire$expr[, ri, drop=T]),                 
        # ranks7. z-score expression
        zExpr    = zscoreExpr(exprMaguire$expr)[, ri, drop=T]              
        )
})
names(ranks) <- reps.i
fitnessUH <- fitness(model         = hmodel.weight,
                     ranks         = ranks,
                     rescue.weight = rescue.weight,
                     step          = step,
                     draw.num      = draw.num,
                     mc.cores      = mc.cores)
submnetsUH <- submnet(model        = hmodel.weight,
                      fn           = fitnessUH,
                      rank.best    = "expr",
                      gene.sets    = list("GO:0006696"=
                          c("euk:Q6C8C2_YARLI","euk:ERG6_YARLI","euk:Q6C231_YARLI",
                            "euk:Q6CFB6_YARLI","euk:Q6C6W3_YARLI","euk:Q6C8J1_YARLI",
                            "euk:Q6CB38_YARLI","euk:Q6CEF6_YARLI","euk:ERG27_YARLI",
                            "euk:Q6CGM4_YARLI","euk:FDFT_YARLI","euk:F2Z6C9_YARLI",
                            "euk:Q6C5R8_YARLI","euk:Q6C704_YARLI","euk:Q6CDK2_YARLI",
                            "euk:Q6CFP4_YARLI","euk:Q6BZW0_YARLI","euk:Q6C2X2_YARLI",
                            "euk:Q6C6U3_YARLI")),
                      mc.cores     = mc.cores)

## ---- eval=FALSE--------------------------------------------------------------
#  ## not run
#  pseudo.rank <- base::rank(rowSums(exprMaguire$expr),
#                            ties.method='first')/nrow(exprMaguire$expr)*1e-6
#  exprMaguire$expr <- exprMaguire$expr + pseudo.rank
#  ##

## ---- eval=TRUE, warnings=FALSE-----------------------------------------------
submnetsUH$condition
knitr::kable(submnetsUH$gene.del)
knitr::kable(submnetsUH$fitness.random, digits=3)
knitr::kable(submnetsUH$fitness.ranks$UH1, digits=3)

## ---- eval=TRUE---------------------------------------------------------------
data(yarliSubmnets)
str(yarliSubmnets$UH$gene.del)
dim(yarliSubmnets$UH$fitness.random)
str(yarliSubmnets$UH$fitness.ranks)

## ---- eval=FALSE--------------------------------------------------------------
#  ## not run
#  simulateSubmnet(sgd=submnetsUH)
#  ##

## ---- eval=FALSE, message=FALSE-----------------------------------------------
#  ## not run
#  GSE <- metaboGSE(yarliSubmnets, method="perm", nperm=1000, nrand=1000,
#                   mc.cores=mc.cores, prefix="/tmp/summary")
#  ##

## ---- eval=TRUE, message=FALSE------------------------------------------------
data(yarliGSE)
GSE <- yarliGSE
str(GSE[["GO:0006696"]], vec.len=2)
GSE[["GO:0006696"]]$res$p.Val

## ---- eval=TRUE, message=FALSE------------------------------------------------
GS.sig.all <- as.data.frame(t(sapply(GSE, function(gsm) {
    c(GS.ID=gsm$res$GS.ID, 
      Description=gsm$res$Description,
      Statistic=gsm$res$statistic,
      p.Cond=if (is.null(gsm$res$p.Cond)) NA else min(gsm$res$p.Cond),
      p.Val=gsm$res$p.Val)
  })), stringsAsFactors=F)
GS.sig.all$FDR  <- p.adjust(as.numeric(GS.sig.all$p.Val), method="BH")
GS.sig.all <- GS.sig.all[!is.na(GS.sig.all$FDR), ]
dim(GS.sig.all)

## ---- eval=TRUE---------------------------------------------------------------
GS.sig <- GS.sig.all[as.numeric(GS.sig.all$FDR) < 0.05, , drop=F]
head(GS.sig[order(as.numeric(GS.sig$Statistic), decreasing=T), ])

