## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(fig.width=12, fig.height=4, fig.path='./', fig.align='center', echo = TRUE)
knitr::opts_knit$set(root.dir = normalizePath('..'))

## ---- eval=TRUE, message=FALSE, warning=FALSE----------------------------
library(metaboGSE)
SYBIL_SETTINGS("SOLVER", "clpAPI")
SYBIL_SETTINGS("METHOD", "inibarrier")
SYBIL_SETTINGS("OPT_DIRECTION", "max")

## ---- eval=TRUE, message=FALSE, warning=FALSE----------------------------
data(iMK735)
iMK735[1]

## ---- eval=TRUE, message=FALSE, warning=FALSE----------------------------
data(exprMaguire)
names(exprMaguire)
str(exprMaguire$expr, vec.len=3)

## ---- eval=TRUE, message=FALSE, warning=FALSE----------------------------
data(yarli2GO)
length(yarli2GO)
str(head(yarli2GO))

## ---- eval=TRUE, message=FALSE, warning=FALSE----------------------------
data(yarliSubmnets)
names(yarliSubmnets)
data(yarliGSE)
length(yarliGSE)

## ---- eval=TRUE----------------------------------------------------------
target.ratio <- 0.2
hmodel <- iMK735$hypoxia$model
hobj <- iMK735$hypoxia$obj
hmodel.rescue <- rescue(model = hmodel,
                        target = c(target.ratio*hobj),
                        react = c(which(obj_coef(hmodel) == 1)))

## ---- eval=TRUE----------------------------------------------------------
hmodel.rescue$rescue
head(hmodel.rescue$coef)

## ---- eval=TRUE----------------------------------------------------------
SYBIL_SETTINGS("TOLERANCE", 1e-08)

## ---- eval=FALSE, message=FALSE, warning=FALSE---------------------------
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

## ---- eval=TRUE----------------------------------------------------------
hmodel.clean <- iMK735$hypoxia$comp
hmodel.clean

## ---- eval=TRUE----------------------------------------------------------
SYBIL_SETTINGS("OPT_DIRECTION", "min")
hmodel.weight <- changeObjFunc(hmodel.clean, react=rownames(hmodel.rescue$coef), obj_coef=hmodel.rescue$coef)
hmodel.weight
optimizeProb(hmodel.weight)

## ---- eval=TRUE, message=FALSE-------------------------------------------
mc.cores <- 1
rescue.weight <- weightReacts(hmodel.weight, mc.cores=mc.cores, gene.num=1)
str(rescue.weight, vec.len=2)

## ---- eval=TRUE, message=FALSE, warning=FALSE----------------------------
library(topGO)
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

## ---- eval=TRUE----------------------------------------------------------
cond <- "UH"
step <- 50
draw.num <- 4
reps.i <- grep(cond, colnames(exprMaguire$expr))
ranks <- mclapply(reps.i, mc.cores=mc.cores, function(ri) {
    data.frame(
         # ranks1. pkm normalized expression
         pkmExpr  = exprMaguire$pkmExpr[, ri, drop=T], 
         # ranks2. relative expression power 1
         relExpr1 = relativeExpr(exprMaguire$expr, power=1)[, ri, drop=T], 
         # ranks3. relative expression power 2
         relExpr2 = relativeExpr(exprMaguire$expr, power=2)[, ri, drop=T],  
         # ranks4. relative expression power 3
         relExpr3 = relativeExpr(exprMaguire$expr, power=3)[, ri, drop=T], 
         # ranks5. reverse expression (the worst)
         revExpr  = 1/(1 + exprMaguire$expr[, ri, drop=T]),                 
         # ranks6. z-score expression
         zExpr    = zscoreExpr(exprMaguire$expr)[, ri, drop=T]              
	 )
    })
submnetsUH <- submnet(model         = hmodel.weight,
                      expr          = exprMaguire$expr[, reps.i, drop=F],
                      ranks         = ranks,
                      rescue.weight = rescue.weight,
                      step          = step,
                      draw.num      = draw.num,
                      gene.sets     = GO2geneID.interest,
                      mc.cores      = mc.cores)

## ---- eval=TRUE, warnings=FALSE------------------------------------------
submnetsUH$condition
knitr::kable(submnetsUH$gene.del)
knitr::kable(submnetsUH$fitness.random, digits=3)
knitr::kable(submnetsUH$fitness.ranked, digits=3)
knitr::kable(submnetsUH$fitness.ranks$UH1, digits=3)

## ---- eval=TRUE----------------------------------------------------------
data(yarliSubmnets)
str(yarliSubmnets$UH$gene.del)
dim(yarliSubmnets$UH$fitness.random)
str(yarliSubmnets$UH$fitness.ranked)
str(yarliSubmnets$UH$fitness.ranks)

## ---- eval=FALSE---------------------------------------------------------
#  simulateSubmnet(model    = hmodel.weight,
#                  sgd      = submnetsUH,
#                  mc.cores = mc.cores)

## ---- eval=FALSE, message=FALSE------------------------------------------
#  ## not run
#  GSE <- metaboGSE(yarliSubmnets, method="perm", nperm=1000, nrand=1000,
#                         mc.cores=mc.cores, prefix="/tmp/summary")
#  ##

## ---- eval=TRUE, message=FALSE-------------------------------------------
data(yarliGSE)
data(yarliGOdata)
GSE <- yarliGSE
GOdata <- yarliGOdata
str(GSE[["GO:0006696"]], vec.len=2)
GSE[["GO:0006696"]]$res$p.Val

## ---- eval=TRUE, message=FALSE-------------------------------------------
GS.sig.all <- as.data.frame(t(sapply(GSE, function(gsm) {
    c(GS.ID=gsm$res$GS.ID, 
      Description=gsm$res$Description,
      Statistic=gsm$res$statistic,
      p.Cond=if (is.null(gsm$res$p.Cond)) NA else min(gsm$res$p.Cond),
      p.Val=gsm$res$p.Val,
      Genes=paste(intersect(genesInTerm(GOdata, gsm$res$GS.ID)[[1]], 
                            sybil::allGenes(hmodel.clean)), collapse=','))
})), stringsAsFactors=F)
GS.sig.all$FDR  <- p.adjust(as.numeric(GS.sig.all$p.Val), method="BH")
GS.sig.all <- GS.sig.all[!is.na(GS.sig.all$FDR), ]
dim(GS.sig.all)

## ---- eval=TRUE----------------------------------------------------------
GS.sig <- GS.sig.all[as.numeric(GS.sig.all$FDR) < 0.05, , drop=F]
GS.sig <- GS.sig[as.numeric(GS.sig$p.Cond) < 0.01, , drop=F]
dim(GS.sig)
