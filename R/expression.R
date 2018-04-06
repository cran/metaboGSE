#' Relative expression data (log2 counts)
#'
#' This function computes the relative gene expression (log2 counts) in a sample w.r.t all samples.
#' @param expr A matrix of expression data (log2 counts), with gene per row and sample per column.
#' @param power The power n of relative expression: (x^n/<x>)^(1/n). Default: 1.
#' @return A matrix of relative log2 counts.
#' @examples 
#' data(exprMaguire)
#' relExpr1 <- relativeExpr(exprMaguire$expr, 1)
#' @export
relativeExpr <- function(expr, power = 1) {
    expr.avg <- apply(expr, 1, mean, na.rm=T)
    res <- apply(expr, 2, function(x) {
        res.1 <- (x^power / expr.avg)^(1/power)
        res.1[is.na(res.1)] <- 0
        
        return (res.1)
    })
    
    return (res)
}


#' Z-score of expression data
#'
#' This function computes the z-score of expression data across different conditions.
#' @param expr A matrix of expression data (log2 counts), with gene per row and experiment per column.
#' @return A matrix of z-score.
#' @examples 
#' data(exprMaguire)
#' zExpr <- zscoreExpr(exprMaguire$expr)
#' @export
zscoreExpr <- function(expr) {
    return (t(apply(expr, 1, function(x) {sapply(x, function(y) {(y-mean(x))/(sd(x))})})))
}
