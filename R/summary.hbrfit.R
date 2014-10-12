# Based on weighted Wilcoxon code developed by Terpstra & McKean (2005) 
# (http://www.jstatsoft.org/v14/i07) under GPL 2.
summary.hbrfit <- function (object, ...) {
   fit<-object
    wald.test = function(est, varcov, amat, true, n) {
        true = as.matrix(true)
        est = as.matrix(est)
        amat = as.matrix(amat)
        p = dim(est)[1] - 1
        q = dim(amat)[1]
        temp1 = as.matrix(amat %*% est - true)
        temp2 = as.matrix(amat %*% varcov %*% t(amat))
        T2 = t(temp1) %*% solve(temp2) %*% temp1
        T2 = T2/q
        pvalue = 1 - pf(T2, q, n - p - 1)
        list(stat = T2, pvalue = pvalue)
    }
    n = length(fit$y)
    p = fit$qrx1$rank - 1
    vhat <- vcov(fit)
    ans = cbind(fit$coef, sqrt(diag(vhat)))
    ans = cbind(ans, ans[, 1]/ans[, 2])
    ans = cbind(ans, 2 * pt(abs(ans[, 3]), n - p - 1, lower.tail = FALSE))
    coef <- ans
    colnames(coef) <- c("Estimate", "Std. Error", "t.value", 
        "p.value")
    A <- matrix(diag(c(0, rep(1, p)))[-1, ], nrow = p)
    wt <- wald.test(fit$coef, vhat, matrix(diag(c(0, rep(1, p)))[-1, 
        ], nrow = p), rep(0, p), n)
    res <- list(coefficients = coef, dropstat = wt$stat, droppval = wt$pval, 
        varcov = vhat)
    res$call <- fit$call
    class(res) <- "summary.hbrfit"
    res
}
