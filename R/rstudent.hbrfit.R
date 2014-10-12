# Based on weighted Wilcoxon code developed by Terpstra & McKean (2005) 
# (http://www.jstatsoft.org/v14/i07) under GPL 2.
rstudent.hbrfit <-
function (model, ...) 
{
	fit<-model
	Q<-qr.Q(fit$qrx1)
	x<-as.matrix(Q[,2:fit$qrx1$rank])
	q1<-Q[,1]
#    x = as.matrix(fit$x)
    bmat = as.matrix(fit$weights)
    res = as.vector(fit$residuals)
    n = dim(x)[1]
    p = dim(x)[2]
    sigma2 = (mad(res))^2
    tau <- fit$tauhat
    tau1 <- fit$taushat
    K1 = (n/(n - p - 1)) * mean(abs(res))
    K2 = 2 * mean((rank(res)/(n + 1) - 0.5) * res)
#    H = x %*% solve(t(x) %*% x) %*% t(x)
#	H <- x%*%t(x)
	H<-tcrossprod(x)
    I = diag(n)
#    J = matrix(1/n, n, n)
	J<-tcrossprod(q1)
    diag(bmat) = rep(0, n)
    w = -1 * bmat
    diag(w) = bmat %*% as.matrix(rep(1, n))
    w = (1/(sqrt(12) * tau)) * w
    cmat = (1/n^2) * t(x) %*% w %*% x
    cinv = solve(cmat)
    u = (1/n) * (bmat - diag(c(bmat %*% cbind(rep(1, n))))) %*% 
        x
    u = u * (1 - 2 * rank(res)/n)
    vmat = var(u)
    v = sigma2 * I + tau1^2 * J + (1/4) * x %*% ((1/n^2) * cinv) %*% 
        vmat %*% ((1/n^2) * cinv) %*% t(x) - 2 * tau1 * K1 * 
        J - sqrt(12) * tau * K2 * (w %*% x %*% ((1/n^2) * cinv) %*% 
        t(x) + x %*% ((1/n^2) * cinv) %*% t(x) %*% w)
    diag(v)[diag(v) <= 0] = sigma2 * diag(I - (1/n + H))[diag(v) <= 
        0]
    as.vector(res/sqrt(diag(v)))
}
