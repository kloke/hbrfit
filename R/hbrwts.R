# Based on weighted Wilcoxon code developed by Terpstra & McKean (2005) 
# (http://www.jstatsoft.org/v14/i07) under GPL 2.
hbrwts <-
function (x, y, percent = 0.95, ehat0 = ltsreg(x, y)$residuals) 
{
    x = as.matrix(x)
    n = dim(x)[1]
    p = dim(x)[2]
    robdis2 <- robdist.hbrfit(x)
    cut = qchisq(percent, p)
    sigma = mad(ehat0)
    m = psi(cut/robdis2)
    a = ehat0/(sigma * m)
    c = (median(a) + 3 * mad(a))^2
    h = sqrt(c)/a
    tmp = pairup(h)
    psi(abs(tmp[, 1] * tmp[, 2]))
}
