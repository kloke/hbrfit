# Based on weighted Wilcoxon code developed by Terpstra & McKean (2005) 
# (http://www.jstatsoft.org/v14/i07) under GPL 2.
psi <-
function (x) 
{
    x[x == -Inf] = -100
    x[x == Inf] = 100
    -1 * (x <= -1) + x * (-1 < x & x < 1) + 1 * (x >= 1)
}
