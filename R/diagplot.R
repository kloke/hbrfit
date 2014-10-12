# Based on weighted Wilcoxon code developed by Terpstra & McKean (2005) 
# (http://www.jstatsoft.org/v14/i07) under GPL 2.
diagplot <-
function (fit, ...) 
{
    par(mfrow = c(2, 2))
    plot(fitted.values(fit), residuals(fit), xlab = "Fit", ylab = "Residual", 
        main = "Residuals vs. Fits")
    hist(residuals(fit), freq = FALSE, , main = "Histogram of Residuals", 
        xlab = "Residual")
    plot(rstudent(fit), xlab = "Case", ylab = "Studentized Residual", 
        main = "Case Plot of\nStudentized Residuals")
    abline(h = c(-2, 2))
    qqnorm(residuals(fit), main = "Normal Q-Q Plot of Residuals")
    qqline(residuals(fit))
}
