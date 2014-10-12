print.summary.hbrfit <-
function (x, digits = max(5, .Options$digits - 2), ...) 
{
    cat("Call:\n")
    print(x$call)
    cat("\nCoefficients:\n")
    printCoefmat(x$coefficients, P.values = TRUE, has.Pvalue = TRUE)
    cat("\nWald Test:", round(x$dropstat, digits = digits), "p-value:", 
        round(x$droppval, digits = digits), "\n")
    cat("\n")
}
