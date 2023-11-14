fitdiag <- function(x,y,est=c("WIL","HBR"),...) {

  centerx <- function(x) scale(x,scale=FALSE)

  x=as.matrix(centerx(x))
  n=dim(x)[1]
  p=dim(x)[2]
  # Wilcoxon estimate
  tempw=rfit(y ~ x,...)
  residw=tempw$residuals
  vcw=vcov(tempw)
  temphbr=NULL
  templs=NULL
  if (any("WIL"==est) & any("HBR"==est)) {
    temphbr=hbrfit(y ~ x)
    diff=tempw$coef - temphbr$coef
  }  
  if (any("WIL"==est) & any("LS"==est)) {
    templs=lm(y~x)
    diff=tempw$coef - templs$coef
  }
  if (any("HBR"==est) & any("LS"==est)) {
    temphbr=hbrfit(y ~ x)
    templs=lm(y~x)
    diff=temphbr$coef - templs$coef
  }
  if (any("WIL"==est) & any("LTS"==est)) {
    templts=ltsreg(x,y)
    diff=tempw$coef - templts$coef
  }
  if (any("HBR"==est) & any("LTS"==est)) {
    temphbr=hbrfit(y~x)
    templts=ltsreg(x,y)
    diff=temphbr$coef - templts$coef
  }
 if (any("LS"==est) & any("LTS"==est)) {
    templs=lm(y~x)
    templts=ltsreg(x,y)
    diff=templs$coef - templts$coef
  }


  tdbeta=t(cbind(diff))%*%solve(vcw)%*%cbind(diff)
  bmtd=(4*(p+1)^2)/n
  xmat=cbind(rep(1,n),x)
  diffc=xmat%*%diff
  diffvc=xmat%*%vcw%*%t(xmat)
  cfit=diffc/(sqrt(diag(diffvc)))
  bmcf=2*sqrt((p+1)/n)
  se=sqrt(diag(vcw))  
  list(tdbeta=tdbeta,bmtd=bmtd,cfit=cfit,bmcf=bmcf,est=est)
}
