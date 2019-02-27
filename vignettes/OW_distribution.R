## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "##"
)

## ---- echo=FALSE, message=F----------------------------------------------
require(RelDists)
## 5^(1/3) = 1.71; 2^(1/3) = 1.26; 0.5^(1/2) = 0.71
curve(dOW(x, 1.71, 3, 0.7), from=0, to=3.5, ylim=c(0,1.6), col="red",
          ylab=expression(f[' '] (x)), xlab="x", las=1, lwd=2)
curve(dOW(x, 1.26, 3, 0.2), from=0, to=3.5, col="blue", 
      add=TRUE, lwd=2, lty=3) 
curve(dOW(x,0.71, 2, 1.3), from=0, to=3.5, col="darkgreen", 
      add=TRUE, lwd=2, lty=2)
curve(dOW(x,1.7,1,0.8), from=0, to=3.5, col="purple", 
      add=TRUE, lwd=2, lty=4)
legend("topright", legend=c(expression(paste(mu,
       " = ", 1.71,"      " ,sigma," = ",3,"   " ,nu," = ", 0.7)),
                            expression(paste(mu, 
       " = ", 1.26,"      " ,sigma," = ",3,"   " ,nu," = ", 0.2)),
                            expression(paste(mu,
       " = ", 0.71,"      " ,sigma," = ",2,"   " ,nu," = ", 1.3)),
                            expression(paste(mu,
       " = ", 1.7,"        " ,sigma," = ",1,"   " ,nu," = ", 0.8))),
       col=c("red","blue","darkgreen","purple"),lty=c(1,3,2,4), 
       bty="n", lwd=2)

## ---- echo=FALSE---------------------------------------------------------
## 5^(1/3) = 1.71; 2^(1/3) = 1.26; 0.5^(1/2) = 0.71
curve(pOW(x, 1.71, 3, 0.7), from=0, to=3.5, col="red",
      ylab=expression(F[' '] (x)), xlab="x", las=1, lwd=2)
curve(pOW(x, 1.26, 3, 0.2), from=0, to=3.5, col="blue", 
      add=TRUE, lwd=2, lty=3) 
curve(pOW(x,0.71, 2, 1.3), from=0, to=3.5, col="darkgreen", 
      add=TRUE, lwd=2, lty=2)
curve(pOW(x,1.7,1,0.8), from=0, to=3.5, col="purple", 
      add=TRUE, lwd=2, lty=4)
legend("bottomright", legend=c(expression(paste(mu, 
       " = ", 1.71,"      " ,sigma," = ",3,"   " ,nu," = ", 0.7)),
                            expression(paste(mu,
       " = ", 1.26,"      " ,sigma," = ",3,"   " ,nu," = ", 0.2)),
                            expression(paste(mu, 
       " = ", 0.71,"      " ,sigma," = ",2,"   " ,nu," = ", 1.3)),
                            expression(paste(mu,
       " = ", 1.7,"        " ,sigma," = ",1,"   " ,nu," = ", 0.8))),
       col=c("red","blue","darkgreen","purple"),lty=c(1,3,2,4), 
       bty="n", lwd=2)

## ---- echo=FALSE---------------------------------------------------------
par(mgp=c(3,0.7,0))
curve(hOW(x, 1/85, 9, 0.7), from=0, to=69, ylim=c(0,0.035), col="red",
      ylab=expression(h[' '] (x)), xlab="x", las=1, lwd=2, lty=1)
curve(hOW(x, 1/100, 0.5, 0.3), from=0, to=70, col="blue", 
      add=TRUE, lwd=2, lty=3) 
curve(hOW(x, 1/50, 1, 1), from=0, to=70, col="darkgreen", 
      add=TRUE, lwd=2, lty=2)
curve(hOW(x, 1/45, 8, 0.01), from=0, to=69, col="purple", 
      add=TRUE, lwd=2, lty=4)
curve(hOW(x, 1/75, -1.5, -0.1), from=0, to=70, col="black", 
      add=TRUE, lwd=2)
legend("topright", legend=c(expression(paste(mu, 
       " = ", 0.012,"     " ,sigma," = ",9,"        " ,nu," = ", 0.7)),
                            expression(paste(mu, 
       " = ", 0.01,"       " ,sigma," = ",0.5,"     " ,nu," = ", 0.3)),
                            expression(paste(mu, 
       " = ", 0.02,"       " ,sigma," = ",1,"        " ,nu," = ", 1)),
                            expression(paste(mu, 
       " = ", 0.0222,"  " ,sigma," = ",8,"         " ,nu," = ", 0.01)),
                            expression(paste(mu, 
       " = ", 0.0133,"  " ,sigma," = ",-1.5,"   " ,nu," = ", -0.1))), 
       col=c("red","blue","darkgreen","purple","black"),lty=c(1,3,2,4,1), 
       bty="n", lwd=2)

