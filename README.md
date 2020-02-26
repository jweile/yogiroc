# yogiroc
Simple ROC and PRC curves

Example:
```R
#install library
devtools::install_github("jweile/yogiroc")

#load library
library(yogiroc)

#create some fake example data
fakeData <- data.frame(
  pathogenic = c(rep(TRUE,10),rep(FALSE,10)),
  prediction1 = c(rnorm(10,mean=2,sd=1),rnorm(10,mean=1,sd=1)),
  prediction2 = c(rnorm(10,mean=1.5,sd=1),rnorm(10,mean=1,sd=1))
)

#create yogiroc objects for the two predictions in the fake data
yr1 <- new.yogiroc(truth=fakeData$pathogenic,scores=fakeData$prediction1)
yr2 <- new.yogiroc(truth=fakeData$pathogenic,scores=fakeData$prediction2)

#draw the PRC curves
yr1$draw.prc(col="black")
yr2$draw.prc(col="red",add=TRUE)

#draw legend with AUPRC values
legend("bottomleft",c(
    sprintf("Prediction 1 (AUPRC=%.02f)",yr1$auprc()),
    sprintf("Prediction 2 (AUPRC=%.02f)",yr2$auprc())
), col=c("black","red"), lwd=2)
```
