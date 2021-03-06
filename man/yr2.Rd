% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/yr2.R
\name{yr2}
\alias{yr2}
\title{YogiRoc2 object constructor}
\usage{
yr2(truth, scores, names = colnames(scores), high = TRUE)
}
\arguments{
\item{truth}{a boolean vector indicating the classes of the reference set}

\item{scores}{a matrix of scores, with rows for each entry in truth, and one column for each predictor}

\item{names}{the names of the predictors}

\item{high}{a boolean vector indicating for each predictor whether its scoring high-to-low (or low-to-high)}
}
\value{
a yogiroc2 object
}
\description{
YogiRoc2 object constructor
}
\examples{
#generate fake data
truth <- c(rep(TRUE,10),rep(FALSE,8))
scores <- cbind(
  pred1=c(rnorm(10,1,0.2),rnorm(8,.9,0.1)),
  pred2=c(rnorm(10,1.1,0.2),rnorm(8,.9,0.2))
)
#create yogiroc2 object
yrobj <- yr2(truth,scores)
#draw PRC curve
draw.prc(yrobj)
#calculate recall at 90\% precision
recall.at.prec(yrobj,0.9)
}
