

#' confidence interval for precision or recall
#'
#' @param i the number of correct calls
#' @param n the number of total calls
#' @param p the confidence interval probability range e.g. 0.025;0.975
#' @param res the resolution at which to sample the call rates
#'
#' @return the confidence interval (numerical vector)
prcCI <- function(i,n,p=c(0.025,0.975),res=0.001) {
  rates <- seq(0,1,res)
  dens <- dbinom(i,n,rates)
  cdf <- c(0,cumsum(dens[-1]*res)/sum(dens*res))
  # plot(rates,cdf,type="l")
  sapply(p,function(.p) rates[max(which(cdf < .p))])
}

#' YogiRoc2 object constructor
#'
#' @param truth a boolean vector indicating the classes of the reference set
#' @param scores a matrix of scores, with rows for each entry in truth, and one column for each predictor
#' @param names the names of the predictors
#' @param high a boolean vector indicating for each predictor whether its scoring high-to-low (or low-to-high)
#'
#' @return a yogiroc2 object
#' @export
#'
#' @examples
#' #generate fake data
#' truth <- c(rep(TRUE,10),rep(FALSE,8))
#' scores <- cbind(
#'   pred1=c(rnorm(10,1,0.2),rnorm(8,.9,0.1)),
#'   pred2=c(rnorm(10,1.1,0.2),rnorm(8,.9,0.2))
#' )
#' #create yogiroc2 object
#' yrobj <- yr2(truth,scores)
#' #draw PRC curve
#' draw.prc(yrobj)
#' #calculate recall at 90% precision
#' recall.at.prec(yrobj,0.9)
yr2 <- function(truth, scores, names=colnames(scores), high=TRUE) {
  
  #make sure all input is of correct datatype
  stopifnot(is.logical(truth), 
            is.data.frame(scores) || is.matrix(scores), 
            is.numeric(scores[1,1]), 
            is.character(names),
            is.logical(high)
  )
  #make sure all input is of correct size
  stopifnot(length(truth) == nrow(scores),
            length(names) == ncol(scores),
            length(high) == 1 || length(high) == ncol(scores)
  )
  
  #apply flipping to any scores that are not high-to-low
  if (length(high) == 1 && !high) {
    scores <- -scores
  } else if (any(!high)) {
    scores[,which(!high)] <- -scores[,which(!high)]
  }
  
  #the sample prior is the share of true cases out of all cases
  prior <- sum(truth)/length(truth)

  #calculate and return the roc/prc tables for each score
  tables <- setNames(lapply(1:ncol(scores), function(coli) {
    #build a table for the ROC/PRC curves by iterating over all possible score thresholds
    ts <- na.omit(c(-Inf,sort(scores[,coli]),Inf))
    data <- do.call(rbind,lapply(ts, function(t) {
      #which scores fall above the current threshold?
      calls <- scores[,coli] >= t
      #calculate True Positives, True Negatives, False Positives and False Negatives
      tp <- sum(calls & truth,na.rm=TRUE)
      tn <- sum(!calls & !truth,na.rm=TRUE)
      fp <- sum(calls & !truth,na.rm=TRUE)
      fn <- sum(!calls & truth,na.rm=TRUE)
      #calculate PPV/precision, TPR/sensitivity/recall, and FPR/fallout
      ppv.prec <- tp/(tp+fp)
      tpr.sens <- tp/(tp+fn)
      fpr.fall <- fp/(tn+fp)
      ppv.prec.balanced <- ppv.prec*(1-prior)/(ppv.prec*(1-prior)+(1-ppv.prec)*prior)
      #assess the confidence intervals
      precCI <- prcCI(tp,tp+fp)
      #return the results
      cbind(thresh=t,ppv.prec=ppv.prec,tpr.sens=tpr.sens,fpr.fall=fpr.fall,
            ppv.prec.balanced=ppv.prec.balanced,precCI5=precCI[[1]],precCI95=precCI[[2]]
      )
    }))
    #set precision at infinite score threshold based on penultimate value
    data[nrow(data),c("ppv.prec","ppv.prec.balanced")] <- data[nrow(data)-1,c("ppv.prec","ppv.prec.balanced")]
    
    monotonize <- function(xs) {
      for (i in 2:length(xs)) {
        if (xs[[i]] < xs[[i-1]]) {
          xs[[i]] <- xs[[i-1]]
        }
      }
      xs
    }
    
    #monotonize precision to correct for slumps in PRC curve
    data <- cbind(data,ppv.prec.mono=monotonize(data[,"ppv.prec"]))
    data <- cbind(data,ppv.prec.balanced.mono=monotonize(data[,"ppv.prec.balanced"]))
    return(data)
    
  }),names)
  
  return(structure(tables,class="yr2"))
}


#' print method for yogiroc2 objects
#'
#' @param yr2 the object
#'
#' @return nothing, just prints a description
#' @export
#'
#' @examples
print.yr2 <- function(yr2) {
  cat("YogiROC object\n")
  cat("Reference set size:",nrow(yr2[[1]]-2),"\n")
  cat("Predictors:",paste(names(yr2),collapse=", "),"\n")
}


#' Draw a ROC curve
#'
#' @param yr2 an underlying yogiroc2 object
#' @param col the colors to use for the predictors
#' @param legend the positioning of the legend (e.g. "bottomright). NA to disable legend.
#' @param ... additional graphical parameters (see \code{par})
#'
#' @return nothing, draws a plot
#' @export
#'
#' @examples
#' #generate fake data
#' truth <- c(rep(TRUE,10),rep(FALSE,8))
#' scores <- cbind(
#'   pred1=c(rnorm(10,1,0.2),rnorm(8,.9,0.1)),
#'   pred2=c(rnorm(10,1.1,0.2),rnorm(8,.9,0.2))
#' )
#' #create yogiroc2 object
#' yrobj <- yr2(truth,scores)
#' #draw PRC curve
#' draw.roc(yrobj)
draw.roc <- function(yr2,col=seq_along(yr2),legend="bottomright",...) {
  stopifnot(inherits(yr2,"yr2"))
  plot(
    100*yr2[[1]][,"fpr.fall"],100*yr2[[1]][,"tpr.sens"],
    type="l",
    xlab="False positive rate (%)\n(= 100%-specificity)", ylab="Sensitivity or True positive rate (%)",
    xlim=c(0,100),ylim=c(0,100),col=col[[1]], ...
  )
  if(length(yr2) > 1) {
    for (i in 2:length(yr2)) {
      lines(
        100*yr2[[i]][,"fpr.fall"],100*yr2[[i]][,"tpr.sens"],
        col=col[[i]], ...
      )
    }
  }
  if (!is.na(legend)) {
    legend(legend,sprintf("%s (AUROC=%.02f)",names(yr2),auroc(yr2)),col=col,lty=1)
  }
}


#' Draw Precision-Recall Curve (PRC)
#'
#' @param yr2 the yogiroc2 object
#' @param col vector of colors to use for the predictors
#' @param monotonized whether or not to monotonized the curve
#' @param balanced whether or not to use prior-balancing
#' @param legend the position of the legend, e.g. "bottomleft". NA disables legend
#' @param ... additional graphical parameters (see \code{par})
#'
#' @return nothing. draws a plot
#' @export
#'
#' @examples
#' #generate fake data
#' truth <- c(rep(TRUE,10),rep(FALSE,8))
#' scores <- cbind(
#'   pred1=c(rnorm(10,1,0.2),rnorm(8,.9,0.1)),
#'   pred2=c(rnorm(10,1.1,0.2),rnorm(8,.9,0.2))
#' )
#' #create yogiroc2 object
#' yrobj <- yr2(truth,scores)
#' #draw PRC curve
#' draw.prc(yrobj)
#' #draw non-monotonized PRC curve
#' draw.prc(yrobj,monotonized=FALSE)
#' #draw balanced PRC curve
#' draw.prc(yrobj,balanced=TRUE)
draw.prc <- function(yr2,col=seq_along(yr2),monotonized=TRUE,balanced=FALSE,legend="bottomleft",...) {
  stopifnot(inherits(yr2,"yr2"))
  ppvcol <- if (balanced) {
    if (monotonized) "ppv.prec.balanced.mono" else "ppv.prec.balanced"
  } else {
    if (monotonized) "ppv.prec.mono" else "ppv.prec"
  }
  plot(
    100*yr2[[1]][,"tpr.sens"],100*yr2[[1]][,ppvcol],
    type="l",
    xlab="Recall (%)", ylab="Precision (%)",
    xlim=c(0,100),ylim=c(0,100),col=col[[1]], ...
  )
  if(length(yr2) > 1) {
    for (i in 2:length(yr2)) {
      lines(
        100*yr2[[i]][,"tpr.sens"],100*yr2[[i]][,ppvcol],
        col=col[[i]], ...
      )
    }
  } 
  if (!is.na(legend)) {
    legend(legend,sprintf("%s (AUPRC=%.02f)",names(yr2),auprc(yr2,monotonized)),col=col,lty=1)
  }
}


#' Calculate area under precision recall curve (AUPRC)
#'
#' @param yr2 the yogiroc2 object
#' @param monotonized whether or not use a monotonized PRC curve
#' @param balanced whether or not to use prior-balancing
#'
#' @return a numerical vector with the AUPRC values for each predictor
#' @export
#'
#' @examples
#' #generate fake data
#' truth <- c(rep(TRUE,10),rep(FALSE,8))
#' scores <- cbind(
#'   pred1=c(rnorm(10,1,0.2),rnorm(8,.9,0.1)),
#'   pred2=c(rnorm(10,1.1,0.2),rnorm(8,.9,0.2))
#' )
#' #create yogiroc2 object
#' yrobj <- yr2(truth,scores)
#' #calculate AUPRC
#' auprc(yrobj)
#' #calculate non-monotonized AUPRC
#' auprc(yrobj,monotonized=FALSE)
#' #calculate balanced AUPRC
#' auprc(yrobj,balanced=TRUE)
auprc <- function(yr2, monotonized=TRUE, balanced=FALSE) {
  stopifnot(inherits(yr2,"yr2"))
  ppvcol <- if (balanced) {
    if (monotonized) "ppv.prec.balanced.mono" else "ppv.prec.balanced"
  } else {
    if (monotonized) "ppv.prec.mono" else "ppv.prec"
  }
  sapply(yr2,function(data) {
    sum(sapply(1:(nrow(data)-1),function(i) {
      delta.x <- data[i,"tpr.sens"]-data[i+1,"tpr.sens"]
      y <- (data[i,ppvcol]+data[i+1,ppvcol])/2
      delta.x * y
    }))
  })
}
  

#' Calculate the area under the ROC curve
#'
#' @param yr2 the yogiroc2 object
#'
#' @return a numerical vector with the AUROC for each predictor
#' @export
#'
#' @examples
#' #generate fake data
#' truth <- c(rep(TRUE,10),rep(FALSE,8))
#' scores <- cbind(
#'   pred1=c(rnorm(10,1,0.2),rnorm(8,.9,0.1)),
#'   pred2=c(rnorm(10,1.1,0.2),rnorm(8,.9,0.2))
#' )
#' #create yogiroc2 object
#' yrobj <- yr2(truth,scores)
#' #calculate AUROC
#' auroc(yrobj)
auroc <- function(yr2) {
  stopifnot(inherits(yr2,"yr2"))
  sapply(yr2,function(data) {
    #iterate across x-axis intervals and sum across areas
    sum(sapply(1:(nrow(data)-1),function(i) {
      #calculate interval width between datapoints on x-axis
      delta.x <- data[i,"fpr.fall"]-data[i+1,"fpr.fall"]
      #calculate the average height of the two points on the y-axis
      y <- (data[i,"tpr.sens"]+data[i+1,"tpr.sens"])/2
      #area = x * y
      delta.x * y
    }))
  })
}

#' Calculate maximum recall at given minimum precision
#'
#' @param yr2 the yogiroc2 object
#' @param x the precision cutoff (default 0.9)
#' @param monotonized whether or not to use monotonized PRC
#' @param balanced whether or not to use prior-balancing
#'
#' @return
#' @export
#'
#' @examples
#' #generate fake data
#' truth <- c(rep(TRUE,10),rep(FALSE,8))
#' scores <- cbind(
#'   pred1=c(rnorm(10,1,0.2),rnorm(8,.9,0.1)),
#'   pred2=c(rnorm(10,1.1,0.2),rnorm(8,.9,0.2))
#' )
#' #create yogiroc2 object
#' yrobj <- yr2(truth,scores)
#' #calculate R90P
#' recall.at.prec(yrobj)
#' #calculate non-monotonized R90P
#' recall.at.prec(yrobj,monotonized=FALSE)
#' #calculate balanced R90P
#' recall.at.prec(yrobj,balanced=TRUE)
recall.at.prec <- function(yr2,x=0.9,monotonized=TRUE,balanced=FALSE) {
  stopifnot(inherits(yr2,"yr2"))
  ppvcol <- if (balanced) {
    if (monotonized) "ppv.prec.balanced.mono" else "ppv.prec.balanced"
  } else {
    if (monotonized) "ppv.prec.mono" else "ppv.prec"
  }
  sapply(yr2,function(data) {
    if (any(data[,"ppv.prec"] > x)) {
      max(data[which(data[,"ppv.prec"] > x),"tpr.sens"])
    } else NA
  })
}

  