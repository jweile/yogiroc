

#' confidence interval for precision or recall
#'
#' @param is the number of correct calls
#' @param ns the number of total calls
#' @param p the confidence interval probability range e.g. 0.025;0.975
#' @param res the resolution at which to sample the call rates
#' 
#' @export
#'
#' @return the confidence interval (numerical vector)
prcCI <- function(is, ns, p=c(0.025,0.975),res=0.001) {
  stopifnot(length(is)==length(ns))
  do.call(rbind,mapply(function(i,n) {
    rates <- seq(0,1,res)
    dens <- dbinom(i,n,rates)
    cdf <- c(0,cumsum(dens[-1]*res)/sum(dens*res))
    setNames(sapply(p,function(.p) rates[max(which(cdf < .p))]),p)
  },is,ns,SIMPLIFY=FALSE))
}

# prcCI <- function(i,n,p=c(0.025,0.975),res=0.001) {
#   rates <- seq(0,1,res)
#   dens <- dbinom(i,n,rates)
#   cdf <- c(0,cumsum(dens[-1]*res)/sum(dens*res))
#   # plot(rates,cdf,type="l")
#   sapply(p,function(.p) rates[max(which(cdf < .p))])
# }

#' Helper function to monotonize precision
#'
#' @param xs numerical input vector, representing precision ordered according to increasing t
#'
#' @return the monotonized equivalent vector
monotonize <- function(xs) {
  for (i in 2:length(xs)) {
    if (xs[[i]] < xs[[i-1]]) {
      xs[[i]] <- xs[[i-1]]
    }
  }
  xs
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
  
  #calculate and return the roc/prc tables for each score
  tables <- setNames(lapply(1:ncol(scores), function(coli) {
    #the sample prior is the share of true cases out of all cases
    appl <- which(!is.na(scores[,coli]))
    prior <- sum(truth[appl])/length(truth[appl])
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
      # precCI <- prcCI(tp,tp+fp)
      #return the results
      c(thresh=t,tp=tp,tn=tn,fp=fp,fn=fn,ppv.prec=ppv.prec,tpr.sens=tpr.sens,fpr.fall=fpr.fall,
            ppv.prec.balanced=ppv.prec.balanced#,precCI5=precCI[[1]],precCI95=precCI[[2]]
      )
    }))
    #set precision at infinite score threshold based on penultimate value
    data[nrow(data),c("ppv.prec","ppv.prec.balanced")] <- data[nrow(data)-1,c("ppv.prec","ppv.prec.balanced")]
    
    #monotonize precision to correct for slumps in PRC curve
    # data <- cbind(data,ppv.prec.mono=monotonize(data[,"ppv.prec"]))
    # data <- cbind(data,ppv.prec.balanced.mono=monotonize(data[,"ppv.prec.balanced"]))
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
  ppv <- function(i) {
    raw <- if (balanced) yr2[[i]][,"ppv.prec.balanced"] else yr2[[i]][,"ppv.prec"]
    if (monotonized) monotonize(raw) else raw
  }
  plot(
    100*yr2[[1]][,"tpr.sens"],100*ppv(1),
    type="l",
    xlab="Recall (%)", ylab="Precision (%)",
    xlim=c(0,100),ylim=c(0,100),col=col[[1]], ...
  )
  if(length(yr2) > 1) {
    for (i in 2:length(yr2)) {
      lines(
        100*yr2[[i]][,"tpr.sens"],100*ppv(i),
        col=col[[i]], ...
      )
    }
  } 
  if (!is.na(legend)) {
    legend(legend,sprintf("%s (AUPRC=%.02f;R90P=%.02f)",
           names(yr2),auprc(yr2,monotonized,balanced),recall.at.prec(yr2,0.9,monotonized,balanced)
    ),col=col,lty=1)
  }
}


#' Draw Precision-Recall Curve (PRC) with confidence intervals
#'
#' @param yr2 the yogiroc2 object
#' @param col vector of colors to use for the predictors
#' @param monotonized whether or not to monotonize the curve
#' @param legend the position of the legend, e.g. "bottomleft". NA disables legend
#' @param ... additional graphical parameters (see \code{par})
#'
#' @return nothing. draws a plot
#' @export
#'
#' @examples
#' #generate fake data
#' N <- 100
#' M <- 80
#' truth <- c(rep(TRUE,N),rep(FALSE,M))
#' scores <- cbind(
#'   pred1=c(rnorm(N,1,0.2),rnorm(M,.9,0.1)),
#'   pred2=c(rnorm(N,1.1,0.2),rnorm(M,.9,0.2))
#' )
#' #create yogiroc2 object
#' yrobj <- yr2(truth,scores)
#' #draw PRC curve
#' draw.prc.CI(yrobj)
#' #draw non-monotonized PRC curve
#' draw.prc.CI(yrobj,monotonized=FALSE)
draw.prc.CI <- function(yr2,col=seq_along(yr2),monotonized=TRUE,legend="bottomleft",...) {
  stopifnot(inherits(yr2,"yr2"))
  mon <- function(xs) if (monotonized) monotonize(xs) else xs
  plot(
    100*yr2[[1]][,"tpr.sens"],100*mon(yr2[[1]][,"ppv.prec"]),
    type="l",
    xlab="Recall (%)", ylab="Precision (%)",
    xlim=c(0,100),ylim=c(0,100),col=col[[1]], #...
  )
  if(length(yr2) > 1) {
    for (i in 2:length(yr2)) {
      lines(
        100*yr2[[i]][,"tpr.sens"],100*mon(yr2[[i]][,"ppv.prec"]),
        col=col[[i]], #...
      )
    }
  } 
  for (i in 1:length(yr2)) {
    x <- 100*yr2[[i]][,"tpr.sens"]
    precCI <- prcCI(yr2[[i]][,"tp"],yr2[[i]][,"tp"]+yr2[[i]][,"fp"])
    polygon(c(x,rev(x)),
            c(100*mon(precCI[,"0.025"]),rev(100*mon(precCI[,"0.975"]))),
            col=yogitools::colAlpha(col[[i]],0.1),border=NA
    )
  }
  if (!is.na(legend)) {
    legend(legend,sprintf("%s (AUPRC=%.02f;R90P=%.02f)",
        names(yr2),auprc(yr2,monotonized),recall.at.prec(yr2,0.9,monotonized)
    ),col=col,lty=1)
  }
}

#' Assess the significance of AUPRC differences
#' 
#' The list returned by this functions contains four elements:
#' \describe{
#' \item{auprc}{is simply the empirical area under the precision recall curve 
#' for each predictor.}
#' \item{ci}{is a matrix listing the lower and upper end of the 95% confidence 
#' interval for the AUPRC of each predictor.}
#' \item{llr}{is a matrix with columns and rows corresponding to each predictor.
#' It lists the log likelihood ratio of how much more (or less) likely the row-wise
#' predictor is to have a greater AUPRC than the column-wise predictor.}
#' \item{pval}{is a matrix with columns and rows corresponding to each predictor.
#' It lists the p-value of how likely it would be to observe the AUPRC of the row-wise
#' predictor under the distribution of the column-wise predictor.}
#' }
#'
#' @param yr2 the yogiroc2 object
#' @param monotonized whether or not to monotonize the curve
#' @param res the resolution at which to sample the probability function 
#' (defaults to 0.001)
#'
#' @return a list containing 4 elements: "auprc" (the empirical area under the 
#' precision recall curve), "ci" (the 95%confidence interval around the auprc), 
#' "llr" (the log likelihood ratio matrix, see details), and "pval" (the p-value 
#' of each auprc against each other)
#' @export
#'
#' @examples
#' #generate fake data
#' N <- 100
#' M <- 80
#' truth <- c(rep(TRUE,N),rep(FALSE,M))
#' scores <- cbind(
#'   pred1=c(rnorm(N,1,0.2),rnorm(M,.9,0.1)),
#'   pred2=c(rnorm(N,1.1,0.2),rnorm(M,.9,0.2))
#' )
#' #create yogiroc2 object
#' yrobj <- yr2(truth,scores)
#' auprc.signif(yrobj)
auprc.signif <- function(yr2,monotonized=TRUE,res=0.001) {
  #probability range
  ps <- seq(res,1-res,res)
  #calculate the AUPRC for each probability
  #i.e. the quantiles corresponding to ps
  auprcs <- do.call(cbind,lapply(1:length(yr2),function(i) {
    precCI <- prcCI(yr2[[i]][,"tp"],yr2[[i]][,"tp"]+yr2[[i]][,"fp"],p=ps)
    apply(precCI,2,function(ppv) {
      if (monotonized) {
        ppv <- monotonize(ppv)
      }
      calc.auc(yr2[[i]][,"tpr.sens"],ppv)
    })
  }))
  
  #empirical AUCs of the predictors
  empAUCs <- auprc(yr2,monotonized=monotonized)
  
  #1. build a reverse-lookup table that returns p for a given auprc
  aucRange <- do.call(seq,as.list(c(round(range(auprcs),digits=2),res)))
  aucPs <- do.call(rbind,lapply(aucRange,function(a){
    apply(auprcs,2,function(ladder){
      c(0,ps)[[sum(ladder < a)+1]]
    })
  }))
  
  confInts <- auprcs[c("0.025","0.975"),]
  colnames(confInts) <- names(yr2)
  
  pvals <- do.call(rbind,lapply(1:ncol(aucPs),function(i) {
    sapply(1:ncol(aucPs),function(j) {
      if (i==j) NA else {
        1-c(0,ps)[[sum(auprcs[,j] < empAUCs[[i]])+1]]
      }
    })
  }))
  dimnames(pvals) <- list(names(yr2), names(yr2))
  
  llrs <- do.call(rbind,lapply(1:ncol(aucPs),function(i) {
    sapply(1:ncol(aucPs),function(j) {
      if (i==j) NA else {
        #2. iterate over range of auprcs and calculate p_A(x) * (1-p_B(x)) 
        #  (i.e. the probability that area A is smaller than x AND area B is greater than x)
        log10(calc.auc(aucRange,aucPs[,j]*(1-aucPs[,i]))/calc.auc(aucRange,aucPs[,i]*(1-aucPs[,j])))
      }
    })
  }))
  dimnames(llrs) <- list(names(yr2), names(yr2))
  
  # plot(NA,type="n",xlim=c(0,1),ylim=c(0,1),xlab="AUPRC",ylab="CDF")
  # for (i in 1:ncol(auprcs)) {
  #   lines(c(0,auprcs[,i],1),c(0,ps,1),col=i)
  # }
  # abline(v=empAUCs,col=1:length(yr2),lty="dashed")
  
  return(list(auprc=empAUCs,ci=confInts,llr=llrs,pval=pvals))
}


#' Assess the significance of AUPRC against random guessing
#' #'
#' @param yr2 the yogiroc2 object
#' @param monotonized whether or not to monotonize the curves
#' @param cycles the sample size of the null distribution to use
#'
#' @return The emprical p-values of the AUPRC against random guessing
#' @export
#'
#' @examples
#' #generate fake data
#' N <- 10
#' M <- 8
#' truth <- c(rep(TRUE,N),rep(FALSE,M))
#' scores <- cbind(
#'   pred1=c(rnorm(N,1,0.2),rnorm(M,.9,0.1)),
#'   pred2=c(rnorm(N,1.1,0.2),rnorm(M,.9,0.2))
#' )
#' #create yogiroc2 object
#' yrobj <- yr2(truth,scores)
#' #call pvrandom function
#' auprc.pvrandom(yrobj)
auprc.pvrandom <- function(yr2,monotonized=TRUE,cycles=10000) {
  
  empAUCs <- auprc(yr2,monotonized=monotonized)
  
  #reconstruct truth table
  real <- yr2[[1]][1,"tp"]
  nreal <- yr2[[1]][1,"fp"]
  truth <- c(rep(TRUE,real),rep(FALSE,nreal))
  
  nullAucs <- replicate(cycles,{
    scores <- runif(real+nreal,0,1)
    ts <- na.omit(c(-Inf,sort(scores),Inf))
    pr <- t(sapply(sort(scores), function(t) {
      calls <- scores >= t
      tp <- sum(calls & truth,na.rm=TRUE)
      # tn <- sum(!calls & !truth,na.rm=TRUE)
      fp <- sum(calls & !truth,na.rm=TRUE)
      fn <- sum(!calls & truth,na.rm=TRUE)
      prec <- tp/(tp+fp)
      recall <- tp/(tp+fn)
      # fpr.fall <- fp/(tn+fp)
      c(prec,recall)
    }))
    if (monotonized){
      pr[,1] <- monotonize(pr[,1])
    }
    calc.auc(pr[,2],pr[,1])
  })
  
  pvals <- sapply(empAUCs, function(eauc) {
    sum(nullAucs >= eauc)/cycles
  })
  
  return(pvals)
}
# auprc.CI <- function(yr2,monotonized=TRUE) {
#   do.call(rbind,lapply(1:length(yr2),function(i) {
#     precCI <- prcCI(yr2[[i]][,"tp"],yr2[[i]][,"tp"]+yr2[[i]][,"fp"])
#     auprcs <- apply(precCI,2,function(ppv) {
#       if (monotonized) {
#         ppv <- monotonize(ppv)
#       }
#       calc.auc(yr2[[i]][,"tpr.sens"],ppv)
#     })
#   }))
# }

#' Helper function to calculate area under curve
#'
#' @param xs the x values of the graph
#' @param ys the corresponding y values of the graph
#'
#' @return the area under the curve
calc.auc <- function(xs,ys) {
  #calculate the sum of the areas of individual x-segments
  sum(sapply(1:(length(xs)-1),function(i) {
    #calculate interval width between datapoints on x-axis
    delta.x <- abs(xs[[i]]-xs[[i+1]])
    #calculate the average height of the two points on the y-axis
    y <- (ys[[i]]+ys[[i+1]])/2
    #area = x * y ; geometrically, this works out to be the same as the area of the polygon
    delta.x * y
  }))
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
  ppv <- function(data) {
    raw <- if (balanced) data[,"ppv.prec.balanced"] else data[,"ppv.prec"]
    if (monotonized) monotonize(raw) else raw
  }
  sapply(yr2,function(data) {
    calc.auc(data[,"tpr.sens"],ppv(data))
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
    calc.auc(data[,"fpr.fall"],data[,"tpr.sens"])
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
  ppv <- function(data) {
    raw <- if (balanced) data[,"ppv.prec.balanced"] else data[,"ppv.prec"]
    if (monotonized) monotonize(raw) else raw
  }
  sapply(yr2,function(data) {
    if (any(ppv(data) > x)) {
      max(data[which(ppv(data) > x),"tpr.sens"])
    } else NA
  })
}

  