

# truth <- c(rep(TRUE,10),rep(FALSE,8))
# scores <- cbind(
#   c(rnorm(10,1,0.2),rnorm(8,0.8,0.1)),
#   c(rnorm(10,1,0.2),rnorm(8,0.8,0.2))
# )

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
      #return the results
      cbind(thresh=t,ppv.prec=ppv.prec,tpr.sens=tpr.sens,fpr.fall=fpr.fall,ppv.prec.balanced=ppv.prec.balanced)
    }))
    #set precision at infinite score threshold to 1
    data[nrow(data),c("ppv.prec","ppv.prec.balanced")] <- 1
    
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

print.yr2 <- function(yr2) {
  cat("YogiROC object\n")
  cat("Reference set size:",nrow(yr2[[1]]-2),"\n")
  cat("Predictors:",paste(names(yr2),collapse=", "),"\n")
}


#plot FPR vs sensitivity
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


#Plot Recall vs Precision
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
  
#calculate AUC for ROC curve
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

#return maximum recall at given minimum precision
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

  