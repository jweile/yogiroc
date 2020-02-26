
#' new yogiroc
#' 
#' Create a new yogiroc object used to draw ROC and PRC curves and calculate the respective AUC.
#' 
#' The object exposes the following methods:
#' \itemize{
#' \item{draw.roc(col,main,lwd,add)} Draws a ROC curve. Parameter "col" controls the line color, 
#' "main" sets the plot title, "lwd" sets the line width, "add" allows the plot to be added as a line
#' to an existing plot.
#' \item{draw.prc(col,main,lwd,add)} Draws a PRC curve, Parameters as above.
#' \item{auroc()} Returns the Area under the ROC curve.
#' \item{auprc()} Returns the Area under the PRC curve.
#' }
#' 
#' @param truth a boolean vector for the ground truth classification
#' @param scores the scores to be evaluated agains the ground truth
#' @return A new yogiroc object
#' @export
#'    
new.yogiroc <- function(truth,scores) {

	#build a table for the ROC/PRC curves by iterating over all possible score thresholds
	data <- do.call(rbind,lapply(c(-Inf,sort(scores),Inf), function(t) {
		#which scores fall above the current threshold?
		calls <- scores >= t
		#calculate True Positives, True Negatives, False Positives and False Negatives
		tp <- sum(calls & truth)
		tn <- sum(!calls & !truth)
		fp <- sum(calls & !truth)
		fn <- sum(!calls & truth)
		#calculate PPV/precision, TPR/sensitivity/recall, and FPR/fallout
		ppv.prec <- tp/(tp+fp)
		tpr.sens <- tp/(tp+fn)
		fpr.fall <- fp/(tn+fp)
		#return the results
		cbind(thresh=t,ppv.prec=ppv.prec,tpr.sens=tpr.sens,fpr.fall=fpr.fall)
	}))

	#monotonize precision to correct for slumps in PRC curve
	for (i in 2:nrow(data)) {
		if (!is.na(data[i,"ppv.prec"]) && data[i,"ppv.prec"] < data[i-1,"ppv.prec"]) {
			data[i,"ppv.prec"] <- data[i-1,"ppv.prec"]
		}
		if (is.nan(data[nrow(data),"ppv.prec"])) {
			data[nrow(data),"ppv.prec"] <- data[nrow(data)-1,"ppv.prec"]
		}
	}
	
	#plot FPR vs sensitivity
	draw.roc <- function(col="steelblue3",main="",lwd=2,add=FALSE) {
		if (add) {
			lines(
				100*data[,"fpr.fall"],100*data[,"tpr.sens"],
				col=col,lwd=lwd
			)
		} else {
			plot(
				100*data[,"fpr.fall"],100*data[,"tpr.sens"],
				type="l",
				xlab="False positive rate (%)\n(= 100%-specificity)", ylab="Sensitivity or True positive rate (%)",
				main=main,
				xlim=c(0,100),ylim=c(0,100),col=col,lwd=lwd
			)
		}
	}

	#Plot Recall vs Precision
	draw.prc <- function(col="firebrick3",main="",lwd=2,add=FALSE) {
		if (add) {
			lines(
				100*data[,"tpr.sens"],100*data[,"ppv.prec"],
				col=col,lwd=lwd
			)
		} else {
			plot(
				100*data[,"tpr.sens"],100*data[,"ppv.prec"],
				type="l",
				xlab="Recall (%)", ylab="Precision (%)",
				main=main,
				xlim=c(0,100),ylim=c(0,100),col=col,lwd=lwd
			)
		}
	}

	#calculate AUC for ROC curve
	auroc <- function() {
		#iterate across x-axis intervals and sum across areas
		sum(sapply(1:(nrow(data)-1),function(i) {
			#calculate interval width between datapoints on x-axis
			delta.x <- data[i,"fpr.fall"]-data[i+1,"fpr.fall"]
			#calculate the average height of the two points on the y-axis
			y <- (data[i,"tpr.sens"]+data[i+1,"tpr.sens"])/2
			#area = x * y
			delta.x * y
		}))
	}

	#calculate AUC rof PRC curve
	auprc <- function() {
		# data <- na.omit(data)
		#same as in auroc, just using the appropriate data
		sum(sapply(1:(nrow(data)-1),function(i) {
			delta.x <- data[i,"tpr.sens"]-data[i+1,"tpr.sens"]
			y <- (data[i,"ppv.prec"]+data[i+1,"ppv.prec"])/2
			delta.x * y
		}))
	}

	#return maximum recall at given minimum precision
	recall.at.prec <- function(x=0.9) {
		if (any(data[,"ppv.prec"] > x)) {
			max(data[which(data[,"ppv.prec"] > x),"tpr.sens"])
		} else NA
	}

	#return object
	structure(list(
		draw.roc=draw.roc,
		draw.prc=draw.prc,
		auroc=auroc,
		auprc=auprc,
		recall.at.prec=recall.at.prec,
		get.table=function() data
	),class="yogiroc")
}

