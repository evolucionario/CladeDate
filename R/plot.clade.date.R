#' @title Plot Empirical Calibration Densities

#' @description
#' This function plots empirical calibration densities generated using [clade.date].

#' @param object an object of class \code{"clade.date"}.
#' @param breaks Specifies the algorithm used to determine the number of histogram cells in the plot (see [hist] for alternatives, if needed).
#' @param ... additional arguments passed to [par].

#' @details
#' This function plots empirical calibration densities from a \code{clade.date} object.

#' @return Creates a plot from a \code{clade.date} object.

#' @import stats grDevices graphics

#' @examples
#'
#' ## Matrix of fossil ages including maximum and minimum bounds for each
#' Fages <- cbind(c(50, 38, 30, 14, 3.5), c(56, 45, 35, 14, 6))
#'
#' ## Clade.date analysis
#' X <- clade.date(Fages)
#'
#' ## Plot
#' plot(X)
#'
#' @author Santiago Claramunt, \email{sclaramunt@@rom.on.ca}

#' @references

#' Claramunt, S. *In review* CladeDate: : empirical calibration density generator for divergence time estimation
#'
#'
#' @seealso [clade.date] for generating empirical calibration densities.

#' @export

plot.clade.date <- function(object, breaks="FD", ...) {
	
		if(!is(object, "clade.date")) { stop("object must be of class clade.date") }
		
		ages <- object$ages # make clade.date return ages
		
		rA <- object$rep.values
		
		# Set x limits
		XX <- c(max(0, min(ages)-min(ages)/8), quantile(rA, prob=0.98))
		
		quartz("CladeDate", 8, 3.5)
		par(lend=1, las=1, mar=c(5,5,1,1), ...)

		HIST <- hist(rA, breaks=breaks, freq=FALSE, xlim=XX, col=rgb(.25,.7,.9,0.5), border="#4DCCFF", xlab="Time", main="", yaxs="i")


		if(ncol(ages) == 2) {
			# Plot fossil age intervals
			rect(xleft=ages[,1], ybottom=0, xright=ages[,2], ytop=max(HIST$density)/40, col=gray(0.5,0.75), border=NA)
			}
		
				
		if(!is.null(object$PDFfit)) {
						
			if(object$PDFfit.model=="lognormal") {
				curve(dlnorm(x-max(ages[,1]), meanlog=object$PDFfit$estimate[1], sdlog=object$PDFfit$estimate[2]), n=301, from=max(ages[,1]), to=XX[2], col="#CC6600", lwd=3, add=TRUE, xpd=TRUE)
				legend("right", legend="log-normal density", text.col="#CC6600", bty="n")
			}

			if(object$PDFfit.model=="gamma") {
				curve(dgamma(x-max(ages[,1]), shape=object$PDFfit$estimate[1], rate=object$PDFfit$estimate[2]), n=301, from=max(ages[,1]), to=XX[2], col="#CC6600", lwd=3, add=TRUE, xpd=TRUE)
				legend("right", legend="gamma density", text.col="#CC6600", bty="n")
			}

			if(object$PDFfit.model=="exponential") {
				curve(dexp(x-max(ages[,1]), rate=object$PDFfit$estimate[1]), n=301, from=max(ages[,1]), to=XX[2], col="#CC6600", lwd=3, add=TRUE, xpd=TRUE)
				legend("right", legend="exponential density", text.col="#CC6600", bty="n")
			}
			if(object$PDFfit.model=="skewnormal") {
				curve(dsnorm(x, mean=object$PDFfit$par["mean"], sd = object$PDFfit$par["sd"], xi = object$PDFfit$par["xi"]), n=301, from=XX[1], to=XX[2], col="#CC6600", lwd=3, add=TRUE, xpd=TRUE)
				legend("right", legend="skew-normal density", text.col="#CC6600", bty="n")
			}
			if(object$PDFfit.model=="skewStudent") {
				curve(dsstd(x, mean=object$PDFfit$estimate["mean"], sd = object$PDFfit$estimate["sd"], nu=object$PDFfit$estimate["nu"], xi = object$PDFfit$estimate["xi"]), n=301, from=XX[1], to=XX[2], col="#CC6600", lwd=3, add=TRUE, xpd=TRUE)
				legend("right", legend="skew-Student density", text.col="#CC6600", bty="n")
			}
		}
				
		mean.ages <- rowMeans(ages)

		points(mean.ages, rep(0, length(mean.ages)), pch=17, xpd=TRUE, cex=1.2)

		points(quantile(rA, prob=c(0.5,0.95)), c(0,0), pch=c(17,18), col="#0078B3", xpd=TRUE, cex=1.2)
		
		legend("topleft", legend=c("fossil ages","MC distribution", "median","95% quantile"), pch=c(17,22, 17,18), col=c("black","#4DCCFF","#0078B3","#0078B3"), pt.bg=c("black",rgb(.25,.7,.9,0.5)), bty="n") # inset=c(0, -0.1), xpd=TRUE
	
}	