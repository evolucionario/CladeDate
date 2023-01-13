#' @title Utility Functions for Dating Clades Based on the Fossil Record

#' @description
#' Point estimate (\code{pdate}), quantile function (\code{qdate}), and random number generator (\code{rdate}) used for estimating the age of a clade and its associated uncertainty.

#' @param ages vector of fossil ages.
#' @param CI vector of two number specifying the confidence interval (degault to the 0 and 95th percentil).
#' @param method a character string specifying the method: \code{"StraussSadler"} (default), \code{"RobsonWhitlock"}, \code{"Beta"}, \code{"NorrisPenGap"}, \code{"NorrisGhostLin"}, or \code{"OLE"}. See Details.
#' @param KStest if TRUE, a Kolmogorov-Smirnov test evaluates the null hypothesis of uniformity of fossil ages.
#' @param k the number of fossil ages to use in the “OLE” method, which only used the k oldest ages. 
#' @param p vector of probabilities.
#' @param n number of samples requested
#' @param ... arguments passed to \code{qdate}

#' @details
#' The "StraussSadler" (Strauss and Sadler 1989) method (default) assumes that fossil ages are uniformly distributed, and a warning is returned if a Kolmogorov-Smirnov test rejects the uniformity hypothesis. The “Beta” method is a different parameterization of the Strauss and Sadler's method (Wang et al. 2009) that uses the \code{qbeta()} function, since the ratio between the observed maximum fossil age and the clade age is Beta-distributed with parameters N and 1 (Wang et al. 2009). The "RobsonWhitlock" method (Robson and Whitlock 1964, Solow 2003) does not assume uniformity but is based on the two oldest ages only. The “OLE” method (for Optimal Linear Estimator) does not assume uniformity and uses a weighted sum of the oldest k ages (Cooke 1980, Roberts and Solow 2003, Solow 2005). "NorrisPenGap" and "NorrisGhostLin" implement the methods of Norris et al.(2015) based on the two oldest fossils only and the log-logistic distribution. "NorrisPenGap" is used when the precise phylogenetic placement of fossils is not known (with results identical to the "RobsonWhitlock" method), whereas "NorrisGhostLin" is used when one fossil from each daughter lineage is used.

#' @return \code{pdate} returns a point estimate, lower and upper bounds (based on the specified quantiles), and results of the Kolmogorov-Smirnov test (optional), \code{qdate} returns quantiles, and \code{rdate} returns a random sample of \code{n} values.

#' @import stats

#' @examples
#' #Point estimates of clade age
#'   pdate(ages=c(54, 30, 25, 14, 5), p=c(0.1, 0.5, 0.9), KStest=TRUE)
#'   pdate(ages=c(54, 30, 25, 14, 5), p=c(0.1, 0.5, 0.9), method="Beta")
#'   pdate(ages=c(54, 30, 25, 14, 5), p=c(0.1, 0.5, 0.9), method="RobsonWhitlock")

#' #Quantiles
#'   qdate(ages=c(54, 30, 25, 14, 5), p=c(0.025, 0.5, 0.975))

#' #Random numbers
#'   rdate(ages=c(54, 30, 25, 14, 5), n=10)

#' @author Santiago Claramunt, \email{s.claramunt@@utoronto.ca}

#' @references

#' Claramunt, S. (2022). CladeDate: calibration information generator for divergence time estimation. Methods in Ecology and Evolution **13**(11):2331-2338. <https://doi.org/10.1111/2041-210X.13977>

#' Cooke, P. (1980). Optimal linear estimation of bounds of random variables. *Biometrika*, **67**, 257--258.

#' Norris, R. W., Strope, C. L., McCandlish, D. M. and Stoltzfus, A. (2015). Bayesian priors for tree calibration: Evaluating two new approaches based on fossil intervals. *bioRxiv*, 014340 doi:<https://doi.org/10.1101/014340>

#' Roberts, D. L. and Solow, A. R. (2003). When did the dodo become extinct?. *Nature*, **426**(6964), 245--245.

#' Robson, D. S., and Whitlock, J. H. (1964). Estimation of a truncation point. Biometrika **51(1)**, 33--39.

#' Solow, A. R. (2003) Estimation of stratigraphic ranges when fossil finds are not randomly distributed. *Paleobiology* **29**(2):181--185.

#' Solow, A. R. (2005) Inferring extinction from a sighting record. Mathematical *Biosciences* **195**:47-–55.

#' Strauss D. and Sadler, P. M.  (1989) Classical confidence intervals and Bayesian probability estimates for ends of local taxon ranges. *Mathematica Geology* **21**,411--427.

#' Wang, S. C., Chudzicki, D. J. and Everson, P. J. (2009) Optimal estimators of the position of a mass extinction when recovery potential is uniform. *Paleobiology* **35**(3), 447–-459.


#' @export

pdate <- function(ages, method="StraussSadler", CI=c(0,0.95), KStest=FALSE, k=min(length(ages),10)) {
 	
 	# Sort ages, just in case
 	
 	ages <- sort(ages, decreasing=TRUE)

	if(method=="StraussSadler") {		

		n <- length(ages)
		
		range <- max(ages)-min(ages)
				
		AGE <- max(ages) + range/n
		
		if(!is.null(CI)) { ci <- range*(1-CI)^-(1/n) + min(ages) }
		
		# Test for uniform distribution
		if(KStest && n > 2) {
					
			# add a small random variance to the inner values to avoid ties		
			ages[2:(n-1)] <- ages[2:(n-1)] + rnorm(n-2,0,AGE/100000)
			
			# Estimate the minimum bound under the null hypothesis
			 Min <- min(ages)-range/n
							
			Pks <- ks.test(ages, 'punif', min=max(c(0, Min)), max=AGE)$p.value
			
			if(Pks >= 0.05) cat("Uniform distribution not rejected (Kolmogorov-Smirnov P = ", signif(Pks, digits=2),")\n", sep="")
			
			if(Pks < 0.05) warning("The data may not be uniformly distributed (Kolmogorov-Smirnov P = ", signif(Pks, digits=2), ")", call.=FALSE)
			}
		
		}

	if(method=="Beta") {
						
		n <- length(ages)

		range <- max(ages) - min(ages)

		X <- range/qbeta(p=0.5, shape1=n, shape2=1, lower.tail=FALSE)

		AGE <- X + min(ages)

		if(!is.null(CI)) { ci <- range/qbeta(p=CI, shape1=n, shape2=1, lower.tail=FALSE) + min(ages) }

		}

	if(method=="RobsonWhitlock") {
						
		AGE <- ages[1] + (ages[1]-ages[2])
		
		if(!is.null(CI)) { ci <- ages[1] + (ages[1]-ages[2])*CI/(1-CI) }

		}

	if(method=="OLE") { 
	
	# Select the k oldest observations (5 by default)
	
		ages <- ages[1:k]

	# Estimate the shape paremter of the joint Weibull distribution
			
		t <- numeric()
	
		for(i in 1:(k-2)) { t[i] <- (ages[1] - ages[k]) / (ages[1] - ages[i+1]) }

		v <- 1/(k-1) * sum(log(t))

	# Construct the e vector
	
		e <- rep.int(1,k)
	
	# Obtain the matrix Lamda
	
		Lambda <- matrix(ncol=k, nrow=k)
		for(i in 1:k) {
			for(j in 1:k) {
				if(j>i) next	# this ensures j =< i and fills the lower triange		
				Lambda[i,j] <- ( gamma(2*v+i) * gamma(v+j)) / (gamma(v+i)*gamma(j)) 
			}
		}

		# Fill the upper triange
		Lambda[upper.tri(Lambda)] <- t(Lambda)[upper.tri(Lambda)]
	
	# Compute the vector of weights
	a <- as.numeric(solve( t(e) %*% solve(Lambda) %*% e)) * solve(Lambda) %*%e
	
	# Compute the point estimate using the weights
	AGE <- sum(a*ages)
	
	if(!is.null(CI)) {
		# Calculate Su (Robert and Solow 2003) or C.alpha (Solow 2005 eq. 18)
		Su <- (-log(1-CI)/k)^-v
		ci <- max(ages) + ( max(ages) - min(ages)) / ( Su - 1 )
		}
	}
	
	if(method=="NorrisPenGap") {

		PenG <- ages[1] - ages[2]
		AGE <- exp(qlogis(p=0.5, location=log(PenG)))  + ages[1]
				
		if(!is.null(CI)) { ci <- exp(qlogis(p=CI, location=log(PenG))) + ages[1] }
	}
	
	if(method=="NorrisGhostLin") {
		
		GLin <- ages[1] - ages[2] 
		AGE <- exp(qlogis(p=0.5, location=log(GLin/2))) + ages[1]
				
		if(!is.null(CI)) { ci <- exp(qlogis(p=CI, location=log(GLin/2))) + ages[1] }
	}
	
	return(c(estimate=AGE, lower=ci[1], upper=ci[2]))
}



#' @rdname pdate
#' @export

qdate <- function(ages, p=0.5, method="StraussSadler", k=min(length(ages),10)) {
 	
 	if(length(ages) < 2) stop("More than one fossil age is needed for inference")
 	 	
 	ages <- sort(ages, decreasing=TRUE)

	if(method=="StraussSadler") {

		n <- length(ages)
		
		range <- max(ages)-min(ages)
		
		AGE <- range*(1-p)^-(1/n) + min(ages)
		}
		
	if(method=="Beta") {
		
		n <- length(ages)

		range <- max(ages) - min(ages)

		X <- range/qbeta(p=p, shape1=n, shape2=1, lower.tail=FALSE)

		AGE <- X + min(ages)
		
		 }
		
	if(method=="RobsonWhitlock") {
				
		AGE <- ages[1] + (ages[1]-ages[2])*p/(1-p) 

		}

	if(method=="NorrisPenGap") {	

		PenG <- ages[1] - ages[2]

		UltG <- exp(qlogis(p=p, location=log(PenG)))
		
		AGE <- UltG + ages[1] }
	
	if(method=="NorrisGhostLin") {
		
		GLin <- ages[1] - ages[2] 

		UltG <- exp(qlogis(p=p, location=log(GLin/2)))
		
		AGE <- UltG + ages[1]}
	
	if(method=="OLE") {
		
	# Select the k oldest observations (5 by default)
	
		ages <- ages[1:k]

	# Estimate the shape paremter of the joint Weibull distribution
			
		t <- numeric()
	
		for(i in 1:(k-2)) { t[i] <- (ages[1] - ages[k]) / (ages[1] - ages[i+1]) }

		v <- 1/(k-1) * sum(log(t))

	# Calculate Su (Robert and Solow 2003, = c(alpha) Solow 2005 eq. 18)
		Su <- (-log(1-p)/k)^-v

	# But Su must be >= 1 (otherwise the estimate may be younger than older fossil) so:
		Su <- pmax(Su, 1) # use pmax instead of max to allow for vectors of p and Su
	
	# Calculate the quantile (Solow 2005 eq. 17) (Something seems wrong)	
	#AGE <- ( ages[1] - Su*ages[k]) / ( 1 - Su)
	
	# Calculate the quantile (Roberts and Solow 2003)
	
		AGE <- ages[1] + ( ages[1] - ages[k]) / ( Su - 1 )

	}
	
	return(AGE)
}


#' @rdname pdate
#' @export

rdate <- function(ages, n, ...) {

	ages <- as.matrix(ages)

	if(ncol(ages) > 2) stop("'ages' must be a vector or a two-column matrix")

	# When fossil ages are known exactly
	if(ncol(ages) == 1) {
		
		rA <- qdate(p=runif(n), ages=ages[,1], ...)
	}

	# When fossil ages are known from time intervals (max-min)
	else if(ncol(ages) == 2) {
		
		# Check that bounds are correct for all fossils
		if(any(apply(ages, 1, diff) < 0)) {
			stop("Detected at least one fossil with older bound < younger bound") }
	
		rA <- numeric(length=n)
	
		for(i in 1:n) {
	
		A <- runif(nrow(ages), min=ages[,1], max=ages[,2])
	
		rA[i] <- qdate(p=runif(1), ages=A, ...)
		
		}
	}
	
	return(rA)	
}
