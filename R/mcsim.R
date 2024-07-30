#' @title Monte-Carlo simulation on aoristic data
#' @description Samples multiple sets of random dates from aoristic weights
#' @param x A ProbMat class object
#' @param nsim Number of Monte-Carlo simulations
#' @details  The function randomly assigns to each event a time-block based on its probability values (i.e. aoristic weight) and computes, for each time-block, the total number of simulated events. This process is repeated \code{nsim} time, allowing to estimate percentile-based intervals on the number of events per time-block (Crema 2012). It should be noted that while this approach accounts for chronological uncertainty, it provides only a description of the sample rather than the underlying population, and can be biased how the underlying archaeological periodisations define the time-spans of each event (see also Crema 2024 for discussion on limitations).
#' @return  An object of class \code{mcsimres} containing relevant metadata and a matrix with the number of events per time-block per Monte-Carlo simulation. 
#' @references 
#
#' Crema, E. R. (2012). Modelling Temporal Uncertainty in Archaeological Analysis. Journal of Archaeological Method and Theory, 19(3), 440–461. doi:10.1007/s10816-011-9122-3
#' Crema, E.R. (2024). A Bayesian alternative to Aoristic analyses in archaeology. Archaeometry. doi:10.1111/arcm.12984
#' 
#' @importFrom stats rmultinom
#' @export

mcsim  <- function(x,nsim=1000)
{
	sims <- sapply(1:x$n,function(x,p,n){rmultinom(n=n,size=1,prob=p[x,])},p=x$pmat,n=nsim,simplify='array')
	sums  <- apply(sims,1:2,sum)
	results <- list(z=x$z,sums=sums,nsim=nsim,tblocks=x$tblocks,n=x$n,resolution=x$resolution,timeRange=x$timeRange)
	class(results) <- c('mcsimres',class(results))
	return(results)
}
