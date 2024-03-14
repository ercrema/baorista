#' @method summary fittedExp
#' @importFrom stats median
#' @import coda
#' @export
summary.fittedExp <- function(object,...)
{
	out <- data.frame('r',median(object$posterior.r$r),paste0(HPDinterval(mcmc(object$posterior.r)),collapse='~'),object$rhat[[1]][c('r'),1],object$ess)
	names(out) <- c('Parameter','Posterior Median','95% HPDI','Rhat','ESS')
	return(out)
}

#' @method summary fittedLogistic
#' @importFrom stats median
#' @import coda
#' @export
summary.fittedLogistic <- function(object,...)
{
	out <- data.frame(c('r','m'),c(median(object$posterior.r$r),median(object$posterior.m$m)),
			  c(paste0(HPDinterval(mcmc(object$posterior.r)),collapse='~'),paste0(HPDinterval(mcmc(object$posterior.m)),collapse='~')),
			  object$rhat[[1]][c('r','m'),1],
			  object$ess)
	names(out) <- c('Parameter','Posterior Median','95% HPDI','Rhat','ESS')
	return(out)
}


#' @method summary fittedICAR
#' @importFrom stats median
#' @import coda
#' @export
summary.fittedICAR <- function(object,...)
{
	out <- data.frame(row.names(object$rhat[[1]]),
			  c(apply(object$posterior.p,2,median),median(object$posterior.sigma)),
			  c(apply(object$posterior.p,2,function(x){paste0(HPDinterval(mcmc(x)),collapse='~')}),paste0(HPDinterval(mcmc(object$posterior.sigma)),collapse='~')),
			  object$rhat[[1]][,1],
			  object$ess)
	names(out) <- c('Parameter','Posterior Median','95% HPDI','Rhat','ESS')
	return(out)
}
