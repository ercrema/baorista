#' @title Estimate Exponential Growth rate from  Aoristic data
#' @description Fits an exponential growth model to \code{ProbMat} calss objects.
#' @param x A ProbMat class object
#' @param niter Number of MCMC iterations. Default is 500,000.
#' @param nburnin Mumber of iterations discarded for burnin. Default is 250,000.
#' @param thin Thinning interval
#' @param nchains Number of MCMC chains
#' @param rPrior A string defining prior for the growth parameter r. Default is 'dnorm(mean=0,sd=0.05)'. 
#' @param rSampler A list containing settings for the MCMC sampler. Default is null and employs nimble's Default sampler (RW sampler).
#' @param parallel Logical specifying whether the chains should be run in parallel or not.
#' @param seeds Random seed for each chain. Default is 1:4.
#' @details (To be completed)
#' @import nimble
#' @import coda
#' @import parallel
#' @export


expfit  <- function(x,niter=100000,nburnin=50000,thin=10,nchains=4,rPrior='dnorm(mean=0,sd=0.05)',rSampler=NULL,parallel=FALSE,seeds=1:4)
{
	require(nimble)
	require(coda)

	# Initial Warnings
	if (nchains==1) {warning('Running MCMC on single chain')}

	# Define Data
	d  <- list(theta=x$pmat)
	# Define Constant
	constants  <- list()
	constants$n.tblocks  <- ncol(d$theta) 
	constants$weights <- icar.struct(constants$n.tblocks)$weight
	constants$num <- icar.struct(constants$n.tblocks)$num
	constants$adj <- icar.struct(constants$n.tblocks)$adj
	constants$L  <- length(constants$adj)

	if (!parallel)
 	{
	expmodel  <- nimbleCode({
		theta[,] ~ dAoristicExponentialGrowth_vector(r=r,z=n.tblocks)
		r ~ dnorm(mean=0,sd=0.05)

	})

	expmodel <- gsub('dnorm\\(mean=0,sd=0.05\\)', rPrior, deparse(expmodel)) |> parse(text=_)

	inits  <- vector('list',length=nchains)
	for (k in 1:nchains)
	{
		set.seed(seeds[k])
		inits[[k]]  <- list(r=rnorm(1,0,0.05))
	}
	print('Compiling nimble model...')
	suppressMessages(model  <- nimbleModel(expmodel,constants=constants,data=d,inits=inits[[1]]))
	suppressMessages(cModel <- compileNimble(model))
	suppressMessages(conf <- configureMCMC(model))
	if (!is.null(rSampler))
	{
		suppressMessages(conf$removeSamplers('r'))
		# 	rSampler=list('sigma',type='slice')
		do.call(conf$addSampler,rSampler)
	}
	suppressMessages(conf$addMonitors('r'))
	MCMC <- buildMCMC(conf)
	suppressMessages(cMCMC <- compileNimble(MCMC))
	results <- runMCMC(cMCMC, niter = niter, thin=thin,nburnin = nburnin,inits=inits,samplesAsCodaMCMC = T,nchains=nchains,progressBar=TRUE,setSeed=seeds)

 	}

	if (parallel)
	{
		require(parallel)
		print('Running in parallel - progress bar will no be visualised')
		runfun  <- function(seed,constants,d,niter,thin,nburnin,rPrior,rSampler)
		{
			require(nimble)
			dAExp=nimbleFunction(
					     run = function(x = double(2),z=integer(0),r=double(0), log = integer(0)) {
						     returnType(double(0))
						     t = 1:z
						     n = numeric(z)
						     for (i in 1:z)
						     {
							     n[i] = (1+r)^t[i]
						     }
						     p = n/sum(n)
						     pg = x %*% p
						     logProb = sum(log(pg))
						     if(log) {
							     return(logProb)
						     } else {
							     return(exp(logProb))
						     }
					     })   
			assign('dAExp',dAExp,envir=.GlobalEnv)

			expmodel  <- nimbleCode({
				theta[,] ~ dAExp(r=r,z=n.tblocks)
				r ~ dnorm(mean=0,sd=0.05)
			})

			expmodel <- gsub('dnorm\\(mean=0,sd=0.05\\)', rPrior, deparse(expmodel)) |> parse(text=_)
			set.seed(seeds)
			inits  <- list(r=rnorm(1,0,0.05))
			model  <- nimbleModel(expmodel,constants=constants,data=d,inits=inits)
			assign('rAExp',rAExp,envir=.GlobalEnv)
			cModel <- compileNimble(model)
			conf <- configureMCMC(model)
			conf$addMonitors('r')
			if (!is.null(rSampler))
			{
				suppressMessages(conf$removeSamplers('sigma'))
				# 	rSampler=list('sigma',type='slice')
				do.call(conf$addSampler,rSampler)
			}
			MCMC <- buildMCMC(conf)
			cMCMC <- compileNimble(MCMC)
			results <- runMCMC(cMCMC, niter = niter, thin=thin,nburnin = nburnin,samplesAsCodaMCMC = T,setSeed=seeds)
		}

		ncores  <- nchains
		cl  <- makeCluster(ncores)
		out  <- parLapply(cl=cl,X=seeds,fun=runfun,d=d,constants=constants,nburnin=nburnin,niter=niter,thin=thin,rPrior,rSampler)
		stopCluster(cl)
		results <- out
	}
	diagnostic  <- gelman.diag(results,multivariate=FALSE) 
	if (any(diagnostic[[1]][,1]>1.01)){warning(paste0('Rhat value above 1.01 (',round(max(diagnostic[[1]][,1]),3),'). Consider rerunning the model with a higher number of iterations'))}
	posterior.r <- do.call(rbind.data.frame,results)
 	posterior.r  <- posterior.r/x$resolution #scale to match resolution
	results  <- list(x=x,posterior.r=posterior.r,diagnostic=diagnostic)
	class(results)  <- c('fittedExp',class(results))
	return(results)
}
