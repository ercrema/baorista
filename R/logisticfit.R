#' @title Fits a Logistic growth model on Aoristic data
#' @description Fits an exponential growth model to \code{ProbMat} calss objects.
#' @param x A ProbMat class object
#' @param niter Number of MCMC iterations. Default is 100,000.
#' @param nburnin Mumber of iterations discarded for burnin. Default is 50,000.
#' @param thin Thinning interval
#' @param nchains Number of MCMC chains
#' @param rPrior A string defining prior for the growth parameter r. Default is 'dexp(1/0.01)'. 
#' @param mPrior A string defining prior for the point of maximum growth rate m. Default is 'dunif(1,z)', where 'z' is the number of time-blocks. 
#' @param rSampler A list containing settings for the MCMC sampler for the parameter 'r'. Default is null and employs nimble's Default sampler (RW sampler).
#' @param mSampler A list containing settings for the MCMC sampler for the parameter 'm'. Default is null and employs nimble's Default sampler (RW sampler).
#' @param parallel Logical specifying whether the chains should be run in parallel or not.
#' @param seeds Random seed for each chain. Default is 1:4.
#' @details (To be completed)
#' @import nimble
#' @import coda
#' @import parallel
#' @export

logisticfit  <- function(x,niter=100000,nburnin=50000,thin=10,nchains=4,rPrior='dexp(1/0.001)',mPrior='dunif(1,z)',rSampler=NULL,mSampler=NULL,parallel=FALSE,seeds=1:4)
{
	require(nimble)
	require(coda)

	# Extract mid points
	mids <- apply(x$tblocks,1,median)

	# Initial Warnings
	if (nchains==1) {warning('Running MCMC on single chain')}

	# Define Data
	d  <- list(theta=x$pmat)
	# Define Constant
	constants  <- list()
	constants$z  <- ncol(d$theta) 

	if (!parallel)
 	{
	logisticmodel  <- nimbleCode({
		theta[,] ~ dAoristicLogisticGrowth_vector(r=r,z=z,m=m)
		r ~ rPrior
		m.raw ~ mPrior
		m  <- round(m.raw)
	})

	assign('dAoristicLogisticGrowth_vector',dAoristicLogisticGrowth_vector,envir=.GlobalEnv)


	logisticmodel <- gsub('rPrior', rPrior, deparse(logisticmodel)) 
	logisticmodel <- gsub('mPrior', mPrior, logisticmodel) |> parse(text=_)

	inits  <- vector('list',length=nchains)
	for (k in 1:nchains)
	{
		set.seed(seeds[k])
		inits[[k]]  <- list(r=rexp(1,1/0.01),m.raw=runif(1,1,constants$z))
	}
	print('Compiling nimble model...')
	suppressMessages(model  <- nimbleModel(logisticmodel,constants=constants,data=d,inits=inits[[1]]))
	assign('rAoristicLogisticGrowth_vector',rAoristicLogisticGrowth_vector,envir=.GlobalEnv)
	suppressMessages(cModel <- compileNimble(model))
	suppressMessages(conf <- configureMCMC(model))
	if (!is.null(rSampler))
	{
		suppressMessages(conf$removeSamplers('r'))
		# 	rSampler=list('sigma',type='slice')
		suppressMessages(do.call(conf$addSampler,rSampler))
	}

	if (!is.null(mSampler))
	{
		suppressMessages(conf$removeSamplers('m.raw'))
		suppressMessages(do.call(conf$addSampler,mSampler))
	}
	suppressMessages(conf$addMonitors('r'))
	suppressMessages(conf$addMonitors('m'))
	suppressMessages(MCMC <- buildMCMC(conf))
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
			dALog=nimbleFunction(
					     run = function(x = double(2),z=integer(0),r=double(0),m=integer(0), log = integer(0)) {
						     returnType(double(0))
						     t = 1:z
						     n = 1/(1+exp(-r*(t-m)))
						     p = n/sum(n)
						     pg = x %*% p
						     logProb = sum(log(pg))
						     if(log) {
							     return(logProb)
						     } else {
							     return(exp(logProb))
						     }
					     })   

			assign('dALog',dALog,envir=.GlobalEnv)

			logisticmodel  <- nimbleCode({
				theta[,] ~ dALog(r=r,z=z,m=m)
				r ~ rPrior
				m.raw ~ mPrior
				m  <- round(m.raw)
			})

			logisticmodel <- gsub('rPrior', rPrior, deparse(logisticmodel)) 
			logisticmodel <- gsub('mPrior', mPrior, logisticmodel) |> parse(text=_)

			set.seed(seeds)
			inits  <- list(r=rexp(1,1/0.01),m=runif(1,1,constants$z))
			model  <- nimbleModel(logisticmodel,constants=constants,data=d,inits=inits)
			assign('rALog',rALog,envir=.GlobalEnv)
			cModel <- compileNimble(model)
			conf <- configureMCMC(model)
			conf$addMonitors('r')
			conf$addMonitors('m')

			if (!is.null(rSampler))
			{
				suppressMessages(conf$removeSamplers('r'))
				# 	rSampler=list('sigma',type='slice')
				suppressMessages(do.call(conf$addSampler,rSampler))
			}

			if (!is.null(mSampler))
			{
				suppressMessages(conf$removeSamplers('m.raw'))
				suppressMessages(do.call(conf$addSampler,mSampler))
			}
			suppressMessages(conf$addMonitors('r'))
			suppressMessages(conf$addMonitors('m'))

			MCMC <- buildMCMC(conf)
			cMCMC <- compileNimble(MCMC)
			results <- runMCMC(cMCMC, niter = niter, thin=thin,nburnin = nburnin,samplesAsCodaMCMC = T,setSeed=seed)
		}

		ncores  <- nchains
		cl  <- makeCluster(ncores)
		out  <- parLapply(cl=cl,X=seeds,fun=runfun,d=d,constants=constants,nburnin=nburnin,niter=niter,thin=thin,rPrior,rSampler)
		stopCluster(cl)
		results <- out
	}

	diagnostic  <- gelman.diag(results,multivariate=FALSE) 
	if (any(diagnostic[[1]][,1]>1.01)){warning(paste0('Rhat value above 1.01 (',round(max(diagnostic[[1]][,1]),3),'). Consider rerunning the model with a higher number of iterations'))}
	posterior <- do.call(rbind.data.frame,results)
 	posterior.r  <- posterior[,'r']/x$resolution #scale to match resolution
	posterior.m <- mids[posterior[,'m']]
	posterior.m.index <- posterior[,'m']


	results  <- list(x=x,posterior.r=posterior.r,posterior.m=posterior.m,posterior.m.index=posterior.m.index,diagnostic=diagnostic)
	class(results)  <- c('fittedLogistic',class(results))
	return(results)
}
