#' @title Estimate Exponential Growth rate from  Aoristic data
#' @description Fits an exponential growth model to \code{ProbMat} class objects.
#' @param x A ProbMat class object
#' @param niter Number of MCMC iterations. Default is 500,000.
#' @param nburnin Number of iterations discarded for burn-in. Default is 250,000.
#' @param thin Thinning interval
#' @param nchains Number of MCMC chains
#' @param rPrior A string defining prior for the growth parameter r. Default is 'dnorm(mean=0,sd=0.05)'. 
#' @param rSampler A list containing settings for the MCMC sampler. Default is null and employs nimble's Default sampler (RW sampler).
#' @param parallel Logical specifying whether the chains should be run in parallel or not.
#' @param seeds Random seed for each chain. Default is 1:4.
#' @details The function fits a discrete bounded exponential growth model on the observed data using MCMC as implemented by the nimble package. The Bayesian model consists of a single growth rate parameter (r), and users can define suitable priors using character strings for the argument \code{rPrior} (for details on how this should be specified please consult the nimble manual). Please note that the function returns posterior of the growth rate normalised by the resolution defined in the \code{ProbMat} class object.  MCMC settings such as the choice the sampler, number of iterations, chains, etc can also be specified.  
#' @return A \code{fittedExp} class object containing the original ProbMat class object, posterior of the growth rate, along with its Gelman Rubin statistic and effective sample sizes. 
#' @import nimble
#' @import coda
#' @import parallel
#' @importFrom stats rnorm
#' @export


expfit  <- function(x,niter=100000,nburnin=50000,thin=10,nchains=4,rPrior='dnorm(mean=0,sd=0.05)',rSampler=NULL,parallel=FALSE,seeds=1:4)
{
	#Handle cleaning of GlobalEnv on exit
	envobj <- ls(envir=.GlobalEnv)
	on.exit(rm(list=ls(envir=.GlobalEnv)[which(!ls(envir=.GlobalEnv)%in%envobj)],envir=.GlobalEnv))
	#Addresses R CMD Check NOTES
	returnType <- m.raw <- nimStop <- nimMatrix <-  dAExp <- rAExp <- runfun <-  NULL
	# Initial Warnings
	if (nchains==1) {warning('Running MCMC on single chain')}

	# Define Data
	d  <- list(theta=x$pmat)
	# Define Constant
	constants  <- list()
	constants$n.tblocks  <- ncol(d$theta) 

	if (!parallel)
 	{
		
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
		rAExp=nimbleFunction(
				     run = function(n=integer(0),z=integer(0),r=double(0)) {
					     returnType(double(2))
					     nimStop("user-defined distribution dAExp provided without random generation function.")
					     x <- nimMatrix()
					     return(x)
				     })
		pos <- 1
		assign('dAExp',dAExp,envir=as.environment(pos))
		assign('rAExp',rAExp,envir=as.environment(pos))
		

		expmodel  <- nimbleCode({
			theta[,] ~ dAExp(r=r,z=n.tblocks)
			r ~ dnorm(mean=0,sd=0.05)
		})

		expmodel <- gsub('dnorm\\(mean=0,sd=0.05\\)', rPrior, deparse(expmodel)) |> parse(text=_)

		inits  <- vector('list',length=nchains)
		for (k in 1:nchains)
		{
			set.seed(seeds[k])
			inits[[k]]  <- list(r=rnorm(1,0,0.05))
		}
		message('Compiling nimble model...')
		suppressMessages(model  <- nimbleModel(expmodel,constants=constants,data=d,inits=inits[[1]]))
		suppressMessages(cModel <- compileNimble(model))
		suppressMessages(conf <- configureMCMC(model))
		if (!is.null(rSampler))
		{
			suppressMessages(conf$removeSamplers('r'))
			# 	rSampler=list('sigma',type='slice')
			suppressMessages(do.call(conf$addSampler,rSampler))
		}
		suppressMessages(conf$addMonitors('r'))
		suppressMessages(MCMC <- buildMCMC(conf))
		suppressMessages(cMCMC <- compileNimble(MCMC))
		results <- runMCMC(cMCMC, niter = niter, thin=thin,nburnin = nburnin,inits=inits,samplesAsCodaMCMC = T,nchains=nchains,progressBar=TRUE,setSeed=seeds)
 	}

	if (parallel)
	{
		message('Running in parallel - progress bar will no be visualised')
		runfun  <- function(seed,constants,d,niter,thin,nburnin,rPrior,rSampler)
		{

			returnType <- nimStop <- nimMatrix <- dAExp <- rAExp <-  NULL
# 			require(nimble)
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

			rAExp=nimbleFunction(
					     run = function(n=integer(0),z=integer(0),r=double(0)) {
						     returnType(double(2))
						     nimStop("user-defined distribution dAExp provided without random generation function.")
						     x <- nimMatrix()
						     return(x)
					     })
			pos <- 1
			assign('dAExp',dAExp,envir=as.environment(pos))
			assign('rAExp',rAExp,envir=as.environment(pos))

			expmodel  <- nimbleCode({
				theta[,] ~ dAExp(r=r,z=n.tblocks)
				r ~ dnorm(mean=0,sd=0.05)
			})

			expmodel <- gsub('dnorm\\(mean=0,sd=0.05\\)', rPrior, deparse(expmodel)) |> parse(text=_)
			set.seed(seed)
			inits  <- list(r=rnorm(1,0,0.05))
			model  <- nimbleModel(expmodel,constants=constants,data=d,inits=inits)
			cModel <- compileNimble(model)
			conf <- configureMCMC(model)
			conf$addMonitors('r')
			if (!is.null(rSampler))
			{
				suppressMessages(conf$removeSamplers('sigma'))
				do.call(conf$addSampler,rSampler)
			}
			MCMC <- buildMCMC(conf)
			cMCMC <- compileNimble(MCMC)
			results <- runMCMC(cMCMC, niter = niter, thin=thin,nburnin = nburnin,samplesAsCodaMCMC = T,setSeed=seed)
		}

		ncores  <- nchains
		cl  <- makeCluster(ncores)
		clusterEvalQ(cl,{library(nimble)})
		out  <- parLapply(cl=cl,X=seeds,fun=runfun,d=d,constants=constants,nburnin=nburnin,niter=niter,thin=thin,rPrior=rPrior,rSampler=rSampler)
		stopCluster(cl)
		results <- out
	}
	rhat  <- gelman.diag(results,multivariate=FALSE) 
	ess  <- effectiveSize(results)
	if (any(rhat[[1]][,1]>1.01)){warning(paste0('Rhat value above 1.01 (',round(max(rhat[[1]][,1]),3),'). Consider rerunning the model with a higher number of iterations'))}
	posterior.r <- do.call(rbind.data.frame,results)
 	posterior.r  <- posterior.r/x$resolution #scale to match resolution
	results  <- list(x=x,posterior.r=posterior.r,rhat=rhat,ess=ess)
	class(results)  <- c('fittedExp',class(results))
	
	envobj.current <- ls(envir=.GlobalEnv)
	toremove  <- envobj.current[which(!envobj.current%in%envobj)]
	rm(list=toremove,envir=.GlobalEnv)
	return(results)
}


