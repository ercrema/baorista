#' @title Fits a Logistic growth model on Aoristic data
#' @description Fits an exponential growth model to \code{ProbMat} class objects.
#' @param x A ProbMat class object
#' @param niter Number of MCMC iterations. Default is 100,000.
#' @param nburnin Number of iterations discarded for burn-in. Default is 50,000.
#' @param thin Thinning interval
#' @param nchains Number of MCMC chains
#' @param rPrior A string defining prior for the growth parameter r. Default is 'dexp(1/0.01)'. 
#' @param mPrior A string defining prior for the point of maximum growth rate m. Default is 'dunif(1,z)', where 'z' is the number of time-blocks. 
#' @param rSampler A list containing settings for the MCMC sampler for the parameter 'r'. Default is null and employs nimble's Default sampler (RW sampler).
#' @param mSampler A list containing settings for the MCMC sampler for the parameter 'm'. Default is null and employs nimble's Default sampler (RW sampler).
#' @param parallel Logical specifying whether the chains should be run in parallel or not.
#' @param seeds Random seed for each chain. Default is 1:4.
#' @details The function fits a discrete bounded logistic growth model on the observed data using MCMC as implemented by the nimble package. The Bayesian model consists of two parameters, a growth rate (r) and a midpoint (m) defining the inflection point of the growth curve. Priors of the two parameters can be defined by the arguments \code{rPrior} and \code{mPrior}. In the latter case the object \code{z} is a placeholder for the number of blocks (e.g. the default 'dunif(1,z)` is a uniform across all blocks). Priors are defined by character strings following the syntax used by nimble. The distribution parameters defined in \code{rPrior} and \code{mPrior} are also used to generate initialisation values for the MCMC. Please note that the function returns posterior of the growth rate normalised by the resolution defined in the \code{ProbMat} class object.  MCMC settings such as the choice the sampler, number of iterations, chains, etc can also be specified.  
#' @return A \code{fittedLogistic} class object containing the original ProbMat class object, posteriors of the growth rate and midpoint and their MCMC diagnostics (i.e. Gelman Rubin statistic and effective sample sizes).
#' @import nimble
#' @import coda
#' @import parallel
#' @importFrom stats runif rexp rchisq rt rlnorm rweibull rnorm rgamma rlogis rbeta rt rlnorm
#' @export

logisticfit  <- function(x,niter=100000,nburnin=50000,thin=10,nchains=4,rPrior='dexp(rate=1/0.001)',mPrior='dunif(min=1,max=z)',rSampler=NULL,mSampler=NULL,parallel=FALSE,seeds=1:4)
{
	#Handle cleaning of GlobalEnv on exit
	envobj <- ls(envir=.GlobalEnv)
	on.exit(rm(list=ls(envir=.GlobalEnv)[which(!ls(envir=.GlobalEnv)%in%envobj)],envir=.GlobalEnv))
	#Addresses R CMD Check NOTES
	returnType <- m.raw <- nimStop <- nimMatrix <- rAoristicLogisticGrowth_vector <- runfun <- dALog <- rALog <-  NULL
	# Extract mid points
	mids <- apply(x$tblocks,1,median)

	# Initial Warnings
	if (nchains==1) {warning('Running MCMC on single chain')}
	# Check prior definitions
	supported.distributions  <- data.frame(d=c('dnorm','dexp','dbeta','dchisq','dgamma','dlogis','dunif','dt','dweib','dlnorm'),r=c('rnorm','rexp','rbeta','rchisq','rgamma','rlogis','runif','rt','rweibull','rlnorm'))
	if(!sub("\\(.*","",rPrior)%in%supported.distributions$d | !sub("\\(.*","",mPrior)%in%supported.distributions$d){stop(paste0('Unsupported distribution in prior definition. Please use one of the following: ',paste(supported.distributions$d,collapse=', ')))}
	rPrior.rand1  <- supported.distributions$r[supported.distributions$d==strsplit(rPrior,"\\(")[[1]][1]]
	rPrior.rand2 <- gsub("\\)","",strsplit(rPrior,"\\(")[[1]][2])
	rPrior.rand <- paste0(rPrior.rand1,"(n=1,",rPrior.rand2,")")

	mPrior.rand1  <- supported.distributions$r[supported.distributions$d==strsplit(mPrior,"\\(")[[1]][1]]
	mPrior.rand2 <- gsub("\\)","",strsplit(mPrior,"\\(")[[1]][2])
	mPrior.rand2 <- gsub('z','constants$z',mPrior.rand2)
	mPrior.rand <- paste0(mPrior.rand1,"(n=1,",mPrior.rand2,")")

	# Define Data
	d  <- list(theta=x$pmat)
	# Define Constant
	constants  <- list()
	constants$z  <- ncol(d$theta) 

	if (!parallel)
 	{

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

		rALog=nimbleFunction(
				     run = function(n=integer(0),z=integer(0),r=double(0),m=integer(0)) {
					     returnType(double(2))
					     nimStop("user-defined distribution dAExp provided without random generation function.")
					     x <- nimMatrix()
					     return(x)
				     })
		pos <- 1
		assign('dALog',dALog,envir=as.environment(pos))
		assign('rALog',rALog,envir=as.environment(pos))

		logisticmodel  <- nimbleCode({
			theta[,] ~ dALog(r=r,z=z,m=m)
			r ~ rPrior
			m.raw ~ mPrior
			m  <- round(m.raw)
		})


		logisticmodel <- gsub('rPrior', rPrior, deparse(logisticmodel)) 
		logisticmodel <- gsub('mPrior', mPrior, logisticmodel) |> parse(text=_)

		inits  <- vector('list',length=nchains)
		for (k in 1:nchains)
		{
			set.seed(seeds[k])
			inits[[k]]  <- list(r=eval(parse(text=rPrior.rand)),m.raw=eval(parse(text=mPrior.rand)))
		}
		message('Compiling nimble model...')
		suppressMessages(model  <- nimbleModel(logisticmodel,constants=constants,data=d,inits=inits[[1]]))
		suppressMessages(cModel <- compileNimble(model))
		suppressMessages(conf <- configureMCMC(model))
		if (!is.null(rSampler))
		{
			suppressMessages(conf$removeSamplers('r'))
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
		rm(dALog,rALog,envir=.GlobalEnv) #clean temporary objects from GlobalEnv

 	}

	if (parallel)
	{
		message('Running in parallel - progress bar will no be visualised')
		runfun  <- function(seed,constants,d,niter,thin,nburnin,rPrior,rPrior.rand,rSampler,mPrior,mPrior.rand,mSampler)
		{
			#Addresses R CMD Check NOTES
			returnType <- m.raw <- nimStop <- nimMatrix <-  dALog <- rALog <- NULL
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
			rALog=nimbleFunction(
					     run = function(n=integer(0),z=integer(0),r=double(0),m=integer(0)) {
						     returnType(double(2))
						     nimStop("user-defined distribution dAExp provided without random generation function.")
						     x <- nimMatrix()
						     return(x)
					     })
			pos <- 1
			assign('dALog',dALog,envir=as.environment(pos))
			assign('rALog',rALog,envir=as.environment(pos))

			logisticmodel  <- nimbleCode({
				theta[,] ~ dALog(r=r,z=z,m=m)
				r ~ rPrior
				m.raw ~ mPrior
				m  <- round(m.raw)
			})

			logisticmodel <- gsub('rPrior', rPrior, deparse(logisticmodel)) 
			logisticmodel <- gsub('mPrior', mPrior, logisticmodel) |> parse(text=_)

			set.seed(seed)
			inits  <- list(r=eval(parse(text=rPrior.rand)),m.raw=eval(parse(text=mPrior.rand)))
			model  <- nimbleModel(logisticmodel,constants=constants,data=d,inits=inits)
			cModel <- compileNimble(model)
			conf <- configureMCMC(model)
			conf$addMonitors('r')
			conf$addMonitors('m')

			if (!is.null(rSampler))
			{
				suppressMessages(conf$removeSamplers('r'))
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
		clusterEvalQ(cl,{library(nimble)})
		out  <- parLapply(cl=cl,X=seeds,fun=runfun,d=d,constants=constants,nburnin=nburnin,niter=niter,thin=thin,rPrior=rPrior,rPrior.rand=rPrior.rand,rSampler=rSampler,mPrior=mPrior,mPrior.rand=mPrior.rand,mSampler=mSampler)
		stopCluster(cl)
		results <- out
	}

	rhat  <- gelman.diag(results,multivariate=FALSE) 
	ess  <- effectiveSize(results)
	if (any(rhat[[1]][,1]>1.01)){warning(paste0('Rhat value above 1.01 (',round(max(rhat[[1]][,1]),3),'). Consider rerunning the model with a higher number of iterations'))}
	posterior <- do.call(rbind.data.frame,results)
 	posterior.r  <- posterior[,'r']/x$resolution #scale to match resolution
	posterior.m <- mids[posterior[,'m']]
	posterior.m.index <- posterior[,'m']


	results  <- list(x=x,posterior.r=posterior.r,posterior.m=posterior.m,posterior.m.index=posterior.m.index,rhat=rhat,ess=ess)
	class(results)  <- c('fittedLogistic',class(results))
	return(results)
}
