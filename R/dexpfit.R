#' @title Estimate Exponential Growth rate from  Aoristic data
#' @description Fits a double exponential growth model to \code{ProbMat} class objects.
#' @param x A ProbMat class object
#' @param niter Number of MCMC iterations. Default is 500,000.
#' @param nburnin Number of iterations discarded for burn-in. Default is 250,000.
#' @param thin Thinning interval
#' @param nchains Number of MCMC chains
#' @param r1Prior A string defining prior for the growth parameter r1. Default is 'dnorm(mean=0,sd=0.05)'. 
#' @param r2Prior A string defining prior for the growth parameter r2. Default is 'dnorm(mean=0,sd=0.05)'. 
#' @param etaPrior A string defining prior for the change point parameter eta. Default is 'dunif(1,z)', where 'z' is the number of time-blocks. 
#' @param r1Sampler A list containing settings for the MCMC sampler. Default is null and employs nimble's Default sampler (RW sampler).
#' @param r2Sampler A list containing settings for the MCMC sampler. Default is null and employs nimble's Default sampler (RW sampler).
#' @param etaSampler A list containing settings for the MCMC sampler. Default is null and employs nimble's Default sampler (RW sampler).
#' @param parallel Logical specifying whether the chains should be run in parallel or not.
#' @param seeds Random seed for each chain. Default is 1:4.
#' @details The function fits a discrete bounded double exponential growth model on the observed data using MCMC as implemented by the nimble package. The Bayesian model consists of a two growth rate parameters (r1 and r2), with the change from r1 and r2 occurring at inferred point in time eta. Users can define suitable priors using character strings for the argument \code{rPrior1},\code{rPrior2}, and \code{cPrior} (for details on how this should be specified please consult the nimble manual). Please note that the function returns posterior of the growth rate normalised by the resolution defined in the \code{ProbMat} class object.  MCMC settings such as the choice the sampler, number of iterations, chains, etc can also be specified.  
#' @return A \code{fitteddoubleExp} class object containing the original ProbMat class object, posterior of the growth rate, along with its Gelman Rubin statistic and effective sample sizes. 
#' @import nimble
#' @import coda
#' @import parallel
#' @importFrom stats rnorm
#' @export


dexpfit  <- function(x,niter=100000,nburnin=50000,thin=10,nchains=4,r1Prior='dnorm(mean=0,sd=0.05)',r2Prior='dnorm(mean=0,sd=0.05)',etaPrior='dunif(min=1,max=z)',r1Sampler=NULL,r2Sampler=NULL,etaSampler=NULL,parallel=FALSE,seeds=1:4)
{
	#Handle cleaning of GlobalEnv on exit
	envobj <- ls(envir=.GlobalEnv)
	on.exit(rm(list=ls(envir=.GlobalEnv)[which(!ls(envir=.GlobalEnv)%in%envobj)],envir=.GlobalEnv))
	#Addresses R CMD Check NOTES
	returnType <- m.raw <- nimStop <- nimMatrix <-  dDAExp <- rDAExp <- runfun <-  NULL
	# Extract mid points
	mids <- apply(x$tblocks,1,median)
	# Initial Warnings
	if (nchains==1) {warning('Running MCMC on single chain')}

	# Check prior definitions
	supported.distributions  <- data.frame(d=c('dnorm','dexp','dbeta','dchisq','dgamma','dlogis','dunif','dt','dweib','dlnorm'),r=c('rnorm','rexp','rbeta','rchisq','rgamma','rlogis','runif','rt','rweibull','rlnorm'))
	if(!sub("\\(.*","",r1Prior)%in%supported.distributions$d|!sub("\\(.*","",r2Prior)%in%supported.distributions$d|!sub("\\(.*","",etaPrior)%in%supported.distributions$d){stop(paste0('Unsupported distribution in prior definition. Please use one of the following: ',paste(supported.distributions$d,collapse=', ')))}
	r1Prior.rand1  <- supported.distributions$r[supported.distributions$d==strsplit(r1Prior,"\\(")[[1]][1]]
	r1Prior.rand2 <- gsub("\\)","",strsplit(r1Prior,"\\(")[[1]][2])
	r1Prior.rand <- paste0(r1Prior.rand1,"(n=1,",r1Prior.rand2,")")
	r2Prior.rand1  <- supported.distributions$r[supported.distributions$d==strsplit(r2Prior,"\\(")[[1]][1]]
	r2Prior.rand2 <- gsub("\\)","",strsplit(r2Prior,"\\(")[[1]][2])
	r2Prior.rand <- paste0(r2Prior.rand1,"(n=1,",r2Prior.rand2,")")
	etaPrior.rand1  <- supported.distributions$r[supported.distributions$d==strsplit(etaPrior,"\\(")[[1]][1]]
	etaPrior.rand2 <- gsub("\\)","",strsplit(etaPrior,"\\(")[[1]][2])
	etaPrior.rand2 <- gsub("z","constants$z",etaPrior.rand2)
	etaPrior.rand <- paste0(etaPrior.rand1,"(n=1,",etaPrior.rand2,")")


	# Define Data
	d  <- list(theta=x$pmat)
	# Define Constant
	constants  <- list()
	constants$z  <- ncol(d$theta) 

	if (!parallel)
 	{
		
		dDAExp=nimbleFunction(
				     run = function(x = double(2),z=integer(0),r1=double(0),r2=double(0),eta=integer(0), log = integer(0)) {
					     returnType(double(0))
					     t1 = 1:eta
					     t2 = 1:(z-eta)
					     t1length = eta
					     t2length = z-eta
					     n1 = numeric(t1length)
					     n2 = numeric(t2length)
					     for (i in 1:t1length)
					     {
						     n1[i] = (1+r1)^t1[i]
					     }
					     for (i in 1:t2length)
					     {
						     n2[i] = ((1+r1)^(eta)) * (1+r2)^t2[i]
					     }
					     n = c(n1,n2)
					     p = n/sum(n)
					     pg = x %*% p
					     logProb = sum(log(pg))
					     if(log) {
						     return(logProb)
					     } else {
						     return(exp(logProb))
					     }
				     })   
		rDAExp=nimbleFunction(
				     run = function(n=integer(0),z=integer(0),r1=double(0),r2=double(0),eta=integer(0)) {
					     returnType(double(2))
					     nimStop("user-defined distribution dDAExp provided without random generation function.")
					     x <- nimMatrix()
					     return(x)
				     })
		pos <- 1
		assign('dDAExp',dDAExp,envir=as.environment(pos))
		assign('rDAExp',rDAExp,envir=as.environment(pos))
		

		double.expmodel  <- nimbleCode({
			theta[,] ~ dDAExp(r1=r1,r2=r2,eta=eta,z=z)
			r1 ~ r1Prior
			r2 ~ r2Prior
			eta ~ etaPrior
		})

		double.expmodel <- gsub('r1Prior', r1Prior, deparse(double.expmodel)) 
		double.expmodel <- gsub('r2Prior', r2Prior, double.expmodel)
		double.expmodel <- gsub('etaPrior', etaPrior, double.expmodel) |> parse(text=_)

		inits  <- vector('list',length=nchains)
		for (k in 1:nchains)
		{
			set.seed(seeds[k])
			inits[[k]]  <- list(r1=eval(parse(text=r1Prior.rand)),r2=eval(parse(text=r2Prior.rand)),eta=eval(parse(text=etaPrior.rand)))
		}
		message('Compiling nimble model...')
		suppressMessages(model  <- nimbleModel(double.expmodel,constants=constants,data=d,inits=inits[[1]]))
		suppressMessages(cModel <- compileNimble(model))
		suppressMessages(conf <- configureMCMC(model))
		if (!is.null(r1Sampler))
		{
			suppressMessages(conf$removeSamplers('r1'))
			suppressMessages(do.call(conf$addSampler,r1Sampler))
		}
		if (!is.null(r2Sampler))
		{
			suppressMessages(conf$removeSamplers('r2'))
			suppressMessages(do.call(conf$addSampler,r2Sampler))
		}
		if (!is.null(etaSampler))
		{
			suppressMessages(conf$removeSamplers('eta'))
			suppressMessages(do.call(conf$addSampler,etaSampler))
		}

		suppressMessages(conf$addMonitors(c('r1','r2','eta')))
		suppressMessages(MCMC <- buildMCMC(conf))
		suppressMessages(cMCMC <- compileNimble(MCMC))
		results <- runMCMC(cMCMC, niter = niter, thin=thin,nburnin = nburnin,inits=inits,samplesAsCodaMCMC = T,nchains=nchains,progressBar=TRUE,setSeed=seeds)
 	}

	if (parallel)
	{
		message('Running in parallel - progress bar will no be visualised')
		runfun  <- function(seed,constants,d,niter,thin,nburnin,r1Prior,r1Prior.rand,r2Prior,r2Prior.rand,etaPrior,etaPrior.rand,r1Sampler,r2Sampler,etaSampler)
		{

			returnType <- nimStop <- nimMatrix <- dDAExp <- rDAExp <-  NULL
# 			require(nimble)

			dDAExp=nimbleFunction(
					      run = function(x = double(2),z=integer(0),r1=double(0),r2=double(0),eta=integer(0), log = integer(0)) {
						      returnType(double(0))
						      t1 = 1:eta
						      t2 = 1:(z-eta)
						      t1length = eta
						      t2length = z-eta
						      n1 = numeric(t1length)
						      n2 = numeric(t2length)
						      for (i in 1:t1length)
						      {
							      n1[i] = (1+r1)^t1[i]
						      }
						      for (i in 1:t2length)
						      {
							      n2[i] = ((1+r1)^(eta)) * (1+r2)^t2[i]
						      }
						      n = c(n1,n2)
						      p = n/sum(n)
						      pg = x %*% p
						      logProb = sum(log(pg))
						      if(log) {
							      return(logProb)
						      } else {
							      return(exp(logProb))
						      }
					      })   
			rDAExp=nimbleFunction(
					      run = function(n=integer(0),z=integer(0),r1=double(0),r2=double(0),eta=integer(0)) {
						      returnType(double(2))
						      nimStop("user-defined distribution dDAExp provided without random generation function.")
						      x <- nimMatrix()
						      return(x)
					      })
			pos <- 1
			assign('dDAExp',dDAExp,envir=as.environment(pos))
			assign('rDAExp',rDAExp,envir=as.environment(pos))


			double.expmodel  <- nimbleCode({
				theta[,] ~ dDAExp(r1=r1,r2=r2,eta=eta,z=z)
				r1 ~ r1Prior
				r2 ~ r2Prior
				eta ~ etaPrior
			})

			double.expmodel <- gsub('r1Prior', r1Prior, deparse(double.expmodel)) 
			double.expmodel <- gsub('r2Prior', r2Prior, double.expmodel)
			double.expmodel <- gsub('etaPrior', etaPrior, double.expmodel) |> parse(text=_)


			set.seed(seed)
			inits  <- list(r1=eval(parse(text=r1Prior.rand)),r2=eval(parse(text=r2Prior.rand)),eta=eval(parse(text=etaPrior.rand)))
			model  <- nimbleModel(double.expmodel,constants=constants,data=d,inits=inits)
			cModel <- compileNimble(model)
			conf <- configureMCMC(model)
			conf$addMonitors(c('r1','r2','eta'))

			if (!is.null(r1Sampler))
			{
				suppressMessages(conf$removeSamplers('r1'))
				suppressMessages(do.call(conf$addSampler,r1Sampler))
			}
			if (!is.null(r2Sampler))
			{
				suppressMessages(conf$removeSamplers('r2'))
				suppressMessages(do.call(conf$addSampler,r2Sampler))
			}
			if (!is.null(etaSampler))
			{
				suppressMessages(conf$removeSamplers('eta'))
				suppressMessages(do.call(conf$addSampler,etaSampler))
			}

			MCMC <- buildMCMC(conf)
			cMCMC <- compileNimble(MCMC)
			results <- runMCMC(cMCMC, niter = niter, thin=thin,nburnin = nburnin,samplesAsCodaMCMC = T,setSeed=seed)
		}

		ncores  <- nchains
		cl  <- makeCluster(ncores)
		clusterEvalQ(cl,{library(nimble)})
		out  <- parLapply(cl=cl,X=seeds,fun=runfun,d=d,constants=constants,nburnin=nburnin,niter=niter,thin=thin,r1Prior=r1Prior,r1Prior.rand=r1Prior.rand,r2Prior=r2Prior,r2Prior.rand=r2Prior.rand,etaPrior=etaPrior,etaPrior.rand=etaPrior.rand,r1Sampler=r1Sampler,r2Sampler=r2Sampler,etaSampler=etaSampler)
		stopCluster(cl)
		results <- out
	}
	rhat  <- gelman.diag(results,multivariate=FALSE) 
	ess  <- effectiveSize(results)
	if (any(rhat[[1]][,1]>1.01)){warning(paste0('Rhat value above 1.01 (',round(max(rhat[[1]][,1]),3),'). Consider rerunning the model with a higher number of iterations'))}
	posterior <- do.call(rbind.data.frame,results)
	posterior.r1 <- posterior[,'r1']/x$resolution
	posterior.r2 <- posterior[,'r2']/x$resolution
	posterior.eta <- mids[posterior[,'eta']]
	posterior.eta.index <- posterior[,'eta']

	results  <- list(x=x,posterior.r1=posterior.r1,posterior.r2=posterior.r2,posterior.eta=posterior.eta,posterior.eta.index=posterior.eta.index,rhat=rhat,ess=ess)
	class(results)  <- c('fittedDExp',class(results))
	
	envobj.current <- ls(envir=.GlobalEnv)
	toremove  <- envobj.current[which(!envobj.current%in%envobj)]
	rm(list=toremove,envir=.GlobalEnv)
	return(results)
}


