#' @title Fits a random walk ICAR model to Aoristic data
#' @description Estimates parameters of a multinomial probability distribution that describes the shape of the of the time-frequency distribution of an observed sets of events with chronological uncertainty. The function is wrapper for fitting a 1D random walk ICAR model via nimble. 
#' @param x A ProbMat class object
#' @param niter Number of MCMC iterations. Default is 500,000.
#' @param nburnin Mumber of iterations discarded for burnin. Default is 250,000.
#' @param thin Thinning interval
#' @param nchains Number of MCMC chains
#' @param sigmaPrior A string defining prior for the sigma parameter. Default is 'dexp(1)'. 
#' @param sigmaSampler A list containing settings for the MCMC sampler. Default is null and employs nimble's Default sampler (RW sampler).
#' @param parallel Logical specifying whether the chains should be run in parallel or not.
#' @param seeds Random seed for each chain. Default is 1:4.
#' @details (To be completed)
#' @import nimble
#' @import coda
#' @import parallel
#' @export


icarfit  <- function(x,niter=100000,nburnin=50000,thin=10,nchains=4,sigmaPrior='dexp(1)',sigmaSampler=NULL,parallel=FALSE,seeds=1:4)
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
	icarmodel  <- nimbleCode({
		theta[,] ~ dAoristicGeneral_vector(p=p[1:n.tblocks])
		for (i in 1:n.tblocks)
		{
# 			p[i]  <- lpseq[i]^2 / lpseqSS
			p[i]  <- exp(lpseq[i]) / lpseqSS
		}
		lpseq[1:n.tblocks] ~ dcar_normal(adj[1:L], weights[1:L], num[1:n.tblocks], tau, zero_mean = 0)
# 		lpseqSS  <- sum(lpseq[1:n.tblocks]^2)
		lpseqSS  <- sum(exp(lpseq[1:n.tblocks]))
		tau <- 1/sigma^2
		sigma  ~ dexp(1)
	})
	assign("dAoristicGeneral_vector",dAoristicGeneral_vector,envir=.GlobalEnv)
	icarmodel <- gsub('dexp\\(1\\)', sigmaPrior, deparse(icarmodel)) |> parse(text=_)

	inits  <- vector('list',length=nchains)
	for (k in 1:nchains)
	{
		inits[[k]]  <- list(sigma=rexp(1),lpseq=rnorm(constants$n.tblocks,0,0.5))
	}
	print('Compiling nimble model...')
	suppressMessages(model  <- nimbleModel(icarmodel,constants=constants,data=d,inits=inits[[1]]))
	assign("rAoristicGeneral_vector",rAoristicGeneral_vector,envir=.GlobalEnv)
	suppressMessages(cModel <- compileNimble(model))
	suppressMessages(conf <- configureMCMC(model))
	if (!is.null(sigmaSampler))
	{
		suppressMessages(conf$removeSamplers('sigma'))
		# 	sigmaSampler=list('sigma',type='slice')
		suppressMessages(do.call(conf$addSampler,sigmaSampler))
	}
	suppressMessages(conf$addMonitors('p'))
	suppressMessages(MCMC <- buildMCMC(conf))
	suppressMessages(cMCMC <- compileNimble(MCMC))
	results <- runMCMC(cMCMC, niter = niter, thin=thin,nburnin = nburnin,inits=inits,samplesAsCodaMCMC = T,nchains=nchains,progressBar=TRUE,setSeed=seeds)

 	}

	if (parallel)
	{
		require(parallel)
		print('Running in parallel - progress bar will no be visualised')
		runfun  <- function(seed,constants,d,niter,thin,nburnin,sigmaPrior,sigmaSampler)
		{
			require(nimble)
			dAOG=nimbleFunction(run = function(x = double(2),p=double(1),log = integer(0))
							       {
								       returnType(double(0))
								       pg = x %*% p
								       logProb = sum(log(pg))
								       if(log) {
									       return(logProb)
								       } else {
									       return(exp(logProb))
								       }
							       })   
			assign('dAOG',dAOG,envir=.GlobalEnv)

			icarmodel  <- nimbleCode({
				theta[,] ~ dAOG(p=p[1:n.tblocks])
				for (i in 1:n.tblocks)
				{
# 					p[i]  <- lpseq[i]^2 / lpseqSS
					p[i]  <- exp(lpseq[i]) / lpseqSS
				}
				lpseq[1:n.tblocks] ~ dcar_normal(adj[1:L], weights[1:L], num[1:n.tblocks], tau, zero_mean = 0)
# 				lpseqSS  <- sum(lpseq[1:n.tblocks]^2)
				lpseqSS  <- sum(exp(lpseq[1:n.tblocks]))
				tau <- 1/sigma^2
				sigma  ~ dexp(1)
			})
			icarmodel <- gsub('dexp\\(1\\)', sigmaPrior, deparse(icarmodel)) |> parse(text=_)
			inits  <- list(sigma=rexp(1),lpseq=rnorm(constants$n.tblocks,0,0.5))
			model  <- nimbleModel(icarmodel,constants=constants,data=d,inits=inits)
			assign('rAOG',rAOG,envir=.GlobalEnv)
			cModel <- compileNimble(model)
			conf <- configureMCMC(model)
			conf$addMonitors('p')
			if (!is.null(sigmaSampler))
			{
				suppressMessages(conf$removeSamplers('sigma'))
				# 	sigmaSampler=list('sigma',type='slice')
				do.call(conf$addSampler,sigmaSampler)
			}
			MCMC <- buildMCMC(conf)
			cMCMC <- compileNimble(MCMC)
			results <- runMCMC(cMCMC, niter = niter, thin=thin,nburnin = nburnin,samplesAsCodaMCMC = T,setSeed=seed)
		}

		ncores  <- nchains
		cl  <- makeCluster(ncores)
		out  <- parLapply(cl=cl,X=seeds,fun=runfun,d=d,constants=constants,nburnin=nburnin,niter=niter,thin=thin,sigmaPrior,sigmaSampler)
		stopCluster(cl)
		results <- out
	}
	rhat  <- gelman.diag(results,multivariate=FALSE) 
	ess  <- effectiveSize(results)
	if (any(rhat[[1]][,1]>1.01)){warning(paste0('Highest Rhat value above 1.01 (',round(max(rhat[[1]][,1]),3),'). Consider rerunning the model with a higher number of iterations'))}
	posteriors <- do.call(rbind.data.frame,results)
	posterior.sigma <- posteriors[,'sigma']
	posterior.p <- posteriors[,grep('p\\[',colnames(posteriors))]
	results  <- list(x=x,posterior.p=posterior.p,posterior.sigma=posterior.sigma,rhat=rhat,ess=ess)
	class(results)  <- c('fittedICAR',class(results))
	return(results)
}


