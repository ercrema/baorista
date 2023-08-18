library(nimble)

# Discrete Uniform Distribution 
dunifdisc<-function(x, min=0, max=1) ifelse(x>=min & x<=max & round(x)==x, 1/(max-min+1), 0)
punifdisc<-function(q, min=0, max=1) ifelse(q<min, 0, ifelse(q>=max, 1, (floor(q)-min+1)/(max-min+1)))
qunifdisc<-function(p, min=0, max=1) floor(p*(max-min+1))
runifdisc<-function(n, min=0, max=1) sample(min:max, n, replace=T)

# Conversion to BCAD/BP
BPtoBCAD <- function(x){
    index <- !is.na(x)
    if (any(x[index] < 0)){ stop("Post-bomb dates (<0 BP) are not currently supported.") }
    res <- matrix(c(x, rep(NA,length(x))), ncol=2)
    res[index & x < 1950,2] <- 1950-res[index & x < 1950,1]
    res[index & x >= 1950,2] <- 1949-res[index & x >= 1950,1]
    return(res[,2])
}

BCADtoBP <- function(x){
    index <- !is.na(x)
    if (any(x[index] == 0)){ stop("0 BC/AD is not a valid year.") }
    if (any(x[index] > 1950)){ stop("Post-bomb dates (> AD 1950) are not currently supported.") }
    res <- matrix(c(x, rep(NA,length(x))), ncol=2)
    res[index & x > 0,2] <- abs(res[index & x > 0,1] - 1950)
    res[index & x < 0,2] <- abs(res[index & x < 0,1] - 1949)
    return(res[,2])
}

icar.struct  <- function(x)
{
	num  <- c(1,rep(2,x-2),1)
	adj  <- numeric()
	for (i in 1:x)
	{
		if(i==1){adj = 2}
		if(i==x){adj = c(adj,x-1)}
		if(i>1&i<x){adj = c(adj,c(i-1,i+1))}
	}
	weight <- rep(1,length(adj))
	return(list(adj=adj,weight=weight,num=num))
}


dAoristicGeneral_vector=nimbleFunction(
  run = function(x = double(2),p=double(1),log = integer(0)) {
    returnType(double(0))
    pg = x %*% p
    logProb = sum(log(pg))
    if(log) {
      return(logProb)
    } else {
      return(exp(logProb))
    }
  })   

suppressMessages(registerDistributions(list(
  dAoristicGeneral_vector = list(
    BUGSdist = "dAoristicGeneral_vector(p)",
    pqAvail = FALSE,
    types = c('value = double(2)', 'p = double(1)')
  ))))


dAoristicExponentialGrowth_vector=nimbleFunction(
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

suppressMessages(registerDistributions(list(
  dAoristicExponentialGrowth_vector = list(
    BUGSdist = "dAoristicExponentialGrowth_vector(z,r)",
    pqAvail = FALSE,
    types = c('value = double(2)', 'z = integer(0)','r=double(0)')
  ))))
