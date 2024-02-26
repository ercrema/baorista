#' @title Creates a probMat class object from user data
#' @description Converts either a data.frame with the start and end date of each event or matrix of probabilities values into a \code{probMat} class object.
#' @param x A data.frame containing the start and end date of the timespan of each event. Dates should be in BP, with the first column defining the start and the second column defining the end of the timespan.
#' @param pmat A matrix of aoristic weights (probabilities), with row representing events and column representing timeblocks.
#' @param timeRange A vector of two numerical values representing the start and end of the window of analysis in BP.
#' @param resolution Resolution of the timeblock. Ignored if \code{pmat} is provided.
#' @return An object of class \code{probMat}.
#' @export 



createProbMat <- function(x=NULL,pmat=NULL,timeRange=NULL,resolution=NULL)
{
	#Warnings and Checks
	if (is.null(x)&is.null(pmat)){stop("Either x or pmat are needed to generate the probMat objects")}
	if (!is.null(x)&!is.null(pmat)){stop("Both x and pmat are provided.")}
	if (is.null(timeRange))	{stop("timeRange not provided.")}
	if (is.null(resolution)) {stop("resolution not provided.")}
	if (length(timeRange)!=2 | timeRange[1]<timeRange[2]){stop("Incorrect format of timeRange argument")}
	if (!is.null(pmat)){if (any(apply(pmat,1,sum)!=1)){stop("Sum of event probabilities are not always equal to 1")}}
	starts  <- seq(from=timeRange[1],to=timeRange[2],-resolution)
	ends  <- seq(from=timeRange[1]-resolution + 1, to=timeRange[2],-resolution) 
	if (length(starts)!=length(ends)){stop("resolution does not break timeRange in equally sized time blocks")}
	tblocks  <-  data.frame(starts = seq(from=timeRange[1],to=timeRange[2],-resolution), ends = seq(from=timeRange[1]-resolution + 1, to=timeRange[2],-resolution))
	if (!is.null(pmat)){if(ncol(pmat)!= nrow(tblocks)){stop("Incorrect number of columns in pmat given the resolution and timeRange provided")}}
	if (!is.null(x)){if(any(x[,1]<=x[,2])){stop("Some events have a start point  of timespan later than its end point")}}
	if (!is.null(x)){if(any(x>timeRange[1]|x<timeRange[2])){stop("Some event have timespan outside the timeRange provided")}}
	if (!is.null(x))
	{
		pmat  <- matrix(0,nrow=nrow(x),ncol=nrow(tblocks))
		for (i in 1:ncol(pmat))
		{
			pmat[,i]  <- apply(x,1,function(x,a,b){sum(dunifdisc(a:b,max=x[1],min=x[2]))},a=tblocks[i,2],b=tblocks[i,1])
		}

	}

	probMat.list  <- list()
	probMat.list$pmat  <- pmat
	probMat.list$n <- nrow(probMat.list$pmat)
	probMat.list$z  <- ncol(probMat.list$pmat)
	probMat.list$tblocks  <- tblocks
	probMat.list$timeRange  <- timeRange
	probMat.list$resolution  <- resolution

	class(probMat.list)  <- c("probMat",class(probMat.list))
	return(probMat.list)
}
