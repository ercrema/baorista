#' @title Plot Aoristic Sums
#'
#' @description Plot summed probabilities of aoristic weights.
#' @param x \code{probMat} class object generated using the \code{generateProbMat()}.
#' @param xlab Label for the x-axis. Default based on \code{calendar}.
#' @param ylab Label for the y-axis. Default is 'Summed Probability' (if \code{type}='asum') or 'Probability Mass' (when \code{type}='dens').
#' @param type Either 'asum' for Aoristic Sum, 'dens' for probability mass. Default is 'asum'.
#' @param calendar Either \code{'BP'} or \code{'BCAD'}. Indicate whether the x-axis should be displayed in BP or BC/AD. Default is  \code{'BP'}.
#' @param lwd Line width. Default is 1.
#' @param col Line col. Default is 'black'
#' @param minortick Interval for minor ticks in the x-axis label. Default is estimated based on timescale
#' @param add if set to \code{TRUE} adds the line and point graph on existing plot.
#' @param ... Additional arguments affecting the plot.
#' @return No return value (plot function)
#' @method plot probMat
#' @export  

plot.probMat <- function(x,xlab=NULL,ylab=NULL,type='asum',calendar='BP',lwd=1,col=1,minortick=NULL,add=FALSE,...)
{
	midPoints  <- apply(x$tblocks,1,median)
	asum  <- apply(x$pmat,2,sum)
	if(!type%in%c('asum','dens')){stop('The argument "type" must be equal to "asum" or "dens"')}
	if (type=='asum')
	{
		yplot <- asum
		ylab <- ifelse(is.null(ylab),'Summed Probability',ylab)
	}

	if (type=='dens')
	{
		yplot <- asum/sum(asum)
		ylab <- ifelse(is.null(ylab),'Probability Mass',ylab)
	}


	scl  <- diff(pretty(midPoints))[1]
	minortick <- ifelse(is.null(minortick),round(scl/5),minortick)

	#Setting calendar and xlim
	if (calendar=="BP"){
		xlabel <- ifelse(is.null(xlab),"Years cal BP",xlab)
		labs  <- labs.pos <- pretty(midPoints)
		minortick.pos  <- seq(max(labs.pos+scl),min(labs.pos-scl),-minortick)
	}
	if (calendar=="BCAD"){
		xlabel <- ifelse(is.null(xlab),"Years BC/AD",xlab)
		labs  <- pretty(BPtoBCAD(midPoints))
		if (all(labs<0)){xlabel <- ifelse(is.null(xlab),"Years BC",xlab)}
		if (all(labs>0)){xlabel <- ifelse(is.null(xlab),"Years AD",xlab)}
		xlabel <- ifelse(is.null(xlab),"Years BC/AD",xlab)
		scl.2  <- diff(labs)
		minortick.pos  <- BCADtoBP(seq(min(labs-scl),max(labs+scl),minortick))
		labs.pos  <- BCADtoBP(labs)
		labs  <- abs(labs)
	}


	if (add==TRUE)
	{
		lines(midPoints,yplot,type='b',col=col,lwd=lwd,pch=20,...)
	}

	if (add==FALSE)
	{
		plot(midPoints,yplot,type='b',col=col,pch=20,lwd=lwd,axes=FALSE,xlab=xlabel,ylab=ylab,xlim=x$timeRange,...)
		axis(1,at=labs.pos,labels=labs)
		axis(1,at=minortick.pos,labels=NA,tck=-0.01)
		axis(2)
		box()
	}
}




