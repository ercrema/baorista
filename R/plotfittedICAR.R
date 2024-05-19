#' @title Plot 1D ICAR model fitted to aoristic data
#' @description Plot posterior estimates of \code{fittedICAR} class objects.
#' @param x An \code{fittedICAR} class object
#' @param hpd A vector with two values defining the highest posterior density interval to display. Default is 0.5 and 0.9.
#' @param minortick Interval for minor ticks in the x-axis label. Default is estimated based on timescale.
#' @param ylim Limits of the y-axis. Default estimated from posterior ranges.
#' @param xlab Label for the x-axis. Default based on \code{calendar}.
#' @param ylab Label for the y-axis. Default is "Probability Mass".
#' @param calendar Either \code{'BP'} or \code{'BCAD'}. Indicate whether the x-axis should be displayed in BP or BC/AD. Default is \code{'BP'}.
#' @param col1 Fill color for the first (inner) HPD interval. Default is 'steelblue'.
#' @param col2 Fill color for the second (outer) HPD interval. Default is 'lightblue'.
#' @param pch Point symbol used to display mean posteriors. Default is 20.
#' @param plot.legend Logical indicating whether to display a legend or not (default is TRUE).
#' @param legend.arg List containing arguments to be directed to the \code{legend()} function.
#' @param ... Additional arguments affecting the plot.
#' @return No return value (plot function)
#' @method plot fittedICAR
#' @import graphics
#' @export


plot.fittedICAR <- function(x,hpd=c(0.5,0.9),minortick=NULL,ylim=NULL,xlab=NULL,ylab='Probability Mass',calendar='BP',col1='steelblue',col2='lightblue',pch=20,plot.legend=TRUE,legend.arg=NULL,...)
{
	oldpar <- par(no.readonly = TRUE)
	on.exit(par(oldpar))
	midPoints  <- apply(x$x$tblocks,1,median)
	scl  <- diff(pretty(midPoints))[1]
	minortick <- ifelse(is.null(minortick),round(scl/5),minortick)

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


	lohi1 <- HPDinterval(as.mcmc(x$posterior.p),prob=hpd[1])
	lohi2 <- HPDinterval(as.mcmc(x$posterior.p),prob=hpd[2])

	if(is.null(ylim)){ylim=c(0,max(lohi2))}

	par(lend=2)
	plot(NULL,axes=FALSE,xlab=xlabel,ylab=ylab,xlim=x$x$timeRange,ylim=ylim,...)

	for (i in 1:nrow(lohi1))
	{
		rect(ybottom=lohi2[i,1],ytop=lohi1[i,1],xleft=midPoints[i]+x$x$resolution/2,xright=midPoints[i]-x$x$resolution/2,border=NA,col=col2)
		rect(ybottom=lohi1[i,1],ytop=lohi1[i,2],xleft=midPoints[i]+x$x$resolution/2,xright=midPoints[i]-x$x$resolution/2,border=NA,col=col1)
		rect(ybottom=lohi1[i,2],ytop=lohi2[i,2],xleft=midPoints[i]+x$x$resolution/2,xright=midPoints[i]-x$x$resolution/2,border=NA,col=col2)
	}

	points(midPoints,apply(x$posterior.p,2,mean),pch=pch)
	axis(1,at=labs.pos,labels=labs)
	axis(1,at=minortick.pos,labels=NA,tck=-0.01)
	axis(2)
	box()
	if (plot.legend)
	{
		tmp.list <- list(legend = c(paste0(round(hpd[1]*100),'% HPDI'),paste0(round(hpd[2]*100),'% HPDI'),'Posterior Mean'),col = c(col1,col2,1),lwd=c(8,8,NA),pch=c(NA,NA,pch))
		if(is.null(legend.arg))
			{
				legend.arg2=c(list(x='topleft'),tmp.list)
				do.call(legend,legend.arg2)
			}
		if(!is.null(legend.arg))
			{
				legend.arg2 <- c(legend.arg,tmp.list)
				if (!'x' %in% names(legend.arg2))
				{legend.arg2 <- c(list(x='topleft'),legend.arg2)}
				do.call(legend,legend.arg2)
			}
	}
}
















