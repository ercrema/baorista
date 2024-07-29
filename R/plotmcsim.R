#' @title Plot Monte-Carlo simulation results on aoristic data
#' @description Plot Monte-Carlo simulation based percentile  intervals on frequency or rate of change of events.
#' @param x A \code{mcsimres} class object generated using the \code{mcsim()} function.
#' @param interval A value between 0 and 1 defining the percentile interval. Default is 0.9.
#' @param minortick Interval for minor ticks in the x-axis label. Default is estimated based on timescale.
#' @param ylim Limits of the y-axis. Default estimated from posterior ranges.
#' @param xlab Label for the x-axis. Default based on \code{calendar}.
#' @param ylab Label for the y-axis. Default is "Probability Mass".
#' @param calendar Either \code{'BP'} or \code{'BCAD'}. Indicate whether the x-axis should be displayed in BP or BC/AD. Default is \code{'BP'}.
#' @param col Color of Monte-Carlo simulation mean. Default is black.
#' @param lwd Line width of Monte-Carlo  mean. Default is 1.
#' @param lty Line type Monte-Carlo  mean. Default is 1.
#' @param col.fill Fill color for the first (inner) percentile  interval. Default is 'lightblue'.
#' @param pch Point symbol used to display mean posteriors. Default is 20.
#' @param type Determine whether to display total number of events (if set to 'sum') or the rate of change ('roc'), computed as (t0/t1)^(1/r)-1, where t0 is the number of events in given time-block t, t1 is the number of events of the next time-block t+1, and r is the size (in years) of the time-blocks. Defaults is 'sum'.
#' @param plot.legend Logical indicating whether to display a legend or not (default is TRUE).
#' @param legend.arg List containing arguments to be directed to the \code{legend()} function.
#' @param ... Additional arguments affecting the plot.
#' @return No return value (plot function)
#' @method plot mcsimres
#' @import graphics
#' @importFrom stats quantile
#' @export


plot.mcsimres <- function(x,interval=0.9,minortick=NULL,ylim=NULL,xlab=NULL,ylab=NULL,calendar='BP',col='black',lwd=1,lty=1,col.fill='lightblue',pch=20,type='sum',plot.legend=TRUE,legend.arg=NULL,...)
{
	oldpar <- par(no.readonly = TRUE)
	on.exit(par(oldpar))

	if (type=='sum'){

		midVals <- apply(x$sums,1,median) 
		lo  <- apply(x$sums,1,quantile,prob=(1-interval)/2)
		hi  <- apply(x$sums,1,quantile,prob= 1-(1-interval)/2)
		midPoints  <- apply(x$tblocks,1,median)
		ylabel <- ifelse(is.null(ylab),"Number of Events",ylab)
	}

	if (type=='roc'){

		roc.mat <- matrix(NA,nrow=x$z-1,ncol=x$nsim)
		r <- x$res
		for (i in 1:(x$z-1))
		{
			t0 <- x$sums[i,]
			t1 <- x$sums[i+1,]
			roc.mat[i,]  <- (t0/t1)^(1/r)-1
		}
		roc.mat[which(roc.mat==Inf|roc.mat==-Inf|is.nan(roc.mat),arr.ind=T)] <- NA
		midVals <- apply(roc.mat,1,median,na.rm=T) 
		lo  <- apply(roc.mat,1,quantile,prob=(1-interval)/2,na.rm=T)
		hi  <- apply(roc.mat,1,quantile,prob= 1-(1-interval)/2,na.rm=T)
		midPoints  <- apply(x$tblocks,1,median) + x$res/2
		midPoints  <- midPoints[-x$z]
		ylabel <- ifelse(is.null(ylab),"Rate of Change",ylab)
	}



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


	if(is.null(ylim) & type=='sum'){ylim=c(0,max(hi))}
	if(is.null(ylim) & type=='roc'){ylim=c(min(lo),max(hi))}

	par(lend=2)
	plot(NULL,axes=FALSE,xlab=xlabel,ylab=ylabel,xlim=x$timeRange,ylim=ylim,...)
	polygon(x=c(midPoints,rev(midPoints)),y=c(lo,rev(hi)),border=NA,col=col.fill)
	lines(midPoints,midVals,lwd=lwd,lty=lty,type='b',pch=pch)

	axis(1,at=labs.pos,labels=labs)
	axis(1,at=minortick.pos,labels=NA,tck=-0.01)
	axis(2)
	if(type=='roc'){abline(h=0,lty=2)}
	box()
	if (plot.legend)
	{
		tmp.list <- list(legend = c(paste0(round(interval*100),'th percentile'),'Mean'),col = c(col.fill,1),lwd=c(8,1),pch=c(NA,pch))
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

