#' @title Sample aoristic data (data.frame)
#' @description Sample datasets to illustrate data formats required for \code{createProbMat()}.
#' @format A data.frame class object with two columns defining the start and the end of each even (\code{sample.df}) 
#' @examples
#' data(sampledf)
#' x  <- createProbMat(x=sampledf,timeRange = c(6500,4001),resolution= 100)
"sampledf"

#' @title Sample aoristic data (matrix)
#' @description Sample datasets to illustrate data formats required for \code{createProbMat()}.
#' @format A matrix class object storing the probability of each event (row) in each time-block (column) 
#' @examples
#' data(samplemat)
#' x <- createProbMat(pmat=samplemat,timeRange = c(5000,3001),resolution=100)
#' plot(x)
"samplemat"
