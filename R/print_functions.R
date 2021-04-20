#' Print function for lgcpRealPlot
#'
#' @param x An lgcpRealPlot object
#' @param ... ...
#' @return Plots the object which is a list of ggplot2 plots.
#' @export
print.lgcpRealPlot <- function(x, ...){
  print(ggpubr::ggarrange(x[[1]],x[[2]],nrow=1))
}

#' Print function for lgcpRealSumm
#'
#' @param x An lgcpRealSumm object
#' @param digits Number of digits to round the summary to
#' @param ... ...
#' @return A summary of the output of an lgcp model.
#' @export
print.lgcpRealSumm <- function(x, ...,digits=3){
  print(round(x$posterior,digits))
}

#' Print function for lgcpReal
#'
#' @param x An lgcpRealPlot object
#' @param ... ...
#' @export
print.lgcpReal <- function(x, ...){
  print(paste0("An lgcpReal model fit with ",nrow(x$lgcpRunInfo$timetaken)," chains."))
}
