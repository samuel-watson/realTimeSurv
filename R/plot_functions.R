#' Extract data from an lgcpReal object
#'
#' Helper function to create plotting data from an lgcpReal object
#'
#'@param samps MCMC samples extracted from lgcpReal object
#'@param grid.data Grid used for sampling and plotting
#'@param lg lgcpReal model fit
#'@param nchains Number of MCMC chains in model fit
#'@param idx.mapping Matrix specifying plot order of grid cells
#'@param cellwidth Cellwidth of the grid
#'@param covariates SpatialPolygonsDataFrame with the original covariate data
#'@param plotlag Laglength for the plot
#'@param data.out Logical indicating whether to return data frame with only
#'plotting data (FALSE) or a list containing each predicted model subcomponent (TRUE)
#'@importFrom stats model.matrix model.frame sd
#'@export
plot_lgcp_dat <- function(samps,
                          grid.data,
                          lg,
                          nchains,
                          idx.mapping=NULL,
                          cellwidth=0.005,
                          covariates,
                          plotlag=0,
                          data.out=FALSE){
  beta <- lg$beta
  ins <- lg$cellInside
  ins <- c(ins)
  grid.data <- grid.data[order(grid.data$x,-grid.data$y),]
  dat1 <- grid.data
  grid.data <- grid.data[idx.mapping,]
  dat1 <- dat1[idx.mapping,]
  dat1pts <- dat1[,c('x','y')]
  sp::coordinates(dat1pts) <- ~x+y
  if(!sp::identicalCRS(dat1pts,covariates)){
    sp::proj4string(dat1pts) <- sp::proj4string(covariates)
  }
  covs <- sp::over(dat1pts,covariates)
  dat1 <- cbind(dat1,covs)

  if(lg$formulae$form.sp=="X ~ 1"){
    lg$formulae$form.sp <- NULL
  }

  #not general enough in case of different naming

  if(!is.null(lg$formulae$form.pop)){
    pop_cols <- which(grepl("pop",colnames(dat1)))
    for(i in pop_cols){
      dat1[is.na(dat1[,i])|dat1[,i]==0,i] <- min(dat1[dat1[,i]>0&!is.na(dat1[,i]),i],na.rm=T)
    }
  }

  dat1 <- dat1[I(ins==1),]
  samps <- samps[I(ins==1),,]
  dat1$X <- 1

  if(!is.null(lg$formulae$form.t)){
    lint <- model.matrix(~dow,data=lg$data.t,na.action="na.pass")
  }

  if(!is.null(lg$formulae$form.sp)){
    linpred <- model.matrix(lg$formulae$form.sp,
                            model.frame(~ ., dat1, na.action="na.pass"))
  }

  if(!is.null(lg$formulae$form.sp)&&
     !is.null(lg$formulae$form.t)&&
     length(which(lint[which(lg$data.t$t==(max(lg$data.t$t)- plotlag))[1],-1]==1))==1){
    linpred <- exp(linpred[,2:ncol(linpred)]%*%t(lg$beta[,2:ncol(linpred)]) +
                     lg$beta[,ncol(linpred)+which(lint[which(lg$data.t$t==(max(lg$data.t$t)- plotlag))[1],-1]==1)])
  } else if(!is.null(lg$formulae$form.sp)){
    linpred <- exp(linpred[,2:ncol(linpred)]%*%t(lg$beta[,2:ncol(linpred)]))
  } else {
    linpred <- matrix(1,nrow=nrow(samps),ncol=ncol(samps))
  }

  if(!is.null(lg$formulae$form.pop)){
    pop <- matrix(dat1[,all.vars(lg$formulae$form.pop)[2]]*cellwidth^2,ncol=1)%*%
      matrix(exp(lg$beta[,1]+mean(samps[,,dim(samps)[3]-plotlag])),nrow=1)
  } else {
    pop <- matrix(cellwidth^2,nrow=nrow(dat1),ncol=1)%*%
      matrix(exp(lg$beta[,1]+mean(samps[,,dim(samps)[3]-plotlag])),nrow=1)
  }

  samps <- exp(samps[,,dim(samps)[3]-plotlag]-mean(samps[,,dim(samps)[3]-plotlag]))

  if(!is.null(lg$formulae$form.pop)){
    dat1$pop <- dat1[,all.vars(lg$formulae$form.pop)[2]]
  } else {
    dat1$pop <- 1
  }

  dat1$poppred <- rowMeans(pop)
  dat1$linpred <- rowMeans(linpred)
  dat1$xpred <- rowMeans(samps)

  v0 <- pop*linpred*samps

  dat1$value <- rowMeans(v0)
  dat1$sd <- log(apply(v0,1,sd))

  if(!data.out){
    return(dat1)
  } else {
    return(list(pop=pop,linpred=linpred,xpred=samps,dat1=dat1))
  }

}

#' Extract samples from \code{lgcp} model
#'
#' Extract MCMC samples from a call to \code{lgcp}
#'
#' @param dirname Rootname of the directory in call to \code{lgcp}.
#' @param nchains Integer, number of chains from call to \code{lgcp}.
#' @return An array of dimension number of gridcells x number of iterations x
#' number of time periods.
#' @importFrom utils flush.console
#' @export
lgcpExtract <- function(dirname, nchains){
  mcmc <- list()
  c1 <- list()
  for(i in 1:nchains){
    mcmc[[i]] <- ncdf4::nc_open(paste0(dirname,".",i,"/simout.nc"))
    c1[[i]] <- ncdf4::ncvar_get(mcmc[[i]],'simrun')
  }
  d1 <- dim(c1[[1]])
  out <- array(NA, dim=c(d1[1]*d1[2],d1[4]*nchains,d1[3]))
  for(t in 1:d1[3]){
    for(i in 1:nchains){
      for(j in 1:d1[4]){
        out[,j+(i-1)*d1[4],t] <- c(c1[[i]][,,t,j])
      }
    }
    cat("\r|",rep("=",floor((t/d1[3])/0.05)),rep(" ",ceiling(((d1[3]-t)/d1[3])/0.05)),"| ",
        floor((t/d1[3])/0.05)*5,"%",sep="");flush.console()
  }

  attr(out,"model") <- "lgcp"
  attr(out, "dirname") <- dirname
  return(out)
}

#' Real-time surveillance plot
#'
#' Plot incidence, model components, and their changes over time.
#'
#' @param x An lgcpReal object, output from a call to \code{lgcp}
#' @param covariates A \code{spatialPolygonsDataFrame} covering the area of interest and containing
#' the covariate and population density data. Typically the same object as specified in the
#' \code{covariates} argument in the call to \code{lgcp}.
#' @param osm A logical value whether a map from OpenStreetMap should be included in the plots.
#' @param per.days Integer, the number of person-days to use for incidence, default is 10,000.
#' @param change.lag If not NULL, then plots are created of the change in outputs compared to
#' this number of periods prior.
#' @param relative A logical value indicating whether the comparisons (if change.lag set) should be relative (default),
#' i.e. incidence rate ratios and ratios of relative risks, or absolute.
#' @param msq Integer, the denominator of the population density, default is hectares (population per
#' 10,000m^2)
#' @param rr_lim Integer, for plotting the relative risk, the maximum value of the colour scale. Useful
#' when comparing multiple plots to put colour gradient on same scale.
#' @param ... ...
#' @return A list of two ggplot objects. The first is the incidence (or change in incidence) the
#' second is a plot of four components: (i) the expected case count in each cell, (ii) the
#' relative risk due to included covariates, (iii) the relative risk associated with the
#' latent Gaussian process, and (iv) the posterior standard deviation of the incidence. An object
#' \code{outl} is exported to the global environment to reduce needing to reload sampling
#' data on further calls to the same \code{lgcpReal} object. This can be removed if needed as
#' it can be large.
#' @seealso aggregator, plot_hotspot, generate_report
#' @importFrom ggplot2 ggplot aes geom_raster theme theme_bw element_blank scale_fill_viridis_c scale_fill_viridis_d
#' @importFrom ggplot2 scale_fill_gradientn ggtitle coord_equal element_rect geom_tile facet_wrap geom_point geom_path
#' @importFrom methods is
#' @importFrom rlang .data
#' @export
plot.lgcpReal <- function(x,
                          covariates,
                          osm=FALSE,
                          per.days=10000,
                          change.lag=NULL,
                          relative=TRUE,
                          msq = 10000,
                          rr_lim=NULL,
                          ...){
  if(missing(covariates))stop("Specify the covariate spatial polygons data frame")
  if(!is(x,"lgcpReal"))stop("x must be of class lgcpReal")
  OW <- lgcp::selectObsWindow(x$xyt, cellwidth = x$cellwidth)

  grid.data <- expand.grid(x=OW$xvals,y=OW$yvals)
  idx.mapping <- matrix(1:nrow(grid.data),nrow=length(OW$yvals),ncol=length(OW$xvals))
  idx.mapping <- c(t(apply(idx.mapping,2,rev)))

  if(file.exists(paste0(tempdir(),"\\outl.RDS"))){
    outl <- readRDS(paste0(tempdir(),"\\outl.RDS"))
  } else {
    outl <- lgcpExtract(x$dirname,nrow(x$lgcpRunInfo$timetaken))
    saveRDS(outl,paste0(tempdir(),"\\outl.RDS"))
  }

  if(attr(outl, "dirname")!=x$dirname){
    outl <- lgcpExtract(x$dirname,nrow(x$lgcpRunInfo$timetaken))
    saveRDS(outl,paste0(tempdir(),"\\outl.RDS"))
  }

  # if(!exists("outl") |(exists("outl")&&attr(outl, "dirname")!=x$dirname)){
  #   print("Extracting posterior samples...")
  #   outl <- lgcpExtract(x$dirname,nrow(x$lgcpRunInfo$timetaken))
  #   assign("outl",outl,.GlobalEnv)
  # }

  res1 <- suppressWarnings( plot_lgcp_dat(outl,
                                          grid.data,
                                          x,
                                          x$nchains,
                                          idx.mapping,
                                          covariates = covariates,
                                          cellwidth = x$cellwidth))

  if(is.null(rr_lim)){
    rr_lim1 <- ceiling(max(c(res1$linpred,res1$xpred),na.rm=T))
    rr_lim2 <- NULL
  } else {
    rr_lim1 <- rr_lim
  }
  col_lim <- c(0.01,rr_lim1)
  col_vals <- c(0,1/(rr_lim1*2),seq(1/rr_lim1,1,length.out=4))

  main_title <- paste0("Incidence per ", per.days," person-days")
  rr_title <- "RR"

  if(!is.null(change.lag)){
    reslag <- suppressWarnings(plot_lgcp_dat(outl,
                                             grid.data,
                                             x,
                                             x$nchains,
                                             idx.mapping,
                                             covariates = covariates,
                                             plotlag = change.lag))

    if(relative){
      res1$poppred <- res1$poppred - reslag$poppred
      res1$linpred <- res1$linpred/reslag$linpred
      res1$xpred <- res1$xpred/reslag$xpred
      res1$value <- res1$value/reslag$value
      res1$sd <- exp(res1$sd + reslag$sd)
      rr_title <- "RRR"
      if(!is.null(rr_lim)){
        rr_lim2 <- rr_lim
      } else {
        rr_lim2 <- ceiling(max(c(res1$value),na.rm=T))
      }

      col_lim2 <- c(0.01,rr_lim2)
      col_vals2 <- c(0,1/(rr_lim2*2),seq(1/rr_lim2,1,length.out=4))
    } else {
      res1$poppred <- res1$poppred - reslag$poppred
      res1$linpred <- res1$linpred- reslag$linpred
      res1$xpred <- res1$xpred - reslag$xpred
      res1$value <- res1$value - reslag$value
      res1$sd <- exp(res1$sd) + exp(reslag$sd)
      col_lim <- c(-1,1)
      col_vals <- c(0,0.2,0.5,0.6,0.8,1)
      rr_title <- "Diff RR"
    }


    main_title <- paste0("Change in incidence vs. ",change.lag," periods ago")
  } else {
    if(!is.null(per.days)){
      cent <- sp::coordinates(rgeos::gCentroid(covariates))
      res1$value <- res1$value*per.days/(res1$pop*cos(cent[2]*pi/180)*
                                           111*1000^2*111.321*x$cellwidth^2/msq)
    }


  }

  ppop <- ggplot(data=res1,aes(x=.data$x,y=.data$y,fill=.data$poppred))+
    geom_raster()+
    scale_fill_viridis_c(name="")+
    theme_bw()+
    theme(panel.grid = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA, size=1))+
    ggtitle("Expected")+
    coord_equal()

  plin <- ggplot(data=res1,aes(x=.data$x,y=.data$y,fill=.data$linpred))+
    geom_raster()+
    scale_fill_gradientn(name=rr_title,colours = c("purple","blue","yellow","orange","red","brown"),
                         values = col_vals,limits=col_lim)+
    theme_bw()+
    theme(panel.grid = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA, size=1))+
    ggtitle("Observed")+
    coord_equal()

  px <- ggplot(data=res1,aes(x=.data$x,y=.data$y,fill=.data$xpred))+
    geom_raster()+
    scale_fill_gradientn(name=rr_title,colours = c("purple","blue","yellow","orange","red","brown"),
                         values = col_vals,limits=col_lim)+
    theme_bw()+
    theme(panel.grid = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA, size=1))+
    ggtitle("Latent")+
    coord_equal()

  psd <- ggplot(data=res1,aes(x=.data$x,y=.data$y,fill=exp(.data$sd)))+
    geom_raster()+
    scale_fill_viridis_c()+
    theme_bw()+
    theme(panel.grid = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA, size=1))+
    ggtitle("Posterior SD")+
    coord_equal()

  ppred <- ggplot(data=res1,aes(x=.data$x,y=.data$y,fill=.data$value))+
    geom_raster()+
    theme_bw()+
    theme(panel.grid = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA, size=1))+
    ggtitle(main_title)+
    coord_equal()

  if(!is.null(change.lag)){
    ppred <- ppred + scale_fill_gradientn(name="IRR",
                                          colours = c("purple","blue","yellow","orange","red","brown"),
                                          values = col_vals2,
                                          limits=col_lim2)
  } else {
    ppred <- ppred + scale_fill_viridis_c()
  }

  if(osm){
    if(requireNamespace("ggmap", quietly=TRUE)){

      xrange <- range(ppred$data$x)
      yrange <- range(ppred$data$y)
      #our background map
      mad_map <- get_map2(c(left=xrange[1],bottom=yrange[1],right=xrange[2],top=yrange[2]),
                          source="stamen",
                          maptype = "toner")

      ppred <- ggmap::ggmap(mad_map) +
        geom_tile(data=ppred$data[!is.na(ppred$data$value),],
                  aes(x=.data$x,y=.data$y,fill=.data$value),alpha=0.4)+
        theme(panel.border = element_rect(colour = "black", fill=NA, size=1))

      if(!is.null(change.lag)){
        ppred <- ppred + scale_fill_gradientn(name="IRR",colours = c("purple","blue","yellow","orange","red","brown"),
                                              values = col_vals2,limits=col_lim2)
      } else {
        ppred <- ppred + scale_fill_viridis_c(name="",option="B")
      }
    } else {
      warning("ggmap package required for osm plotting.")
    }

  }

  prow <- ggpubr::ggarrange(ppop,plin,px,psd,nrow=2,ncol=2)
  #print(ggpubr::ggarrange(ppred,prow,nrow=1,ncol=2))
  #prow$scales$scales[[3]] <- c(rr_lim,rr_lim2)
  out <- list(ppred,prow)
  class(out) <- "lgcpRealPlot"
  attr(out,"type") <- "main"
  attr(out,"rr_lim") <- c(rr_lim1,rr_lim2)
  #return(invisible(out))
  out
}

#' Hotspot mapping and visualisation
#'
#' A function for mapping hotspots according to user defined criteria.
#'
#' A ``hotspot'' is defined as an area that exceeds a user-defined criterion with
#' probability of at least p. The criterion can be a function of one or two variables
#' derived from the model; where two variables are used then there are four possible
#' hotspot classifications, where only one is used then there are two classifications
#' (above or below the threshold).
#'
#' The log-linear model can be divided into a set of multiplicative components:
#'
#' (A) population density x (B) size of the area x (C) average disease rate x
#'          (D) RR observed covariates x (E) RR latent process
#'
#' A threshold can be any combination of these factors, or their difference over time.
#' The user can specify the combination using the labels
#' (A)x(C) \code{poppp}
#' (A)x(B)x(C) \code{pop}
#' (D) \code{obs}
#' (E) \code{latent}
#' in the argument to \code{threshold.var} as an additive sum. For example, to specify
#' the incidence (in person-days) as the variable 'poppp+obs+latent', or to specify
#' the overall relative risk of an area 'obs+latent'. To difference the variable with
#' respect to t time periods prior, add '+lag(t)'. So to use the incidence rate ratio
#' relative to 7 days prior, we can specify 'poppp+obs+latent+lag(7)'. The 'hotspot' is
#' an area where Pr(variable > threshold) > p.
#'
#' Hotspots are labelled in the following way. For a single variable definition, the labels are given
#' as \code{c(a,b)} where
#'
#' a = Pr(variable > threshold) <= p
#'
#' b = Pr(variable > threshold) > p
#'
#' For a two variable definition the labels are \code{c(a,b,c,d)} where
#'
#' a = Pr(variable 1 > threshold 1) <= p1 & Pr(variable 2 > threshold 2) <= p2
#'
#' b = Pr(variable 1 > threshold 1) > p1 & Pr(variable 2 > threshold 2) <= p2
#'
#' c = Pr(variable 1 > threshold 1) <= p1 & Pr(variable 2 > threshold 2) > p2
#'
#' d = Pr(variable 1 > threshold 1) > p1 & Pr(variable 2 > threshold 2) > p2
#'
#' The labels do not need to be unique.
#'
#' @param lg Output from a call to \code{lgcp}
#' @param covariates A \code{spatialPolygonsDataFrame} covering the area of interest and containing
#' the covariate and population density data. Typically the same object as specified in the
#' \code{covariates} argument in the call to \code{lgcp}.
#' @param threshold.var A vector of one or two strings specifying the variables to define the hotspots,
#' see Details for how to specify.
#' @param threshold.value A vector or one or two values indicating the threshold(s) for determining
#' a hotspot. Given in the same order as threshold.var.
#' @param labels A vector of two or four labels for the hotspots, see Details.
#' @param threshold.prob A vector of one or two values specifying the exceedence probabilities.
#' @param relative A logical value. If one or both of the variable is with respect to a previous time period, whether the comparison
#' should be relative (TRUE) or absolute (FALSE)
#' @param osm A logical value indicating Whether to include a Open Street Map map under the plot.
#' @param per.days If one or both of the variables is incidence, the denominator number of person-days.
#' @param msq The denominator for the population density in m^2. Default is hectares (per 10,000m^2)
#' @return  An lgcpRealPlot object comprising a list of two ggplot objects.
#' The first is the hotspot classifications, the second the exceedence probabilities. An object
#' \code{outl} is exported to the global environment to reduce needing to reload sampling
#' data on further calls to the same \code{lgcpReal} object. This can be removed if needed as
#' it can be large.
#' @importFrom methods as is
#' @importFrom rlang .data
#' @examples
#' \dontrun{
#' p1 <- plot_hotspot(lg,
#'                    covariates=lsoa,
#'                    threshold.var=c('poppp+obs+latent'),
#'                    threshold.value = c(1),
#'                    labels = c('low','high'),
#'                    osm = TRUE)
#'
#' p2 <- plot_hotspot(lg,
#'                    covariates=lsoa,
#'                    threshold.var=c('poppp+obs+latent','poppp+obs+latent+lag(7)'),
#'                    threshold.value = c(1,1.5),
#'                    labels = c('low','high','rising','both'),
#'                    threshold.prob=0.5,
#'                    osm = TRUE)
#' }
#' @export
plot_hotspot <- function(lg,
                         covariates,
                         threshold.var=NULL,
                         threshold.value=NULL,
                         labels,
                         threshold.prob=0.8,
                         relative=TRUE,
                         osm=FALSE,
                         per.days=10000,
                         msq=10000){

  if(!length(threshold.var)%in%1:2)stop("Name only one or two threshold variables.")
  if((length(labels)!=2&length(threshold.value)==1)|
     (length(labels)!=4&length(threshold.value)==2))stop("For 1/2 thresholds, 2/4 labels are required.")
  if(length(labels)==4&length(relative)==1){
    relative <- c(relative,relative)
  }

  OW <- lgcp::selectObsWindow(lg$xyt, cellwidth = lg$cellwidth)
  grid.data <- expand.grid(x=OW$xvals,y=OW$yvals)
  idx.mapping <- matrix(1:nrow(grid.data),nrow=length(OW$yvals),ncol=length(OW$xvals))
  idx.mapping <- c(t(apply(idx.mapping,2,rev)))

  if(file.exists(paste0(tempdir(),"\\outl.RDS"))){
    outl <- readRDS(paste0(tempdir(),"\\outl.RDS"))
  } else {
    outl <- lgcpExtract(lg$dirname,nrow(lg$lgcpRunInfo$timetaken))
    saveRDS(outl,paste0(tempdir(),"\\outl.RDS"))
  }

  if(attr(outl, "dirname")!=lg$dirname){
    outl <- lgcpExtract(lg$dirname,nrow(lg$lgcpRunInfo$timetaken))
    saveRDS(outl,paste0(tempdir(),"\\outl.RDS"))
  }

  # if(!exists("outl") |(exists("outl")&&attr(outl, "dirname")!=lg$dirname)){
  #   print("Extracting posterior samples...")
  #   outl <- lgcpExtract(lg$dirname,nrow(lg$lgcpRunInfo$timetaken))
  #   assign("outl",outl,.GlobalEnv)
  # }

  str1 <- unlist(strsplit(threshold.var[1],"\\+"))

  if(any(grepl("lag",str1))){
    lag1 <- str1[which(grepl("lag",str1))]
    lag1 <- as.numeric(gsub("\\D", "", lag1))
  } else {
    lag1 <- 0
  }

  res1 <- plot_lgcp_dat(outl,
                        grid.data,
                        lg,
                        lg$nchains,
                        idx.mapping,
                        covariates = covariates,
                        data.out = TRUE)

  v0 <- 1
  if(any(str1=="pop") & any(str1=="poppp")) stop("Cannot have pop and poppp in the same string.")
  if(any(str1=="pop")){
    v0 <- v0*res1$pop
  }
  if(any(str1=="poppp")){
    cent <- sp::coordinates(rgeos::gCentroid(covariates))
    area <- cos(cent[2]*pi/180)*
      111*1000^2*111.321*lg$cellwidth^2/msq
    mult <- per.days/(res1$dat1$pop*area)
    v0x <- sapply(1:nrow(res1$pop),function(i) return(res1$pop[i,]*mult[i]))
    v0 <- v0*t(v0x)
  }

  if(any(grepl("obs",str1))){
    v0 <- v0*res1$linpred
  }
  if(any(grepl("latent",str1))){
    v0 <- v0*res1$xpred
  }

  if(lag1>0){
    res2 <- plot_lgcp_dat(outl,
                          grid.data,
                          lg,
                          lg$nchains,
                          idx.mapping,
                          covariates = covariates,
                          plotlag = lag1,
                          data.out = TRUE)
    v0b <- 1
    if(any(str1=="pop")){
      v0b <- v0b*res2$pop
    }
    if(any(str1=="poppp")){
      v0bx <- sapply(1:nrow(res2$pop),function(i) return(res2$pop[i,]*mult[i]))
      v0b <- v0b*t(v0bx)
    }
    if(any(grepl("obs",str1))){
      v0b <- v0b*res2$linpred
    }
    if(any(grepl("latent",str1))){
      v0b <- v0b*res2$xpred
    }
    rm(res2)
    if(relative[1]){
      v0 <- v0/v0b
    } else {
      v0 <- v0 - v0b
    }

  }


  if(length(threshold.var)==2){
    str2 <- unlist(strsplit(threshold.var[2],"\\+"))

    if(any(grepl("lag",str2))){
      lag2 <- str2[which(grepl("lag",str2))]
      lag2 <- as.numeric(gsub("\\D", "", lag2))
    } else {
      lag2 <- 0
    }

    res3 <- plot_lgcp_dat(outl,
                          grid.data,
                          lg,
                          lg$nchains,
                          idx.mapping,
                          covariates = covariates,
                          data.out = TRUE)

    v1 <- 1
    if(any(str2=="pop") & any(str2=="poppp")) stop("Cannot have pop and poppp in the same string.")
    if(any(str2=="pop")){
      v1 <- v1*res3$pop
    }
    if(any(str2=="poppp")){
      cent <- sp::coordinates(rgeos::gCentroid(covariates))
      area <- cos(cent[2]*pi/180)*
        111*1000^2*111.321*lg$cellwidth^2/msq
      mult <- per.days/(res3$dat1$pop*area)
      v1x <- sapply(1:nrow(res3$pop),function(i) return(res3$pop[i,]*mult[i]))
      v1 <- v1*t(v1x)
    }
    if(any(grepl("obs",str2))){
      v1 <- v1*res3$linpred
    }
    if(any(grepl("latent",str2))){
      v1 <- v1*res3$xpred
    }

    if(lag2>0){
      res4 <- plot_lgcp_dat(outl,
                            grid.data,
                            lg,
                            lg$nchains,
                            idx.mapping,
                            covariates = covariates,
                            plotlag = lag2,
                            data.out = TRUE)
      v1b <- 1
      if(any(str2=="pop")){
        v1b <- v1b*res4$pop
      }
      if(any(str2=="poppp")){
        v1bx <- sapply(1:nrow(res4$pop),function(i) return(res4$pop[i,]*mult[i]))
        v1b <- v1b*t(v1bx)
      }
      if(any(grepl("obs",str2))){
        v1b <- v1b*res4$linpred
      }
      if(any(grepl("latent",str2))){
        v1b <- v1b*res4$xpred
      }
      rm(res4)
      if(relative[2]){
        v1 <- v1/v1b
      } else {
        v1 <- v1 - v1b
      }

    }


    datprobs <- data.frame(a=rep(NA,nrow(v1)),
                           b=rep(NA,nrow(v1)),
                           c=rep(NA,nrow(v1)),
                           d=rep(NA,nrow(v1)))


    datprobs$b <- round(sapply(1:nrow(v1),function(i)length(v0[i,][v0[i,] > threshold.value[1] &
                                                                     v1[i,] <= threshold.value[2]])/
                                 length(v0[i,])),3)

    datprobs$c <- round(sapply(1:nrow(v1),function(i)length(v0[i,][v0[i,] <= threshold.value[1] &
                                                                     v1[i,] > threshold.value[2]])/
                                 length(v0[i,])),3)

    datprobs$d <- round(sapply(1:nrow(v1),function(i)length(v0[i,][v0[i,] > threshold.value[1] &
                                                                     v1[i,] > threshold.value[2]])/
                                 length(v0[i,])),3)

    datprobs$a <- round(1 - datprobs$b - datprobs$c - datprobs$d,3)

    catid <- unique(labels)
    ncats <- length(catid)
    datprobs2 <- matrix(NA,ncol=ncats,nrow=nrow(v0))
    colnames(datprobs2) <- catid

    for(i in 1:ncats){
      if(length(which(labels==catid[i]))>1){
        datprobs2[,i] <- rowSums(datprobs[,which(labels==catid[i])])
      } else {
        datprobs2[,i] <- datprobs[,which(labels==catid[i])]
      }
    }

    res1$dat1$class <- labels[apply(datprobs2,1,which.max)]

    resp <- cbind(res1$dat1[,c('x','y')],datprobs2)
    res1$dat1 <- cbind(res1$dat1,datprobs2)
    resp <- reshape2::melt(resp,id.vars=1:2)
    resp <- resp[!is.na(res1$dat1$value),]
    res1$dat1$class <- factor(res1$dat1$class,levels=labels,ordered=TRUE,labels=labels)


    pclass <- ggplot(data=res1$dat1[!is.na(res1$dat1$value),],aes(x=.data$x,y=.data$y,fill=.data$class))+
      geom_raster()+
      scale_fill_viridis_d(drop=FALSE)+
      theme_bw()+
      theme(panel.grid = element_blank(),
            panel.border = element_rect(colour = "black", fill=NA, size=1))+
      ggtitle("Classification based on modal probability")+
      coord_equal()


    pclass_prop <- ggplot(data=resp,aes(x=.data$x,y=.data$y,fill=.data$value))+
      geom_raster()+
      scale_fill_gradientn(name="Prob",colours = c("purple","blue","green","yellow","red","brown"),
                           values = seq(0,1,length.out=6),limits=c(0,1))+
      facet_wrap(~variable)+
      theme_bw()+
      theme(panel.grid = element_blank(),
            panel.border = element_rect(colour = "black", fill=NA, size=1))+
      ggtitle("Probabilities")+
      coord_equal()


  } else {
    res1$dat1$class_prop <- apply(v0,1,function(i)return(length(i[i > threshold.value[1]])/length(i)))
    res1$dat1$class <- I(res1$dat1$class_prop > threshold.prob)*1
    res1$dat1$class <- factor(res1$dat1$class,levels=c(0,1),labels=labels)

    pclass <- ggplot(data=res1$dat1[!is.na(res1$dat1$value),],aes(x=.data$x,y=.data$y,fill=.data$class))+
      geom_raster()+
      scale_fill_viridis_d(drop=FALSE)+
      theme_bw()+
      theme(panel.grid = element_blank(),
            panel.border = element_rect(colour = "black", fill=NA, size=1))+
      ggtitle(paste0("Classification if Pr(",threshold.var[1],">",threshold.value[1],") > ",round(threshold.prob*100,0),"%"))+
      coord_equal()

    pclass_prop <- ggplot(data=res1$dat1[!is.na(res1$dat1$value),],aes(x=.data$x,y=.data$y,fill=.data$class_prop))+
      geom_raster()+
      scale_fill_gradientn(name="Prob",colours = c("purple","blue","green","yellow","red","brown"),
                           values = seq(0,1,length.out=6),limits=c(0,1))+
      theme_bw()+
      theme(panel.grid = element_blank(),
            panel.border = element_rect(colour = "black", fill=NA, size=1))+
      ggtitle("Probability")+
      coord_equal()

  }

  if(osm){
    if(requireNamespace("ggmap",quietly = TRUE)){
      xrange <- range(res1$dat1$x)
      yrange <- range(res1$dat1$y)

      #our background map
      mad_map <- get_map2(c(left=xrange[1],bottom=yrange[1],right=xrange[2],top=yrange[2]),
                          source="stamen",
                          maptype = "toner")
      pclass <- ggmap::ggmap(mad_map) +
        geom_tile(data=res1$dat1[!is.na(res1$dat1$value),],aes(x=.data$x,y=.data$y,fill=.data$class),alpha=0.4)+
        scale_fill_viridis_d(name="",option="B")+
        theme(panel.border = element_rect(colour = "black", fill=NA, size=1))
    } else {
      stop("ggmap is required for osm plotting")
    }


  }

  #print(ggpubr::ggarrange(pclass,pclass_prop,nrow=1))

  out <- list(pclass,pclass_prop)
  attr(out,"str") <- threshold.var
  attr(out,"threshold") <- threshold.prob
  attr(out,"vals") <- threshold.value
  attr(out,"labs") <- labels
  attr(out,"type") <- "hotspot"
  class(out) <- "lgcpRealPlot"
  #return(invisible(out))
  out
}

#' Summarise lgcp output
#'
#' Summarise output from a call to \code{lgcp}
#'
#' @param object An \code{lgcpReal} object from a call to \code{lgcp}
#' @param linear A logical value indicating whether results should be reported on linear or exponential scales
#' @param plot A logical value indicating whether to produce plots of the prior and posterior distributions of model parameters
#' @param verbose A logical value indicating whether to print running time and chains
#' @param ... ...
#' @return A table with posterior mean, SD, and quantiles of posterior and prior distributions
#' @importFrom stats terms sd quantile qnorm rnorm
#' @importFrom ggplot2 geom_density scale_linetype_discrete scale_color_discrete
#' @importFrom rlang .data
#' @export
summary.lgcpReal <- function(object,
                             linear=TRUE,
                             plot=TRUE,
                             verbose=TRUE,
                             ...){
  rown <- c("Mean","SD","2.5%","10%","25%","50%","75%","90%","97.5%")
  labs <- c("(Intercept)")
  if(!is.null(object$formulae$form.sp)){
    labs <- c(labs,attr(terms(object$formulae$form.sp),"term.labels"))
  }
  if(!is.null(object$formulae$form.t)){
    labs.t <- attr(terms(object$formulae$form.t),"term.labels")
    labs <- c(labs,levels(factor(object$data.t[,labs.t]))[-1])
  }
  #Posteriors for bet
  ans <- do.call(data.frame,
                 list(
                   mean = colMeans(object$beta),
                   SD = apply(object$beta,2,sd),
                   q = t(apply(object$beta,2,function(i)quantile(i,c(0.025,0.1,0.25,0.5,0.75,0.9,0.975))))
                 ))
  rownames(ans) <- labs
  colnames(ans) <- rown
  if(verbose){
    cat("Summary for model\n")
    hrs <- floor(max(object$lgcpRunInfo$timetaken))
    mins <- round((max(object$lgcpRunInfo$timetaken)%%1)*60,0)
    cat("Running time: ",hrs," hours ",mins," minutes\n")
    cat("Number of chains: ",nrow(object$lgcpRunInfo$timetaken),"\n")
  }

  # cat("\n------------------------------------------------------------------\n
  #        Posterior samples:\n")

  eta.labs <- c("Sigma^2","Spatial range","Temporal range")
  ans.eta <- do.call(data.frame,
                     list(
                       mean = colMeans(object$eta),
                       SD = apply(object$eta,2,sd),
                       q = t(apply(object$eta,2,function(i)quantile(i,c(0.025,0.1,0.25,0.5,0.75,0.9,0.975))))
                     ))
  rownames(ans.eta) <- eta.labs
  colnames(ans.eta) <- rown
  #ans <- round(ans,3)
  #ans[nrow(ans)+1,] <- NA
  #rownames(ans)[nrow(ans)] <- ""
  if(linear){
    #print(knitr::kable(rbind(ans,ans.eta),"simple",digits=3,options=list(knitr.kable.NA="")))
    output <- rbind(ans,ans.eta)
  } else {
    #print(knitr::kable(exp(rbind(ans,ans.eta)),"simple",digits=3,options=list(knitr.kable.NA="")))
    output <- exp(rbind(ans,ans.eta))
  }

  conv.res <- convergence(object,plots=FALSE)
  output <- cbind(output,conv.res)

    ans.prior <- do.call(data.frame,
                         list(
                           mean = object$lgcpRunInfo$priors$betaprior$mean,
                           SD = sqrt(diag(object$lgcpRunInfo$priors$betaprior$variance)),
                           q = t(sapply(1:length(object$lgcpRunInfo$priors$betaprior$mean),
                                        function(i)
                                          qnorm(c(0.025,0.1,0.25,0.5,0.75,0.9,0.975),
                                                object$lgcpRunInfo$priors$betaprior$mean[i],
                                                sqrt(diag(object$lgcpRunInfo$priors$betaprior$variance))[i])))
                         ))

    rownames(ans.prior) <- labs
    colnames(ans.prior) <- rown

    ans.prior.eta <- do.call(data.frame,
                             list(
                               mean = object$lgcpRunInfo$priors$etaprior$mean,
                               SD = sqrt(diag(object$lgcpRunInfo$priors$etaprior$variance)),
                               q = t(sapply(1:length(object$lgcpRunInfo$priors$etaprior$mean),
                                            function(i)
                                              qnorm(c(0.025,0.1,0.25,0.5,0.75,0.9,0.975),
                                                    object$lgcpRunInfo$priors$etaprior$mean[i],
                                                    sqrt(diag(object$lgcpRunInfo$priors$etaprior$variance))[i])))
                             ))
    rownames(ans.prior.eta) <- eta.labs
    colnames(ans.prior.eta) <- rown
    #ans.prior <- round(ans.prior,3)
    #ans.prior[nrow(ans.prior)+1,] <- NA
    #rownames(ans.prior)[nrow(ans.prior)] <- ""

    # cat("\n-------------------------------------------------------------\n
    #    Prior distributions:\n")

    if(linear){
      #print(knitr::kable(rbind(ans.prior,ans.prior.eta),"simple",digits=3,options=list(knitr.kable.NA="")))
      output.prior <- rbind(ans.prior,ans.prior.eta)

    } else {
      #print(knitr::kable(exp(rbind(ans.prior,ans.prior.eta)),"simple",digits=3,options=list(knitr.kable.NA="")))
      output.prior <- exp(rbind(ans.prior,ans.prior.eta))
    }



  if(plot){
    niter <- nrow(object$beta)
    betadf <- reshape2::melt(object$beta)
    betadf$post <- "Posterior"
    betadf2 <- reshape2::melt(sapply(1:ncol(object$beta),
                                     function(i)rnorm(niter,
                                                      object$lgcpRunInfo$priors$betaprior$mean[i],
                                                      sqrt(diag(object$lgcpRunInfo$priors$betaprior$variance))[i])))
    betadf2$post <- "Prior"
    betadf <- suppressWarnings(rbind(betadf,betadf2))
    betadf$varname <- labs[betadf$Var2]
    if(!linear){
      betadf$value <- exp(betadf$value)
    }

    p.beta <- ggplot(data=betadf,aes(x=.data$value,lty=.data$post,color=.data$post))+
      geom_density()+
      facet_wrap(~varname,scales ="free")+
      theme_bw()+
      theme(panel.grid=element_blank())+
      scale_linetype_discrete(name="")+
      scale_color_discrete(name="")+
      ggtitle("Beta parameters")

    etadf <- reshape2::melt(object$eta)
    etadf$post <- "Posterior"
    etadf2 <- reshape2::melt(sapply(1:ncol(object$eta),
                                    function(i)rnorm(niter,
                                                     object$lgcpRunInfo$priors$etaprior$mean[i],
                                                     sqrt(diag(object$lgcpRunInfo$priors$etaprior$variance))[i])))
    etadf2$post <- "Prior"
    etadf <- suppressWarnings(rbind(etadf,etadf2))
    etadf$varname <- eta.labs[etadf$Var2]
    if(!linear){
      etadf$value <- exp(etadf$value)
    }

    p.eta <- ggplot(data=etadf,aes(x=.data$value,lty=.data$post,color=.data$post))+
      geom_density()+
      facet_wrap(~varname,scales="free")+
      theme_bw()+
      theme(panel.grid=element_blank())+
      scale_linetype_discrete(name="")+
      scale_color_discrete(name="")+
      ggtitle("Eta parameters")

    print(ggpubr::ggarrange(p.beta,p.eta,ncol=1))
  }

    res <- list(posterior = output,
                prior = output.prior)
    class(res) <- "lgcpRealSumm"
    res

}


#' Aggregate lgcp output to larger geography
#'
#' Take a lgcpRealPlot object and aggregates the output to a larger geography specified by a \code{spatialPolygons} object.
#'
#' This function provides a way of producing aggregated model output for larger geographies. The model fitting takes place on
#' a fine regular lattice, this function provides a way to aggregate this to non-regular polygons such as administrative or
#' political boundaries.
#'
#' @param obj An lgcpRealPlot produced by \code{plot} or \code{plot_hotspot}. NOTE: the call \code{plot} or \code{plot_hotspot} must have
#' had \code{osm=FALSE} set to work with this function.
#' @param aggpoly A \code{spatialPolygons} or \code{spatialPolygonsDataFrame} object specifying the geography to aggregate to.
#' @param osm A logical value indicating whether to overlay the plot on an OpenStreetMap map
#' @return An lgcpRealPlot object comprising a list of two ggplot objects.
#' @importFrom methods as is
#' @importFrom utils flush.console
#' @importFrom rlang .data
#' @export
aggregator <- function(obj,
                       aggpoly,
                       osm=FALSE){

  if(!is(obj,"lgcpRealPlot"))stop("obj should be an lgcpRealPlot")
  if((!is(aggpoly,"SpatialPolygons")|!is(aggpoly,"SpatialPolygonsDataFrame")))
    stop("aggpoly must be of class SpatialPolygons or SpatialPolygonsDataFrame")
  if(nrow(obj[[1]]$data)<10)stop("Please replot without the osm option and add the OSM option to
                                 this function.")

  dat1 <- obj[[1]]$data
  dat2 <- dat1[,c('x','y')]
  sp::coordinates(dat2) <- ~ x+y
  dat2 <- as(dat2,"SpatialPixels")
  sp::proj4string(aggpoly) <- sp::CRS("+init=epsg:4326")
  sp::proj4string(dat2) <- sp::proj4string(aggpoly)
  dat2 <- as(dat2,"SpatialPolygons")
  dat_area <- dat2@polygons[[1]]@area


  if(is(aggpoly,"SpatialPolygonsDataFrame")){
    map <- as(aggpoly,"SpatialPolygons")
  } else {
    map <- aggpoly
  }

  dat1$ID <- unlist(lapply(dat2@polygons,function(i)return(i@ID)))
  dataagg <- data.frame(ID=unlist(lapply(map@polygons,function(i)return(i@ID))),
                        area_hect = geosphere::areaPolygon(aggpoly)/10000,
                        value = NA,pop=NA,poppred=NA,linpred=NA,xpred=NA,sd=NA,value=NA)

  if(attr(obj,"type")=="main"){
    extr_vars <- c('poppred','linpred','xpred','sd','value')
  } else if(attr(obj,"type")=="hotspot"){
    if(length(attr(obj,"str"))==1){
      extr_vars <- c("class_prop")
    } else {
      extr_vars <- attr(obj,"labs")
    }
  }

  cat("\nAggregating results:\n")
  for(i in 1:nrow(dataagg)){
    map_int <- rgeos::gIntersection(dat2,map[i],byid = TRUE, drop_not_poly = TRUE)
    if(!is.null(map_int)&length(map_int)>0){
      map_df <- data.frame(
        ID=unlist(lapply(map_int@polygons,function(i)return(i@ID))),
        area = unlist(lapply(map_int@polygons,function(i)return(i@area))),
        stringsAsFactors = FALSE
      )
      map_df_id <- do.call(rbind,strsplit(map_df$ID,split = " "))
      map_df$dat_id <- map_df_id[,1]
      map_df$map_id <- map_df_id[,2]
      map_df$area_prop <- map_df$area/dat_area
      map_df$pop <- dat1[match(map_df$dat_id,dat1$ID),'pop']
      for(var in extr_vars){
        map_df$tmp <- dat1[match(map_df$dat_id,dat1$ID),var]
        if(var%in%c('linpred','xpred')){
          dataagg[i,var] <- exp(with(map_df[!is.na(map_df$tmp),],
                                     weighted.mean(log(tmp),pop*area_prop)))
        } else if(var=="poppred"){
          dataagg[i,var] <- sum(map_df$tmp*map_df$area_prop,na.rm = TRUE)
        } else {
          dataagg[i,var] <- with(map_df[!is.na(map_df$tmp),],
                                 weighted.mean(tmp,pop*area_prop))
        }
      }
    }

    cat("\r|",rep("=",floor((i/nrow(dataagg))/0.05)),
        rep(" ",ceiling(((nrow(dataagg)-i)/nrow(dataagg))/0.05)),"| ",
        floor((i/nrow(dataagg))/0.05)*5,"%",sep="");flush.console()
  }
  rownames(dataagg) <- dataagg$ID
  res <- sp::SpatialPolygonsDataFrame(map,dataagg)

  if(attr(obj,"type")=="main"){
    rr_max <- attr(obj,"rr_lim")[1]
    col_lim <- c(0.01,rr_max)
    col_vals <- c(0,1/(rr_max*2),seq(1/rr_max,1,length.out = 4))

    ppred <- ggplot()+
      ggspatial::layer_spatial(data=res,aes(fill=.data$value))+
      theme_bw()+
      theme(panel.grid = element_blank())+
      ggtitle(obj[[1]]$labels$title)

    px <- ggplot()+
      ggspatial::layer_spatial(data=res,aes(fill=.data$xpred))+
      scale_fill_gradientn(name="RR",colours = c("purple","blue","yellow","orange","red","brown"),
                           values = col_vals,limits=col_lim)+
      theme_bw()+
      theme(panel.grid = element_blank())+
      ggtitle("Latent")

    plin <- ggplot()+
      ggspatial::layer_spatial(data=res,aes(fill=.data$linpred))+
      scale_fill_gradientn(name="RR",colours = c("purple","blue","yellow","orange","red","brown"),
                           values = col_vals,limits=col_lim)+
      theme_bw()+
      theme(panel.grid = element_blank())+
      ggtitle("Observed")

    ppop <- ggplot()+
      ggspatial::layer_spatial(data=res,aes(fill=.data$poppred))+
      scale_fill_viridis_c()+
      theme_bw()+
      theme(panel.grid = element_blank())+
      ggtitle("Expected")

    psd <- ggplot()+
      ggspatial::layer_spatial(data=res,aes(fill=.data$xpred))+
      scale_fill_viridis_c()+
      theme_bw()+
      theme(panel.grid = element_blank())+
      ggtitle("Posterior SD")

    if(osm){
      if(requireNamespace("ggmap",quietly=TRUE)){
        xrange <- range(obj[[1]]$data$x)
        yrange <- range(obj[[1]]$data$y)
        #our background map
        mad_map <- get_map2(c(left=xrange[1],bottom=yrange[1],right=xrange[2],top=yrange[2]),
                            source="stamen",
                            maptype = "toner")
        ppred <- ggmap::ggmap(mad_map) +
          ggspatial::layer_spatial(data=res,aes(fill=.data$value),alpha=0.4)+
          theme_bw()+
          theme(panel.grid = element_blank())+
          ggtitle(obj[[1]]$labels$title)

      } else {
        stop("ggmap required for osm plotting.")
      }

    }

    if(grepl("Change",obj[[1]]$labels$title)){
      rr_max <- attr(obj,"rr_lim")[2]
      col_lim <- c(0.01,rr_max)
      col_vals <- c(0,1/(rr_max*2),seq(1/rr_max,1,length.out=4))
      ppred <- ppred + scale_fill_gradientn(name="IRR",
                                            colours = c("purple","blue","yellow","orange","red","brown"),
                                            values = col_vals,limits=col_lim)
    } else {
      ppred <- ppred + scale_fill_viridis_c(name="")
    }

    prow <- ggpubr::ggarrange(ppop,plin,px,psd,nrow=2,ncol=2)
    out <- list(ppred,prow)
    #print(ggpubr::ggarrange(ppred,prow,nrow=1,ncol=2))
  } else if(attr(obj,"type")=="hotspot"){
    if(length(attr(obj,"str"))==1){

      res@data$class <- ifelse(res@data$class_prop>attr(obj,"threshold"),attr(obj,"labs")[2],attr(obj,"labs")[1])
      res@data$class <- factor(res@data$class,levels=attr(obj,"labs"),ordered = TRUE)

      pclass <- ggplot()+
        ggspatial::layer_spatial(data=res,aes(fill=.data$class))+
        scale_fill_viridis_d()+
        theme_bw()+
        theme(panel.grid = element_blank())+
        ggtitle(paste0("Classification if Pr(",attr(obj,"str"),">",attr(obj,"vals"),") > ",round(attr(obj,"threshold")*100,0),"%"))

      ppred <- ggplot()+
        ggspatial::layer_spatial(data=res,aes(fill=.data$class_prop))+
        scale_fill_gradientn(name="Prob",colours = c("purple","blue","green","yellow","red","brown"),
                             values = seq(0,1,length.out=6),limits=c(0,1))+
        theme_bw()+
        theme(panel.grid = element_blank())+
        ggtitle("Probability")

    } else {
      plist <- list()
      nlabs <- length(attr(obj,"labs"))
      res@data$class <- NA
      res@data$class[!is.na(res@data[,attr(obj,"labs")[1]])] <- attr(obj,"labs")[unname(unlist(apply(res@data[,attr(obj,"labs")],1,which.max)))]
      res@data$class <- factor(res@data$class,levels=attr(obj,"labs"),ordered = TRUE)


      pclass <- ggplot()+
        ggspatial::layer_spatial(data=res,aes(fill=.data$class))+
        scale_fill_viridis_d()+
        theme_bw()+
        theme(panel.grid = element_blank())+
        ggtitle("Classification based on modal probability.")


      for(i in 1:nlabs){
        res@data$tmp <- res@data[,attr(obj,"labs")[i]]
        plist[[i]] <- ggplot()+
          ggspatial::layer_spatial(data=res,aes(fill=.data$tmp))+
          scale_fill_gradientn(name="Prob",colours = c("purple","blue","yellow","orange","red","brown"),
                               values = seq(0,1,length.out=6),limits=c(0,1))+
          theme_bw()+
          theme(panel.grid = element_blank())+
          ggtitle(attr(obj,"labs")[i])
      }
      if(nlabs==2){
        ppred <- ggpubr::ggarrange(plist[[1]],plist[[2]])
      } else if(nlabs==3){
        ppred <- ggpubr::ggarrange(plist[[1]],plist[[2]],plist[[3]])
      } else {
        ppred <- ggpubr::ggarrange(plist[[1]],plist[[2]],plist[[3]],plist[[4]])
      }

    }

    if(osm){
      if(requireNamespace("ggmap",quietly = TRUE)){
        xrange <- range(obj[[1]]$data$x)
        yrange <- range(obj[[1]]$data$y)
        #our background map
        mad_map <- get_map2(c(left=xrange[1],bottom=yrange[1],right=xrange[2],top=yrange[2]),
                            source="stamen",
                            maptype = "toner")
        pclass <- ggmap::ggmap(mad_map) +
          ggspatial::layer_spatial(data=res,aes(fill=.data$class),alpha=0.4)+
          scale_fill_viridis_d()+
          theme_bw()+
          theme(panel.grid = element_blank())+
          ggtitle(paste0("Classification if Pr(",attr(obj,"str"),">",attr(obj,"vals"),") > ",round(attr(obj,"threshold")*100,0),"%"))

      } else {
        stop("ggmap required for osm plotting.")
      }

    }

    out <- list(pclass,ppred)
    #print(ggpubr::ggarrange(pclass,ppred,nrow=1))
  }

  class(out) <- "lgcpRealPlot"
  out
  #return(invisible(out))

}

#' Alternative \code{get_map} function
#'
#' An alternative to \code{ggmap}'s \code{get_map} function
#'
#' An error in the CRAN available version of get_map means it will only plot Google Maps objects rather than Stamen or OSM maps. When ggmap is
#' updated this function will be removed. See \code{help(get_map)} for more details.
#'
#' @param location an address, longitude/latitude pair (in that order), or
#'   left/bottom/right/top bounding box
#' @param zoom map zoom, an integer from 3 (continent) to 21 (building), default
#'   value 10 (city).  openstreetmaps limits a zoom of 18, and the limit on
#'   stamen maps depends on the maptype.  "auto" automatically determines the
#'   zoom for bounding box specifications, and is defaulted to 10 with
#'   center/zoom specifications.  maps of the whole world currently not
#'   supported.
#' @param scale scale argument of get_googlemap() or get_openstreetmap()
#' @param maptype character string providing map theme. options available are
#'   "terrain", "terrain-background", "satellite", "roadmap", and "hybrid"
#'   (google maps), "terrain", "watercolor", and "toner" (stamen maps)
#' @param source Google Maps ("google"), OpenStreetMap ("osm"), Stamen Maps
#'   ("stamen")
#' @param force force new map (don't use archived version)
#' @param messaging turn messaging on/off
#' @param urlonly return url only
#' @param filename destination file for download (file extension added according
#'   to format). Default \code{NULL} means a random tempfile().
#' @param crop (stamen and cloudmade maps) crop tiles to bounding box
#' @param color color ("color") or black-and-white ("bw")
#' @param language language for google maps
#' @param ... ...
#' @return a ggmap object (a classed raster object with a bounding box
#'   attribute)
#' @export
get_map2 <- function (location = c(lon = -95.3632715, lat = 29.7632836),
                      zoom = "auto",
                      scale = "auto",
                      maptype = c("terrain","terrain-background", "satellite", "roadmap",
                                  "hybrid", "toner", "watercolor", "terrain-labels",
                                  "terrain-lines", "toner-2010", "toner-2011",
                                  "toner-background", "toner-hybrid", "toner-labels",
                                  "toner-lines", "toner-lite"),
                      source = c("google","osm", "stamen"),
                      force = ifelse(source == "google", TRUE, FALSE),
                      messaging = FALSE,
                      urlonly = FALSE,
                      filename = NULL,
                      crop = TRUE,
                      color = c("color", "bw"),
                      language = "en-EN", ...)
{
  if(requireNamespace("ggmap",quietly = TRUE)){
    args <- as.list(match.call(expand.dots = TRUE)[-1])
    if ("verbose" %in% names(args)) {
      .Deprecated(msg = "verbose argument deprecated, use messaging.")
      messaging <- eval(args$verbose)
    }
    if ("center" %in% names(args)) {
      .Deprecated(msg = "center argument deprecated, use location.")
      location <- eval(args$center)
    }
    source <- match.arg(source)
    color <- match.arg(color)
    if (missing(maptype)) {
      if (source != "cloudmade") {
        maptype <- "terrain"
      }
      else {
        maptype <- 1
      }
    }
    if (source == "stamen") {
      if (!(maptype %in% c("terrain", "terrain-background",
                           "terrain-labels", "terrain-lines", "toner",
                           "toner-2010", "toner-2011", "toner-background",
                           "toner-hybrid", "toner-labels", "toner-lines",
                           "toner-lite", "watercolor"))) {
        stop("invalid stamen maptype, see ?get_stamenmap",
             call. = FALSE)
      }
    }
    if (source == "google" & (maptype %in% c("terrain-background",
                                             "terrain-labels", "terrain-lines", "toner",
                                             "toner-2010", "toner-2011", "toner-background",
                                             "toner-hybrid", "toner-labels", "toner-lines",
                                             "toner-lite", "watercolor"))) {
      message(paste0("maptype = \"", maptype, "\" is only available with source = \"stamen\"."))
      message(paste0("resetting to source = \"stamen\"..."))
      source <- "stamen"
    }
    location_stop <- TRUE
    if (is.character(location) && length(location) == 1) {
      location_type <- "address"
      location_stop <- FALSE
    }
    if (is.data.frame(location) && ncol(location) == 2) {
      location <- colMeans(location)
    }
    if (is.numeric(location) && length(location) == 2) {
      location_type <- "lonlat"
      location_stop <- FALSE
      if (!is.null(names(location))) {
        loc_names <- names(location)
        if (all(loc_names == c("long", "lat"))) {
          names(location) <- c("lon", "lat")
        }
        else if (all(loc_names == c("lat", "lon"))) {
          message("note : locations should be specified in the lon/lat format, not lat/lon.")
          location <- location[c("lon", "lat")]
        }
        else if (all(loc_names == c("lat", "long"))) {
          message("note : locations should be specified in the lon/lat format, not lat/lon.")
          location <- location[c("long", "lat")]
          names(location) <- c("lon", "lat")
        }
      }
      else {
        names(location) <- c("lon", "lat")
      }
    }
    if (is.numeric(location) && length(location) == 4) {
      location_type <- "bbox"
      location_stop <- FALSE
      #source <- "stamen"
      #maptype <- "terrain"
      if (length(names(location)) > 0) {
        if (!all(names(location) %in% c("left", "bottom",
                                        "right", "top"))) {
          stop("bounding boxes should have name left, bottom, right, top)",
               call. = FALSE)
        }
        location <- location[c("left", "bottom",
                               "right", "top")]
      }
      else {
        names(location) <- c("left", "bottom",
                             "right", "top")
      }
    }
    if (location_stop) {
      stop("improper location specification, see ?get_map.",
           call. = F)
    }
    if (zoom == "auto" && location_type == "bbox") {
      if (zoom == "auto") {
        lon_range <- location[c("left", "right")]
        lat_range <- location[c("bottom", "top")]
        if (missing(zoom)) {
          lonlength <- diff(lon_range)
          latlength <- diff(lat_range)
          zoomlon <- ceiling(log2(360 * 2/lonlength))
          zoomlat <- ceiling(log2(180 * 2/latlength))
          zoom <- max(zoomlon, zoomlat)
        }
      }
    }
    else if (zoom == "auto" && location_type != "bbox") {
      zoom = 10
    }
    if (scale == "auto") {
      if (source == "google")
        scale <- 2
      if (source == "osm")
        scale <- ggmap::OSM_scale_lookup(zoom)
    }
    if (source == "google") {
      if (location_type == "bbox") {
        warning("bounding box given to google - spatial extent only approximate.",
                call. = FALSE, immediate. = TRUE)
        message("converting bounding box to center/zoom specification. (experimental)")
        user_bbox <- location
        location <- c(lon = mean(location[c("left",
                                            "right")]), lat = mean(location[c("bottom",
                                                                              "top")]))
      }
      map <- ggmap::get_googlemap(center = location, zoom = zoom,
                                  maptype = maptype, scale = scale, messaging = messaging,
                                  urlonly = urlonly, force = force, filename = filename,
                                  color = color, language = language)
      if (FALSE) {
        bb <- attr(map, "bb")
        mbbox <- c(left = bb$ll.lon, bottom = bb$ll.lat,
                   right = bb$ur.lon, top = bb$ur.lat)
        size <- dim(map)
        if (location_type == "bbox") {
          slon <- seq(mbbox["left"], mbbox["right"],
                      length.out = size[1])
          slat <- seq(mbbox["top"], mbbox["bottom"],
                      length.out = size[2])
          keep_x_ndcs <- which(user_bbox["left"] <=
                                 slon & slon <= user_bbox["right"])
          keep_y_ndcs <- which(user_bbox["bottom"] <=
                                 slat & slat <= user_bbox["top"])
          map <- map[keep_y_ndcs, keep_x_ndcs]
          class(map) <- c("ggmap", "raster")
          attr(map, "bb") <- data.frame(ll.lat = user_bbox["bottom"],
                                        ll.lon = user_bbox["left"], ur.lat = user_bbox["top"],
                                        ur.lon = user_bbox["right"])
        }
      }
      return(map)
    }
    if (source == "osm") {
      if (location_type != "bbox") {
        gm <- ggmap::get_googlemap(center = location, zoom = zoom,
                                   filename = filename)
        location <- as.numeric(attr(gm, "bb"))[c(2,
                                                 1, 4, 3)]
      }
      return(ggmap::get_openstreetmap(bbox = location, scale = scale,
                                      messaging = messaging, urlonly = urlonly, filename = filename,
                                      color = color))
    }
    if (source == "stamen") {
      if (location_type != "bbox") {
        gm <- ggmap::get_googlemap(center = location, zoom = zoom,
                                   filename = filename)
        location <- as.numeric(attr(gm, "bb"))[c(2,
                                                 1, 4, 3)]
      }
      return(ggmap::get_stamenmap(bbox = location, zoom = zoom, maptype = maptype,
                                  crop = crop, messaging = messaging, urlonly = urlonly,
                                  filename = filename, force = force, color = color))
    }
  } else {
    stop("ggmap package required for osm plotting.")
  }

}
