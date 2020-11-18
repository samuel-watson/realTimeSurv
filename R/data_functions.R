#' Create spatial polygons data frame
#'
#' Aggregate to larger geography and output data
#'
#' This function aggregates the output of plot or plot_hotspot to a
#' larger geography and outputs a spatial polygons data frame that
#' can be used to create other plots.
#'
#' @param obj A lgcpReal object from the lgcp function
#' @param aggpoly A \code{spatialPolygons} or \code{spatialPolygonsDataFrame} describing the
#' geography to aggregate to.
#' @return A \code{spatialPolygonsDataFrame}
#' @export
aggregator_data <- function(obj,aggpoly,osm=FALSE){
  if(!(class(aggpoly)=="SpatialPolygons"|class(aggpoly)=="SpatialPolygonsDataFrame"))
    stop("aggpoly must be of class ''SpatialPolygons or SpatialPolygonsDataFrame")
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


  if(class(aggpoly)=="SpatialPolygonsDataFrame"){
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





    if(grepl("Change",obj[[1]]$labels$title)){
      rr_max <- attr(obj,"rr_lim")[2]
      col_lim <- c(0.01,rr_max)
      col_vals <- c(0,1/(rr_max*2),seq(1/rr_max,1,length.out=4))

    }


  }
  if(attr(obj,"type")=="hotspot"){
    if(length(attr(obj,"str"))==1){

      res@data$class <- ifelse(res@data$class_prop>attr(obj,"threshold"),attr(obj,"labs")[2],attr(obj,"labs")[1])
      res@data$class <- factor(res@data$class,levels=attr(obj,"labs"),ordered = TRUE)



    } else {
      plist <- list()
      nlabs <- length(attr(obj,"labs"))
      res@data$class <- NA
      res@data$class[!is.na(res@data[,attr(obj,"labs")[1]])] <- attr(obj,"labs")[unname(unlist(apply(res@data[,attr(obj,"labs")],1,which.max)))]
      res@data$class <- factor(res@data$class,levels=attr(obj,"labs"),ordered = TRUE)





      for(i in 1:nlabs){
        res@data$tmp <- res@data[,attr(obj,"labs")[i]]
        names(res)[length(names(res))] <- attr(obj,"labs")[i]

      }


    }



  }


  if("class"%in%names(res)){
    if("class_prop"%in%names(res)){
      res <- res[,c("class","class_prop")]
      names(res) <- c("Classification","Probabilities")
    } else {
      colns <- which(names(res)=="class")
      res <- res[,names(res)[colns:length(names(res))]]
    }

  } else {
    res <- res[,c("value","poppred","linpred","xpred")]
    names(res) <- c("Incidence","Expected","Observed_rr","Latent_rr")
  }

  return(res)

}

#' Create spatial pixels data frame
#'
#' Output spatial pixels data frame of lgcp predictions
#'
#' This function generates a spatial pixels data frame with the output
#' from lgcp.
#'
#' @param lg An lgcpReal object, output from a call to \code{lgcp}
#' @param covariates A \code{spatialPolygonsDataFrame} covering the area of interest and containing
#' the covariate and population density data. Typically the same object as specified in the
#' \code{covariates} argument in the call to \code{lgcp}.
#' @param per.days Integer, the number of person-days to use for incidence, default is 10,000.
#' @param change.lag If not NULL, then plots are created of the change in outputs compared to
#' this number of periods prior.
#' @param relative A logical value indicating whether the comparisons (if change.lag set) should be relative (default),
#' i.e. incidence rate ratios and ratios of relative risks, or absolute.
#' @param msq Integer, the denominator of the population density, default is hectares (population per
#' 10,000m^2)
#' @param rr_lim Integer, for plotting the relative risk, the maximum value of the colour scale. Useful
#' when comparing multiple plots to put colour gradient on same scale.
#' @return A code{spatialPolygonsDataFrame}
#' @seealso aggregator_data, plot_hotspot
#' @export
plot_data <- function(lg,
                      covariates,
                      per.days=10000,
                      change.lag=NULL,
                      relative=TRUE,
                      msq = 10000,
                      rr_lim=NULL){
  if(missing(covariates))stop("Specify the covariate spatial polygons data frame")
  if(class(lg)!="lgcpReal")stop("lg must be of class lgcpReal")
  OW <- lgcp::selectObsWindow(lg$xyt, cellwidth = lg$cellwidth)

  grid.data <- expand.grid(x=OW$xvals,y=OW$yvals)
  idx.mapping <- matrix(1:nrow(grid.data),nrow=length(OW$yvals),ncol=length(OW$xvals))
  idx.mapping <- c(t(apply(idx.mapping,2,rev)))

  if(!exists("outl") |(exists("outl")&&attr(outl, "dirname")!=lg$dirname)){
    print("Extracting posterior samples...")
    outl <- lgcpExtract(lg$dirname,nrow(lg$lgcpRunInfo$timetaken))
    assign("outl",outl,.GlobalEnv)
  }

  res1 <- suppressWarnings( realTimeSurv:::.plot_lgcp_dat(outl,
                                                          grid.data,
                                                          lg,
                                                          lg$nchains,
                                                          idx.mapping,
                                                          covariates = covariates,
                                                          cellwidth = lg$cellwidth))

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
    reslag <- suppressWarnings(.plot_lgcp_dat(outl,grid.data,lg,lg$nchains,
                                              idx.mapping,covariates = covariates,
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
                                           111*1000^2*111.321*lg$cellwidth^2/msq)
    }


  }
  sp::coordinates(res1) <- ~x+y
  res1 <- as(res1,"SpatialPixelsDataFrame")
  sp::proj4string(res1) <- sp::CRS("+init=epsg:4326")
  res1 <- res1[!is.na(res1@data$value),]
  res1 <- res1[,c("value","poppred","linpred","xpred")]
  names(res1) <- c("Incidence","Expected","Observed_rr","Latent_rr")
  return(res1)

}

#' Produce hotspot spatial pixels data
#'
#' Generate a spatial pixels data frame classifying hotspots and probabilities
#'
#' @param lg Output from a call to \code{lgcp}
#' @param covariates A \code{spatialPolygonsDataFrame} covering the area of interest and containing
#' the covariate and population density data. Typically the same object as specified in the
#' \code{covariates} argument in the call to \code{lgcp}.
#' @param threshold.var A vector of one or two strings specifying the variables to define the hotspots,
#' see Details for how to specify.
#' @param threshold.value A vector or one or two values indicating the threshold(s) for determining
#' a hotspot. Given in the same order as threshold.var.
#' @param lables A vector of two or four labels for the hotspots, see Details.
#' @param threshold.prob A vector of one or two values specifying the exceedence probabilities.
#' @param relative A logical value. If one or both of the variable is with respect to a previous time period, whether the comparison
#' should be relative (TRUE) or absolute (FALSE)
#' @param per.days If one or both of the variables is incidence, the denominator number of person-days.
#' @param msq The denominator for the population density in m^2. Default is hectares (per 10,000m^2)
#' @return  A \code{spatialPixelsDataFrame}
#' @export
plot_hotspot_data <- function(lg,
                              covariates,
                              threshold.var=NULL,
                              threshold.value=NULL,
                              labels,
                              threshold.prob=0.8,
                              relative=TRUE,
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

  if(!exists("outl") |(exists("outl")&&attr(outl, "dirname")!=lg$dirname)){
    print("Extracting posterior samples...")
    outl <- lgcpExtract(lg$dirname,nrow(lg$lgcpRunInfo$timetaken))
    assign("outl",outl,.GlobalEnv)
  }

  str1 <- unlist(strsplit(threshold.var[1],"\\+"))

  if(any(grepl("lag",str1))){
    lag1 <- str1[which(grepl("lag",str1))]
    lag1 <- as.numeric(gsub("\\D", "", lag1))
  } else {
    lag1 <- 0
  }

  res1 <- realTimeSurv:::.plot_lgcp_dat(outl,
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
    res2 <- realTimeSurv:::.plot_lgcp_dat(outl,
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

    res3 <- realTimeSurv:::.plot_lgcp_dat(outl,
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
      res4 <- .plot_lgcp_dat(outl,
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




  } else {
    res1$dat1$class_prop <- apply(v0,1,function(i)return(length(i[i > threshold.value[1]])/length(i)))
    res1$dat1$class <- I(res1$dat1$class_prop > threshold.prob)*1
    res1$dat1$class <- factor(res1$dat1$class,levels=c(0,1),labels=labels)



  }
  res1 <- res1$dat1
  sp::coordinates(res1) <- ~x+y
  res1 <- as(res1,"SpatialPixelsDataFrame")
  sp::proj4string(res1) <- sp::CRS("+init=epsg:4326")
  res1 <- res1[!is.na(res1@data$value),]
  res1 <- res1[,names(res1)[(which(names(res1)=="class")):length(names(res1))]]
  return(res1)
}
