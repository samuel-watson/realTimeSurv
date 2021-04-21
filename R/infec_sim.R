#' Infection simulation
#'
#' Simulate data from an infectious disease process modelled as a Poisson process
#'
#' 1.	Background incidence is a Poisson process with intensity lambda x (pop density)
#' 2.	A small fraction p of background events trigger a temporary increase in local
#' intensity to lambda x (pop density) x (1+rho) in a disc of  radius delta for the next k time-periods.
#'
#' @param region A \code{spatialPolygon} defining the area to simulate data for
#' @param t.win A vector indicating the time window to simulate data for, eg. c(1, 30)
#' @param covariates A \code{spatialPolygonsDataFrame} containing population density information for the area of interest. The population
#' density variable should be called \code{popdens}.
#' @param mean.val Integer, the mean number of cases per time period
#' @param p Probability a case generates a cluster
#' @param delta A vector of two values: the spatial range and temporal range parameters
#' @param rho Multiplicative factor of effect of a cluster
#' @param beta Parameters of latent field covariates
#' @param t.off Temporal offset parameter - number of time periods to displace the peak of the cluster intensity.
#' @param cov.pars Covariance parameters of latent Gaussian field covariates: partial sill, range, and nugget
#' @param grid.size Size of the computational grid to simulate the infectious process
#' @return A list: (1) simulated data for each computational grid cell, (2) simulated case locations and time, (3) plot
#' of Poisson intensity, (4) plot of simulated case locations, (5) stpp object of case locations for use with stpp functions.
#' @export
infecSim <- function(region,
                     t.win,
                     covariates,
                     mean.val,
                     p = 1/mean.val,
                     delta,
                     rho,
                     beta=c(1,1),
                     t.off=0,
                     cov.pars=c(1,0.0075,0.5),
                     grid.size=64^2){

  region.df = ggplot2::fortify(region)

  grid <- sp::makegrid(region,n=grid.size)
  x.size <- unique(diff(unique(grid[,1])))[1]
  spgrd <- sp::SpatialPoints(grid, proj4string = sp::CRS(sp::proj4string(region)))
  spgrd <- spgrd[region,]

  sp::proj4string(spgrd) <- sp::proj4string(covariates)
  pop <- sp::over(spgrd,covariates)
  pop <- pop$popdens
  pop[is.na(pop)] <- 0

  spgrd <- matrix(c(spgrd$x1,spgrd$x2),ncol=2)

  cov1 <- geoR::grf(n=nrow(spgrd)*2,
              grid=spgrd,
              borders = as.matrix(region.df[,c('long','lat')],ncol=2),
              nx=100,
              ny=100,
              cov.model = "exponential",
              cov.pars = cov.pars[1:2],
              nugget = cov.pars[3])

  cov2 <- geoR::grf(n=nrow(spgrd)*2,
              grid=spgrd,
              borders = as.matrix(region.df[,c('long','lat')],ncol=2),
              nx=100,
              ny=100,
              cov.model = "exponential",
              cov.pars = cov.pars[1:2],
              nugget = cov.pars[3])

  data <- cbind(as.data.frame(spgrd), cov=cov1$data, cov2=cov2$data, pop=pop)
  colnames(data)[1:2] <- c('x','y')
  #intensity in each cell
  imat <- matrix(data$pop*exp(beta[1]*data$cov - beta[1]^2*(cov.pars[1]+cov.pars[3])/2 +
                   beta[2]*data$cov2 - beta[2]^2*(cov.pars[1]+cov.pars[3])/2)*
                   mean.val*(1-p)/sum(data$pop),
                 nrow=nrow(data),ncol=length(t.win[1]:t.win[2])+1)

  imat_in <- matrix(data$pop*exp(beta[1]*data$cov - beta[1]^2*(cov.pars[1]+cov.pars[3])/2 +
                                   beta[2]*data$cov2 - beta[2]^2*(cov.pars[1]+cov.pars[3])/2)*mean.val*p/sum(data$pop),
                 nrow=nrow(data),ncol=length(t.win[1]:t.win[2])+1)

  # for each time period,1. generate cases, 2. with probability p case raises intensity

  print("Generating simulation")
  for(t in 1:length(t.win[1]:t.win[2])){
    dt1 <- data
    dt1$t <- t
    dt1$n <- rpois(nrow(data),imat[,t])
    dt1$n_in <- rpois(nrow(data),imat_in[,t]) #infectious parents
    # probability of spreading
    if(any(dt1$n_in>0)){
        idx <- which(dt1$n_in>0)
        for(i in idx){
          dis1 <- sqrt((data$x - data$x[i])^2 + (data$y - data$y[i])^2)
          #print(cat("T: ",t))
          for(s in (t+1):length(t.win[1]:t.win[2])){
            imat[,s] <- imat[,s]*(1+rho*exp(-dis1/delta[1])*exp(-((s-t.off-t)^2)/delta[2]))
          }
      }
    }
    dt1$intens <- imat[,t] + imat_in[,t]
    dt1$n_tot <- dt1$n+dt1$n_in
    if(t==1){
      dat <- dt1
    } else {
      dat <- rbind(dat,dt1)
    }
  }

  #points locations
  idx <- which(dat$n>0)
  pts <- dat[rep(which(dat$n>0),dat[dat$n>0,'n']),]
  pts$x <- pts$x+runif(nrow(pts),-x.size/2,x.size/2)
  pts$y <- pts$y+runif(nrow(pts),-x.size/2,x.size/2)

  # pl1 <- ggplot(data=dat,aes(x=x,y=y,fill=log(intens)))+
  #   geom_raster()+
  #   scale_fill_viridis_c()+
  #   facet_wrap(~t)+
  #   theme_bw()+
  #   theme(panel.grid = element_blank())
  #
  # pl2 <- ggplot(data=pts,aes(x=x,y=y))+
  #   geom_point(color="red",size=0.1)+
  #   geom_path(data=region.df,aes(x=long,y=lat))+
  #   scale_fill_viridis_c(name="Log\nintensity")+
  #   facet_wrap(~t)+
  #   theme_bw()+
  #   theme(panel.grid = element_blank())

  # datxyt <- as.matrix(pts[,c('x','y','t')])
  # attr(datxyt,"class") <- "stpp"

  #return(list(dat,pts[,c('x','y','t')],pl1,pl2,datxyt))

  return(pts[,c('x','y','t')])
}


