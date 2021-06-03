#' Generate geocoded and time-stamped data
#'
#' Use Google Maps API to produce longitude and latitude data
#'
#' A wrapper to functions from package \code{googleway} that checks for errors and
#' fetches longitude and latitude coordinates using the Google Maps API. Dates are converted
#' to integer values with t=1 being the earliest date provided.
#' @param df A data frame containing columns with the address (either "Address" or "address") and date (either "Date" or "date")
#' of the cases.
#' @param api_key A string. A valid Google Developers Geocode API key.
#' @return A data frame with the columns for longitude, latitude, and time.
#' @examples
#' \dontrun{
#' #requires working API to run
#' tmp <- data.frame(address=c("Buckingham palace","Big ben, Westminster","Marble arch, London"),
#'                   date = c("01/01/2020","02/01/2020","03/01/2020"))
#' geocode_st(tmp, api_key = "ENTER_KEY")
#' }
#' @importFrom utils flush.console
#' @export
geocode_st <- function(df,api_key){
  if(!"address"%in%tolower(colnames(df)))stop("Column name for addresses should be 'Address' or 'address'")
  if(!"date"%in%tolower(colnames(df)))stop("Column name for dates should be 'Date' or 'date'")

  addrs <- df[,which(tolower(colnames(df))=="address")]
  addrs <- as.character(addrs)
  dts <- df[,which(tolower(colnames(df))=="date")]
  dts <- lubridate::dmy(dts)
  dts <- as.numeric(dts)
  dts <- dts - (min(dts)-1)

  data <- data.frame(lat=rep(NA,nrow(df)),
                     long=rep(NA,nrow(df)),
                     t=dts)

  message("\nProcessing addresses...\n")
  nmat <- 0
  nnomat <- 0
  for(i in 1:length(addrs)){

    if(addrs[i]!="NULL"){
      add <- tryCatch(googleway::geocode_coordinates(googleway::google_geocode(addrs[i], key=api_key)),
                      error = function(e)NA)
      if(length(add)>1){
        data$lat[i] <- add$lat[1]
        data$long[i] <- add$lng[1]
        nmat <- nmat + 1
      } else {
        nnomat <- nnomat + 1
      }
    }
    message("\rAddresses: ",nmat," matched; ",nnomat," unable to match, of ",length(addrs)," total.");flush.console()
  }

  return(data)
}

#' Get day of the week
#'
#' Returns a data frame with day of the week covariates
#' @param df A data frame containing columns with the address (either "Address" or "address") and date (either "Date" or "date")
#' of the cases. Typically the same data frame used in a call to \code{geocode_st}
#' @return A data frame with value for time period and day of the week.
#' @examples
#' tmp <- data.frame(address=c("Buckingham palace","Big ben, Westminster","Marble arch, London"),
#'                   date = c("01/01/2020","02/01/2020","03/01/2020"))
#' get_day(tmp)
#' @export
get_day <- function(df){
  if(!"date"%in%tolower(colnames(df)))stop("Column name for dates should be 'Date' or 'date'")

  dts <- df[,which(tolower(colnames(df))=="date")]
  dts <- lubridate::dmy(dts)
  dts2 <- as.numeric(dts)
  dts2 <- dts2 - (min(dts2)-1)

  data_t <- data.frame(t = min(dts2):max(dts2),
                       day = lubridate::wday(lubridate::ymd(seq(min(dts),max(dts),by=1)),label=TRUE))

  return(data_t)
}



#' Alternative to lgcpPredictSpatioTemporalPlusPars
#'
#' Internal function and alternative to \code{lgcp::lgcpPredictSpatioTemporalPlusPars}.
#'
#' A copy of \code{lgcp::lgcpPredictSpatioTemporalPlusPars} that parallelises the sampler and
#' produces a reduced output and lgcpReal object. See \code{help(lgcpPredictSpatioTemporalPlusPars)}
#' for more information.
#'
#' @param formula a formula object of the form X ~ var1 + var2 etc. The name of the dependent variable must be "X". Only accepts 'simple' formulae, such as the example given.
#' @param xyt An object of class stppp
#' @param T the time point of interest
#' @param laglength the number of previous time points to include in the analysis
#' @param ZmatList A list of design matrices Z constructed with getZmat and possibly addTemporalCovariates see the details below and Bayesian_lgcp vignette for details on how to construct this.
#' @param model.priors model priors, set using lgcpPrior
#' @param model.inits model initial values. The default is NULL, in which case lgcp will use the prior mean to initialise eta and beta will be initialised from an oversispersed glm fit to the data. Otherwise use lgcpInits to specify.
#' @param spatial.covmodel choice of spatial covariance function. See ?CovFunction
#' @param cellwidth the width of computational cells
#' @param poisson.offset A list of SpatialAtRisk objects (of length the number of types) defining lambda_k (see below)
#' @param mcmc.control MCMC paramters, see ?mcmcpars
#' @param output.control output choice, see ?setoutput
#' @param gradtrunc truncation for gradient vector equal to H parameter Moller et al 1998 pp 473. Default is Inf, which means no gradient truncation, which seems to work in most settings.
#' @param ext integer multiple by which grid should be extended, default is 2. Generally this will not need to be altered, but if the spatial correlation decays slowly, increasing 'ext' may be necessary.
#' @param inclusion criterion for cells being included into observation window. Either 'touching' or 'centroid'. The former, the default, includes all cells that touch the observation window, the latter includes all cells whose centroids are inside the observation window.
#' @importFrom lgcp selectObsWindow spatialAtRisk genFFTgrid fftinterpolate cov.interp.fft
#' @importFrom lgcp mcmcLoop getCounts nullFunction MALAlgcpSpatioTemporal.PlusPars lgcpInits
#' @importFrom lgcp setoutput mcmcProgressTextBar
#' @importFrom spatstat.geom as.polygonal
#' @importFrom utils object.size menu
#' @importFrom stats runif
#' @return an object of class lgcpPredictSpatioTemporalPlusParameters
#' @export
lgcpST <- function (formula, xyt, T, laglength, ZmatList = NULL, model.priors,
                    model.inits = lgcpInits(), spatial.covmodel, cellwidth = NULL,
                    poisson.offset = NULL, mcmc.control, output.control = setoutput(),
                    gradtrunc = Inf, ext = 2, inclusion = "touching")
{
  regionalcovariates <- NULL
  pixelcovariates <- NULL
  nprobe <- 1e+05
  gridsize <- NULL
  starttime <- Sys.time()
  if (!inherits(T, "integer")) {
    warning("Converting T into integer value, see ?as.integer",
            immediate. = TRUE)
    T <- as.integer(T)
  }
  if (!inherits(laglength, "integer")) {
    warning("Converting laglength into integer values, see ?as.integer",
            immediate. = TRUE)
    laglength <- as.integer(laglength)
  }
  if (!inherits(xyt$tlim, "integer")) {
    warning("Converting xyt$tlim into integer values, see ?as.integer",
            immediate. = TRUE)
    xyt$tlim <- as.integer(xyt$tlim)
  }
  if (!inherits(xyt$t, "integer")) {
    warning("Converting xyt$t into integer values, see ?as.integer",
            immediate. = TRUE)
    xyt$t <- as.integer(xyt$t)
  }
  if (xyt$window$type == "rectangle") {
    xyt$window <- as.polygonal(xyt$window)
  }
  if (class(model.priors) != "lgcpPrior") {
    stop("Argument model.priors must be of class lgcpPrior, see ?lgcpPrior")
  }
  if (is.null(cellwidth) & is.null(gridsize)) {
    stop("Either cell width OR grid size must be specified")
  }
  if (!is.null(cellwidth) & !is.null(gridsize)) {
    stop("Either cell width OR grid size must be specified")
  }
  # if (!all(sapply(gridsize, is.pow2))) {
  #   stop("All elements of gridsize must be a power of 2")
  # }
  if (!is.null(gridsize)) {
    approxcw <- diff(xyt$window$xrange)/gridsize[1]
    cwseq <- seq(approxcw/2, 2 * approxcw, length.out = 500)
    cwfun <- function(cw) {
      ow <- lgcp::selectObsWindow(xyt, cw)
      return(c(ow$M, ow$N))
    }
    gsmat <- t(sapply(cwseq, cwfun))
    tf <- apply(gsmat, 1, function(x) {
      return(all(x == gridsize))
    })
    if (sum(tf) == 0) {
      stop("No sensible observation window found: either change gridsize, or specify cellwidth instead")
    }
    else {
      cellwidth <- cwseq[min(which(tf))]
    }
  }
  if (!is.null(gradtrunc)) {
    if (gradtrunc < 0) {
      stop("gradtrunc must be non-negative")
    }
  }
  if (mcmc.control$burnin > mcmc.control$mala.length) {
    stop("Number of burnin iterations must be less than the total number of iterations")
  }
  aggtimes <- T - laglength:0
  nobser <- 0
  for (i in 1:(laglength + 1)) {
    nobser <- nobser + sum(xyt$t == aggtimes[i])
  }
  if (nobser == 0) {
    cat("NOTE: time data should be integer-valued.\n")
    stop("No data in chosen time interval")
  }
  tdiff <- c(Inf, diff(aggtimes))
  numt <- length(tdiff)
  ow <- lgcp::selectObsWindow(xyt, cellwidth)
  xyt <- ow$xyt
  M <- ow$M
  N <- ow$N
  if (M * N >= (256^2)) {
    Sys.sleep(1)
    cat("\n")
    warning("USING LARGE FFT GRID: COMPUTATION MAY BE SLOW ON SOME MACHINES ...",
            .immediate = TRUE)
    cat("\n")
  }
  cat(paste("FFT Grid size: [", ext * M, " , ",
            ext * N, "]\n", sep = ""))
  Sys.sleep(1)
  rm(ow)
  if (is.null(poisson.offset)) {
    poisson.offset <- list()
    for (i in 1:numt) {
      poisson.offset[[i]] <- list(X = seq(xyt$window$xrange[1],
                                          xyt$window$xrange[2], length.out = 100), Y = seq(xyt$window$yrange[1],
                                                                                           xyt$window$yrange[2], length.out = 100), Zm = matrix(1,
                                                                                                                                                100, 100))
    }
  }
  else {
    if (!is.list(poisson.offset)) {
      po <- list()
      for (i in 1:numt) {
        po[[i]] <- poisson.offset
      }
      poisson.offset <- po
      rm(po)
      gc()
    }
    else {
      if (length(poisson.offset) != numt) {
        stop(paste("Poisson offset should have length",
                   numt))
      }
    }
  }
  spatial <- list()
  for (i in 1:numt) {
    if (!inherits(poisson.offset[[i]], "spatialAtRisk")) {
      spatial[[i]] <- lgcp::spatialAtRisk(poisson.offset[[i]])
    }
    else {
      spatial[[i]] <- poisson.offset[[i]]
    }
    if (inherits(spatial[[i]], "fromXYZ")) {
      spatial[[i]]$Zm <- spatial[[i]]$Zm * attr(spatial[[i]],
                                                "NC")
    }
    if (inherits(spatial[[i]], "fromSPDF")) {
      spatial[[i]]$atrisk <- spatial[[i]]$atrisk * attr(spatial[[i]],
                                                        "NC")
      spatial[[i]]$spdf$atrisk <- spatial[[i]]$atrisk
    }
  }
  study.region <- xyt$window
  if (!is.null(attr(ZmatList[[1]], "gridobj"))) {
    gridobj <- attr(ZmatList[[1]], "gridobj")
  }
  else {
    gridobj <- lgcp::genFFTgrid(study.region = study.region, M = M,
                          N = N, ext = ext, inclusion = inclusion)
  }
  del1 <- gridobj$del1
  del2 <- gridobj$del2
  Mext <- gridobj$Mext
  Next <- gridobj$Next
  mcens <- gridobj$mcens
  ncens <- gridobj$ncens
  cellarea <- gridobj$cellarea
  cellInside <- gridobj$cellInside
  x <- gridobj$mcens
  y <- gridobj$ncens
  xidx <- rep(1:Mext, Next)
  yidx <- rep(1:Next, each = Mext)
  dxidx <- pmin(abs(xidx - xidx[1]), Mext - abs(xidx - xidx[1]))
  dyidx <- pmin(abs(yidx - yidx[1]), Next - abs(yidx - yidx[1]))
  d <- sqrt(((x[2] - x[1]) * dxidx)^2 + ((y[2] - y[1]) * dyidx)^2)
  spatial.offset <- list()
  for (i in 1:numt) {
    spatial.offset[[i]] <- lgcp::fftinterpolate(spatial[[i]], mcens,
                                          ncens, ext = ext)
    spatial.offset[[i]] <- spatial.offset[[i]] * cellInside
  }
  spatialOnlyCovariates <- FALSE
  if (is.null(ZmatList)) {
    ZmatList <- list()
    if (!inherits(regionalcovariates, "list") & !inherits(pixelcovariates,
                                                          "list")) {
      ZmatList <- lgcp::cov.interp.fft(formula = formula, W = study.region,
                                 regionalcovariates = regionalcovariates, pixelcovariates = pixelcovariates,
                                 mcens = mcens[1:M], ncens = ncens[1:N], cellInside = cellInside[1:M,
                                                                                                 1:N])
      spatialOnlyCovariates <- TRUE
    }
    else {
      if ((!inherits(regionalcovariates, "list") |
           !inherits(pixelcovariates, "list"))) {
        stop("regionalcovariates and pixelcovariates must EITHER both be list objects OR SpatialPolygonsDataFrame and SpatialPixelsDataFrame objects respectively.")
      }
      else {
        ZmatList[[i]] <- lgcp::cov.interp.fft(formula = formula,
                                        W = study.region, regionalcovariates = regionalcovariates[[i]],
                                        pixelcovariates = pixelcovariates[[i]], mcens = mcens[1:M],
                                        ncens = ncens[1:N], cellInside = cellInside[1:M,
                                                                                    1:N])
      }
    }
  }
  else {
    for (i in 1:numt) {
      if (inherits(ZmatList, "matrix")) {
        if (!isTRUE(all.equal(mcens[1:M], attr(ZmatList,
                                               "mcens"))) | !isTRUE(all.equal(ncens[1:N],
                                                                              attr(ZmatList, "ncens")))) {
          stop(paste("FFT grid and ZmatList[[",
                     i, "]] are on different grids. Please recompute ZmatList using 'getZmat'.",
                     sep = ""))
        }
        spatialOnlyCovariates <- TRUE
      }
      else {
        if (!isTRUE(all.equal(mcens[1:M], attr(ZmatList[[i]],
                                               "mcens"))) | !isTRUE(all.equal(ncens[1:N],
                                                                              attr(ZmatList[[i]], "ncens")))) {
          stop(paste("FFT grid and ZmatList[[",
                     i, "]] are on different grids. Please recompute ZmatList using 'getZmat'.",
                     sep = ""))
        }
      }
    }
  }
  mLoop = lgcp::mcmcLoop(N = mcmc.control$mala.length, burnin = mcmc.control$burnin,
                   thin = mcmc.control$retain, progressor = lgcp::mcmcProgressTextBar)
  nsamp <- floor((mLoop$N - mLoop$burnin)/mLoop$thin)
  if (!is.null(output.control$gridfunction) & class(output.control$gridfunction)[1] ==
      "dump2dir") {
    cat("WARNING: disk space required for saving is approximately ",
        round(nsamp * object.size(array(runif((M) * (N) *
                                                (length(aggtimes))), dim = c((M), (N), (length(aggtimes)))))/1024^2,
              2), " Mb, ", sep = "")
    if (!output.control$gridfunction$forceSave) {
      m <- menu(c("yes", "no"), title = "continue?")
      if (m == 1) {
        cat("Note: to bypass this menu, set forceSave=TRUE in dump2dir\n")
        Sys.sleep(2)
      }
      else {
        stop("Stopped")
      }
    }
  }
  nis <- list()
  for (i in 1:numt) {
    if (sum(xyt$t == aggtimes[i]) > 0) {
      nis[[i]] <- lgcp::getCounts(xyt = xyt, subset = (xyt$t ==
                                                   aggtimes[i]), M = M, N = N, ext = ext)
    }
    else {
      nis[[i]] <- matrix(0, ext * M, ext * N)
    }
    ct1 <- sum(nis[[i]])
    nis[[i]] <- nis[[i]] * (spatial.offset[[i]] > 0)
    ct2 <- sum(nis[[i]])
    if (ct2 < ct1) {
      warning(paste("Time ", aggtimes[i], ": ",
                    ct1 - ct2, " data points lost due to discretisation.",
                    sep = ""), immediate. = TRUE)
    }
  }
  gridfun <- output.control$gridfunction
  if (is.null(gridfun)) {
    gridfun <- lgcp::nullFunction()
  }
  gridav <- output.control$gridmeans
  if (is.null(gridav)) {
    gridav <- lgcp::nullAverage()
  }
  lg <- lgcp::MALAlgcpSpatioTemporal.PlusPars(mcmcloop = mLoop, inits = mcmc.control$inits,
                                        adaptivescheme = mcmc.control$adaptivescheme, M = M,
                                        N = N, Mext = Mext, Next = Next, mcens = mcens, ncens = ncens,
                                        formula = formula, ZmatList = ZmatList, model.priors = model.priors,
                                        model.inits = model.inits, fftgrid = gridobj, spatial.covmodel = spatial.covmodel,
                                        tdiff = tdiff, nis = nis, cellarea = cellarea, spatialvals = spatial.offset,
                                        cellInside = cellInside, MCMCdiag = mcmc.control$MCMCdiag,
                                        gradtrunc = gradtrunc, gridfun = gridfun, gridav = gridav,
                                        d = d, aggtimes = aggtimes, spatialOnlyCovariates = spatialOnlyCovariates)
  endtime <- Sys.time()
  timetaken <- endtime - starttime
  lg$xyt <- xyt
  lg$M <- M
  lg$N <- N
  lg$timetaken <- timetaken
  class(lg) <- c("lgcpPredictSpatioTemporalPlusParameters",
                 "lgcpPredict", "lgcpobject")
  return(lg)
}

#' Minimum constrast spatio-temporal estimation
#'
#' A wrapper to the function \code{lgcp::minimum.contrast.spatiotemporal}
#'
#' \code{mincontrast_st} provides a simplified wrapper to the
#' function \code{lgcp::minimum.contrast.spatiotemporal}. For more details see
#' \code{help(lgcp::minimum.contrast.spatiotemporal).}
#'
#' @param data A data frame consisting of columns \code{x}, \code{y}, and \code{t}, which are two
#' spatial coordinates and time point respectively.
#' @param covariates A \code{spatialPolygonsDataFrame} covering the area of interest with the
#' population density column named \code{popdens}
#' @param boundary A \code{spatialPolygonsDataFrame} of the boundary of the area of interest.
#' @param pop.var Name of the population density variable
#' @param verbose Logical indicating whether to provide progress bars
#' @return Returned values are the minimum contrast estimates of phi, sigma^2 and theta,
#' as well as the overall squared discrepancy between the parametric and nonparametric forms
#' of the spatial function used corresponding to these estimates.
#' @importFrom methods is
#' @export
mincontrast_st <- function(data,
                           covariates,
                           boundary,
                           pop.var = "popdens",
                           verbose=TRUE){
  if(!is(data,"data.frame")|any(!colnames(data)%in%c('x','y','t')))stop("Data needs to be a data frame with columns x,y, and t")
  if(!is(boundary,"SpatialPolygonsDataFrame"))stop("Boundary needs to be of class SpatialPolygonsDataFrame")
  if(!is.null(covariates)&!is(covariates,"SpatialPolygonsDataFrame"))stop("Covariates needs to be of class SpatialPolygonsDataFrame")

  win <- maptools::as.owin.SpatialPolygons(boundary)

  if(!is.null(covariates)){
    requireNamespace("spatstat")
    popVal <- function(x,y){
      spp <- sp::SpatialPoints(data.frame(x=x,y=y))
      crsN <- sp::CRS("+init=epsg:4326")
      sp::proj4string(spp) <- crsN
      sp::proj4string(covariates) <- crsN
      val <- sp::over(spp,covariates)
      return(val[,pop.var])
    }

    pop <- spatstat.geom::as.im(popVal,win)
    xyt <- lgcp::stppp(list(data = data, tlim = range(data$t), window = win))
    vars.est <- suppressWarnings(lgcp::minimum.contrast.spatiotemporal(data=xyt,
                                                      model="exponential",
                                                      spatial.dens = pop,
                                                      temporal.intens = lgcp::muEst(xyt),
                                                      verbose = verbose))
  } else {
    xyt <- lgcp::stppp(list(data = data, tlim = range(data$t), window = win))
    vars.est <- suppressWarnings(lgcp::minimum.contrast.spatiotemporal(data=xyt,
                                                      model="exponential",
                                                      temporal.intens = lgcp::muEst(xyt),
                                                      verbose=verbose))
  }


  return(vars.est$estimates)
}

#' Spatio-temporal Log-Gaussian Cox Process Model
#'
#' Bayesian inference for a spatio-temporal LGCP model with or without covariates.
#'
#' The \code{lgcp} function provides a wrapper to several functions from the \code{lgcp} package.
#' It simplifies the workflow described in the vignette for that package, providing a single
#' function to generate the appropriate grid, covariate matrices and lists, and perform inference
#' with the function \code{lgcp::lgcpPredictSpatioTemporalPlusPars}. See the vignette for this
#' package for a description of the model. The implementation here allows for spatially and/or
#' temporally varying covariates but not spatio-temporally varying covariates, as in the time-scales
#' relevant to real-time surveillance applications these are not generally available. For users
#' requiring additional functionality, please refer to the \code{lgcp} package documentation.
#' @param data A data frame consisting of columns \code{x}, \code{y}, and \code{t}, which are two
#' spatial coordinates and time point respectively.
#' @param data.t A data frame containing any temporal covariates with a column \code{t} with the
#' time period and subsequent columns describing the value(s) of the covariates.
#' @param sp.covs A vector with the names of spatially-varying covariate to use in the model
#' (can be \code{NULL}). These must match the column names in \code{covariates}.
#' @param t.covs A vector with the names of temporally-varying covariates (can be \code{NULL}).
#' These names must match the names of columns in \code{data.t}.
#' @param pop.var The name of the population density variable to be used for the population
#' offset (can be \code{NULL}). This must match the name of a column in \code{covariates}.
#' @param boundary A \code{spatialPolygonsDataFrame} of the boundary of the area of interest.
#' @param covariates A \code{spatialPolygonsDataFrame} covering the area of interest and containing
#' the covariate and population density data.
#' @param cellwidth The width of cells of the computational grid.
#' @param laglength The number of time periods to include. The maximum value of \code{t} in \code{data}
#' is used as the present period, and time periods are counted back from this value.
#' @param dirname The directory root name to save model output. A directory is created for each
#' MCMC chain as \code{dirname.1}, \code{dirname.2}, etc. If NULL then a temporary directory is used,
#' this will result in the data being lost after the session in closed though.
#' @param prevRun Used to set prior distributions. Either output from a previous call to \code{lgcp}
#' to use posterior distributions from previous period, or a call to \code{lgcp::lgcpPrior}.
#' @param mala.pars Parameters for the MCMC sampler. A vector of three numbers: the total number
#' of iterations, the number of warmup iterations, and the number to thin.
#' @param nchains The number of MCMC chains, default is \code{parallel::detectCores()}
#' @param lib Library location if not the default, otherwise NULL
#' @return An object of class lgcpReal
#' @importFrom methods is
#' @importFrom stats as.formula var
#' @examples
#' \donttest{
#' data(dat,square,square_pop)
#' lg1 <- lgcp(data=dat,
#'             pop.var = c("popdens"),
#'             boundary=square,
#'             covariates=square_pop,
#'             cellwidth=0.1,
#'             laglength = 7,
#'             mala.pars=c(200,100,1),
#'             nchains=2)
#' }
#' @export
lgcp <- function(data,
                 data.t=NULL,
                 sp.covs=NULL,
                 t.covs=NULL,
                 pop.var=NULL,
                 boundary,
                 covariates=NULL,
                 cellwidth,
                 laglength,
                 dirname=NULL,
                 prevRun=NULL,
                 mala.pars=c(26250,20000,50),
                 nchains=parallel::detectCores(),
                 lib=NULL){

  if(!is(data,"data.frame")|any(!colnames(data)%in%c('x','y','t')))stop("Data needs to be a data frame with columns x,y, and t")
  if(!is(boundary,"SpatialPolygonsDataFrame"))stop("Boundary needs to be of class SpatialPolygonsDataFrame")
  if(!is.null(covariates)&!is(covariates,"SpatialPolygonsDataFrame"))stop("Covariates needs to be of class SpatialPolygonsDataFrame")
  if(any(is.na(data$x)|is.na(data$y)))warning(paste0(sum(is.na(data$x)|is.na(data$y))," rows have NA values and will be removed\n"))
  if(is.null(dirname))warning('Dirname is NULL so any model fits will be lost once the session is closed\n')

  data <- data[!is.na(data$x)&!is.na(data$y),]
  data <- data[,c('x','y','t')]
  data <- data[order(data$t),]

  tlim <- c(1,laglength)
  data <- data[data$t %in% c(max(data$t):(max(data$t)-laglength)) ,]
  data$t <- data$t - (max(data$t)-laglength)
  data <- data[data$t%in%c(tlim[1]:tlim[2]),]

  win <- maptools::as.owin.SpatialPolygons(boundary)

  xyt <- lgcp::stppp(list(data = data, tlim = tlim, window = win))
  Owin <- lgcp::selectObsWindow(xyt,cellwidth)

  if(!is.null(sp.covs)){
    form.sp <- "X ~ "
    for(i in 1:length(sp.covs)){
      if(i==1){
        form.sp <- paste0(form.sp," ",sp.covs[i])
      } else {
        form.sp <- paste0(form.sp," + ",sp.covs[i])
      }
    }
  } else {
    form.sp <- "X ~ 1"
  }


  if(!is.null(t.covs)){
    form.t <- "t ~ -1 "
    for(i in 1:length(t.covs)){
      form.t <- paste0(form.t," + ",t.covs[i])
      form <- paste0(form.sp," + ",t.covs[i])
    }
    form <- as.formula(form)
    form.t <- as.formula(form.t)
  } else {
    form <- as.formula(form.sp)
    form.t <- NULL
  }

  T <- tlim[2]

  if(!is.null(pop.var)|!is.null(sp.covs)){
    covariates@data <- lgcp::guessinterp(covariates@data)

  }

  if(!is.null(sp.covs)){
    form.sp <- as.formula(form.sp)

    polyolay <- lgcp::getpolyol(data = xyt,
                          cellwidth = cellwidth,
                          regionalcovariates = covariates,
                          ext = 2)

    #covariates@data <- lgcp::guessinterp(covariates@data)

    Zmat <- lgcp::getZmat(
      formula = form.sp,
      data = xyt,
      cellwidth = cellwidth,
      regionalcovariates = covariates,
      ext = 2,
      overl = polyolay
    )
  } else {

    form.sp <- X ~ 1
    attr(form.sp, ".Environment") <- .GlobalEnv
    polyolay <- lgcp::getpolyol(data = xyt,
                          cellwidth = cellwidth,
                          regionalcovariates = covariates,
                          ext = 2)


    Zmat <- lgcp::getZmat(
      formula = form.sp,
      data = xyt,
      cellwidth = cellwidth,
      ext = 2,
      overl = polyolay
    )
  }

  if(!is.null(pop.var)){
    form.pop <- as.formula(paste0("X ~ -1 + ",pop.var))

    Zmat_pop <- lgcp::getZmat(
      formula = form.pop,
      data = xyt,
      cellwidth = cellwidth,
      regionalcovariates = covariates,
      ext = 2,
      overl = polyolay
    )
    mm <- length(attr(Zmat_pop, "mcens"))
    nn <- length(attr(Zmat_pop, "ncens"))
    offset <- lgcp::spatialAtRisk(list(X = attr(Zmat_pop, "mcens"),
                                 Y = attr(Zmat_pop, "ncens"),
                                 Zm = matrix(Zmat_pop, mm, nn)))
    offsetList <- list()
    for(i in 1:(laglength))offsetList[[i]] <- offset
  } else {
    form.pop <- NULL
    offset <- NULL
    offsetList <- NULL
  }





  message("\nAdding temporal covariates\n")
  if(!is.null(t.covs)){
    Zmat <- lgcp::addTemporalCovariates(temporal.formula = form.t,
                                  T = T,
                                  laglength = laglength-1,
                                  tdata = data.t,
                                  Zmat = Zmat)

  }


  if(is.list(Zmat)){
    nvar <- ncol(Zmat[[1]])
  }else{
    nvar <- ncol(Zmat)
  }

  if(is.null(prevRun)){
    gprior <- lgcp::PriorSpec(lgcp::GaussianPrior(mean = rep(0,nvar),variance = diag(rep(25,nvar),nvar)))
    if(!is.null(pop.var)){
      popVal <- function(x,y){
        spp <- sp::SpatialPoints(data.frame(x=x,y=y))
        crsN <- sp::CRS("+init=epsg:4326")
        sp::proj4string(spp) <- crsN
        sp::proj4string(covariates) <- crsN
        val <- sp::over(spp,covariates)
        return(val[,pop.var])
      }

      pop <- spatstat.geom::as.im(popVal,win)
      vars.est <- suppressWarnings(lgcp::minimum.contrast.spatiotemporal(data=xyt,
                                                  model="exponential",
                                                  spatial.dens = pop,
                                                  temporal.intens = lgcp::muEst(xyt)))
    } else {
      vars.est <- suppressWarnings(lgcp::minimum.contrast.spatiotemporal(data=xyt,
                                                  model="exponential",
                                                  temporal.intens = lgcp::muEst(xyt)))
    }


    INITS <- lgcp::lgcpInits(etainit = log(c(vars.est$estimates[2],
                                       vars.est$estimates[1],
                                       vars.est$estimates[3])),
                       betainit = NULL )

    lgprior <- lgcp::PriorSpec(lgcp::LogGaussianPrior(mean = log(c(vars.est$estimates[2],
                                                       vars.est$estimates[1],
                                                       vars.est$estimates[3])),
                                          variance = diag(rep(log(5),3), 3)))

    priors <- lgcp::lgcpPrior(etaprior = lgprior, betaprior = gprior)
  } else if(is(prevRun,"lgcpReal")){
    priors <- lgcp::lgcpPrior(etaprior = prevRun$lgprior, betaprior = prevRun$gprior)
    INITS <- prevRun$INITS
  } else if(is(prevRun,"lgcpPrior")){
    priors <- prevRun
  }

  CF <- lgcp::CovFunction(lgcp::exponentialCovFct)

  ## parellise
  message("\nStarting sampling... This may take a long time.\n")

  if(is.null(dirname)){
    dir1 <- tempdir()
    dir1 <- paste0(dir1,"\\",as.numeric(Sys.time()))
  } else {
    dir1 <- dirname
  }

  pbapply::pboptions(type="none")

  if(nchains > 1){
    cl <- parallel::makeCluster(nchains)
    if(!is.null(lib)){
      parallel::clusterExport(cl,c('lib'),envir = environment())
      parallel::clusterEvalQ(cl,.libPaths(lib))
    }
    #parallel::clusterEvalQ(cl,library(realTimeSurv))
    parallel::clusterEvalQ(cl,library(lgcp))
    parallel::clusterCall(cl, assign, "lgcpST", lgcpST, envir = .GlobalEnv)
    #parallel::clusterCall(cl, requireNamespace("lgcp"))

    parallel::clusterExport(cl,c('form','xyt','T','laglength','Zmat','priors','INITS',
                                 'CF','cellwidth','dir1','mala.pars',"offsetList"),
                            envir = environment())

    lg.out <- pbapply::pbsapply(1:nchains,function(i)lgcpST(formula = form,
                                                            xyt = xyt,
                                                            T = T,
                                                            laglength = laglength-1,
                                                            ZmatList = Zmat,
                                                            poisson.offset = offsetList,
                                                            model.priors = priors,
                                                            model.inits = INITS,
                                                            spatial.covmodel = CF,
                                                            cellwidth = cellwidth,
                                                            mcmc.control = lgcp::mcmcpars(mala.length = mala.pars[1],
                                                                                          burnin = mala.pars[2],
                                                                                          retain = mala.pars[3],
                                                                                          adaptivescheme = lgcp::andrieuthomsh(inith = 1,
                                                                                                                               alpha = 0.5,
                                                                                                                               C = 1,
                                                                                                                               targetacceptance = 0.574)),
                                                            output.control = lgcp::setoutput(gridfunction =
                                                                                               lgcp::dump2dir(dirname = file.path(paste0(dir1,".",i)),
                                                                                                              lastonly = F,
                                                                                                              forceSave = TRUE)),
                                                            ext = 2),
                                cl = cl)
    parallel::stopCluster(cl)
  } else {
    lg.out <- pbapply::pbsapply(1:nchains,function(i)lgcpST(formula = form,
                                                            xyt = xyt,
                                                            T = T,
                                                            laglength = laglength-1,
                                                            ZmatList = Zmat,
                                                            poisson.offset = offsetList,
                                                            model.priors = priors,
                                                            model.inits = INITS,
                                                            spatial.covmodel = CF,
                                                            cellwidth = cellwidth,
                                                            mcmc.control = lgcp::mcmcpars(mala.length = mala.pars[1],
                                                                                          burnin = mala.pars[2],
                                                                                          retain = mala.pars[3],
                                                                                          adaptivescheme = lgcp::andrieuthomsh(inith = 1,
                                                                                                                               alpha = 0.5,
                                                                                                                               C = 1,
                                                                                                                               targetacceptance = 0.574)),
                                                            output.control = lgcp::setoutput(gridfunction =
                                                                                               lgcp::dump2dir(dirname = file.path(paste0(dir1,".",i)),
                                                                                                              lastonly = F,
                                                                                                              forceSave = TRUE)),
                                                            ext = 2))
  }




  message("\nSampling complete at: ",Sys.time())


  if(is(lg.out,"matrix")){
    eta <- do.call(rbind,lapply(1:nchains,function(i)return(lg.out[,i]$etarec)))
    beta <- do.call(rbind,lapply(1:nchains,function(i)return(lg.out[,i]$betarec)))
    timetaken <- do.call(rbind,lapply(1:nchains,function(i)return(lg.out[,i]$timetaken)))
    lasth <- do.call(rbind,lapply(1:nchains,function(i)return(lg.out[,i]$lasth)))
  } else if(is(lg.out,"list")){
    eta <- do.call(rbind,lapply(1:nchains,function(i)return(lg.out[[i]]$etarec)))
    beta <- do.call(rbind,lapply(1:nchains,function(i)return(lg.out[[i]]$betarec)))
    timetaken <- do.call(rbind,lapply(1:nchains,function(i)return(lg.out[[i]]$timetaken)))
    lasth <- do.call(rbind,lapply(1:nchains,function(i)return(lg.out[[i]]$lasth)))
  } else {
    eta <- do.call(rbind,lapply(lg.out,function(i)return(i$etarec)))
    beta <- do.call(rbind,lapply(lg.out,function(i)return(i$betarec)))
    timetaken <- do.call(rbind,lapply(1:nchains,function(i)return(i$timetaken)))
    lasth <- do.call(rbind,lapply(1:nchains,function(i)return(i$lasth)))
  }

  environment(form.sp) <- NULL
  environment(form) <- NULL
  environment(form.pop) <- NULL
  environment(form.t) <- NULL

  if(is.list(Zmat)){
    ins <- attr(Zmat[[1]], "gridobj")$cellInside[1:Owin$M,1:Owin$N]
  } else {
    ins <- attr(Zmat, "cellInside")[1:Owin$M,1:Owin$N]
  }

  out <- list(
    eta = eta,
    beta = beta,
    cellInside = ins,
    lgcpRunInfo = list(
      timetaken=unlist(timetaken),
      lasth = unlist(lasth),
      priors=priors,
      M = Owin$M,
      N = Owin$N),
    lgprior = lgcp::PriorSpec(lgcp::LogGaussianPrior(mean = colMeans(eta),variance = diag(apply(eta,2,var), 3))),
    INITS = lgcp::lgcpInits(etainit = colMeans(eta),betainit = NULL ),
    gprior = lgcp::PriorSpec(lgcp::GaussianPrior(mean = colMeans(beta),variance = diag(apply(beta,2,var),nvar))),
    formulae = list(form.sp=form.sp,
                    form=form,
                    form.pop=form.pop,
                    form.t=form.t),
    dirname=dir1,
    boundary = boundary,
    cellwidth=cellwidth,
    data=data,
    data.t=data.t,
    xyt=xyt,
    nchains=nchains
  )

  class(out) <- "lgcpReal"

  return(out)

}

#' MCMC convergence diagnostics
#'
#' Diagnostics for convergence of the MCMC chains from a call to \code{lgcp}.
#'
#' Produces a traceplot of the model parameters and prints R-hat and ESS statistics.
#'
#' @param lg Output from a call to \code{lgcp}
#' @param plots Logical indicating whether to plot MCMC traceplots
#' @return Traceplot of the MCMC chains of parameters from the linear predictor and
#' covariance function is plotted, and R-hat and ESS statistics are printed.
#' @importFrom methods is
#' @importFrom stats terms
#' @export
convergence <- function(lg,
                        plots=TRUE){
  if(!is(lg,"lgcpReal"))stop("lg must be of class lgcpReal")
  nchains <- nrow(lg$lgcpRunInfo$timetaken)
  iter <- nrow(lg$beta)
  iter2 <- iter/nchains
  labs <- c("(Intercept)")
  if(!is.null(lg$formulae$form.sp)){
    labs <- c(labs,attr(terms(lg$formulae$form.sp),"term.labels"))
  }
  if(!is.null(lg$formulae$form.t)){
    labs.t <- attr(terms(lg$formulae$form.t),"term.labels")
    labs <- c(labs,levels(factor(lg$data.t[,labs.t]))[-1])
  }
  eta.labs <- c("Sigma^2","Spatial range","Temporal range")

  betalist <- lapply(1:nchains,function(i){
    out <- lg$beta[(1+(i-1)*iter2):(iter2+(i-1)*iter2),]
    if(is.null(dim(out))){
      out <- matrix(out,ncol=1)
    }
    colnames(out) <- labs
    return(out)
  })
  betalist_mcmc <- coda::mcmc.list(lapply(betalist,function(i)coda::mcmc(i)))
  etalist <- lapply(1:nchains,function(i){
    out <- lg$eta[(1+(i-1)*iter2):(iter2+(i-1)*iter2),]
    colnames(out) <- eta.labs
    return(out)
  })
  etalist_mcmc <- coda::mcmc.list(lapply(etalist,function(i)coda::mcmc(i)))
  if(plots){
    p1 <- bayesplot::mcmc_trace(betalist)
    #p1 <- p1 + scale_colour_brewer(palette = "Set1")
    p2 <- bayesplot::mcmc_trace(etalist)
    #p2 <- p2 + scale_colour_brewer(palette = "Set1")
    print(ggpubr::ggarrange(p1,p2,ncol=1,heights = c(2,1)))
  }


  beta_r <- coda::gelman.diag(betalist_mcmc)
  beta_ess <- coda::effectiveSize(betalist_mcmc)
  eta_r <- coda::gelman.diag(etalist_mcmc)
  eta_ess <- coda::effectiveSize(etalist_mcmc)
  res <- cbind(as.data.frame(beta_r$psrf[,1]),as.data.frame(beta_ess))
  res_e <- cbind(as.data.frame(eta_r$psrf[,1]),as.data.frame(eta_ess))
  colnames(res) <- c("R-hat","ESS")
  colnames(res_e) <- c("R-hat","ESS")

  return(rbind(res,res_e))
}

#' Variance partition coefficient
#'
#' Reports the proportion of variance attributable to the latent Gaussian process.
#'
#' The total variance in either the log incidence or log relatve risk
#' can be partitioned into observed and latent components. It is assumed that
#' the log incidence and log relative risk are normally distributed so that the VPC is
#' equal to the variance of the latent Gaussian field over the total variance.
#' The variance depends on the intensity of the Poisson process so \code{vpc}
#' reports quantiles of the distribution of the variance partition coeffient across
#' the predicted values.
#' @param lg Output from a call to \code{lgcp}
#' @param covariates A \code{spatialPolygonsDataFrame} covering the area of interest and containing
#' the covariate and population density data. Typically the same object as specified in the
#' \code{covariates} argument in the call to \code{lgcp}.
#' @param rr Whether to report the VPC for the log relative risk (TRUE) or for the log
#' incidence (FALSE)
#' @return Prints the quantiles of the VPC. An object called \code{outl} is exported to the
#' global environment which contains the samples from the model. Used to reduce loading time
#' of samples. This object is large so remove if no further analysis required.
#' @importFrom stats quantile var
#' @export
vpc <- function(lg,
                covariates,
                rr=FALSE){
  if(is.null(covariates))stop("Please specify covariates.")
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

  tmp <- suppressWarnings( plot_lgcp_dat(outl,
                                         grid.data,
                                         lg,
                                         lg$nchains,
                                         idx.mapping,
                                         covariates = covariates,
                                         data.out = TRUE))

  tmp$dat1$vpc <- NA
  for(i in 1:nrow(tmp$xpred)){
    re_samp <- sample(tmp$xpred,ncol(tmp$xpred))
    if(!rr){
      tmp$dat1$vpc[i] <- var(tmp$pop[i,]*tmp$linpred[i,]*re_samp,na.rm=TRUE)/(
        var(tmp$pop[i,]*tmp$linpred[i,]*re_samp,na.rm=TRUE)+
          mean(tmp$pop[i,]*tmp$linpred[i,]*re_samp,na.rm=TRUE))
    } else {
      tmp$dat1$vpc[i] <- var(log(tmp$linpred[i,]))/(var(log(tmp$linpred[i,]))+var(log(re_samp)))
    }
  }

  res <- t(as.data.frame(quantile(tmp$dat1$vpc,c(0.025,0.25,0.5,0.75,0.975),na.rm=TRUE)))
  vals <- res[,1:ncol(res)]
  res[,1:ncol(res)] <- paste0(round(res[,1:ncol(res)]*100,0),"%")
  rownames(res) <- "VPC"

  res[,1:ncol(res)] <- round(as.numeric(vals),3)
  print(res)
  return(invisible(res))
}
