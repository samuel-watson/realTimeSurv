#' Alternative to lgcpPredictSpatioTemporalPlusPars
#'
#' Internal function and alternative to \code{lgcp::lgcpPredictSpatioTemporalPlusPars}.
#'
#' A copy of \code{lgcp::lgcpPredictSpatioTemporalPlusPars} that parallelises the sampler and
#' produces a reduced output and lgcpReal object
.lgcpST <- function (formula, xyt, T, laglength, ZmatList = NULL, model.priors,
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
  if (!all(sapply(gridsize, is.pow2))) {
    stop("All elements of gridsize must be a power of 2")
  }
  if (!is.null(gridsize)) {
    approxcw <- diff(xyt$window$xrange)/gridsize[1]
    cwseq <- seq(approxcw/2, 2 * approxcw, length.out = 500)
    cwfun <- function(cw) {
      ow <- selectObsWindow(xyt, cw)
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
  ow <- selectObsWindow(xyt, cellwidth)
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
      spatial[[i]] <- spatialAtRisk(poisson.offset[[i]])
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
    gridobj <- genFFTgrid(study.region = study.region, M = M,
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
    spatial.offset[[i]] <- fftinterpolate(spatial[[i]], mcens,
                                          ncens, ext = ext)
    spatial.offset[[i]] <- spatial.offset[[i]] * cellInside
  }
  spatialOnlyCovariates <- FALSE
  if (is.null(ZmatList)) {
    ZmatList <- list()
    if (!inherits(regionalcovariates, "list") & !inherits(pixelcovariates,
                                                          "list")) {
      ZmatList <- cov.interp.fft(formula = formula, W = study.region,
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
        ZmatList[[i]] <- cov.interp.fft(formula = formula,
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
  mLoop = mcmcLoop(N = mcmc.control$mala.length, burnin = mcmc.control$burnin,
                   thin = mcmc.control$retain, progressor = mcmcProgressTextBar)
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
      nis[[i]] <- getCounts(xyt = xyt, subset = (xyt$t ==
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
    gridfun <- nullFunction()
  }
  gridav <- output.control$gridmeans
  if (is.null(gridav)) {
    gridav <- nullAverage()
  }
  lg <- MALAlgcpSpatioTemporal.PlusPars(mcmcloop = mLoop, inits = mcmc.control$inits,
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
#' @return Returned values are the minimum contrast estimates of phi, sigma^2 and theta,
#' as well as the overall squared discrepancy between the parametric and nonparametric forms
#' of the spatial function used corresponding to these estimates.
#' @export
mincontrast_st <- function(data,covariates,boundary){
  requireNamespace("spatstat")
  popVal <- function(x,y){
    spp <- sp::SpatialPoints(data.frame(x=x,y=y))
    sp::proj4string(spp) <- sp::proj4string(covariates)
    val <- sp::over(spp,covariates)
    return(val[,"popdens"])
  }
  win <- as.owin.SpatialPolygons(boundary)
  pop <- as.im(popVal,win)


  xyt <- lgcp::stppp(list(data = data, tlim = range(data$t), window = win))
  vars.est <- lgcp::minimum.contrast.spatiotemporal(data=xyt,
                                              model="exponential",
                                              spatial.dens = pop,
                                              temporal.intens = muEst(xyt))
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
#' MCMC chain as \code{dirname.1}, \code{dirname.2}, etc.
#' @param prevRun Used to set prior distributions. Either output from a previous call to \code{lgcp}
#' to use posterior distributions from previous period, or a call to \code{lgcp::lgcpPrior}.
#' @param mala.pars Parameters for the MCMC sampler. A vector of three numbers: the total number
#' of iterations, the number of warmup iterations, and the number to thin.
#' @param nchains The number of MCMC chains, default is \code{parallel::detectCores()}
#' @return An object of class lgcpReal
#' @export
lgcp <- function(data,
                 data.t=NULL,
                 sp.covs=NULL,
                 t.covs=NULL,
                 pop.var=NULL,
                 boundary,
                 covariates,
                 cellwidth,
                 laglength,
                 dirname,
                 prevRun=NULL,
                 mala.pars=c(26250,20000,50),
                 nchains=parallel::detectCores()){

  if(class(data)!="data.frame"|any(!colnames(data)%in%c('x','y','t')))stop("Data needs to be a data frame with columns x,y, and t")
  if(class(boundary)!="SpatialPolygonsDataFrame")stop("Boundary needs to be of class SpatialPolygonsDataFrame")
  if(class(covariates)!="SpatialPolygonsDataFrame")stop("Covariates needs to be of class SpatialPolygonsDataFrame")
  #if(length(pop.var)!=1)stop("Name one population variable.")
  requireNamespace(spatstat)
  tlim <- c(1,laglength)
  data <- data[data$t %in% c(max(data$t):(max(data$t)-laglength)) ,]
  data$t <- data$t - (max(data$t)-laglength)
  data <- data[data$t%in%c(tlim[1]:tlim[2]),]

  win <- as.owin.SpatialPolygons(boundary)

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

  covariates@data <- lgcp::guessinterp(covariates@data)

  if(!is.null(sp.covs)){
    form.sp <- as.formula(form.sp)

    polyolay <- lgcp::getpolyol(data = xyt,
                          cellwidth = cellwidth,
                          regionalcovariates = covariates,
                          ext = 2)

    covariates@data <- lgcp::guessinterp(covariates@data)

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





  cat("\nAdding temporal covariates\n")
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
    gprior <- lgcp::PriorSpec(GaussianPrior(mean = rep(0,nvar),variance = diag(rep(25,nvar),nvar)))
    if(!is.null(pop.var)){
      popVal <- function(x,y){
        spp <- sp::SpatialPoints(data.frame(x=x,y=y))
        sp::proj4string(spp) <- sp::proj4string(covariates)
        val <- sp::over(spp,covariates)
        return(val[,pop.var])
      }

      pop <- as.im(popVal,win)
      vars.est <- lgcp::minimum.contrast.spatiotemporal(data=xyt,
                                                  model="exponential",
                                                  spatial.dens = pop,
                                                  temporal.intens = muEst(xyt))
    } else {
      vars.est <- lgcp::minimum.contrast.spatiotemporal(data=xyt,
                                                  model="exponential",
                                                  temporal.intens = muEst(xyt))
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
  } else if(class(prevRun)=="lgcpReal"){
    priors <- lgcp::lgcpPrior(etaprior = prevRun$lgprior, betaprior = prevRun$gprior)
    INITS <- prevRun$INITS
  } else if(class(prevRun)[1]=="lgcpPrior"){
    priors <- prevRun
  }

  CF <- lgcp::CovFunction(exponentialCovFct)

  ## parellise
  cat("\nStarting sampling... This may take a long time.\n")
  cl <- parallel::makeCluster(nchains)
  parallel::clusterEvalQ(cl,require(lgcp))
  parallel::clusterExport(cl,c('form','xyt','T','laglength','Zmat','priors','INITS',
                     'CF','cellwidth','dirname','mala.pars',"offsetList",".lgcpST"),
                envir = environment())

  pbapply::pboptions(type="none")
  lg.out <- pbapply::pbsapply(1:8,function(i).lgcpST(formula = form,
                                           xyt = xyt,
                                           T = T,
                                           laglength = laglength-1,
                                           ZmatList = Zmat,
                                           poisson.offset = offsetList,
                                           model.priors = priors,
                                           model.inits = INITS,
                                           spatial.covmodel = CF,
                                           cellwidth = cellwidth,
                                           mcmc.control = mcmcpars(mala.length = mala.pars[1],
                                                                   burnin = mala.pars[2],
                                                                   retain = mala.pars[3],
                                                                   adaptivescheme = andrieuthomsh(inith = 1,
                                                                                                  alpha = 0.5,
                                                                                                  C = 1,
                                                                                                  targetacceptance = 0.574)),
                                           output.control = setoutput(gridfunction =
                                                                        dump2dir(dirname = file.path(paste0(dirname,".",i)),
                                                                                 lastonly = F,
                                                                                 forceSave = TRUE)),
                                           ext = 2),
                     cl = cl)
  cat("\nSampling complete at: ",Sys.time())
  stopCluster(cl)

  if(class(lg.out)[1]=="matrix"){
    eta <- do.call(rbind,lapply(1:8,function(i)return(lg.out[,i]$etarec)))
    beta <- do.call(rbind,lapply(1:8,function(i)return(lg.out[,i]$betarec)))
    timetaken <- do.call(rbind,lapply(1:8,function(i)return(lg.out[,i]$timetaken)))
    lasth <- do.call(rbind,lapply(1:8,function(i)return(lg.out[,i]$lasth)))
  } else if(class(lg.out)[1]=="list"){
    eta <- do.call(rbind,lapply(1:8,function(i)return(lg.out[[i]]$etarec)))
    beta <- do.call(rbind,lapply(1:8,function(i)return(lg.out[[i]]$betarec)))
    timetaken <- do.call(rbind,lapply(1:8,function(i)return(lg.out[[i]]$timetaken)))
    lasth <- do.call(rbind,lapply(1:8,function(i)return(lg.out[[i]]$lasth)))
  } else {
    eta <- do.call(rbind,lapply(lg.out,function(i)return(i$etarec)))
    beta <- do.call(rbind,lapply(lg.out,function(i)return(i$betarec)))
    timetaken <- do.call(rbind,lapply(1:8,function(i)return(i$timetaken)))
    lasth <- do.call(rbind,lapply(1:8,function(i)return(i$lasth)))
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
    gprior = lgcp::PriorSpec(GaussianPrior(mean = colMeans(beta),variance = diag(apply(beta,2,var),nvar))),
    formulae = list(form.sp=form.sp,
                    form=form,
                    form.pop=form.pop,
                    form.t=form.t),
    dirname=dirname,
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

#' Report generating function
#' @export
generateReport <- function(lg,
                           plot.opts=list(),
                           dirname=getwd(),
                           change.lag=NULL,
                           par.summ = FALSE,
                           hotspot.opts=NULL,
                           breakdown.plots=TRUE,
                           aggregate=NULL,
                           repdate=format(Sys.time(), '%d %B, %Y'),
                           add3D=FALSE,
                           pdf=FALSE) {

  if(is.null(plot.opts)|length(plot.opts)==0)stop("No plot options provided.\n")
  if(!is.null(aggregate)&!(class(aggregate)=="SpatialPolygons"|
                           class(aggregate)=="SpatialPolygonsDataFrame"))stop("Aggregate must be of class SpatialPolygons or SpatialPolygonsDataFrame")
  set.plot.osm1 <- FALSE
  if(!is.null(plot.opts$osm)){
    if(plot.opts$osm==TRUE){
      plot.opts$osm <- FALSE
      set.plot.osm1 <- TRUE
    }
  }

  set.plot.osm2 <- FALSE
  if(!is.null(hotspot.opts$osm)){
    if(hotspot.opts$osm==TRUE){
      hotspot.opts$osm <- FALSE
      set.plot.osm2 <- TRUE
    }
  }

  rmd_file_name <- paste0("surveillance_",format(Sys.time(), '%d_%m_%Y'),".Rmd")
  nm <-deparse(substitute(lg))
  nm_dat <-deparse(substitute(cov))
  datetoday <- repdate

  if(!pdf){
    content <- paste0(
      "---",
      "\n",
      "title: Surveillance report\n",
      "output:  html_document\n",
      "date: ",
      datetoday,
      "\n",
      "---"
    )
  } else {
    content <- paste0(
      "---",
      "\n",
      "title: Surveillance report\n",
      "output:  pdf_document\n",
      "date: ",
      datetoday,
      "\n",
      "---"
    )
  }


  #   opts <- "```{r setup, include=FALSE}
  # knitr::opts_chunk$set(echo = TRUE)
  # ```"

  tmp <- do.call(plot.lgcpReal,append(plot.opts,list(lg),after=0))
  if(!is.null(aggregate)){
    tmp <- aggregator(tmp,aggregate,osm=set.plot.osm1)
  }

  content_mainplot <- paste0(
    "## Incidence\n\n",
    "The incidence of cases is shown in the plot below as the number of cases per 10,000 person days.\n",
    "```{r, echo=FALSE}\ntmp[[1]]\n",
    "```\n\n"
  )

  content <- paste0(content,
                    "\n\n",
                    content_mainplot)

  if(breakdown.plots){
    content_secondplot <- paste0(
      "### Components\n\n",
      "The breakdown of the different components of incidence is shown below:\n\n
    1. 'Expected': the expected number of cases from each location,\n
    2. 'Observed': differences in risk of a case associated with observed factors (as a relative risk),\n
    3. 'latent': unexplained differences in the risk of cases (as a relative risk),\n
    4. 'posterior SD': the standard deviation of the prediction of incidence.\n\n",
      "```{r, echo=FALSE}\ntmp[[2]]\n",
      "```\n\n"
    )


    content <- paste0(content,
                      "\n\n",
                      content_secondplot)
  }



  if(!is.null(change.lag)){

    chang.opts <- append(plot.opts,c(change.lag=change.lag),after=0)
    chang.opts <- append(chang.opts,list(lg),after=0)

    tmp2 <- do.call(plot.lgcpReal,chang.opts)
    if(!is.null(aggregate)){
      tmp2 <- aggregator(tmp2,aggregate)
    }

    content_addplot <- paste0(
      "## Change to incidence\n\n",
      "The change to incidence compared to ",change.lag," days ago. The values represent _incidence
      rate ratio_, for example a value of 1 indicates no change, a value of 2 is that on average the incidence
      has doubled from what it was ",change.lag," days previous, and a value of 0.5 is a halving.
      The yellow to red colours indicate values greater than one for
       relative risks and incidence rate ratios, and blue to purple is for values below 1.\n\n\n",
      "```{r, echo=FALSE}\ntmp2[[1]]\n",
      "```\n\n\n")

    if(breakdown.plots){
      content_addplot <- paste0(content_addplot,
                                "And comparisons of the different components over this period.\n\n\n",
                                "```{r, echo=FALSE}\ntmp2[[2]]\n",
                                "```\n\n"
      )
    }


    content <- paste0(content,content_addplot)
  }

  if(!is.null(hotspot.opts)){

    tmp3 <- do.call(plot_hotspot,append(hotspot.opts,list(lg),after=0))
    tstr <- attr(tmp3,'str')
    tval <- attr(tmp3,'vals')
    tlab <- attr(tmp3,'labs')

    if(!is.null(aggregate)){
      tmp3 <- aggregator(tmp3,aggregate,osm=set.plot.osm2)
    }

    if(length(tstr)==1){

      if(grepl("pop",tstr)&grepl("obs",tstr)&grepl("latent",tstr)&!grepl("lag",tstr)){
        tstr <- "Incidence"
      } else if(grepl("pop",tstr)&grepl("obs",tstr)&grepl("latent",tstr)&grepl("lag",tstr)){
        tstr <- "Incidence rate ratio"
      }

      hotspot_desc <- paste0(
        "The hotspot is defined as:\n\n",
        "1. ",tstr," < ",tval,": ",tlab[1],"\n\n",
        "2. ",tstr," >= ",tval,": ",tlab[2],"\n\n"
      )
    } else{
      hotspot_desc <- paste0(
        "The hotspot is defined as:\n\n",
        "1. ",tstr[1]," < ",tval[1]," **and** ",
        tstr[2]," < ",tval[2],
        ": ",tlab[1],"\n\n",
        "2. ",tstr[1]," >= ",tval[1]," **and** ",
        tstr[2]," < ",tval[2],
        ": ",tlab[2],"\n\n",
        "3. ",tstr[1]," < ",tval[1]," **and** ",
        tstr[2]," >= ",tval[2],
        ": ",tlab[3],"\n\n",
        "4. ",tstr[1]," >= ",tval[1]," **and** ",
        tstr[2]," >= ",tval[2],
        ": ",tlab[4],"\n\n"
      )
    }

    # hotplot_func <- "```{r, echo=FALSE, results='hide', warning=FALSE, message=FALSE,fig.show='hide'}
    #   \ntmp3 <- do.call(plot_hotspot,hotspot.opts)\ntvals <- attr(tmp3,'vals')\n
    # tprob <- attr(tmp3,'threshold')\n tstr <- attr(tmp3,'str')\n```\n\n\n"

    # hotspot_desc2 <- "The terms above describe different combinations of terms from the model.
    # 'pop' indicates the population density has been included, 'obs' is the observed covariates, and
    # 'latent' is the unexplained differences. Where 'pop' is not included the results are relative risks
    # that multiply the expected risk in each area, where it is included, then the results are
    # an incidence. For example 'pop+obs+latent' is the complete model a predicts absolute incidence,
    # whereas 'obs+latent' is the relative risk accounted for by both observed and unobserved
    # components. If a term 'lag(t)' is included, this indicates that the preceding term is differenced
    # with *t* periods prior.\n\n"

    hotspot_desc2 <- "The following hotspot classification indicates that if there is a greater than
    80% probability that the incidence has increased by 50% or more in the past 7 days. The second plot
    shows the raw probability of being a 'hotspot'."

    content_hotplot <- paste0(
      "## Hotspots\n\n",
      "The area below is classified into hotspots based on the definition provided.\n\n",hotspot_desc,hotspot_desc2,"\n\n",
      "```{r, echo=FALSE}\ntmp3[[1]]\n",
      "```\n\n\n\n",
      "```{r, echo=FALSE}\ntmp3[[2]]\n",
      "```\n\n\n"
    )

    content <- paste0(content,content_hotplot)
  }

  if(par.summ){
    summtab <- summary_html(lg)

    content_summary <- paste0(
      "## Model parameters\n\n",summtab[[3]],summtab[[1]],"\n\n\n### Plot of model parameters\n\n\n",
      "```{r, echo=FALSE}\nsummtab[[2]]\n",
      "```\n\n\n"
    )

    content <- paste0(content,content_summary)
  }

  if(add3D){
    content_3d_head <- "## 3D renders\n\n### Incidence\n\n"
    content_3d <- paste0(content_3d_head,"```{r, echo=FALSE}\nrayshader::plot_gg(tmp[[1]],
            multicore = TRUE,width=5,height=5,scale=125,zoom=0.5, windowsize = c(1400,866),
            phi=30,theta=35)\nrayshader::render_snapshot(clear=T)\n",
                         "```\n\n")
    # if(osm){
    #   content_3d<- paste0(content_3d_head,
    #                       "```{r, echo=FALSE, results='hide', warning=FALSE, message=FALSE,fig.show='hide'}\n
    # tmpb <- plot(",nm,")\n",
    #                       "```{r, echo=FALSE}\nrayshader::plot_gg(tmpb[[1]],
    #         multicore = TRUE,width=5,height=5,scale=125,zoom=0.5, windowsize = c(1400,866),
    #         phi=30,theta=35)\nrayshader::render_snapshot(clear=T)\n",
    #                       "```\n\n")
    # } else {
    #   content_3d <- paste0(content_3d_head,"```{r, echo=FALSE}\nrayshader::plot_gg(tmp[[1]],
    #         multicore = TRUE,width=5,height=5,scale=125,zoom=0.5, windowsize = c(1400,866),
    #         phi=30,theta=35)\nrayshader::render_snapshot(clear=T)\n",
    #                        "```\n\n")
    # }


    content <- paste0(content,content_3d)

    if(!is.null(hotspot.opts)){
      content_3d2 <- paste0("### Hotspot probabilities\n\n",
                            "```{r, echo=FALSE}\nrayshader::plot_gg(tmp3[[2]],
            multicore = TRUE,width=5,height=5,scale=100,zoom=0.5, windowsize = c(1400,866),
            phi=30,theta=35)\nrayshader::render_snapshot(clear=T)\n",
                            "```\n\n")

      content <- paste0(content,content_3d2)
    }
  }

  write(content, paste0(dirname,"/",rmd_file_name))

  rmarkdown::render(paste0(dirname,"/",rmd_file_name))
}

#' MCMC convergence diagnostics
#'
#' Diagnostics for convergence of the MCMC chains from a call to \code{lgcp}.
#'
#' Produces a traceplot of the model parameters and prints R-hat and ESS statistics.
#'
#' @param lg Output from a call to \code{lgcp}
#' @return Traceplot of the MCMC chains of parameters from the linear predictor and
#' covariance function is plotted, and R-hat and ESS statistics are printed.
#' @export
convergence <- function(lg){
  if(class(lg)!="lgcpReal")stop("lg must be of class lgcpReal")
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
  p1 <- bayesplot::mcmc_trace(betalist)
  #p1 <- p1 + scale_colour_brewer(palette = "Set1")
  p2 <- bayesplot::mcmc_trace(etalist)
  #p2 <- p2 + scale_colour_brewer(palette = "Set1")

  print(ggpubr::ggarrange(p1,p2,ncol=1,heights = c(2,1)))

  #p3 <- bayesplot::mcmc_acf(betalist)
  #p4 <- bayesplot::mcmc_acf(etalist)

  # print(p3)
  # print(p4)

  beta_r <- coda::gelman.diag(betalist_mcmc)
  beta_ess <- coda::effectiveSize(betalist_mcmc)
  eta_r <- coda::gelman.diag(etalist_mcmc)
  eta_ess <- coda::effectiveSize(etalist_mcmc)
  res <- cbind(as.data.frame(beta_r$psrf[,1]),as.data.frame(beta_ess))
  res_e <- cbind(as.data.frame(eta_r$psrf[,1]),as.data.frame(eta_ess))
  colnames(res) <- c("R-hat","ESS")
  colnames(res_e) <- c("R-hat","ESS")


  print(knitr::kable(rbind(res,res_e),"simple",digits=3,options=list(knitr.kable.NA="")))
  return(invisible(list(knitr::kable(rbind(res,res_e),"simple",digits=3,options=list(knitr.kable.NA="")),
                        ggpubr::ggarrange(p1,p2,ncol=1,heights = c(2,1)))))
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
#' @export
vpc <- function(lg,covariates,rr=FALSE){
  OW <- lgcp::selectObsWindow(lg$xyt, cellwidth = lg$cellwidth)
  grid.data <- expand.grid(x=OW$xvals,y=OW$yvals)
  idx.mapping <- matrix(1:nrow(grid.data),nrow=length(OW$yvals),ncol=length(OW$xvals))
  idx.mapping <- c(t(apply(idx.mapping,2,rev)))

  if(!exists("outl") |(exists("outl")&&attr(outl, "dirname")!=lg$dirname)){
    print("Extracting posterior samples...")
    outl <- lgcpExtract(lg$dirname,nrow(lg$lgcpRunInfo$timetaken))
    assign("outl",outl,.GlobalEnv)
  }

  tmp <- suppressWarnings( .plot_lgcp_dat(outl,
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
      tmp$dat1$vpc[i] <- var(tmp$pop[i,]*tmp$linpred[i,]*re_samp,na.rm=T)/(
        var(tmp$pop[i,]*tmp$linpred[i,]*re_samp,na.rm=T)+
          mean(tmp$pop[i,]*tmp$linpred[i,]*re_samp,na.rm=T))
    } else {
      tmp$dat1$vpc[i] <- var(log(tmp$linpred[i,]))/(var(log(tmp$linpred[i,]))+var(log(re_samp)))
    }
  }

  res <- t(as.data.frame(quantile(tmp$dat1$vpc,c(0.025,0.25,0.5,0.75,0.975),na.rm=T)))
  vals <- res[,1:ncol(res)]
  res[,1:ncol(res)] <- paste0(round(res[,1:ncol(res)]*100,0),"%")
  rownames(res) <- "VPC"

  print(knitr::kable(res,"simple",digits=3,options=list(knitr.kable.NA="")))
  res[,1:ncol(res)] <- round(as.numeric(vals),3)
  return(invisible(res))
}
