test_that("data sim works.",{
  #.libPaths("C:/docs/R")
  #suppressWarnings(require(sp))
  square <- data.frame(x=c(0,0,1,1),y=c(0,1,1,0))
  square <- sp::SpatialPolygons(list(sp::Polygons(list(sp::Polygon(square)),"s1")))
  square <- as(square,"SpatialPolygonsDataFrame")
  square_pop <- expand.grid(x=seq(0,1,by=0.1),y=seq(0,1,by=0.1))
  square_pop <- sp::SpatialPointsDataFrame(coords= square_pop, data=square_pop)
  sp::gridded(square_pop) <- TRUE
  square_pop <- as(square_pop,"SpatialPolygons")
  square_pop <- as(square_pop,"SpatialPolygonsDataFrame")
  square_pop@data$popdens <- abs(rnorm(nrow(square_pop@data)))

  dat <- infecSim(region = square,
                  t.win = c(1,7),
                  covariates = square_pop,
                  mean.val= 10,
                  p =1/10,
                  delta = c(0.01,4),
                  rho=3,
                  t.off = 4,
                  cov.pars = c(0.9,0.03,0.1),
                  grid.size = 10^2)

  expect_s3_class(dat,"data.frame")

  mintest <- mincontrast_st(dat,
                            square_pop,
                            square)

  expect_type(mintest,"double")

  lg1 <- lgcp(data=dat,
              pop.var = c("popdens"),
              boundary=square,
              covariates=square_pop,
              cellwidth=0.1,
              laglength = 7,
              mala.pars=c(200,100,1),
              nchains=2,
              dirname="~/test")

  expect_s3_class(lg1,"lgcpReal")

  s1 <- summary(lg1, plot=FALSE)

  expect_type(s1,"list")
  expect_s3_class(s1,"lgcpRealSumm")

  p1 <- plot(lg1,
             covariates=square_pop)
  expect_s3_class(p1,"lgcpRealPlot")
  expect_s3_class(p1[[1]],"gg")

  p2 <- plot_hotspot(lg1,
                     covariates = square_pop,
                     threshold.var = c("poppp+obs+latent",
                                       "poppp+obs+latent+lag(3)"),
                     threshold.value = c(0.1,1),
                     threshold.prob=0.8,
                     labels=c('low','high incidence',
                              'rising incidence','both'))

  expect_s3_class(p2,"lgcpRealPlot")
  expect_s3_class(p2[[1]],"gg")

})

