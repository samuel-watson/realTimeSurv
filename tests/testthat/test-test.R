test_that("data sim works.",{

  data(dat,square,square_pop)

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
              nchains=2)

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

  unlink(paste0(lg1$dirname,".1"),recursive = TRUE)
  unlink(paste0(lg1$dirname,".2"),recursive = TRUE)
  file.remove(paste0(tempdir(),"\\outl.RDS"))

})

