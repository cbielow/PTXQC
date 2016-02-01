library(PTXQC)
#library(testthat)

context("qcMetric.R")

test_that("qcMetric", {
  
  dd = data.frame(x=1:10, y=11:20)
  
  a = qcMetric$new(helpText="small help text", 
                   workerFcn=function(.self, data, gtit)
                   {
                     # usually, plots are produced, but they are hard to check. So we generate simple values
                     pl = lapply(1:2, function(xx) ggplot(data) + geom_point(aes(x=x*xx,y=y)) + ggtitle(paste(gtit, xx)))
                     qcScores = data.frame(r=1, z=9:10)
                     return(list(plots = pl, qcScores = qcScores))
                   }, 
                   qcCat="LC", 
                   qcName="MS/MS Peak shape", 
                   orderNr = 30)
  
  a$setData(dd, "title assigned by worker")
  s = a$qcScores
  expect_equivalent(s, data.frame(r=1, z=9:10))

  p = a$plots
  expect_equal(length(p), 2)
  a$getPlots(TRUE) ## just check if it works
  a$getPlots(FALSE) ## just check if it works
  expect_equal(paste("title assigned by worker", 1:2), a$getTitles())
  
  expect_equal("small help text", a$helpText)
  expect_equal("MS/MS Peak shape", a$qcName)
  
  expect_equal(30, a$orderNr)
})

