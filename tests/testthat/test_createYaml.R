library(PTXQC)
context("createYaml.R")

test_that("createYaml", {

  ##
  ##test empty yc object, no parameter or metric specifications
  ##
  yc <- YAMLClass$new(list())
  expect_equal(length(createYaml(yc)$param), 19)

  ##
  ##test invalid parameter input
  ##
  parameter <- list()
  parameter$nonsense <- c(1,4)
  expect_equivalent(createYaml(yc, param = parameter)$param$nonsense, NULL)
  
  ##
  ##test valid parameter input
  ##
  yc <- YAMLClass$new(list())
  parameter$param_PG_intThresh <- 30
  parameter$param_OutputFormats <- "txt"
  expect_equivalent(createYaml(yc, param = parameter)$param$param_PG_intThresh, 30)
  expect_equivalent(createYaml(yc, param = parameter)$param$param_OutputFormats, "txt")
  expect_equivalent(createYaml(yc, param = parameter)$param$add_fs_col, 14)
  
  ##
  ##test metrics deactivation (all except qcMetric_PAR)
  ##
  yc <- YAMLClass$new(list())
  mets <- "qcMetric_PAR" 
  expect_equal(createYaml(yc, param = parameter, metrics = mets)$yc$yamlObj$order$qcMetric_EVD_UpSet, -1)
  expect_equal(createYaml(yc, param = parameter, metrics = mets)$yc$yamlObj$order$qcMetric_PAR, 1)
  
  
  
})
