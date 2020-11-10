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
  
  yc <- YAMLClass$new(list())
  
  parameter <- list()
  parameter$nonsense1 <- c(1,4)
  parameter$nonsense2 <- "test"
  parameter$param_PG_intThresh <- 30
  
  expect_null(createYaml(yc, param = parameter)$param$nonsense1)
  expect_null(createYaml(yc, param = parameter)$param$nonsense2)
  expect_equivalent(createYaml(yc, param = parameter)$param$param_PG_intThresh, 30)
  expect_equal(length(createYaml(yc)$param), 19)
  
  ##
  ##test valid parameter input
  ##
  
  yc <- YAMLClass$new(list())
  
  parameter$add_fs_col <- 14
  parameter$param_OutputFormats <- "txt"
  
  expect_equivalent(createYaml(yc, param = parameter)$param$add_fs_col, 14)
  expect_equivalent(createYaml(yc, param = parameter)$param$param_OutputFormats, "txt")
  
  ##test default values
  yc <- YAMLClass$new(list())
  
  expect_equivalent(createYaml(yc)$param$param_PG_intThresh, 25)
  expect_equivalent(createYaml(yc)$param$param_OutputFormats, c("html", "plainPDF"))

  
  ##
  ##test metrics deactivation (all except qcMetric_PAR)
  ##
  yc <- YAMLClass$new(list())
  
  mets <- "qcMetric_PAR" 
  expect_equal(createYaml(yc, param = parameter, metrics = mets)$yc$yamlObj$order$qcMetric_EVD_UpSet, -1)
  expect_equal(createYaml(yc, param = parameter, metrics = mets)$yc$yamlObj$order$qcMetric_PAR, 1)
  
  ##test no deactivation 
  yc <- YAMLClass$new(list())
  expect_false(any(createYaml(yc)$yc$yamlObj$order < 1))
  
})


