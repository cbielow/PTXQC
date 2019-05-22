library(PTXQC)
context("fcn_qualities.R")

test_that("qualLinThresh", {
  expect_error(qualLinThresh(-3))
  expect_error(qualLinThresh(c(1,2,-3)))
  expect_equal(qualLinThresh(3, 10), 0.3)
  expect_equal(qualLinThresh(3), 1)
  expect_equal(qualLinThresh(3, 3), 1)
  expect_equal(qualLinThresh(c(0,1,3,NA), 3), c(0, 1/3, 1, 0))
})

test_that("qualCentered", {
  expect_equal(qualCentered(3), 0)
  ## median is 0 -- perfect 
  expect_equal(qualCentered(-5:5), 1)
  ## some outliers, but still median is 0 -- perfect
  expect_equal(qualCentered(c(-5:5, -9999, 10)), 1)
  ## median 2, max=3: 1 - 2/3
  expect_equal(qualCentered( c(1,2,3)), 1 - 2/3)
  expect_equal(qualCentered(-c(1,2,3)), 1 - 2/3)
  ## further off ...
  expect_equal(qualCentered( c(5,6,7)), 1 - 6/7)
  expect_equal(qualCentered(-c(5,6,7)), 1 - 6/7)
})

test_that("qualCenteredRef", {
  expect_error(qualCenteredRef(3, 0))
  expect_error(qualCenteredRef(3, -1))
  
  ## median is 13
  expect_equal(qualCenteredRef(c(10:16), 13), 0)
  expect_equal(qualCenteredRef(c(10:16), 13*2), 0.5)
  expect_equal(qualCenteredRef(c(10:16), 13*4), 0.75)
  
  ## median is 13, ref-interval is smaller -- worst score
  expect_equal(qualCenteredRef(c(10:16), 12), 0)
  ## median is perfectly 0 -- ref-interval does not really matter
  expect_equal(qualCenteredRef(-5:5, 1), 1)
  expect_equal(qualCenteredRef(-5:5, 100), 1)
})

test_that("qualMedianDist", {
  ## values need be in [0,1]
  expect_error(qualMedianDist(-3))
  expect_error(qualMedianDist(c(0.1,0.3,-0.3)))
  
  ## median: 0.5
  x = (1:9)/10
  expect_equal(qualMedianDist(x), 1 - abs(0.5 - x))
})

test_that("qualUniform", {
 expect_equal(qualUniform(c(3,3,3)), 1)
 expect_equal(qualUniform(c(4,0,0)), 0)         

 ## how 'uniform' is a vector where only a single index has weight?-- answer: very
 expect_equal(qualUniform(c(4,0,0), c(1,0,0)), 1)
 expect_equal(qualUniform(c(4,0,0), c(0,1,0)), 1)
 expect_equal(qualUniform(c(0,4,0)), 0)
 expect_lt(abs(qualUniform(c(3,2,1))-0.58578), 0.0001)
 expect_lt(abs(qualUniform(c(1,2,3))-0.58578), 0.0001)
 expect_equal(qualUniform(c(1,2,3), c(0,1,0)), 1)
 expect_lt(abs(qualUniform(c(1,2,3))-0.58578), 0.0001)
 expect_lt(abs(qualUniform(c(1,2,3), c(0,1,1))- 0.590316), 0.0001)
 expect_lt(abs(qualUniform(c(2,3), c(1,1))-0.552786), 0.0001)
 expect_lt(abs(qualUniform(1:120)-0.38661), 0.0001)
})


test_that("qualGaussDev", {
  expect_error(qualGaussDev(-3))
  expect_equal(qualGaussDev(0, -4), NaN)
  
  ## perfect
  expect_equal(qualGaussDev(0, 1), 1)
  expect_equal(qualGaussDev(0, 3), 1)
  
  ## shifted slighty
  expect_equal(qualGaussDev(1, 1), 0.6065307, tolerance = .001)
  ## shifted more...
  expect_equal(qualGaussDev(2, 1), 0.1353353, tolerance = .001)
  
})

test_that("qualHighest", {
  expect_equal(qualHighest(c(0,0,0,16), 4), 1)
  expect_equal(qualHighest(c(16,0,0,0), 4), 0)
  expect_equal(qualHighest(c(1,1,1,1), 4) , 0.5)
  expect_equal(qualHighest(c(0,16,0,0), 4), 1/3)
})


test_that("qualBestKS", {
  ## sample some Gaussians
  g1 = rnorm(100)
  g2 = rnorm(100)
  g3 = rnorm(100)
  go = rnorm(100, mean=4) ## this will be our outlier distribution
  x = list(g1,g2,g3,go)
  
  r = qualBestKS(x)
    
  expect_equal(max(r$ks_best[1:3]), 1)
  expect_true(all(r$ks_best[1:3] > 0.7)) ## the three similar ones, should score good (one of them the the reference)
  expect_lt(r$ks_best[4], 0.2) ## outlier should score badly
})


