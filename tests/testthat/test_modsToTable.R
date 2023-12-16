context("modstoTable.R")

test_that("modstoTable", {
  r = modsToTable(c("A", "Unmodified", "A,B", "C", "Unmodified", "Unmodified"))
 

  exp = data.frame(modification_names = factor(c("A", "B", "C", "Unmodified")), 
                   Freq = c(200/6, 100/6, 100/6, 300/6))
  testthat::expect_equal(r, exp)
})