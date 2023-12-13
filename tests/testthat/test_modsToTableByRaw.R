context("modstoTableByRaw.R")

test_that("modstoTableByRaw", {
  dt_in = data.frame(  fc.raw.file = rep(c("file A", "file B"), each = 3),
                     modifications = c("A", "Unmodified", "A,B", "C", "Unmodified", "Unmodified"))
  r = modsToTableByRaw(dt_in)
 
  testthat::expect_identical(colnames(r), c('fc.raw.file', 'modification_names', 'Freq'))
  
  exp = data.frame(fc.raw.file = rep(c("file A",                     "file B"), times = c(3, 2)),
                   modification_names = factor(c("A", "B", "Unmodified",    "C", "Unmodified"), levels = levels(r$modification_names)),
                   Freq = c(200/3, 100/3, 100/3,                     100/3, 200/3))
  testthat::expect_equal(r, exp)
  
  dt_in = data.frame(  fc.raw.file = rep(c("file A", "file B"), each = 3),
                       modifications = c("C", "Unmodified", "Unmodified",    "Unmodified", "A", "Unmodified"))
  r_rename = modsToTableByRaw(dt_in, name_unmod = "Unmodified", name_unmod_inverse = "MOD_TOTAL")
  exp_rename = data.frame(fc.raw.file = rep(c("file A",                     "file B"), times = c(2, 2)),
                          modification_names = factor(c("C", "MOD_TOTAL",    "A", "MOD_TOTAL"), levels = levels(r_rename$modification_names)),
                          Freq = c(100/3, 100/3,                             100/3, 100/3))
  testthat::expect_equal(r_rename, exp_rename)

})
