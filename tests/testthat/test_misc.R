library(PTXQC)
context("fcn_misc.R")

test_that("delLCP", {
  ## using expect_equivalent(), since delLCP returns a named vector, but we only check the values
  expect_equivalent(delLCP(c("TK12345_H1"), min_out_length=0), "")
  expect_equivalent(delLCP(c("TK12345_H1"), min_out_length=4), "5_H1")
  expect_equivalent(delLCP(c("TK12345_H1"), min_out_length=4, add_dots = TRUE), "..5_H1")
  expect_equivalent(delLCP(c("TK12345_H1", "TK12345_H2"), min_out_length=4), c("5_H1", "5_H2"))
  expect_equivalent(delLCP(c("TK12345_H1", "TK12345_H2"), min_out_length=4, add_dots = TRUE), c("..5_H1", "..5_H2"))
  expect_equivalent(delLCP(c("TK12345_H1", "TK12345_H2"), min_out_length=8), c("12345_H1", "12345_H2"))
  
  ## (unchanged, even though its 10 chars, since '..' would add another two)
  expect_equivalent(delLCP(c("TK12345_H1", "TK12345_H2"), min_out_length=8, add_dots = TRUE), c("TK12345_H1", "TK12345_H2"))
  ## (unchanged)
  expect_equivalent(delLCP(c("TK12345_H1", "TK12345_H2"), min_out_length=60), c("TK12345_H1", "TK12345_H2"))
  ## (unchanged)
  expect_equivalent(delLCP(c("TK12345_H1", "TK12345_H2"), min_out_length=60, add_dots = TRUE), c("TK12345_H1", "TK12345_H2"))

})

test_that("delLCS", {
  ## using expect_equivalent(), since delLCP returns a named vector, but we only check the values
  expect_equivalent(delLCS(c("TK12345_H1")), "")
  expect_equivalent(delLCS(c("TK12345_H1", "TK12345_H2")), c("TK12345_H1", "TK12345_H2"))
  expect_equivalent(delLCS(c("TK12345_H1", "TK12!45_H1")), c("TK123", "TK12!"))
})

test_that("lcpCount", {
  expect_equal(lcpCount(c("TK12345_H1")), 10)
  expect_equal(lcpCount(c("TK12345_H1", "TK12345!_H2")), 7)
  expect_equal(lcpCount(c("TK12345_H1", "TK12345_H2")), 9)
})

test_that("lcsCount", {
  expect_equal(lcsCount(c("TK12345_H1")), 10)
  expect_equal(lcsCount(c("TK12345_H1", "TK12345!_H2")), 0)
  expect_equal(lcsCount(c("TK12345!_H2", "TK12345_H2")), 3)
})

test_that("LCS", {
  expect_equal(LCS("TK12345_H1", "TK12345!_H2"), "TK12345")
  expect_equal(LCS("PRE1_TK12345!_H2", "PRE2_TK12345_H2"), "_TK12345")
  expect_equal(LCS("PRE1_TK12345_H2", "PRE2_TK12345_H2"), "_TK12345_H2")
  expect_equal(LCS("TK12345_H1", ""), "")
  expect_equal(LCS("", "TK12345_H1"), "")
})

test_that("LCSn", {
  expect_equal(LCSn(c("1_abcde...")), "1_abcde...")
  expect_equal(LCSn(c("123")), "123")
  expect_equal(LCSn(c("123"), min_LCS_length= 7), "")
  expect_equal(LCSn(c("1_abcde...", "2_abcd...", "x_abc...")), "_abc")
  expect_equal(LCSn(c("1_abcde...", "2_abcd...", "x_abc..."), min_LCS_length= 7), "")
})


test_that("simplifyNames", {
  expect_equal(simplifyNames(c('TK20130501_H2M1_010_IMU008_CISPLA_E3_R1.raw')),
               c("TK20130501_H2M1_010_IMU008_CISPLA_E3_R1.raw"))
  expect_equal(simplifyNames(c('TK20130501_H2M1_010_IMU008_CISPLA_E3_R1.raw', 'TK20130501_H2M1_026_IMU008_CISPLA_E7_R2.raw'), infix_iterations = 2),
               c("TK.._010_I.._E3_R1.raw","TK.._026_I.._E7_R2.raw"))
  expect_equal(simplifyNames(c('TK20130501_H2M1_010_IMU008_CISPLA_E3_R1.raw', 'TK20130501_H2M1_026_IMU008_CISPLA_E7_R2.raw'), infix_iterations = 1),
               c("TK.._010_IMU008_CISPLA_E3_R1.raw","TK.._026_IMU008_CISPLA_E7_R2.raw"))
  
  ## could be shortened to c("1ab..gh", "2ab..gh"), but thats too short...
  expect_equal(simplifyNames(c("1abcdefgh", "2abcdefgh"), min_out_length=99), c("1abcdefgh", "2abcdefgh"))
  expect_equal(simplifyNames(c("1abcdefgh", "2abcdefgh"), min_out_length=3), c("1ab..gh", "2ab..gh"))
  
  inp = unlist(strsplit("3265_Lumik3_Shakya_150IS_HCD_11A_1 3265_Lumik3_Shakya_150IS_HCD_11B_1 3265_Lumik3_Shakya_150IS_HCD_11C_1 3265_Lumik3_Shakya_150IS_HCD_11D_1 3265_Lumik3_Shakya_150IS_HCD_12A_1 3265_Lumik3_Shakya_150IS_HCD_12B_1 3265_Lumik3_Shakya_150IS_HCD_12C_1 3265_Lumik3_Shakya_150IS_HCD_12D_1 3265_Lumik3_Shakya_150IS_HCD_13A_1 3265_Lumik3_Shakya_150IS_HCD_13B_1 3265_Lumik3_Shakya_150IS_HCD_13C_1 3265_Lumik3_Shakya_150IS_HCD_13D_1 3265_Lumik3_Shakya_150IS_HCD_14A_1 3265_Lumik3_Shakya_150IS_HCD_14B_1 3265_Lumik3_Shakya_150IS_HCD_14C_1 3265_Lumik3_Shakya_150IS_HCD_14D_1", " "))
  expc = unlist(strsplit("32.._11A_1,32.._11B_1,32.._11C_1,32.._11D_1,32.._12A_1,32.._12B_1,32.._12C_1,32.._12D_1,32.._13A_1,32.._13B_1,32.._13C_1,32.._13D_1,32.._14A_1,32.._14B_1,32.._14C_1,32.._14D_1", ","))
  expect_equal(simplifyNames(inp, infix_iterations = 2), expc)

  
  expect_equal(simplifyNames(c("bla", "foo"), min_LCS_length=99), c("bla", "foo"))
  ## internally, min_LCS_length>=6 is required
  expect_error(simplifyNames(c("bla", "foo"), min_LCS_length=5))
  
})

test_that("shortenStrings", {
  ## no not allow duplicates (since they are already ambiguous, even without shortening)
  dups = c("someString", "aDuplicate", "aDuplicate")
  expect_error(shortenStrings(dups))
  expect_equivalent(shortenStrings(dups, allow_duplicates = TRUE), dups)
  
  expect_equivalent(shortenStrings(c("gamg_101", "gamg_101230100451", "jurkat_06_100731121305", "jurkat_06_1"), max_len = 20),
                    c("gamg_101", "gamg_101230100451", "jurkat_06_10073112..", "jurkat_06_1"))

})

test_that("supCount", {
  expect_equal(supCount(c("abcde...", "abcd...", "abc...")), 5)

  x = c("doubled", "doubled", "aLongDummyString")
  ## due to duplicated entries there is no prefix which makes them unique
  expect_equal(substr(x, 1, supCount(x)), x)
})



test_that("correctSetSize", {
  expect_equal(correctSetSize(8, 5), 4)
  expect_equal(correctSetSize(101, 25), 26)
})

test_that("byX", {
  expect_equivalent(byX(data.frame(d=1:10), 1:10, 1, sum, sort_indices = FALSE), as.list(1:10))
  expect_equivalent(byX(data.frame(d=1:10), 1:10, 1, sum, sort_indices = TRUE), as.list(1:10))
  expect_equivalent(byX(data.frame(d=1:10), as.character(1:10), 1, sum, sort_indices = TRUE), as.list(c(1, 10, 2:9)))
  expect_equivalent(byX(data.frame(d=1:10), 1:10, 2, sum, sort_indices = FALSE), as.list(c(1+2,3+4,5+6,7+8,9+10)))
  expect_equivalent(byX(data.frame(d=1:10), 1:10, 7, sum, sort_indices = FALSE), as.list(c(sum(1:7), sum(8:10))))
})

test_that("byXflex", {
  expect_equivalent(byXflex(data.frame(d=1:10), 1:10, 1, sum, sort_indices = FALSE), as.list(1:10))
  expect_equivalent(byXflex(data.frame(d=1:10), 1:10, 1, sum, sort_indices = TRUE), as.list(1:10))
  expect_equivalent(byXflex(data.frame(d=1:10), as.character(1:10), 1, sum, sort_indices = TRUE), as.list(c(1,10, 2:9)))
  expect_equivalent(byXflex(data.frame(d=1:10), 1:10, 2, sum, sort_indices = FALSE), as.list(c(1+2,3+4,5+6,7+8,9+10)))
  ## within 150% (7*1.5), just create one group
  expect_equivalent(byXflex(data.frame(d=1:10), 1:10, 7, sum, sort_indices = FALSE), as.list(c(sum(1:10))))
  ## two groups are maintained, but instead of 6, 4 we use 5,5 for more balanced group size
  expect_equivalent(byXflex(data.frame(d=1:10), 1:10, 6, sum, sort_indices = FALSE), as.list(c(sum(1:5), sum(6:10))))
})

test_that("getMaxima", {
  expect_equivalent(getMaxima(c(1,0,3,4,5,0)), c(TRUE,FALSE,FALSE,FALSE,TRUE,FALSE))
  expect_equivalent(getMaxima(c(1,0,3,4,5,0), thresh_rel=0.5), c(FALSE,FALSE,FALSE,FALSE,TRUE,FALSE))
})

thinOut = function(data, filterColname, binsize)
test_that("thinOut", {
  data = data.frame(x = c(1:100), y = 1:10)
  expect_equal(dim(thinOut(data, "x", 5)), c(21,2))
  # duplicates in 'x' column are irrelevant
  data = data.frame(x = c(1:100, 1:10), y = 1:110)
  expect_equal(dim(thinOut(data, "x", 5)), c(21,2))
})
