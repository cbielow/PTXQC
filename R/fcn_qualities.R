
#'
#' Quality metric with linear response to input, reaching the maximum score at the given threshold.
#' 
#' Ranges between 0 (worst score) and 1 (best score).
#' Useful for performance measures where reaching a certain reference threshold 't'
#' will be enough to reach 100\%.
#' The input range from [0, t] is scored from 0-100\%.
#' 
#' @param x Numeric value(s) between [0, inf]
#' @param t Threshold value, which indicates 100\%
#' @return Value between [0, 1]
#' 
#' @export
#' 
qualLinThresh = function(x, t = 1)
{
  if (length(x) == 0) stop("Error: qualLinThresh() received empty data!")
  
  x[is.na(x)] = 0; ## replace NA with 0
  if (any(x < 0)) stop("qualLinThresh(): negative input 'x' not allowed!")
  if (t < 0) stop("qualLinThresh(): negative threshold 't' not allowed!")
  
  q = pmin(1, x / t)
  return (q)
}


#'
#' Quality metric for 'centeredness' of a distribution around zero.
#' 
#' Ranges between 0 (worst score) and 1 (best score).
#' A median of zero gives the best score of 1.
#' The closer the median is to the most extreme value of the distribution, the smaller the score (until reaching 0).
#' Can be used for calibrated mass errors, as a measure of how well they are centered around 0.
#' E.g. if the median is 0.1, while the range is [-0.5,0.5], the score will be 0.8 (punishing the 20% deviation).
#' If the range of data is asymmetric, e.g. [-1.5,-0.5] and does not include zero, the score cannot reach 1, 
#' since the median can never be zero.
#' 
#' @param x Numeric values (e.g. ppm errors)
#' @return Value between [0, 1]
#' 
#' @export
#' 
qualCentered = function(x)
{
  if (length(x) == 0) stop("Error: qualCentered() received empty data!")
  if (any(is.na(x))) stop("Error: qualCentered() received NA data!")
  
  q = 1 - (abs(median(x)) / max(abs(x)))
  return (q)
}

#'
#' Quality metric for 'centeredness' of a distribution around zero with a user-supplied range threshold.
#' 
#' Ranges between 0 (worst score) and 1 (best score).
#' The best score is achieved when the median of 'x' is close to the center of the interval [-tol, tol].
#' If median of 'x' is close to the border (on either side), the score decreases linearly to zero.
#' Can be used for uncalibrated mass errors, as a measure of how well they are centered around 0.
#' 
#' NA's are removed for all computations.
#' 
#' @param x  Vector of values (hopefully in interval [-tol, tol])
#' @param tol Border of interval (must be positive)
#' @return Value between [0, 1]
#'
#' @export
#' 
qualCenteredRef = function(x, tol)
{
  x = na.omit(x)
  if (length(x) == 0) stop("qualCenteredRef(): input is empty or consists of NA only!")
  if (tol <= 0) stop("qualCenteredRef(): negative or zero interval border not allowed!")
  m = abs(median(x, na.rm = TRUE))
  if (m > tol) warning("qualCenteredRef(): Median of x is outside of interval. Score will be set to 0.")
  q = 1 - (m / tol)
  q = max(0, q) ## avoid negative scores if abs(x)>tol
  return (q)
}

#'
#' Scales a vector of values linearly to [0, 1]
#' If all input values are equal, returned values are all 0
#' @param X Vector of values
#' @return Scaled vector
scale01linear = function(X)
{
  rng = (max(X) - min(X));
  if (rng == 0) rng = 1; ## no range, so divide by 1
  X2 = (X - min(X)) / rng
  #summary(X2)
  return(X2)
}


#'
#' Quality metric which measures the absolute distance from median.
#' 
#' Ranges between 0 (worst score) and 1 (best score).
#' Input must be between [0,1].
#' Deviations from the median of the sample represent the score for each sample point.
#' 
#' @param x A vector numeric values between [0,1]
#' @return A vector of the same size as x, with quality values between [0, 1]
#'
#' @export
#'
qualMedianDist = function(x)
{
  if (length(x) == 0) stop("Error: qualMedianDist() received empty data!")
  if (any(is.na(x))) stop("Error: qualMedianDist() received NA data!")
  if (any(x < 0 | x > 1)) stop("qualMedianDist(): x values out of range [0,1]!")
  q = 1 - abs(x - median(x, na.rm = TRUE))
  return (q)
}


#'
#' Compute deviation from uniform distribution
#' 
#' Ranges between 0 (worst score) and 1 (best score).
#' Input 'x' is a vector of counts (or probabilities) for equally spaced bins in a histogram.
#' A uniform distribution (e.g. c(3,3,3) will get a score of 1. The worst possible case (e.g. c(4,0,0)), will get a score of 0,
#' and a linear increasing function (e.g. c(1,2,3)) will get something in between (0.585 here)
#' 
#' In addition, bin values can be weighted (e.g. by their confidence). The total sum of weights is normalized to 1 internally.
#' 
#' The distance function used is the square root of the absolute difference between a uniform distribution and the input 'x'
#' (summed for each element of 'x').
#' This distance is normalized to the worst possible input (e.g. one bin with 100% counts, all other bins being empty).
#' 
#' @param x Vector of numeric intensity/count values (e.g. ID's per RT bin); bins are assumed to have equal widths
#' @param weight Vector of weights for values in 'x' (same length as 'x').
#' @return Value between [0, 1]
#' 
#' @export
#' 
#' @examples 
#'  stopifnot(qualUniform(c(3,3,3))==1)
#'  stopifnot(qualUniform(c(4,0,0))==0)         
#'
#'  ## how 'uniform' is a vector where only a single index has weight?-- answer: very
#'  stopifnot(qualUniform(c(4,0,0), c(1,0,0))==1)   
#'  stopifnot(qualUniform(c(4,0,0), c(0,1,0))==1)     
#'  stopifnot(qualUniform(c(0,4,0))==0)              
#'  stopifnot(abs(qualUniform(c(3,2,1))-0.58578) < 0.0001)
#'  stopifnot(abs(qualUniform(c(1,2,3))-0.58578) < 0.0001)
#'  stopifnot(qualUniform(c(1,2,3), c(0,1,0))==1)   
#'  stopifnot(abs(qualUniform(c(1,2,3))-0.58578) < 0.0001)
#'  stopifnot(abs(qualUniform(c(1,2,3), c(0,1,1))- 0.590316) < 0.0001)
#'  stopifnot(abs(qualUniform(c(2,3), c(1,1))-0.552786) < 0.0001)
#'  stopifnot(abs(qualUniform(1:120)-0.38661) < 0.0001)
#'  
qualUniform = function(x, weight=vector())
{
  if (length(x) == 0 || any(is.na(x) | (x<0))) stop("Error: qualUniform() received negative or missing data!")
  
  if(length(weight)==0) {
    weight=rep(1, length(x)) 
  } else if (length(weight)!=length(x)) {
    stop("Error: qualUniform() received weights which do not match length of input!")
  }
  weight = weight/sum(weight)  ## normalize to 1
  x = x/sum(x)

  n = length(x)
  #y = 1/n ## expected intensity (no weights)
  y = sum(x*weight) ## expected weighted intensity
  y
  p = 1/2
  ## for worst possible score, assume bin with highest weight has all the counts
  ## i.e. error is 1-expected == 1-y for 100% bin, and abs(0-y)=y for all other bins of zero intensity
  idx_maxWeight = which(weight == max(weight))[1]
  worst = ((1-y)^p)*weight[idx_maxWeight] + (y^p)*sum(weight[-idx_maxWeight])
  ## worst ==0, iff weight-vector has only one non-zero entry, e.g. c(0,1,0)
  #worst = ((1-y)^p) + (n-1)*((1/n)^p)
  worst  
  ## score of current distribution
  sc = sum( (abs(x-y)^p) * weight)
  sc
  q = ifelse(worst==0, 1, (worst - sc) / worst)
  return (q)    
}


#' 
#' Compute probability of Gaussian (mu=m, sd=s) at a position 0, with reference 
#' to the max obtainable probability of that Gaussian at its center.
#' 
#' Measure for centeredness around 0.
#' Highest score is 1, worst score is 0.
#' 
#' @param mu Center of Gaussian
#' @param sd SD of Gaussian
#' @return quality, ranging from 0 (bad agreement) to 1 (perfect, i.e. centered at 0)
#' 
#' @export
#' 
qualGaussDev = function(mu, sd)
{
  q = dnorm(0, mean = mu, sd = sd) / dnorm(mu, mean = mu, sd = sd)
  return (q)
}

#'
#' Score an empirical density distribution of values, where the best possible distribution is right-skewed.
#' 
#' The score is computed according to
#' 
#' q = ((N-1) - sum_i( ((N-i-1)*x_i) ) / (N-1)
#' 
#' Scores range from 0 (worst), to 1 (best).
#' E.g. c(0,0,0,16) would yield a score of 1.
#' c(16,0,0,0,0) gives a score of 0.
#' 
#' @param x Vector of numeric values (e.g. height of histogram bins)
#' @param N Length of x (just a precaution currently)
#' @return Quality score in the range of [0,1]
#' 
#' @export
#' 
#' @examples
#'  qualHighest(c(0,0,0,16), 4)   ## 1
#'  qualHighest(c(16,0,0,0), 4)   ## 0 
#'  qualHighest(c(1,1,1,1), 4)    ## 0.5
#'  qualHighest(c(0,16,0,0), 4)   ## 1/3
#'  
qualHighest = function(x, N)
{
  if (length(x) == 0) stop("Error: qualHighest() received empty data!")
  if (any(is.na(x))) stop("Error: qualHighest() received NA data!")
  if (length(x)!=N) stop("Error in qualHighest(): length of x not equal to N")
  
  # normalize
  x = x / sum(x)
  p = 1 ## exponent needs to be >=1, using larger values will only mildly penalize intermediate values
  worst = (N-1)^p
  penalty = sum( (x*(N:1 - 1))^p )
  q = (worst - penalty) / worst
  return (q)
}


#'
#' From a list of vectors, compute all vs. all Kolmogorov-Smirnoff distance statistics (D)
#' 
#' ... and report the row of the matrix which has maximum sum (i.e the best "reference" distribution).
#' The returned data.frame has as many rows as distributions given and two columns.
#' The first column 'name' gives the name of the list element, the second column 'ks_best' gives '1-statistic' of the
#' Kolmogorov-Smirnoff test to the "reference" distribution (which was picked by maximising the sum of 'ks_best').
#' Thus, the row with a 'ks_best' of 1 is the reference distribution.
#' 
#' @param x List of vectors, where each vector holds a distribution
#' @return A data.frame with ks-test values of the "reference" to all other distributions (see Details)
#'
#' @import stats
#'
#' @export
#'
qualBestKS = function(x) {
  
  if (class(x) != "list") stop("Function bestKS() expects a list!")
  if (length(x) == 0) stop("Error: qualBestKS() received empty data!")

  if (is.null(names(x))) names(x) = 1:length(x)
  
  ## result matrix
  rr = matrix(0, nrow=length(x), ncol=length(x))
  ## compute upper diagonal
  for (i in 1:length(x))
  {
    for (j in (i+1):length(x))
    {
      if (j>length(x)) next;
      rr[i,j] = 1 - suppressWarnings(  ## brags about '-value will be approximate in the presence of ties'
                      ks.test(x[[i]], x[[j]])$statistic
                    )    
    }
  }
  ## fill diagonal with 1's
  diag(rr) = 1
  ## add up all values for each reference distribution (equivalent to just rowSums on full matrix)
  rr_sums = rowSums(rr, na.rm = TRUE) + colSums(rr, na.rm = TRUE)
  ## pick best
  r_max = which.max(rr_sums)
  ## get values
  v = pmax(rr[r_max,], rr[, r_max])
  
  return (data.frame(name = names(x), ks_best = v))
         
}



