#'
#' Quality metric for 'centeredness' of a distribution
#' 
#' Ranges between 0 (worst score) and 1 (best score).
#' Can be used for calibrated mass errors, as a measure of how well they are centered around 0.
#' 
#' @param x Numeric values (e.g. ppm errors)
#' @return Value between [0, 1]
#' 
#' 
qualCentered = function(x)
{
  q = 1 - (abs(median(x)) / max(abs(x)))
  return (q)
}

#'
#' Quality metric for 'centeredness' of a distribution with maximum reference
#' 
#' Ranges between 0 (worst score) and 1 (best score).
#' The best score is achieved when 'x' is close to the center of the interval (given by 'tol').
#' If 'x' is close to the border (on either side), the score decreases linearly to zero.
#' Can be used for uncalibrated mass errors, as a measure of how well they are centered around 0.
#' 
#' @param x  Position in interval [-tol, tol]
#' @param tol Border of interval
#' @return Value between [0, 1]
#' 
#' 
qualCenteredRef = function(x, tol)
{
  if (tol <=0) stop("qualCenteredRef(): negative or zero interval border not allowed!")
  if (abs(x) > tol) warning("qualCenteredRef(): Position is outside of interval. Score will be set to 0.")
  q = 1 - (abs(x) / tol)
  q = max(0, q) ## avoid negative scores if abs(x)>tol
  return (q)
}

#'
#' Compute deviation from uniform distribution
#' 
#' Ranges between 0 (worst score) and 1 (best score).
#' A uniform distribution (e.g. c(3,3,3) will get a score of 1. The worst possible case (e.g. c(4,0,0)), will get a score of 0,
#' and a linear increasing function (e.g. c(1,2,3)) will get something in between (0.46 here)
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
#'  stopifnot(qualUniform(c(4,0,0), c(1,0,0))==1) 
#'  stopifnot(qualUniform(c(4,0,0), c(0,1,0))==1)     
#'  stopifnot(qualUniform(c(0,4,0))==0)              
#'  stopifnot(abs(qualUniform(c(3,2,1))-0.58578) < 0.0001)
#'  stopifnot(abs(qualUniform(c(1,2,3))-0.58578) < 0.0001)
#'  stopifnot(qualUniform(c(1,2,3), c(0,1,0))==1)   
#'  stopifnot(abs(qualUniform(c(1,2,3))-0.58578) < 0.0001)
#'  stopifnot(abs(qualUniform(c(1,2,3), c(0,1,1))- 0.56949) < 0.0001)
#'  stopifnot(abs(qualUniform(c(2,3), c(1,1))-0.552786) < 0.0001)
#'  stopifnot(abs(qualUniform(1:120)-0.38661) < 0.0001)
#'  
qualUniform = function(x, weight=vector())
{
  if (any(is.na(x) | (x<0))) stop("Error: qualUniform() received negative or missing data!")
  
  if(length(weight)==0) {
    weight=rep(1, length(x)) 
  } else if (length(weight)!=length(x)) {
    stop("Error: qualUniform() received weights which do not match length of input!")
  }
  weight = weight/sum(weight)  ## normalize to 1
  x = x/sum(x)

  n = length(x)
  #y = 1/n ## expected intensity (no weights)
  y = sum(x*weight) ## weighted intensity (expected)
  y
  p = 1/2
  ## bin with highest weight has 100% intensity
  idx_maxWeight = which(weight == max(weight))[1]
  worst = ((1-y)^p)*weight[idx_maxWeight] + ((1/n)^p)*sum(weight[-idx_maxWeight])
  #worst = ((1-y)^p) + (n-1)*((1/n)^p)
  worst  
  ## score of current distribution
  sc = sum( (abs(x-y)^p) * weight)
  sc
  q = ifelse(worst==0, 1, (worst - sc) / worst)
  return (q)    
}


#'
#' Test for uniform distribution using Kolmogorov-Smirnoff
#'
#' If only 'x' is given, a one sample KS test is done (min and max are taken from 'x').
#' If both 'x' and 'y' are given, a two-sample KS test is conducted.
#' 
#' Thoughts: this test has multiple problems:
#'      - is much too stringent, leading to p-values  ~ 0 very quickly
#'      - looking at the 'D' statistic instead is also not good:
#'            - using a Gaussian, centered at the middle of the data range, gives
#'              D~0.17, i.e. q=1-D=0.83, which seems to high
#'            - if only one bin dominates, its position strongly influences D:
#'              e.g. 
#'                (all data at the right (or left)):
#'                ks.test(c(100,100,rep(0,100)), y="punif", min=min(x), max=max(x)) ==> d=0.98 (bad fit, thus good metric)
#'                (all data at the center):
#'                ks.test(c(0,100,rep(50,100)), y="punif", min=min(x), max=max(x))  ==> d=0.57 (similarly bad fit, but metric just does not reflect it)
#'
#' @param x Vector of values from distribution 1
#' @param y Vector of values from distribution 2 (or assumed uniform if omitted)
#' @return p.value of KS test
#' 
qualUnifKS = function(x, y = NULL)
{
  if (is.null(y)) return (ks.test(x, y="punif", min=min(x), max=max(x))$p.value)
  else return (ks.test(x, y)$p.value)
}

#'
#' Discete ChiSquare-test on raw (unbinned) data
#' 
#' Data are first assigned to bins of equal width using 'hist()'.
#' 
#' If only 'x' is given, its distribution of compared to a uniform distribution.
#' If 'y' is given as well, 'x' vs. 'y' is compared (after binning).
#' If counts of 'x' and 'y' differ, the binned counts of the vector with higher counts 
#' are rescaled by N/M to match the smaller one (this might lead to bias).
#' 
#' 
#'
#'
# qualUniform_C2 = function(x, y = NULL, bin_count = 30)
# {
#   
#   if (!is.null(y)) {
#     xy = c(x,y)
#     h = hist(c(x,y), breaks = bin_count, plot = FALSE)
#     hx = hist(x, breaks = h$breaks, plot = FALSE)
#     hy = hist(y, breaks = h$breaks, plot = FALSE)
#     if (length(x) != length(y)) { ## rescale
#       fac = length(x) / length(y)
#     }
#   }
#   else xy = x;
#   
#   
#   
#   chisq.test
#   
# }

#' 
#' Compute position of each element in a vector, assuming the rest of the data is a Gaussian
#' 
#' For each index i, compute a Gaussian G from x \\ x_i and assess the prob of x_i in G.
#' The best result is 1, i.e. x_i is in the center of all other values.
#' The worst results is 0 (for very far outliers), i.e. the value is very far from the median.
#' 
#' @param x Numeric vector of values (e.g. counts for charge state of 2)
#' @return Numeric vector of qualities, same length of 'x', ranging from 0 (bad agreement) to 1 (perfect)
#' 
qualGauss = function(x)
{
  r = rep(NA, length(x))
  for (i in 1:length(x))
  {
    med = median(x[-i])
    gsd = 0.05 * med; # constant 5%, otherwise some datapoint is always bad;  old: "mad(x[-i], center = med)"
    p_best = dnorm(med, mean = med, sd = gsd)
    p_obs =  dnorm(x[i], mean = med, sd = gsd)
    r[i]  = p_obs/p_best  ## ratio of current position divided by best possible position (=center)
  }
  return(r)
}


#'
#' Score a distribution of values, where the best possible distribution is right-skewed.
#' 
#' The score is computed according to
#' 
#' q = ((N-1) - sum_i( ((N-i-1)*x_i) ) / (N-1)
#' 
#' Scores range from 0 (worst), to 1 (best).
#' E.g. c(0,0,0,16) would yield a score of 1.
#' c(16,0,0,0,0) gives a score of 0.
#' 
#' @param x Vector of numeric values (height of histogram bins)
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
#' Given a vector of values sampled from a mixture of two Gaussians, compute the area contribution
#' of the higher/more abundant (=target) Gaussian.
#' 
#' @note This is pre-tuned for Match-time-differences
#'
#' @param x Vector of numeric values
#' @return Ratio of area of the target distribution vs. whole area under the two Gaussians
#'
#' @importFrom mixtools normalmixEM2comp
#' 
#' 
gauss2Mix = function(x)
{
  d = na.omit(x)
  fit2 = normalmixEM2comp(d, c(0.3,0.6), c(0,0), c(1.5,0.1), eps= 1e-8, maxit = 1000, verb=FALSE)
  #   fit2$mu
  #   fit2$sigma
  #   fit2$lambda
  #   h=hist(d, 100, plot=F)
  #   x = seq(-4,4,by=0.1)
  #   plot(h$mids, h$density)
  #   lines(x, dnorm(x, fit2$mu[1], fit2$sigma[1]) *fit2$lambda[1], col="red")
  #   lines(x, dnorm(x, fit2$mu[2], fit2$sigma[2]) *fit2$lambda[2], col="green")
  #   lines(x, dnorm(x, fit2$mu[1], fit2$sigma[1]) *fit2$lambda[1] + dnorm(x, fit2$mu[2], fit2$sigma[2]) *fit2$lambda[2])
  ## look at scaling of the one with better SD (order defined above)
  q = fit2$lambda[2] / sum(fit2$lambda)
  return (q)
}



#'
#' From a list of vectors, compute all vs. all Kolmogorov-Smirnoff statistics (D)
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
#'
bestKS = function(x) {
  
  if (class(x) != "list") stop("Function bestKS() expects a list!")
  
  ## result matrix
  rr = matrix(0, nrow=length(x), ncol=length(x))
  ## compute upper diagonal
  for (i in 1:length(x))
  {
    for (j in (i+1):length(x))
    {
      if (j>length(x)) next;
      rr[i,j] = 1 - ks.test(x[[i]], x[[j]])$statistic
    }
  }
  ## fill diagonal with 1's
  diag(rr) = 1
  ## add up all values for each reference distribution (equivalent to just rowSums on full matrix)
  rr_sums = rowSums(rr, na.rm=T) + colSums(rr, na.rm=T)
  ## pick best
  r_max = which.max(rr_sums)
  ## get values
  v = pmax(rr[r_max,], rr[, r_max])
  
  return(data.frame(name = names(x), ks_best = v))
         
}




