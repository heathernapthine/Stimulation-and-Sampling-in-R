#' Heather Napthine, s2065896
#' Add your own function definitions on this file.

#' Log-Exponential density
#'
#' Compute the density or log-density for a Log-Exponential (LogExp)
#' distribution
#'
#' @param x vector of quantiles
#' @param rate vector of rates
#' @param log logical; if TRUE, the log-density is returned

dlogexp <- function(x, rate = 1, log = FALSE) {
  result <- log(rate) + x - rate * exp(x)
  if (!log) {
    exp(result)
  }
  result
}

#' Log-Sum-Exp
#'
#' Convenience function for computing log(sum(exp(x))) in a
#' numerically stable manner
#'
#' @param x numerical vector

log_sum_exp <- function(x) {
  max_x <- max(x, na.rm = TRUE)
  max_x + log(sum(exp(x - max_x)))
}
#' CI
#'
#' Compute the (1-alpha)*100% confidence interval for theta
#' given a sample, y.
#'
#' @param y vector of data
#' @param alpha specified significance level
#'

CI <- function(y, alpha) {
  n <- length(y)
  lambda_hat <- mean(y)
  theta_interval <-
    lambda_hat - sqrt(lambda_hat / n) * qnorm(c(1 - alpha / 2, alpha / 2))
  pmax(theta_interval, 0)
}

#' Estimate Coverage
#'
#' Function to perform interval coverage estimation 
#' 
#' @param CI  function object for confidence interval construction
#' @param alpha 1-alpha is the intended coverage probability
#' @param N the number of simulation replications to use for the coverage estimate
#' @param n the sample size
#' @param lambda the true lambda values for the Poisson model.
#'
estimate_coverage <- function(CI, alpha = 0.1, N = 10000, n, lambda){
  vector <- c()
  j=0
  for (k in (1:N)){
    y <- rpois(n, lambda)
    interval <- CI(y, alpha)
    
    if (lambda <= interval[2] && lambda >= interval[1]){
      j = j+1
    }
  }
  return(j/N)
}

#' Log-Prior density
#'
#' Compute the log-density for the prior distribution
#'
#' @param theta vector of parameters
#' @param params vector of positive parameters
#'
log_prior_density <- function(theta, params){
  density1 <- dnorm(theta[1], 0, params[1]**(1/2), log = TRUE)
  density2 <- dnorm(theta[2], 1, params[2]**(1/2), log = TRUE)
  density3 <- dlogexp(theta[3], rate = params[3], log = TRUE)
  density4 <- dlogexp(theta[4], rate = params[4], log = TRUE)
  logjointdensity <- density1 + density2 + density3 + density4
  return(logjointdensity)
}

#' Log-Likelihood
#'
#' Compute the log-likelihood for our data
#'
#' @param theta vector of parameters
#' @param x vector of data
#' @param y vector of data
#' 
loglik <- function(theta, x, y){
  sum((dnorm(y, theta[1]+ x*theta[2], (exp(theta[3]) + (exp(theta[4]))*(x**2))**(1/2), log = TRUE)))
}

#' Log-posterior Density
#'
#' Compute the log-posterior density for our data
#'
#' @param theta vector of parameters
#' @param x vector of data
#' @param y vector of data
#' @param params vector of positive parameters
#'
log_posterior_density <- function(theta, x, y, params){
  logpriordensity <- log_prior_density(theta, params)
  observationloglik <- loglik(theta, x, y)
  return(logpriordensity + observationloglik)
}

#' Posterior Mode
#'
#' Evaluate the posterior mode for our data
#'
#' @param theta_start vector of initial parameters
#' @param x vector of data
#' @param y vector of data
#' @param params vector of positive parameters
#'
posterior_mode <- function(theta_start, x, y, params){
  opt <- optim(par = theta_start, fn = log_posterior_density, x = x, y = y, params = params, control = list(fnscale = -1), hessian = TRUE)
  thetamodes <- opt$par
  hessian <- opt$hessian
  S <- -solve(hessian)
  returnlist <- list(thetamodes, hessian,  S)
  return(returnlist)
}

#' Gaussian
#'
#' Obtain a Gaussian approximation to the posterior distribution for theta. 
#'
#' @param theta_start vector of initial parameters
#' @param x vector of data
#' @param y vector of data
#' @param params vector of positive parameters
#' @param N the number of samples from which we make the approximation
#'
gaussian <- function(theta_start, x, y, params, N){
  posteriormodelist <- posterior_mode(theta_start, x, y, params)
  muvector <- posteriormodelist[[1]]
  sigmamatrix <- posteriormodelist[[3]]
  normalapproximation <- rmvnorm(N, muvector, sigmamatrix)
  return(normalapproximation)
}

#' Do Importance
#'
#' Use a multivariate Normal approximation as an importance
#' sampling distribution. 
#'
#' @param N the number of samples to generate
#' @param mu the mean vector for the importance distribution
#' @param S the covariance matrix
#' @param params vector of positive parameters
#'
do_importance <- function(theta_start, N, mu, S, x, y, params){
  
  beta1 <- c()
  beta2 <- c()
  beta3 <- c()
  beta4 <- c()
  log_weights <- c()
  
  xsample <- gaussian(theta_start, x, y, params, N)
  beta1 <- xsample[,1]
  beta2 <- xsample[,2]
  beta3 <- exp(xsample[,3])
  beta4 <- exp(xsample[,4])
  
  for (i in (1:N)){ 
    
    log_weight <- log_posterior_density(xsample[i,], x, y, params) -
      dmvnorm(xsample[i,], mean = mu, sigma = S, log = TRUE)
    log_weights <- append(log_weights, log_weight)
  }
  
  normalisingconstant <- log_sum_exp(log_weights)
  normalisedlogweights <- log_weights - normalisingconstant
  
  df <- data.frame(beta1 = beta1, beta2 = beta2, beta3 = beta3, beta4= beta4, NormalisedLogWeights = normalisedlogweights)
  
  return(df)
}




