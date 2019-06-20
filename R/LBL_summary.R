#' Posterior Inference for LBL
#'
#' \code{LBL_summary} Provides inferences based on the posterior samples. Specifically, this function will provide posterior means,
#' credible intervals, and Bayes Factors (BF) estimates for the haplotypic effect coefficients.
#'
#' @param output An object returned by LBL with `summary=FALSE`.
#'
#' @param a  First hyperparameter of the prior for regression coefficients,
#'    \eqn{\beta}. The prior variance of \eqn{\beta} is 2/\eqn{\lambda^2} and \eqn{\lambda} has Gamma(a,b)
#'    prior. The Gamma prior parameters a and b are such that the mean and variance of
#'    the Gamma distribution are \eqn{a/b} and \eqn{a/b^2}.
#'    The value will be transferred from LBL functions.
#'
#' @param b  Second hyperparameter of the Gamma(a,b) distribution described above.
#'  The value will be transferred from LBL functions.
#'
#' @param e A (small) number \eqn{\epsilon} in the null hypothesis of no association,
#'    \eqn{H_0: |\beta| \le \epsilon}. The default is 0.1. Changing e from the default of 0.1 may necessitate choosing a
#'    different threshold for Bayes Factor (one of the outputs) to infer
#'    association.
#'
#' @param ci.level Credible probability.. The probability that the true value of \eqn{beta} will
#'    be within the credible interval. Default is 0.95 which corresponds to a 95\% posterior credible interval. Only used if `summary = TRUE`.
#'
#' @return A list with the following components:
#' \describe{
#' \item{haplotypes}{The list of haplotypes used in the analysis.}
#'
#' \item{OR}{Posterior mean of odds ratio.}
#'
#' \item{OR.CI}{95\% posterior credible sets for the ORs.}
#'
#' \item{BF}{Bayes Factor estimates based on posterior samples.
#'  if the posterior samples are all greater than e, then BF is set to be 999.}
#'
#'}
#'
#' @export
#'

####Function to provide a summary of analysis based on posterior samples.
##including: posterior mean, CI, BF
LBL_summary <- function(output, a, b, e=0.1, ci.level=0.95){
  #OR estimates and CIs
  beta.out <- output$beta
  lambda.out <- output$lambda
  freq.out <- output$freq
  x.length <- ncol(beta.out)
  
  
  beta.ci <- exp(apply(beta.out, 2, stats::quantile, c((1-ci.level)/2, 1-(1-ci.level)/2)))
  beta.est <- exp(colMeans(beta.out))
  
  freq.est <- colMeans(freq.out)
  
  nreps <- nrow(beta.out)
  
  #prior probability of abs(beta_j) > e
  prior.prob<-(b/(e+b))^a
  #prob.alt<-numeric(x.length)
  #BF<-rep(999, x.length)
  # for (i in 1:x.length)
  # {
  #   prob.alt[i]<-length(beta.out[,i][abs(beta.out[,i]) > e])/(num.it-burn.in)
  #   if (prob.alt[i] < 1) BF[i]<-(prob.alt[i]/(1-prob.alt[i]))/(prior.prob/(1-prior.prob))
  # }
  #this is SLIGHTLY faster than the older implementation
  BF <- apply(beta.out,2, function(x)  {
    temp <- sum(abs(x) > e)
    bf <- ifelse(temp < nreps, temp/(nreps-temp)/(prior.prob/(1-prior.prob)),999)
    return(bf)
  }
  )
  #XF: added frequency
  result<-list(haplotypes=output$haplo.names, freq=freq.est,
               OR=beta.est, OR.CI=t(beta.ci), BF=BF)
  return(result)
}
