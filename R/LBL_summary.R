#' Posterior Inference for LBL
#'
#' \code{LBL_summary} Provides inferences based on the posterior samples. Speficially, this function will provide posterior means,
#' confidence intervals and Bayes Factors (BF) estimates for the haplotypic effect coefficients.
#'
#'
#' @param output an obect returned by LBL with `summary=FALSE`.
#'
#' @param a  first hyperparameter of the prior for regression coefficients,
#'    \eqn{\beta}. The prior variance of \eqn{\beta} is 2/\eqn{\lambda^2} and \eqn{\lambda} has Gamma(a,b)
#'    prior. The Gamma prior parameters a and b are such that the mean and variance of
#'    the Gamma distribution are \eqn{a/b} and \eqn{a/b^2}. The default value of a is 15.
#'
#' @param b  b parameter of the Gamma(a,b) distribution described above; default
#'    is 15.
#'
#' @param e a (small) number \eqn{\epsilon} in the null hypothesis of no association,
#'    \eqn{H_0: |\beta| \le \epsilon}. The default is 0.1. Changing e from default of 0.1 may need choosing a
#'    different threshold for Bayes Factor (one of the outputs) to infer
#'    association.
#'
#' @param ci.level Confidence level. The probability that the true value of \eqn{beta} will
#'    be within the confidence interval. Default is 0.95 which corresponds to a 95\% posterior confidence interval. Only used if `summary = TRUE`.
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
#' @references
#' Biswas S, Lin S (2012). Logistic Bayesian LASSO for Identifying
#'   Association with Rare Haplotypes and Application to Age-related Macular
#'   Degeneration. Biometrics, 68(2): 587-97.
#' @export
#'

####Function to provide a summary of analysis based on posterior samples.
##including: posterior mean, CI, BF
LBL_summary <- function(output, a=15,b=15,e=0.1, ci.level=0.95){
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

  #Changed by Xiaofei
  result1<-data.frame(output$haplo.names,freq.est,
                      stringsAsFactors = F)
  result2<-data.frame(beta.est,t(beta.ci),BF,
                      stringsAsFactors = F)
  result<-merge(result1,result2,by=0,all=T)
  result<-result[,-1]
  colnames(result)<-c("Hap","Freq","OR","OR Lower", "OR Upper", "BF")

  #BF[BF==999] <-">100"
  #result<-list(haplotypes=output$haplo.names, OR=beta.est, OR.CI=t(beta.ci), BF=BF)
  return(result)
}

