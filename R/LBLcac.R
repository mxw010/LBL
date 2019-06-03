#' Logistic Bayesian Lasso for detecting Rare (and Common) Haplotyptic Association with
#' Population Based Designs
#'
#' \code{LBLcac} is a Bayesian LASSO method developed to detect association between
#' common/rare haplotypes and dichotomous disease phenotype, based on MCMC algorithm.
#' The "cac" in the title stands for "case and control", as this function will handle
#' independent cases and controls study design. For other types of study designs, see \code{\link{LBLfam}} and
#' \code{\link{LBLcombined}}. This function takes standard pedigree format as input with an individual's
#' genotypes, phenotype and familiar relationships, uses
#' an algorithm developed by Burkett \emph{et al}. to impute possible haplotype configurations
#' (for more detail, see \href{https://www.rdocumentation.org/packages/hapassoc/versions/0.5-1/topics/pre.hapassoc}{pre.hapassoc} from \pkg{hapassoc} package).
#' The input does not allow missing observations, and therefore subjects with
#' missing data are removed. This function returns an object containing summary statistics
#' for all haplotypes that are consistent with the input phenotype data.
#'
#' @param data.cac Input data. data.cac should be either a data frame or a matrix,
#'  consisting of "n" rows and 6+2*p
#'   columns, where n is the number of cases and controls, and p
#'   is the number of SNPs.  The data should be in standard pedigree format, with
#'   the first 6 columns representing the family ID, individual ID, father ID,
#'   mother ID, sex, and affection status. The other 2*p columns are genotype
#'   data in allelic format, with each allele of a SNP taking up one column.
#'   An example can be found in this package under the name \code{"cac"}. For
#'   more information about the format, type \code{"?cac"} into R, or see "Linkage Format"
#'   section of \url{https://www.broadinstitute.org/haploview/input-file-formats}.
#'   Note that since these are independent case-control data, the father ID and mother ID
#'    are missing (coded as 0) and each individual has an unique family ID.
#'
#' @param baseline  Haplotype to be used for baseline coding; default is the most
#'   frequent haplotype according to the initial haplotype frequency estimates.
#'   This argument should be a character, starting with an h and followed by the
#'   SNPs at each marker locus, for example, if the desired baseline haplotype is 0 1 1 0 0,
#'   then baseline should be coded as "h01100".
#'
#' @param a  First hyperparameter of the prior for regression coefficients,
#'    \strong{\eqn{\beta}}. The prior variance of \eqn{\beta} is 2/\eqn{\lambda^2} and \eqn{\lambda} has Gamma(a,b)
#'    prior. The Gamma prior parameters a and b are such that the mean and variance of
#'    the Gamma distribution are \eqn{a/b} and \eqn{a/b^2}. The default value of a is 15.
#'
#' @param b  Second hyperparameter of the Gamma(a,b) distribution described above; default
#'    is 15.
#'
#' @param start.beta Starting value of all regression coefficients, \strong{\eqn{\beta}};
#'    default is 0.01.
#'
#' @param lambda Starting value of the \eqn{\lambda} parameter described above;
#'    default is 1.
#'
#' @param D Starting value of the D parameter, which is the within-population
#'    inbreeding coefficient; default is 0.
#'
#' @param seed Seed to be used for the MCMC in Bayesian Lasso; default is a
#'    random seed. If exact same results need to be reproduced, seed should be
#'    fixed to the same number.
#'
#' @param burn.in Burn-in period of the MCMC sampling scheme; default is 10000.
#'
#' @param num.it Total number of MCMC iterations including burn-in; default is
#'    40000.
#'
#' @param summary Logical; if \code{TRUE}, LBLcac will return a summary of the analysis in the form of a list. If \code{FALSE}, LBLcac will return the posterior samples for all parameters.
#'  Default is set to be TRUE.
#'
#' @param e A (small) number \eqn{\epsilon} in the null hypothesis of no association,
#'    \eqn{H_0: |\beta| \le \epsilon}. The default is 0.1. Changing e from the default of 0.1 may necessitate choosing a
#'    different threshold for Bayes Factor (one of the outputs) to infer
#'    association. Only used if \code{summary = TRUE}.
#'
#' @param ci.level Credible probability. The probability that the true value of \eqn{beta} will
#'    be within the credible interval. Default is 0.95, which corresponds to a 95\% posterior credible interval. Only used if \code{summary = TRUE}.
#'
#'@return If \code{summary = FALSE}, return a list with the following components:
#' \describe{
#' \item{haplotypes}{The list of haplotypes used in the analysis. The last column is the reference haplotype.}
#'
#' \item{beta}{Posterior samples of betas stored in a matrix.}
#'
#' \item{lambda}{a vector of (num.it-burn.in) posterior samples of lambda.}
#'
#' \item{freq}{Posterior samples of the frequencies of haplotypes stored in a matrix format, in the same order as haplotypes.}
#'
#' \item{init.freq}{The haplotype distribution used to initiate the MCMC.}
#'
#'}
#'
#' If \code{summary = TRUE}, return the result of LBL_summary.
#'  For details, see the description of the \code{?LBL_summary} function.
#'
#' @references
#' Biswas S, Lin S (2012). Logistic Bayesian LASSO for Identifying
#'   Association with Rare Haplotypes and Application to Age-related Macular
#'   Degeneration. Biometrics, 68(2): 587-97.
#'
#'
#' @seealso
#' \code{\link{LBLfam}}, \code{\link{LBLcombined}}, \code{\link{LBL_summary}}.
#'
#' @examples
#'  data(cac)
#'  cac.obj<-LBLcac(cac)
#'  cac.obj
#'  print_LBL_summary(cac.obj)
#'
#'
#' @export
#'
#' @useDynLib LBL LBLmcmc
#'
LBLcac <- function(data.cac, baseline="missing", a = 15, b = 15, start.beta = 0.01, lambda = 1, D = 0, seed = NULL, burn.in = 10000,
                   num.it = 40000,summary=T, e = 0.1, ci.level=0.95)
{
  #still problem with I/O between R and C, disable monitor for now.
  #monitor: if monitor == F, do not monitor,
  #                          otherwise, monitor progress every other monitor = n iterations

  #does not check for independency of case control data
  #user should run programs like pedstats ahead of time to make sure there is no inconsistency in the pedigree files
  #also make it so that only one affected offspring per family

  #baseline="missing"; a = 15; b = 15; start.beta = 0.01; lambda = 1; D = 0; seed = NULL; e = 0.1; burn.in = 10000; num.it = 50000; verbose=F

  #already taken care of this in dependency
  # if (!requireNamespace("hapassoc", quietly = TRUE)) {
  #   stop("Package \"hapassoc\" needed for this function to work. Please install it.",
  #     call. = FALSE)
  # }

  #if fam is missing, type can't be fam or combined

  if (missing(data.cac) ) {
    stop(paste("Must provide case-control data!\n\n"))
  }

  if (!is.null(seed)) set.seed(seed)

  data.cac.re <- data.cac

  #reformat case and control
  #case=1, everything else =0
  if (!is.matrix(data.cac)) data.cac <- matrix(unlist(data.cac,use.names=F),nrow=nrow(data.cac))
  status <- ifelse(data.cac[,6] ==2, 1, 0)
  data.cac[,6] <- status

  data.new <- data.cac
  p <- (ncol(data.new)-6)/2
  haplos.list <- pre.hapassoc(data.new, p, pooling.tol = 0, allelic = T, verbose=F)
  haplos.names <- names(haplos.list$initFreq)

  #if (missing(input.freq)) {
  #  freq <- haplos.list$initFreq
  #} else {
  #  freq <- input.freq
  #}
  freq <- haplos.list$initFreq

  #set baseline haplotype
  if (!(baseline %in% colnames(haplos.list$haploDM))& (baseline != "missing")) {
    warning("Baseline haplotype does not exist!\nSetting baseline haplotype as missing...")
    baseline<-"missing"
  }

  if (baseline=="missing") {
    baseline <- haplos.names[which.max(freq)]
  }

  column.subset <- colnames(haplos.list$haploDM) != baseline
  freq.new<-freq[column.subset]
  freq.new[length(freq.new)+1]<-freq[!column.subset]
  freq.new<-as.vector(freq.new)     #numerical frequencies with the baseline freq in the end

  #hdat <- cbind(haplos.list$nonHaploDM[,1], haplos.list$haploDM[,column.subset])
  ID <- haplos.list$ID
  N <- sum(haplos.list$wt)
  y <- as.numeric(haplos.list$nonHaploDM[, 6])
  x <- data.matrix(haplos.list$haploDM[,column.subset])
  colnames(x)<-NULL
  num.haplo.id<-as.vector(table(ID))
  x.length<-as.integer(dim(x)[2])
  unique.x<-unique(x)

  #===========end here=================#

  # XF: these are part of input for MCMC
  beta=rep(start.beta, x.length)
  beta.out<-numeric((num.it-burn.in)*x.length)
  lambda.out<-numeric(num.it-burn.in)
  freq.out<-numeric((num.it-burn.in)*(x.length+1))     #XF: include baseline
  D.out<-numeric(num.it-burn.in)

  out<-.C("LBLmcmc", x=as.integer(x), n=as.integer(dim(x)[1]), as.integer(y), as.integer(N), as.integer(num.haplo.id), as.integer(x.length),
          as.double(freq.new), as.double(D), as.double(beta), as.double(a), as.double(b), as.integer(t(unique.x)), as.integer(dim(unique.x)[1]),
          as.double(lambda), as.integer(num.it), as.integer(burn.in), beta.out=as.double(beta.out), lambda.out=as.double(lambda.out),
          freq.out=as.double(freq.out), D.out=as.double(D.out), monitor=as.integer(FALSE))


  beta.out<-matrix(out$beta.out,nrow=num.it-burn.in, byrow=TRUE)
  #beta.out <- rbind(as.vector(beta),beta.out,deparse.level=0)
  lambda.out<-matrix(out$lambda.out,nrow=num.it-burn.in, byrow=TRUE)
  #lambda.out <- c(lambda,lambda.out)
  freq.out<-matrix(out$freq.out,nrow=num.it-burn.in, byrow=TRUE)
  #freq.out <- rbind(as.vector(freq.new), freq.out)
  haplo.names <- names(haplos.list$initFreq)[column.subset]
  haplo.names <- c(haplo.names, names(haplos.list$initFreq)[!column.subset])   #XF: uncommented this line

  raw.output <- list(haplo.names=haplo.names, beta=beta.out, lambda=lambda.out, freq=freq.out, init.freq=freq.new)

  if (summary==FALSE) {
    output <- raw.output
  } else {
    #calculating the Bayes Factors
    output <- LBL_summary(raw.output, a=a,b=b,e=e, ci.level=ci.level)
  }

  return(output)
}
