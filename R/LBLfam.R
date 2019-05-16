#' Bayesian Lasso for detecting Rare (or Common) Haplotype Association in Case-Parent triad Designs
#'
#' \code{LBLfam} is an MCMC algorithm that generates a posterior samples for family trio data. This function
#' takes standard pedigree format as input, and uses an algorithm based on
#' function from \pkg{hapassoc} package to impute the (phased) haplotypes of individuals.
#' The input does not allow missing observations and subjects with
#' missing data are removed. The function returns an object containing posterior
#' samples based on Markov Chain Monte Carlo algorithm after the burn-in period.
#'
#' @param data.fam The input data. It should consist of "n" rows and 6+2*p columns,
#'   where n is the number of individuals of case-parent trios, and p is the
#'   number of SNPs.  The data should be in standard pedigree format, with the
#'   first 6 columns representing the family ID, individual ID, father ID,
#'   mother ID, sex, and affection status. The other 2*p columns are genotype
#'   data in allelic format, with each allele of a SNP taking up one column.
#'   An example can be found in this package under the name \code{"fam"}.For
#'   more information about the format, type \code{"?fam"} into R see "Linkage Format"
#'   section of \url{https://www.broadinstitute.org/haploview/input-file-formats}.
#'
# @param input.freq  optional. Specify frequency distribution of haplotypes. If
#'   not provided, the algorithm will use the frequency estimated by hapassoc.
#'
#' @param baseline  haplotype to be used for baseline coding; default is the most
#'   frequent haplotype according to the initial haplotype frequency estimates
#'   returned by \link{pre.hapassoc}. This argument should be a character, starting
#'   with an h and followed by the baseline haplotype.
#'
#' @param a  first hyperparameter of the prior for regression coefficients,
#'    \eqn{\beta}. The prior variance of \eqn{\beta} is 2/\eqn{\lambda^2} and \eqn{\lambda} has Gamma(a,b)
#'    prior. The Gamma prior parameters a and b are such that the mean and variance of
#'    the Gamma distribution are \eqn{a/b} and \eqn{a/b^2}. The default value of a is 15.
#'
#' @param b  b parameter of the Gamma(a,b) distribution described above; default
#'    is 15.
#'
#' @param start.beta starting value of all regression coefficients, \eqn{\beta};
#'    default is 0.01.
#'
#' @param lambda starting value of the \eqn{\lambda} parameter described above;
#'    default is 1.
#'
#' @param D starting value of the D parameter, which is the within-population
#'    inbreeding coefficient; default is 0.
#'
#' @param seed the seed to be used for the MCMC in Bayesian Lasso; default is a
#'    random seed. If exactly same results need to be reproduced, seed should be
#'    fixed to the same number.
#'
#' @param burn.in burn-in period of the MCMC sampling scheme; default is 10000.
#'
#' @param num.it total number of MCMC iterations including burn-in; default is
#'    40000.
#'
# @param verbose should the output from \code{\link[hapassoc]{pre.hapassoc}} be printed. Default is
#'    \code{FALSE}.
#'
# @param monitor if true, will monitor the progress of the Markov Chain by
#'    reporting progress every 5,000 iterations.
#' @param summary logical; if TRUE, LBLfam will return a summary of the analysis. If FALSE, LBLfam will return the posterior samples of MCMC.
#'
#' @param e a (small) number \eqn{\epsilon} in the null hypothesis of no association,
#'    \eqn{H_0: |\beta| \le \epsilon}. The default is 0.1. Changing e from default of 0.1 may need choosing a
#'    different threshold for Bayes Factor (one of the outputs) to infer
#'    association. Only used if \code{summary = TRUE}.
#'
#' @param ci.level Credible probability. The probability that the true value of \eqn{beta} will
#'    be within the credible interval. Default is 0.95 which corresponds to a 95\% posterior credible interval. Only used if \code{summary = TRUE}.
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
#' If \code{summary = TRUE}, return the result of LBL_summary. For details, see \code{?LBL_summary}.
#'
#' @references
#'
#' Wang, M., & Lin, S. (2014). FamLBL: detecting rare haplotype disease association
#' based on common SNPs using case-parent triads. Bioinformatics, 30(18), 2611-2618.
#'
#' Burkett, K., Graham, k. and McNeney, B. (2006). "hapassoc: Software for likelihood
#' inference of trait associations with SNP haplotypes and other attributes."
#' Journal of Statistical Software 16(2): 1-19.
#'
#' @author {Meng Wang, Swati Biswas, Shili Lin}
#' @seealso
#' \code{\link{LBL_summary}}.
#'
#' @examples
#'  data(fam)
#'  LBLfam(fam)
#'   fam.obj<-LBLfam(fam)
#'   fam.obj
#'   print_LBL_summary(fam.obj)
#'
#' @export
#'
#' @useDynLib LBL famLBLmcmc
#'
#'
LBLfam <- function(data.fam, baseline="missing", start.beta = 0.01, lambda = 1, D = 0, seed = NULL, a = 15, b = 15, burn.in = 10000,
                   num.it = 40000, summary=TRUE, e = 0.1, ci.level=0.95)
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
  if (missing(data.fam)) {
    stop(paste("Must provide family data!\n\n"))
  }

  if (!is.null(seed)) set.seed(seed)


  if (!is.matrix(data.fam)) data.fam <- matrix(unlist(data.fam,use.names=F),nrow=nrow(data.fam))
  n.fam <- length(unique(data.fam[,1]))
  #only works for trios!
  aff <- mo <- fa <- matrix(NA,ncol=ncol(data.fam),nrow=n.fam)
  t <- 0
  for (fam in unique(data.fam[,1])) {
    temp <- data.fam[data.fam[,1]==fam,]
    #offspring is the individual whose parents are in the pedigee
    #get affected offspring
    #if no affected individuals, skip
    affect.offspring <- which(temp[,6]==2 & rowSums(temp[,3:4]==0)==0 )
    if (length(affect.offspring) == 0) next
    t <- t + 1
    #write father, mother and affected offspring into seperate files
    aff.fam <- aff[t,] <- unlist(temp[affect.offspring,])
    fa[t,] <- unlist(temp[match(aff.fam[3],temp[,2]),])
    mo[t,] <- unlist(temp[match(aff.fam[4],temp[,2]),])
  }
  data.fam <- data.frame(rbind(fa[1:t,], mo[1:t,],aff[1:t,]))
  #names(data.fam.re) <- names(data.cac.re)


  data.new <- data.fam


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


  hap.ID <- haplos.list$ID
  ID <- hap.ID[hap.ID<=(n.fam*3)]
  N <- n.fam*2
  n <- n.fam
  x <- data.matrix(haplos.list$haploDM)[hap.ID<=(n.fam*3),]

  fam.wt<-haplos.list$wt[hap.ID<=(n.fam*3)]
  nonunique <- unique(ID[fam.wt!= 1])

  nonuni.fa <- nonunique[nonunique <= n]
  nonuni.mo <- nonunique[nonunique > n & nonunique <= 2*n]
  nonuni.kid <- nonunique[nonunique > 2*n]
  nonuni.family <- sort(unique(c(nonuni.fa, (nonuni.mo - n), (nonuni.kid - 2*n))))

  xc <- xu <- NULL
  count <- 0
  num.haplo.id <- rep(1, n)
  for ( i in 1:n) {
    if (i %in% nonuni.family){
      for (j in grep(T, ID==i)) {
        for (k in grep(T, ID==(i+n))) {
          for (m in grep(T, ID == (i + 2*n))) {
            child.hap <- grep(FALSE, x[m,] == 0)
            if (length(child.hap) ==1) child.hap = rep(child.hap, 2)
            judge.father <- sapply(child.hap, '%in%', grep(FALSE, x[j,] == 0))
            if (sum(judge.father) == 2) {
              judge <- sum(sapply(child.hap, '%in%', grep(FALSE, x[k,] == 0))) > 0
            } else if (sum(judge.father) == 1) {
              judge <- child.hap[grep(FALSE, judge.father)] %in% grep(FALSE, x[k,] == 0)
            } else if (sum(judge.father) == 0) {
              judge <- 0
            } else stop("error in configuring family haplotype!")

            if (judge > 0) {
              tmp <- x[j,] + x[k,] - x[m,]
              xc <- rbind(xc, x[m,])
              xu <- rbind(xu, tmp)
              count <- count + 1
            }
          }
        }
      }
      num.haplo.id[i] <- count
      count <- 0
    } else {
      xc <- rbind(xc, x[ID == (i + 2*n),])
      xu <- rbind(xu, (x[ID == i,] + x[ID == (i + n),] - x[ID == (i + 2*n),]))
    }
  }

  if (any(num.haplo.id==0)) {
    cat("Incompatile trios have been detected, and will be removed from analysis!\n")
    num.haplo.id <- num.haplo.id[num.haplo.id !=0]
    N <- length(num.haplo.id)*2
    cat("A total of", length(num.haplo.id), "families remain\n")
  } else
    cat("A total of", length(num.haplo.id), "families are in the study\n")

  #XF: changed position
  num.haplo.id.fam <- rep(num.haplo.id,2)
  #doesn't have this in famLBL


  if (any(xc) < 0 | any(xu) <0) stop("Error with seperating haplotypes!")
  #error check here


  x <- rbind(xc, xu)
  rownames(x) <- NULL
  x <- data.matrix(x[,column.subset])
  y <- rep(1:0, c(nrow(xc), nrow(xu)))

  x.length<-as.integer(dim(x)[2])
  unique.x<-unique(x)


  #===========end here=================#

  # XF: these are part of input for MCMC
  beta=rep(start.beta, x.length)
  beta.out<-numeric((num.it-burn.in)*x.length)
  lambda.out<-numeric(num.it-burn.in)
  freq.out<-numeric((num.it-burn.in)*(x.length+1))     #XF: include baseline
  D.out<-numeric(num.it-burn.in)


  out<-.C("famLBLmcmc", x=as.integer(x), n=as.integer(dim(x)[1]), as.integer(y), as.integer(N), as.integer(num.haplo.id.fam),
          as.integer(x.length), as.double(freq.new), as.double(D), as.double(beta), as.double(a), as.double(b), as.integer(t(unique.x)),
          as.integer(dim(unique.x)[1]), as.double(lambda), as.integer(num.it), as.integer(burn.in), beta.out=as.double(beta.out),
          lambda.out=as.double(lambda.out), freq.out=as.double(freq.out), D.out=as.double(D.out), monitor=as.integer(FALSE))


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
