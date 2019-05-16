#' LBL: Bayesian Lasso for detecting Rare (or Common) Haplotype Association
#'
#' LBL features a Bayesian LASSO model that detects association between
#' a phenotype and haplotypes given the (unphased) genotypes of individuals.
#'
#' LBL imposes a double exponential prior on the regression coefficients, effectively penalizing
#' unassociated hapltoypes while enhancing the stability and accuracy of the rare
#' haplotypes. LBL is capable of handling different study designs: this version of the software is capable
#' of handling independent cases and controls, case-parent trios and a mixture of both (provided
#' that the family data is independent from the case-control data). LBL uses algorithms from
#' \pkg{\link[hapassoc:pre.hapassoc]{hapassoc}} to obtain all possible haplotype pairs compatible with an
#' indivdiual's genotype.
#'
#' The posterior samples are obtained via Markov Chain Monte Carlo (MCMC)
#' algorithm and inference on the parameters of interest can be carried out
#' (Bayes Factor, Credible Interval, etc.) based on the posterior samples.
#'
#' @section Functions:
#' \code{\link{LBLcac}}  MCMC algorithm to obtain posterior samples for independent case-control data.
#' \code{\link{LBLfam}}  MCMC algorithm to obtain posterior samples for case-parent trio data.
#' \code{\link{LBLcombined}}  MCMC algorithm to obtain posterior samples for combined data.
#'
#' \code{\link{LBL_summary}} Provide model summary (in the form of list) based on posterior samples.
#' \code{\link{print_LBL_summary}} Print model summary in a nice fashion from the list result of \code{\link{LBL_summary}}.
#'
#' @docType package
#' @name LBL
NULL
