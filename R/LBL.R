#' LBL: Logistic Bayesian Lasso for detecting Rare (or Common) Haplotype Association
#'
#' LBL uses the Bayesian LASSO framework to detect association between
#' a phenotype and haplotypes given the (unphased) genotypes of individuals.
#'
#' LBL uses Bayesian Lasso to detect rare haplotypes that are associated with common diseases. 
#' The current implementation considers dichotomous traits.
#' A future release will include quantitative and survival traits.
#' LBL is capable of handling different study designs: this version of the software is capable of handling independent cases and controls, case-parent trios and a mixture of both (provided
#' that the family data is independent of the case-control data). LBL uses algorithms from
#' \pkg{\link[hapassoc:pre.hapassoc]{hapassoc}} to obtain all possible haplotype pairs compatible with an
#' indivdiual's set of genotypes.
#'
#' The posterior samples are obtained via Markov Chain Monte Carlo (MCMC)
#' algorithm and inference on the parameters of interest can be carried out
#' (Bayes Factor, Credible Interval, etc.) based on the posterior samples.
#'
#' @section Functions:
#' \code{\link{LBL}}:  MCMC algorithm to obtain posterior samples for independent case-control data.
#' \code{\link{famLBL}}:  MCMC algorithm to obtain posterior samples for case-parent trio data.
#' \code{\link{cLBL}}:  MCMC algorithm to obtain posterior samples for combined data.
#'
#' \code{\link{LBL_summary}} provides model summary (in the form of list) based on posterior samples.
#' \code{\link{print_LBL_summary}} prints model summary in a user-friendly format from the list result of \code{\link{LBL_summary}}.
#'
#' @docType package
#' @name LBL
#' 
#' 
#' @references
#'
#' Biswas S, Lin S (2012). Logistic Bayesian LASSO for Identifying
#'   Association with Rare Haplotypes and Application to Age-related Macular
#'   Degeneration. Biometrics, 68(2): 587-97.
#'
#' Wang, M., & Lin, S. (2014). FamLBL: detecting rare haplotype disease association
#' based on common SNPs using case-parent triads. Bioinformatics, 30(18), 2611-2618.
#'
#' Burkett, K., Graham, k. and McNeney, B. (2006). "hapassoc: Software for likelihood
#' inference of trait associations with SNP haplotypes and other attributes."
#' Journal of Statistical Software 16(2): 1-19.
NULL