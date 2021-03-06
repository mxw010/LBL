#' An example file of case-parent trios
#'
#'A dataset containing the pedigree information, phenotypes and genotypes
#' of 250 case-parent trios in pedigree format.
#'
#' @usage data(fam)
#'
#' @docType data
#'
#' @format A data frame with 750 rows (250 trios) and
#'  16 ( = 6 fixed columns + 2 * # of SNPs) columns:
#' \describe{
#'   \item{column 1}{family ID, integers}
#'   \item{column 2}{individual ID, integers}
#'   \item{column 3}{father's ID, 0 for unknown father.}
#'   \item{column 4}{mother's ID, 0 for unknown mother.}
#'   \item{column 5}{individual's sex. 1 = male and 2 = female.}
#'   \item{column 6}{affection status. 0 = unknown, 1 = unaffected and 2 = affected.}
#'   \item{column 7}{}
#'   \item{column 8}{}
#'   \item{column 9}{}
#'   \item{column 10}{}
#'   \item{column 11}{}
#'   \item{column 12}{}
#'   \item{column 13}{}
#'   \item{column 14}{}
#'   \item{column 15}{}
#'   \item{column 16}{marker genotypes. Each marker is represented by two columns: one for each allele.}
#' }
"fam"


#' An example file consisting of independent cases/controls
#'
#' A dataset containing phenotypes and genotypes of 5 SNPs
#' for 500 independent cases and controls, and is presented in
#' a pedigree format for easier integration with family data.
#'
#' @usage data(cac)
#'
#' @docType data
#'
#' @format A data frame with 500 rows (250 cases and 250 controls)
#'  and 16 ( = 6 fixed columns + 2 * # of SNPs) columns:
#' \describe{
#'   \item{column 1}{family ID, integers}
#'   \item{column 2}{individual ID, integers}
#'   \item{column 3}{father's ID, 0 for unknown father.}
#'   \item{column 4}{mother's ID, 0 for unknown mother.}
#'   \item{column 5}{individual's sex: 1 = male and 2 = female.}
#'   \item{column 6}{affection status: 0 = unknown, 1 = unaffected, and 2 = affected.}
#'   \item{column 7}{}
#'   \item{column 8}{}
#'   \item{column 9}{}
#'   \item{column 10}{}
#'   \item{column 11}{}
#'   \item{column 12}{}
#'   \item{column 13}{}
#'   \item{column 14}{}
#'   \item{column 15}{}
#'   \item{column 16}{marker genotypes. Each SNP is represented by two columns: one for each allele.}
#' }
#' @docType data
"cac"



