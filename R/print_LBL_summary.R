#' Posterior Inference for LBL Functions
#'
#' \code{print_LBL_summary} prints list result from LBL_summary function in a more legible style.
#'
#' @param LBL_summary_object Summary object from \code{\link{LBL_summary}}.
#' @export
#'


print_LBL_summary <- function(LBL_summary_object){
  
  result1 <- data.frame(LBL_summary_object$haplotypes, LBL_summary_object$freq,
                        stringsAsFactors = F)
  result2 <- data.frame(LBL_summary_object$OR, LBL_summary_object$OR.CI,
                        LBL_summary_object$BF,
                        stringsAsFactors = F)
  result2$Significance <- rep("", nrow(result2))
  result2$Significance[result2[, 2] > 1 & result2[, 4] > 2] <- "*+"
  result2$Significance[result2[, 3] < 1 & result2[, 4] > 2] <- "*-"
  result <- merge(result1, result2, by = 0, all = T)
  result <- result[, -1]
  
  #Xiaofei: changes made
  result$Significance[is.na(result$Significance)]<-""
  colnames(result) <- c("Hap", "Freq", "OR", "OR Lower", "OR Upper",
                        "BF", "")
  
  print(result)
  cat("---\n")
  cat("Signif.codes: Risk '*+' Protective '*-' Not significant ' ' \n\n")
  
}

