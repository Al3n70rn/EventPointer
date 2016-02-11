#' Wilcoxon test for enrichment analysis
#' 
#' Wilcoxon test for enrichment analysis
#' @param ExprT Expression matrix Used
#' @param t1 t-statistic associated to Pv1
#' @param GeneGO Matrix relating Genes/Transcripts/Events to GO/Domains
#' @param alternative Type of test to be used
#' @return z Resulting Z values for enrichment analysis

Wilcoxon.z <- function(ExprT, GeneGO) {
  
  Transp_ENSTxGO <- Matrix::t(GeneGO)
  RExprT <- rank(ExprT)
  
  Prod <- as.vector(Transp_ENSTxGO%*%RExprT)
  
  # Trick to deal with the three classes analysis of the domains
  Transp_ENSTxGO@x[Transp_ENSTxGO@x<0.5] <- 1
  
  # Calculate the number of elements "0"(ny) and "1"(nx)
  nx <-  as.vector(rep(1,nrow(GeneGO)) %*% GeneGO)
  ny <- nrow(GeneGO) -nx
  
  # Calculate the estimated variance
  Var <- (nx*ny*(nx+ny+1)/12)
  
  # Calculate the estimated mean
  media <- nx*(nx + ny + 1)/2
  
  # Calculate the standard desviation
  std <- sqrt(Var)
  
  z <- (Prod-media)/std
  
  return(z)
  
}
