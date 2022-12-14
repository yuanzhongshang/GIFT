#' @title The main function for conditional fine-mapping for in transcriptome-wide association studies with summary-level data
#' @description GIFT_individual applies a likelihood-based approach, accounting for the correlated cis-SNPs and genes in a region 
#' @param Zscore_1 the Zscore matrix of the cis-SNP effect size matrix, each column for one specific gene in eQTL data
#' @param Zscore_2 the Zscore vector of the cis-SNP effect size vector for one specific trait in GWAS data
#' @param Sigma1 the LD matrix in eQTL data
#' @param Sigma2 the LD matrix in GWAS data,both Sigma2sin and sigma1sin are often the same from the reference panel
#' @param R the estimated correlation matrix of gene expressions
#' @param n1 the sample size of eQTL data
#' @param n2 the sample size of GWAS data
#' @param pindex a vector with each element represents the number of cis-SNPs for each gene.
#' @param max_iterin The maximum iteration, which can be determined by users.
#' @param epsin The convergence tolerance of the absolute value of the difference  between the nth and (n+1)th log likelihood, which can be determined by users.
#' @param Cores The number of cores used in analysis. If the number of cores is greater than 1, analysis will perform with fast parallel computing. The function mclapply() depends on another R package "parallel" in Linux.
#' @return A list of estimated parameters including the p values for the gene-based test. 
#' \item{causal_effect}{The estimates of causal effect for each gene in a specific region}
#' \item{gene_based_test_pvalue}{The p values for each gene by the gene-based test}

GIFT_summary<-function(Zscore1, Zscore2, Sigma1, Sigma2, R, n1, n2, pindex, max_iterin =1000,epsin=1e-4, Cores=1){
  
  betax<-Zscore1/sqrt(n1-1)
  betay<-Zscore2/sqrt(n2-1)
  betax<-as.matrix(betax)
  betay<-as.vector(betay)
  k<-length(pindex)  
  Sigma1 <- as.matrix(Sigma1)
  Sigma2 <- as.matrix(Sigma2)
  constrFactor <- numeric(k)
  H1=GIFT_summarycpp(betax, betay, Sigma1, Sigma2, R, constrFactor, pindex, k, n1, n2, max_iterin, epsin)
  loglik0 <- H1$loglik[length(H1$loglik)]

  Cores=min(c(k,Cores))
  gene_specific_test_pvalue <- NULL 
  if(Cores == 1){
    for (i in 1:k) {
      constrFactor <- numeric(k)
      constrFactor[i] <- 1
      fit <- GIFT_summarycpp(betax, betay, Sigma1, Sigma2, R, constrFactor, pindex, k, n1, n2, max_iterin, epsin)
      gene_specific_test_pvalue[i] <- pchisq(2*(loglik0-fit$loglik[length(fit$loglik)]), 1, lower.tail=F)
    }
  }else{
    perform_func <- function(ti){
      constrFactor <- numeric(k)
      constrFactor[ti] <- 1
      fit <- GIFT_summarycpp(betax, betay, Sigma1, Sigma2, R, constrFactor, pindex, k, n1, n2, max_iterin, epsin)
      pvalue <- pchisq(2*(loglik0-fit$loglik[length(fit$loglik)]), 1, lower.tail=F)
      return(pvalue)
    }
    ti_list=c(1:k)
    gene_specific_test_pvalue <- unlist(mclapply(ti_list, perform_func, mc.cores = Cores))
  }
  
  result=list()
  result$causal_effect=H1$alpha
  result$gene_based_test_pvalue=gene_specific_test_pvalue
  #result$sigma_cisSNP=H1$sigmaZ
  #result$sigma_error_1=H1$sigmaX
  #result$sigma_error_2=H1$sigmaY
  return(result)
  
}
