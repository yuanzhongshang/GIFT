#' @title The main function for conditional fine-mapping for in transcriptome-wide association studies with summary-level data
#' @description GIFT_summary applies a likelihood-based approach, accounting for the correlated cis-SNPs and genes in a region 
#' @param Zscore_1 Zscore matrix of the cis-SNP effect size matrix, each column for one specific gene from eQTL data.
#' @param Zscore_2 Zscore vector of the cis-SNP effect size vector from GWAS data.
#' @param Sigma1 LD matrix from eQTL data.
#' @param Sigma2 LD matrix from GWAS data, both Sigma1 and Sigma2 are often the same from the reference panel.
#' @param R Estimated correlation matrix of gene expressions.
#' @param n1 Sample size of eQTL data.
#' @param n2 Sample size of GWAS data.
#' @param gene The gene names vector.
#' @param pindex A vector with each element represents the number of cis-SNPs for each gene.
#' @param maxiter The maximum iteration, which can be determined by users.
#' @param tol The convergence tolerance of the absolute value of the difference  between the nth and (n+1)th log likelihood, which can be determined by users.
#' @param ncores The number of cores used in analysis. If the number of cores is greater than 1, analysis will perform with fast parallel computing. The function mclapply() depends on another R package "parallel" in Linux.
#' @param in_sample_LD A logical value, which represents whether in-sample LD was used. If in-sample LD was not used, the LD matrix is regularized to be (1-s1)*Sigma1+s1*E and (1-s2)*Sigma2+s2*E, where s1 and s2 are estimated by estimate_s_rss in susieR, and E is an identity matrix. A grid search is performed over the range from 0.1 to 1 if the estimation does not work. The function estimate_s_rss() depends on another R package "susieR".
#' @return A data frame including the causal effect estimates and p values for the gene-based test. 

GIFT_summary<-function(Zscore1, Zscore2, Sigma1, Sigma2, R = NULL, n1, n2, gene, pindex, maxiter = 1000, tol = 1e-4, ncores = 1, in_sample_LD = F){
  
  betax<-Zscore1/sqrt(n1-1)
  betay<-Zscore2/sqrt(n2-1)
  betax<-as.matrix(betax)
  betay<-as.vector(betay)
  k<-length(pindex)  
  Sigma1 <- as.matrix(Sigma1)
  Sigma2 <- as.matrix(Sigma2)
  constrFactor <- numeric(k)
  if(is.null(R)){
	  R=diag(k)
  }
  if(in_sample_LD==T){
	result <- try({
	  H1=GIFT_summarycpp(betax, betay, Sigma1, Sigma2, R, constrFactor, pindex, k, n1, n2, maxiter, tol)
	  loglik0 <- H1$loglik[length(H1$loglik)]

	  Cores=min(c(k,ncores))
	  gene_specific_test_pvalue <- NULL 
	  
	  if(Cores == 1){
		for (i in 1:k) {
		  constrFactor <- numeric(k)
		  constrFactor[i] <- 1
		  fit <- GIFT_summarycpp(betax, betay, Sigma1, Sigma2, R, constrFactor, pindex, k, n1, n2, maxiter, tol)
		  gene_specific_test_pvalue[i] <- pchisq(2*(loglik0-fit$loglik[length(fit$loglik)]), 1, lower.tail=F)
		}
	  }else{
		library(parallel)
		perform_func <- function(ti){
		  constrFactor <- numeric(k)
		  constrFactor[ti] <- 1
		  fit <- GIFT_summarycpp(betax, betay, Sigma1, Sigma2, R, constrFactor, pindex, k, n1, n2, maxiter, tol)
		  pvalue <- pchisq(2*(loglik0-fit$loglik[length(fit$loglik)]), 1, lower.tail=F)
		  return(pvalue)
		}
		ti_list=c(1:k)
		gene_specific_test_pvalue <- unlist(mclapply(ti_list, perform_func, mc.cores = Cores))
	  }
  
	  result <- data.frame(gene = gene, causal_effect = H1$alpha, p = gene_specific_test_pvalue)
	  return(result)
	}, silent = TRUE)
  
    if(sum(grep("decomposition failed",result[1]))==1){
      cat("LD reference is not representative of the population in the analysis\n")
      cat("Regularize the LD matrix as (1-s1)*Sigma1+s1*E and (1-s2)*Sigma2+s2*E \n")
      result <- try({
        library(susieR)
        s2 <- estimate_s_rss(Zscore2, Sigma2, n2, r_tol = 1e-08, method = "null-mle")
        E <- diag(dim(Zscore1)[1])
        s1 <- NULL
        for(i in 1:length(gene)){
          s1[i] <- estimate_s_rss(Zscore1[,i], Sigma1, n1, r_tol = 1e-08, method = "null-mle")
        }
         s1 <- max(s1)
         Sigma1 <- (1-s1)*Sigma1+s1*diag(sum(pindex))
         Sigma2 <- (1-s2)*Sigma2+s2*diag(sum(pindex))
         H1=GIFT_summarycpp(betax, betay, Sigma1, Sigma2, R, constrFactor, pindex, k, n1, n2, maxiter, tol)
         loglik0 <- H1$loglik[length(H1$loglik)]
        
         Cores=min(c(k,ncores))
         gene_specific_test_pvalue <- NULL 
         if(Cores == 1){
           for (i in 1:k) {
             constrFactor <- numeric(k)
             constrFactor[i] <- 1
             fit <- GIFT_summarycpp(betax, betay, Sigma1, Sigma2, R, constrFactor, pindex, k, n1, n2, maxiter, tol)
             gene_specific_test_pvalue[i] <- pchisq(2*(loglik0-fit$loglik[length(fit$loglik)]), 1, lower.tail=F)
           }
         }else{
           perform_func <- function(ti){
             constrFactor <- numeric(k)
             constrFactor[ti] <- 1
             fit <- GIFT_summarycpp(betax, betay, Sigma1, Sigma2, R, constrFactor, pindex, k, n1, n2, maxiter, tol)
             pvalue <- pchisq(2*(loglik0-fit$loglik[length(fit$loglik)]), 1, lower.tail=F)
             return(pvalue)
           }
           ti_list=c(1:k)
           gene_specific_test_pvalue <- unlist(mclapply(ti_list, perform_func, mc.cores = Cores))
         }

         result <- data.frame(gene = gene, causal_effect = H1$alpha, p = gene_specific_test_pvalue)
         return(result)
      }, silent = TRUE)
	
      if(sum(grep("decomposition failed",result[1]))==1){
        for(s in seq(0.1, 1, by = 0.1)){
          result <- try({
            Sigma1 <- (1-s)*Sigma1+s*diag(sum(pindex))
            Sigma2 <- (1-s)*Sigma2+s*diag(sum(pindex))
            H1=GIFT_summarycpp(betax, betay, Sigma1, Sigma2, R, constrFactor, pindex, k, n1, n2, maxiter, tol)
            loglik0 <- H1$loglik[length(H1$loglik)]
          
            Cores=min(c(k,ncores))
            gene_specific_test_pvalue <- NULL 
            if(Cores == 1){
              for (i in 1:k) {
                constrFactor <- numeric(k)
                constrFactor[i] <- 1
                fit <- GIFT_summarycpp(betax, betay, Sigma1, Sigma2, R, constrFactor, pindex, k, n1, n2, maxiter, tol)
                gene_specific_test_pvalue[i] <- pchisq(2*(loglik0-fit$loglik[length(fit$loglik)]), 1, lower.tail=F)
              }
            }else{
              perform_func <- function(ti){
                constrFactor <- numeric(k)
                constrFactor[ti] <- 1
                fit <- GIFT_summarycpp(betax, betay, Sigma1, Sigma2, R, constrFactor, pindex, k, n1, n2, maxiter, tol)
                pvalue <- pchisq(2*(loglik0-fit$loglik[length(fit$loglik)]), 1, lower.tail=F)
                return(pvalue)
              }
              ti_list=c(1:k)
              gene_specific_test_pvalue <- unlist(mclapply(ti_list, perform_func, mc.cores = Cores))
            }
          
            result <- data.frame(gene = gene, causal_effect = H1$alpha, p = gene_specific_test_pvalue)
            return(result)
          }, silent = TRUE)
      
	      if(sum(grep("decomposition failed",result[1]))==0){
            break
          }
        }
		
        if(sum(grep("decomposition failed",result[1]))==1){
          cat("Failed! Please check your inputs carefully. \n")
        }
      }
    }
  }
  
  if(in_sample_LD==F){
    cat("Regularize the LD matrix as (1-s1)*Sigma1+s1*E and (1-s2)*Sigma2+s2*E \n")
    result <- try({
      library(susieR)
      s2 <- estimate_s_rss(Zscore2, Sigma2, n2, r_tol = 1e-08, method = "null-mle")
      E <- diag(dim(Zscore1)[1])
      s1 <- NULL
      for(i in 1:length(gene)){
        s1[i] <- estimate_s_rss(Zscore1[,i], Sigma1, n1, r_tol = 1e-08, method = "null-mle")
      }
      s1 <- max(s1)
      Sigma1 <- (1-s1)*Sigma1+s1*diag(sum(pindex))
      Sigma2 <- (1-s2)*Sigma2+s2*diag(sum(pindex))
      H1=GIFT_summarycpp(betax, betay, Sigma1, Sigma2, R, constrFactor, pindex, k, n1, n2, maxiter, tol)
      loglik0 <- H1$loglik[length(H1$loglik)]
        
      Cores=min(c(k,ncores))
      gene_specific_test_pvalue <- NULL 
      if(Cores == 1){
        for (i in 1:k) {
          constrFactor <- numeric(k)
          constrFactor[i] <- 1
          fit <- GIFT_summarycpp(betax, betay, Sigma1, Sigma2, R, constrFactor, pindex, k, n1, n2, maxiter, tol)
          gene_specific_test_pvalue[i] <- pchisq(2*(loglik0-fit$loglik[length(fit$loglik)]), 1, lower.tail=F)
        }
      }else{
        perform_func <- function(ti){
          constrFactor <- numeric(k)
          constrFactor[ti] <- 1
          fit <- GIFT_summarycpp(betax, betay, Sigma1, Sigma2, R, constrFactor, pindex, k, n1, n2, maxiter, tol)
          pvalue <- pchisq(2*(loglik0-fit$loglik[length(fit$loglik)]), 1, lower.tail=F)
          return(pvalue)
        }
        ti_list=c(1:k)
        gene_specific_test_pvalue <- unlist(mclapply(ti_list, perform_func, mc.cores = Cores))
      }
        
      result <- data.frame(gene = gene, causal_effect = H1$alpha, p = gene_specific_test_pvalue)
      return(result)
    }, silent = TRUE)
    if(sum(grep("decomposition failed",result[1]))==1){
      for(s in seq(0.1, 1, by = 0.1)){
        result <- try({
          Sigma1 <- (1-s)*Sigma1+s*diag(sum(pindex))
          Sigma2 <- (1-s)*Sigma2+s*diag(sum(pindex))
          H1=GIFT_summarycpp(betax, betay, Sigma1, Sigma2, R, constrFactor, pindex, k, n1, n2, maxiter, tol)
          loglik0 <- H1$loglik[length(H1$loglik)]
            
          Cores=min(c(k,ncores))
          gene_specific_test_pvalue <- NULL 
          if(Cores == 1){
            for (i in 1:k) {
              constrFactor <- numeric(k)
              constrFactor[i] <- 1
              fit <- GIFT_summarycpp(betax, betay, Sigma1, Sigma2, R, constrFactor, pindex, k, n1, n2, maxiter, tol)
              gene_specific_test_pvalue[i] <- pchisq(2*(loglik0-fit$loglik[length(fit$loglik)]), 1, lower.tail=F)
            }
          }else{
            perform_func <- function(ti){
              constrFactor <- numeric(k)
              constrFactor[ti] <- 1
              fit <- GIFT_summarycpp(betax, betay, Sigma1, Sigma2, R, constrFactor, pindex, k, n1, n2, maxiter, tol)
              pvalue <- pchisq(2*(loglik0-fit$loglik[length(fit$loglik)]), 1, lower.tail=F)
              return(pvalue)
            }
            ti_list=c(1:k)
            gene_specific_test_pvalue <- unlist(mclapply(ti_list, perform_func, mc.cores = Cores))
          }
            
          result <- data.frame(gene = gene, causal_effect = H1$alpha, p = gene_specific_test_pvalue)
          return(result)
        }, silent = TRUE)
        if(sum(grep("decomposition failed",result[1]))==0){
          break
        }
      }
      if(sum(grep("decomposition failed",result[1]))==1){
        cat("Failed! Please check your inputs carefully. \n")
      }
    }
  }
}

