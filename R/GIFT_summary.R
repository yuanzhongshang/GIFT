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
#' @param pleio The option of controlling the pleiotropy, which can be determined by users. If pleio is set to 0, analysis will perform without controlling any SNP; If pleio is set to 1, analysis will perform  controlling the top SNP; If pleio is set to 2, analysis will perform controlling the top two SNPs.
#' @param ncores The number of cores used in analysis. If the number of cores is greater than 1, analysis will perform with fast parallel computing. The function mclapply() depends on another R package "parallel" in Linux.
#' @param in_sample_LD A logical value, which represents whether in-sample LD was used. If in-sample LD was not used, the LD matrix is regularized to be (1-s1)*Sigma1+s1*E and (1-s2)*Sigma2+s2*E, where s1 and s2 are estimated by estimate_s_rss in susieR, and E is an identity matrix. A grid search is performed over the range from 0.1 to 1 if the estimation does not work. The function estimate_s_rss() depends on another R package "susieR".
#' @param filter A logical value, which represents whether filter the SNP with GWAS P>0.05. This step will improve the computational speed.
#' @return A data frame including the causal effect estimates and p values for the gene-based test. 

GIFT_summary<-function(Zscore1, Zscore2, Sigma1, Sigma2, n1, n2, gene, pindex, R = NULL, maxiter = 100, tol = 1e-3, pleio = 0, ncores = 1, in_sample_LD = F, filter = T){
  
  if(filter==T & n2>100000){
    Zscore2_p<-pchisq(Zscore2^2, df = 1, lower.tail = FALSE) 
    genes<-rep(gene, pindex)
    inclu<-which(Zscore2_p<0.05)
    for (g in unique(gene)) {
      gene_snp_indices <- which(genes == g) 
      if (length(intersect(gene_snp_indices, inclu)) == 0) {
        inclu <- c(inclu, gene_snp_indices[which.min(Zscore2_p[gene_snp_indices])])  
      }
    }
    inclu=sort(inclu)
    Zscore1<-Zscore1[inclu,]
    Zscore2<-Zscore2[inclu]
    Sigma1<-Sigma1[inclu,inclu]
    Sigma2<-Sigma2[inclu,inclu]
    pindex<-as.numeric(table(genes[inclu]))
    pindex_ordered <- names(table(genes[inclu]))
    pindex<-pindex[match(gene, pindex_ordered)]
  }
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
	  if(pleio == 0){
	    H1<-GIFT_summarycpp(betax, betay, Sigma1, Sigma2, R, constrFactor, pindex, k, n1, n2, maxiter, tol)
	    loglik0 <- H1$loglik[length(H1$loglik)]

	    Cores=min(c(k,ncores))
	    gene_specific_test_pvalue <- NULL 
	  
	    if(Cores == 1){
		    for(i in 1:k){
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
      }
	  
	  if(pleio != 0){
	    zscore=abs(Zscore2)
        pleioindex=min(which(zscore==max(zscore)))
        if(pleio == 2){
          zscore1=zscore
          for(i in 1:(length(unique(zscore))-1)){
            zscore1=zscore1[-which(zscore1==max(zscore1))]
            if(abs(Sigma2[pleioindex,min(which(zscore==max(zscore1)))])<0.5){
              pleioindex=c(pleioindex,min(which(zscore==max(zscore1))))
              break
            }
          }
        }
	    H1=GIFT_summarycpppleio(betax, betay, Sigma1, Sigma2, R, constrFactor, pindex, k, n1, n2, maxiter, tol, pleioindex-1)
	    loglik0 <- H1$loglik[length(H1$loglik)]

	    Cores=min(c(k,ncores))
	    gene_specific_test_pvalue <- NULL 
	  
	    if(Cores == 1){
		  for(i in 1:k){
		    constrFactor <- numeric(k)
		    constrFactor[i] <- 1
		    fit <- GIFT_summarycpppleio(betax, betay, Sigma1, Sigma2, R, constrFactor, pindex, k, n1, n2, maxiter, tol, pleioindex-1)
		    gene_specific_test_pvalue[i] <- pchisq(2*(loglik0-fit$loglik[length(fit$loglik)]), 1, lower.tail=F)
		  }
	    }else{
		  library(parallel)
		  perform_func <- function(ti){
		    constrFactor <- numeric(k)
		    constrFactor[ti] <- 1
		    fit <- GIFT_summarycpppleio(betax, betay, Sigma1, Sigma2, R, constrFactor, pindex, k, n1, n2, maxiter, tol, pleioindex-1)
		    pvalue <- pchisq(2*(loglik0-fit$loglik[length(fit$loglik)]), 1, lower.tail=F)
		    return(pvalue)
		  }
		  ti_list=c(1:k)
		  gene_specific_test_pvalue <- unlist(mclapply(ti_list, perform_func, mc.cores = Cores))
	    }
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
		constrFactor <- numeric(k)
		
        if(pleio == 0){
		  H1=GIFT_summarycpp(betax, betay, Sigma1, Sigma2, R, constrFactor, pindex, k, n1, n2, maxiter, tol)
	      loglik0 <- H1$loglik[length(H1$loglik)]

	      Cores=min(c(k,ncores))
	      gene_specific_test_pvalue <- NULL 
	  
	      if(Cores == 1){
		    for(i in 1:k){
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
        }
		
	   if(pleio != 0){
         pleioindex=min(which(Zscore2==max(Zscore2)))
         if(pleio == 2){
          zscore=Zscore2
          for(i in 1:(length(unique(zscore))-1)){
            zscore=zscore[-which(zscore==max(zscore))]
            if(abs(Sigma2[pleioindex,min(which(Zscore2==max(zscore)))])<0.5){
              pleioindex=c(pleioindex,min(which(Zscore2==max(zscore))))
              break
            }
          }
         }
		 H1=GIFT_summarycpppleio(betax, betay, Sigma1, Sigma2, R, constrFactor, pindex, k, n1, n2, maxiter, tol, pleioindex-1)
	     loglik0 <- H1$loglik[length(H1$loglik)]

	     Cores=min(c(k,ncores))
	     gene_specific_test_pvalue <- NULL 
	  
	     if(Cores == 1){
		   for(i in 1:k){
		     constrFactor <- numeric(k)
		     constrFactor[i] <- 1
		     fit <- GIFT_summarycpppleio(betax, betay, Sigma1, Sigma2, R, constrFactor, pindex, k, n1, n2, maxiter, tol, pleioindex-1)
		     gene_specific_test_pvalue[i] <- pchisq(2*(loglik0-fit$loglik[length(fit$loglik)]), 1, lower.tail=F)
		   }
	     }else{
		   perform_func <- function(ti){
		     constrFactor <- numeric(k)
		     constrFactor[ti] <- 1
		     fit <- GIFT_summarycpppleio(betax, betay, Sigma1, Sigma2, R, constrFactor, pindex, k, n1, n2, maxiter, tol, pleioindex-1)
		     pvalue <- pchisq(2*(loglik0-fit$loglik[length(fit$loglik)]), 1, lower.tail=F)
		     return(pvalue)
	       }
		   ti_list=c(1:k)
		   gene_specific_test_pvalue <- unlist(mclapply(ti_list, perform_func, mc.cores = Cores))
		   }
	    }
        result <- data.frame(gene = gene, causal_effect = H1$alpha, p = gene_specific_test_pvalue)
        return(result)
      }, silent = TRUE)
	
      if(sum(grep("decomposition failed",result[1]))==1){
        for(s in seq(0.1, 1, by = 0.1)){
          result <- try({
            Sigma1 <- (1-s)*Sigma1+s*diag(sum(pindex))
            Sigma2 <- (1-s)*Sigma2+s*diag(sum(pindex))
			constrFactor <- numeric(k)
			
            if(pleio == 0){
			  H1=GIFT_summarycpp(betax, betay, Sigma1, Sigma2, R, constrFactor, pindex, k, n1, n2, maxiter, tol)
			  loglik0 <- H1$loglik[length(H1$loglik)]

			  Cores=min(c(k,ncores))
			  gene_specific_test_pvalue <- NULL 
	  
			  if(Cores == 1){
			    for(i in 1:k){
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
            }
			
            if(pleio != 0){
			  pleioindex=min(which(Zscore2==max(Zscore2)))
			  if(pleio == 2){
				zscore=Zscore2
				for(i in 1:(length(unique(zscore))-1)){
				  zscore=zscore[-which(zscore==max(zscore))]
				  if(abs(Sigma2[pleioindex,min(which(Zscore2==max(zscore)))])<0.5){
					pleioindex=c(pleioindex,min(which(Zscore2==max(zscore))))
					break
				  }
				}
			  }
              H1=GIFT_summarycpppleio(betax, betay, Sigma1, Sigma2, R, constrFactor, pindex, k, n1, n2, maxiter, tol, pleioindex-1)
              loglik0 <- H1$loglik[length(H1$loglik)]

              Cores=min(c(k,ncores))
              gene_specific_test_pvalue <- NULL 
	  
              if(Cores == 1){
			    for(i in 1:k){
			      constrFactor <- numeric(k)
			  	  constrFactor[i] <- 1
				  fit <- GIFT_summarycpppleio(betax, betay, Sigma1, Sigma2, R, constrFactor, pindex, k, n1, n2, maxiter, tol, pleioindex-1)
				  gene_specific_test_pvalue[i] <- pchisq(2*(loglik0-fit$loglik[length(fit$loglik)]), 1, lower.tail=F)
			    }
              }else{
			    perform_func <- function(ti){
			      constrFactor <- numeric(k)
			      constrFactor[ti] <- 1
			      fit <- GIFT_summarycpppleio(betax, betay, Sigma1, Sigma2, R, constrFactor, pindex, k, n1, n2, maxiter, tol, pleioindex-1)
			      pvalue <- pchisq(2*(loglik0-fit$loglik[length(fit$loglik)]), 1, lower.tail=F)
			      return(pvalue)
			    }
			    ti_list=c(1:k)
			    gene_specific_test_pvalue <- unlist(mclapply(ti_list, perform_func, mc.cores = Cores))
              }
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
	  
      if(pleio == 0){
	    H1=GIFT_summarycpp(betax, betay, Sigma1, Sigma2, R, constrFactor, pindex, k, n1, n2, maxiter, tol)
	    loglik0 <- H1$loglik[length(H1$loglik)]

	    Cores=min(c(k,ncores))
	    gene_specific_test_pvalue <- NULL 
	   
	    if(Cores == 1){
		  for(i in 1:k){
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
      }
	  
	  if(pleio != 0){
        pleioindex=min(which(Zscore2==max(Zscore2)))
        if(pleio == 2){
          zscore=Zscore2
          for(i in 1:(length(unique(zscore))-1)){
            zscore=zscore[-which(zscore==max(zscore))]
            if(abs(Sigma2[pleioindex,min(which(Zscore2==max(zscore)))])<0.5){
              pleioindex=c(pleioindex,min(which(Zscore2==max(zscore))))
              break
            }
          }
        }
	    H1=GIFT_summarycpppleio(betax, betay, Sigma1, Sigma2, R, constrFactor, pindex, k, n1, n2, maxiter, tol, pleioindex-1)
	    loglik0 <- H1$loglik[length(H1$loglik)]

	    Cores=min(c(k,ncores))
	    gene_specific_test_pvalue <- NULL 
	  
	    if(Cores == 1){
		  for(i in 1:k){
		    constrFactor <- numeric(k)
		    constrFactor[i] <- 1
		    fit <- GIFT_summarycpppleio(betax, betay, Sigma1, Sigma2, R, constrFactor, pindex, k, n1, n2, maxiter, tol, pleioindex-1)
		    gene_specific_test_pvalue[i] <- pchisq(2*(loglik0-fit$loglik[length(fit$loglik)]), 1, lower.tail=F)
		  }
	    }else{
		  library(parallel)
		  perform_func <- function(ti){
		    constrFactor <- numeric(k)
		    constrFactor[ti] <- 1
		    fit <- GIFT_summarycpppleio(betax, betay, Sigma1, Sigma2, R, constrFactor, pindex, k, n1, n2, maxiter, tol, pleioindex-1)
		    pvalue <- pchisq(2*(loglik0-fit$loglik[length(fit$loglik)]), 1, lower.tail=F)
		    return(pvalue)
		  }
		  ti_list=c(1:k)
		  gene_specific_test_pvalue <- unlist(mclapply(ti_list, perform_func, mc.cores = Cores))
	    }
	  }
      result <- data.frame(gene = gene, causal_effect = H1$alpha, p = gene_specific_test_pvalue)
      return(result)
    }, silent = TRUE)
	
    if(sum(grep("decomposition failed",result[1]))==1){
      for(s in seq(0.1, 1, by = 0.1)){
        result <- try({
          Sigma1 <- (1-s)*Sigma1+s*diag(sum(pindex))
          Sigma2 <- (1-s)*Sigma2+s*diag(sum(pindex))
		  constrFactor <- numeric(k)
		  
          if(pleio == 0){
			H1=GIFT_summarycpp(betax, betay, Sigma1, Sigma2, R, constrFactor, pindex, k, n1, n2, maxiter, tol)
			loglik0 <- H1$loglik[length(H1$loglik)]

			Cores=min(c(k,ncores))
			gene_specific_test_pvalue <- NULL 
	  
			if(Cores == 1){
			  for(i in 1:k){
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
          }
		  
          if(pleio != 0){
			pleioindex=min(which(Zscore2==max(Zscore2)))
			if(pleio == 2){
			  zscore=Zscore2
			  for(i in 1:(length(unique(zscore))-1)){
				zscore=zscore[-which(zscore==max(zscore))]
				if(abs(Sigma2[pleioindex,min(which(Zscore2==max(zscore)))])<0.5){
				  pleioindex=c(pleioindex,min(which(Zscore2==max(zscore))))
				break
				}
			  }
			}
			H1=GIFT_summarycpppleio(betax, betay, Sigma1, Sigma2, R, constrFactor, pindex, k, n1, n2, maxiter, tol, pleioindex-1)
			loglik0 <- H1$loglik[length(H1$loglik)]

			Cores=min(c(k,ncores))
			gene_specific_test_pvalue <- NULL 
			
			if(Cores == 1){
			  for(i in 1:k){
				constrFactor <- numeric(k)
				constrFactor[i] <- 1
				fit <- GIFT_summarycpppleio(betax, betay, Sigma1, Sigma2, R, constrFactor, pindex, k, n1, n2, maxiter, tol, pleioindex-1)
				gene_specific_test_pvalue[i] <- pchisq(2*(loglik0-fit$loglik[length(fit$loglik)]), 1, lower.tail=F)
			  }
			}else{
			  perform_func <- function(ti){
			    constrFactor <- numeric(k)
			    constrFactor[ti] <- 1
			    fit <- GIFT_summarycpppleio(betax, betay, Sigma1, Sigma2, R, constrFactor, pindex, k, n1, n2, maxiter, tol, pleioindex-1)
			    pvalue <- pchisq(2*(loglik0-fit$loglik[length(fit$loglik)]), 1, lower.tail=F)
			    return(pvalue)
			  }
			  ti_list=c(1:k)
			  gene_specific_test_pvalue <- unlist(mclapply(ti_list, perform_func, mc.cores = Cores))
			}
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

