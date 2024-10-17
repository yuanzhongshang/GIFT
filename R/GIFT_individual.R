#' @title The main function for conditional fine-mapping for in transcriptome-wide association studies with individual-level data
#' @description GIFT_individual applies a likelihood-based approach, accounting for the correlated cis-SNPs and genes in a region 
#' @param X Complex gene expression matrix, each column is the standardized specific expression.
#' @param Y Standarized trait vector.
#' @param Zx Standardized cis-genotype matrix from eQTL data.
#' @param Zy Standardized cis-genotype matrix from GWAS data.
#' @param gene The gene names vector.
#' @param pindex A vector with each element represents the number of cis-SNPs for each gene.
#' @param maxiter The maximum iteration, which can be determined by users.
#' @param tol The convergence tolerance of the absolute value of the difference between the nth and (n+1)th log likelihood, which can be determined by users.
#' @param pleio The option of controlling the pleiotropy, which can be determined by users. If pleio is set to 0, analysis will perform without controlling any SNP; If pleio is set to 1, analysis will perform  controlling the top SNP; If pleio is set to 2, analysis will perform controlling the top two SNPs.
#' @param ncores The number of cores used in analysis. If the number of cores is greater than 1, analysis will perform with fast parallel computing. The function mclapply() depends on another R package "parallel" in Linux.
#' @param filter A logical value, which represents whether filter the SNP with GWAS P>0.05. This step will improve the computational speed.
#' @return A data frame including the causal effect estimates and p values for the gene-based test. 

GIFT_individual<-function(X, Y, Zx, Zy, gene, pindex, maxiter =1000, tol=1e-4, pleio=0, ncores=1, filter = T){
  
  k<-length(pindex)
  eQTLdata<-na.omit(cbind(X,Zx))
  xin_NA<-c(eQTLdata[,1:k])
  R=cor(eQTLdata[,1:k])
  Zxin_NA<-eQTLdata[,-c(1:k)]

  GWASdata<-na.omit(cbind(Y,Zy))
  Yin_NA<-as.matrix(GWASdata[,1])
  Zyin_NA<-GWASdata[,-1]
  
  Zxin_NA_sd<-as.vector(apply(Zxin_NA,2,sd))
  
  Zyin_NA_sd<-as.vector(apply(Zyin_NA,2,sd))
  
  if(sd(Zxin_NA)!=0 & sd(Zyin_NA)!=0 & sum(Zxin_NA_sd==0)==0 & sum(Zyin_NA_sd==0)==0){
    
    X<-as.vector(xin_NA)
    Y<-as.vector(scale(Yin_NA))
    Zx<-scale(Zxin_NA)
    Zy<-scale(Zyin_NA)
    if(filter==T & sum(pindex)>500){
    Zscore2_p<-pchisq(((t(Zy) %*% Y)/sqrt(dim(Zy)[1]-1))^2, df = 1, lower.tail = FALSE) 
    genes<-rep(gene, pindex)
    inclu<-which(Zscore2_p<0.05)
    for (g in unique(gene)) {
      gene_snp_indices <- which(genes == g) 
      if (length(intersect(gene_snp_indices, inclu)) == 0) {
        inclu <- c(inclu, gene_snp_indices[which.min(Zscore2_p[gene_snp_indices])])  
      }
    }
    inclu=sort(inclu)
    Zx<-Zx[,inclu]
    Zy<-Zy[,inclu]
    pindex<-as.numeric(table(genes[inclu]))
    pindex_ordered <- names(table(genes[inclu]))
    pindex<-pindex[match(gene, pindex_ordered)]
  }
    #################################################################################
    if(pleio == 0){
      constrFactor <- numeric(k)
      H1 = GIFT_individualcpp(X,Y,Zx,Zy,R,constrFactor,pindex,maxiter,tol)
      loglik0 <- H1$loglik[length(H1$loglik)]
      
      Cores=min(c(k,ncores))
      gene_specific_test_pvalue <- NULL 
      if(Cores == 1){
        for(i in 1:k){
          constrFactor <- numeric(k)
          constrFactor[i] <- 1
          fit <- GIFT_individualcpp(X,Y,Zx,Zy,R,constrFactor,pindex,maxiter,tol)
          gene_specific_test_pvalue[i] <- pchisq(2*(loglik0-fit$loglik[length(fit$loglik)]), 1, lower.tail=F)
        }
      }else{
        perform_func <- function(ti){
          constrFactor <- numeric(k)
          constrFactor[ti] <- 1
          fit <- GIFT_individualcpp(X,Y,Zx,Zy,R,constrFactor,pindex,maxiter,tol)
          pvalue <- pchisq(2*(loglik0-fit$loglik[length(fit$loglik)]), 1, lower.tail=F)
          return(pvalue)
        }
        ti_list=c(1:k)
        gene_specific_test_pvalue <- unlist(mclapply(ti_list, perform_func, mc.cores = Cores))
    }
    }
    if(pleio != 0){
      zscore <- as.numeric(abs(t(Zy) %*% Y / sqrt(length(Y))))
      pleioindex=min(which(zscore==max(zscore)))
      if(pleio == 2){
        zscore1=zscore
        for(i in 1:(length(unique(zscore))-1)){
          zscore1=zscore1[-which(zscore1==max(zscore1))]
          if(abs(cor(Zy[,c(pleioindex,min(which(zscore==max(zscore1))))])[1,2])<0.5){
            pleioindex=c(pleioindex,min(which(zscore==max(zscore1))))
            break
          }
        }
      }
      constrFactor <- numeric(k)
      H1 = GIFT_individualcpppleio(X,Y,Zx,Zy,R,constrFactor,pindex,maxiter,tol,pleioindex-1)
      loglik0 <- H1$loglik[length(H1$loglik)]
      
      Cores=min(c(k,ncores))
      gene_specific_test_pvalue <- NULL 
      if(Cores == 1){
        for(i in 1:k){
          constrFactor <- numeric(k)
          constrFactor[i] <- 1
          fit <- GIFT_individualcpppleio(X,Y,Zx,Zy,R,constrFactor,pindex,maxiter,tol,pleioindex-1)
          gene_specific_test_pvalue[i] <- pchisq(2*(loglik0-fit$loglik[length(fit$loglik)]), 1, lower.tail=F)
        }
      }else{
        perform_func <- function(ti){
          constrFactor <- numeric(k)
          constrFactor[ti] <- 1
          fit <- GIFT_individualcpppleio(X,Y,Zx,Zy,R,constrFactor,pindex,maxiter,tol,pleioindex-1)
          pvalue <- pchisq(2*(loglik0-fit$loglik[length(fit$loglik)]), 1, lower.tail=F)
          return(pvalue)
        }
        ti_list=c(1:k)
        gene_specific_test_pvalue <- unlist(mclapply(ti_list, perform_func, mc.cores = Cores))
    }
    }
	
    #result=list()
    #result$causal_effect=H1$alpha
    #result$gene_based_test_pvalue=gene_specific_test_pvalue
    #result$sigma_cisSNP=H1$sigmaZ
    #result$sigma_error_1=H1$sigmaX
    #result$sigma_error_2=H1$sigmaY
    result <- data.frame(gene = gene, causal_effect = H1$alpha, p = gene_specific_test_pvalue)
    return(result)
    
    ######################
  }
}
