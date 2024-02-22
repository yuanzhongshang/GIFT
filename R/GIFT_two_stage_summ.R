#' @title The two-stage version of GIFT for conditional fine-mapping in transcriptome-wide association studies with existing gene expression prediction models
#' @description GIFT_two_stage_summ applies a two-stage approach, accounting for the correlated cis-SNPs and genes in a region. 
#' @param betax Weight matrix for the cis-SNP effect size.
#' @param betay Beta vector for the cis-SNP effect size from GWAS data.
#' @param se_betay Corresponding se vector for the cis-SNP effect size vector from GWAS data.
#' @param Sigma LD matrix from GWAS data.
#' @param n Sample size of GWAS data.
#' @param gene The gene names vector.
#' @param in_sample_LD A logical value, which represents whether in-sample LD was used. If in-sample LD was not used, the LD matrix is regularized to be (1-s)*Sigma+s*E and (1-s2)*Sigma2+s2*E, where s is estimated by estimate_s_rss in susieR, and E is an identity matrix.
#' @return A data frame including the z scores and p values for the gene-based test. 

GIFT_two_stage_summ <- function(betax, betay, se_betay, Sigma, n, gene, in_sample_LD = F){
  
  if(sum(is.na(Sigma)) != 0){
    cat("Correlation Matrix cannot contain NA values.", "\n")
  }else if(var(c(nrow(betax), nrow(betay), nrow(se_betay), nrow(Sigma), ncol(Sigma))) != 0){
    cat("The numbers of rows in betax, betay, se_betay, and Sigma (rows and columns) are not matched.", "\n")
  }else{

  if(in_sample_LD==F){
     library(susieR)
     s <- estimate_s_rss(betay/se_betay, Sigma, n, r_tol = 1e-08, method = "null-mle")
     Sigma <- (1-s)*Sigma + s*diag(length(betay))
  }
    nsnps <- nrow(betay)
    
    ZTY <- matrix(diag(Sigma) * betay, ncol = 1);
    YTY_list <- list()
    for(SNP in 1:nrow(Sigma)){
      YTY_list[[length(YTY_list) + 1]] <- (n)*Sigma[SNP, SNP]*(se_betay^2)[SNP, 1] + betay[SNP,1]*ZTY[SNP,]
    }
    
    YTY <- median((do.call("rbind", YTY_list))[,1])
    
    beta <- matrix(solve(t(betax) %*% Sigma %*% betax) %*% t(betax) %*% (ZTY), ncol = 1)
    sigma0 <- as.numeric(unlist((YTY - t(ZTY)%*%betax%*%solve(t(betax)%*%Sigma%*%betax)%*%t(betax)%*%ZTY)/(n - ncol(betax))))
    se <- matrix(sqrt(diag(solve(t(betax) %*% Sigma %*% betax) * c(sigma0))), ncol = 1)
    Z <- matrix(beta/se, ncol = 1)
    p <- matrix(2*pnorm(-abs(Z), 0, 1), ncol = 1)
    result <- data.frame(gene = gene, z = Z, p = p)
    
    return(result)
  }
}

