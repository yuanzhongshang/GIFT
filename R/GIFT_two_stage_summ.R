#' @title The two-stage version of GIFT for conditional fine-mapping in transcriptome-wide association studies with existing gene expression prediction models
#' @description GIFT_two_stage_summ applies a two-stage approach, accounting for the correlated cis-SNPs and genes in a region. 
#' @param betax Weight matrix for the cis-SNP effect size.
#' @param betay Beta vector for the cis-SNP effect size from GWAS data.
#' @param se_betay Corresponding se vector for the cis-SNP effect size vector from GWAS data.
#' @param Sigma LD matrix from GWAS data.
#' @param n Sample size of GWAS data.
#' @param gene The gene names vector.
#' @return A data frame including the z scores, p values, number of SNPs analyzed for the gene-based test. 

GIFT_two_stage_summ <- function(betax, betay, se_betay, Sigma, n, gene){
  
  if(sum(is.na(Sigma)) != 0){
    cat("Correlation Matrix cannot contain NA values.", "\n")
  }else if(var(c(nrow(betax), nrow(betay), nrow(se_betay), nrow(Sigma), ncol(Sigma))) != 0){
    cat("The numbers of rows in betax, betay, se_betay, and Sigma (rows and columns) are not matched.", "\n")
  }else{
    
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
    result <- data.frame(gene = gene, Z = Z, P = p, NSNP = nsnps)
    
    return(result)
  }
}

