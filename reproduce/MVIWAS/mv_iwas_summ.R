#' MV-IWAS using Summary Statistics
#'
#' This function produces causal estimates for phenotype-disease associations using MV-IWAS and MV-IWAS-Egger using GWAS summary statistics for a disease and exposure phenotypes (Knutson, Deng, and Pan, 2020)
#' @param betaZY A one column matrix with estimated SNP effects on the disease from GWAS summary data. There should be a separate row for each SNP. The order of SNPs in rows should match across betaZY, se_betaZY, betaZX, se_betaZX, and corr_mat.
#' @param se_betaZY A one column matrix with estimated SE of SNP effects on the disease from GWAS summary data. There should be a separate row for each SNP. The order of SNPs in rows should match across betaZY, se_betaZY, betaZX, se_betaZX, and corr_mat.
#' @param betaZX A matrix with estimated SNP effects on each exposure phenotype of interest, taken from GWAS summary data. There should be a separate column for each exposure phenotype of interest and a separate row for each SNP. The order of SNPs in rows should match across betaZY, se_betaZY, betaZX, se_betaZX, and corr_mat. For SNPs which are not to be included in Stage 1 for a given exposure phenotype, set the corresponding SNP row to 0 for that phenotype's column. Do not use NA. If the number of columns is equal to 1, this test will be equivalent to the so-called UV-IWAS test.
#' @param se_betaZX A matrix with estimated SE of SNP effects on each exposure phenotype of interest, taken from GWAS summary data. There should be a separate column for each exposure phenotype of interest and a separate row for each SNP. The order of SNPs in rows should match across betaZY, se_betaZY, betaZX, se_betaZX, and corr_mat. For SNPs which are not to be included in Stage 1 for a given exposure phenotype, set the corresponding SNP row to 0 for that phenotype's column. Do not use NA.
#' @param corr_mat The covariance or correlation matrix of SNPs, estimated from a reference panel of similar ancestry to the GWAS summary data. The order of SNPs in the columns and rows should match each other and the row orders of betaZY, se_betaZY, betaZX, se_betaZX.
#' @param n The sample size used to estimate the disease GWAS summary data
#' @param trait_type Either "Continuous" or "Binary", for continous or binary disease traits
#' @param n_case The number of cases used in the disease GWAS sample. Required for trait_type = "Binary", default for Continuous trait is NULL
#' @param n_control The number of controls used in the disease GWAS sample. Required for trait_type = "Binary", default for Continuous trait is NULL
#' @import dplyr
#' @export
#' @examples
#' mv_iwas_summ()

mv_iwas_summ <- function(betaZY, se_betaZY, betaZX, corr_mat, n, trait_type, n_case = NULL, n_control = NULL){
  
  
  cat("Have you made sure that all betaZY, se_betaZY, and betaZX are the same order by SNP ID ?", "\n")
  
  if(sum(c(missing(betaZY), missing(betaZX), missing(se_betaZY), missing(corr_mat), missing(n), missing(trait_type))) > 0){
    cat("At least one of the arguments betaZY, se_betaZY, betaZX, corr_mat, n, or trait_type is missing \n")
  }else if(sum(is.na(corr_mat)) != 0){
    cat("Correlation Matrix cannot contain NA values.", "\n")
  }else if(!(trait_type %in% c("Continuous", "Binary"))){
    cat("Trait Type must be specified as Continuous or Binary \n")
  }else if((trait_type == "Binary" && sum(c(is.null(n_case), is.null(n_control))) > 0)){
    cat("Trait is Binary and either n_case and/or n_control are missing \n")
  }else if((trait_type == "Binary" && n_case + n_control != n)){
    cat("n_case + n_control not equal to n! \n")
  }else if(sum(is.na(betaZY)) != 0){
    cat("betaZY cannot contain NA values.", "\n")
  }else if(sum(is.na(betaZX)) != 0){
    cat("betaZX cannot contain NA values.", "\n")
  }else if(sum(is.na(se_betaZY)) != 0){
    cat("se_betaZY cannot contain NA values.", "\n")
  }else if(ncol(betaZY) != 1){
    cat("betaZY must be a matrix with only 1 column", "\n")
  }else if(ncol(se_betaZY) != 1){
    cat("se_betaZY must be a matrix with only 1 column", "\n")
  }else if(var(c(nrow(betaZX), nrow(betaZY), nrow(se_betaZY), nrow(corr_mat), ncol(corr_mat))) != 0){
    cat("the number of rows in betaZX, betaZY, se_betaZY, and corr_mat (rows and columns) must be equal.", "\n")
  }else{
    
    if(trait_type == "Binary"){
      cat("Converting OR to Linear Coefficients \n")
      exp_b0 <- n_control/n_case
      fac <- exp_b0/(1 + exp_b0)^2
      betaZY <- betaZY*fac
      se_betaZY <- se_betaZY*fac
    }
    nsnps <- nrow(betaZY)
    
    ZTY <- matrix(diag(corr_mat) * betaZY, ncol = 1);
    YTY_list <- list()
    for(SNP in 1:nrow(corr_mat)){
      YTY_list[[length(YTY_list) + 1]] <- (n)*corr_mat[SNP, SNP]*(se_betaZY^2)[SNP, 1] + betaZY[SNP,1]*ZTY[SNP,]
    }
    
    YTY <- median((do.call("rbind", YTY_list))[,1])
    
    beta <- matrix(solve(t(betaZX) %*% corr_mat %*% betaZX) %*% t(betaZX) %*% (ZTY), ncol = 1)
    sigma <- as.numeric(unlist((YTY - t(ZTY)%*%betaZX%*%solve(t(betaZX)%*%corr_mat%*%betaZX)%*%t(betaZX)%*%ZTY)/(n - ncol(betaZX))))
    se <- matrix(sqrt(diag(solve(t(betaZX) %*% corr_mat %*% betaZX) * c(sigma))), ncol = 1)
    Z <- matrix(beta/se, ncol = 1)
    p <- matrix(2*pnorm(-abs(Z), 0, 1), ncol = 1)
    result <- data.frame(MODEL = "MV-IWAS", TERM = paste0("Phenotype", 1:ncol(betaZX)), BETA = beta, SE = se, Z = Z, P = p, NSNP = nsnps, DiseaseType = trait_type)
    
    cat("Finished with MV-IWAS \n");    
    return(result)
    
    
  }
}
