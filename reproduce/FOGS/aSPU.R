#' Weighted Adaptive Sum of Powered Score tests (waSPU) test for single trait with single or multiple weights for each genetic marker
#'
#' It returns p-values.
#'
#' @param U Score vector for the genetic marker set (p by 1 matrix).
#'
#' @param V Corresponding covariance matrix for the score vector (p by p matrix).
#'
#' @param weight Single weight for each genetic markers (p by 1 matrix).
#'
#' @param pow Gamma sets. Different gamma corresponds to different existing test. For example, gamma = 1 equals Sum test; gamma = 2 equals SSU or SKAT with linear kernel; gamma = Inf is similar to minimum p value test. We recommend use pow = c(1:6, Inf) to maintain high power under various scenarios, although other choice may be slightly powerful in some specific situations.
#'
#' @param n.perm number of permutations.
#'
#' @export
#'
#' @return P-values for aSPU test.
#'
#' @author Chong Wu, Zhiyuan Xu, and Wei Pan

aSPU <- function(U,V,weight, pow = c(1:8,Inf),n.perm = 1000) {
    
    weight = as.matrix(weight)
    
    # remove SNPs corresponding to zero weight
    weight.tmp = abs(weight)
    index = rowSums(weight.tmp) > 0
    U = as.matrix(U)
    U = U[index,]
    V = V[index,]
    V = V[,index]
    weight = weight[index,]
    
    weight = as.matrix(weight)
    
    Ts <- rep(0, length(pow))
    for(j in 1:length(pow)){
        if (pow[j] < Inf) {
            Ts[j] <- sum((U * weight)^pow[j])
        } else {
            Ts[j] <- max(abs(weight * U))
        }
    }
    
    eV <- eigen(V)
    eV$values[eV$values<0] = 0
    
    CovSsqrt <- t(eV$vectors %*% (t(eV$vectors) * sqrt(eV$values)))
    pow[pow==Inf] <- 0 # pass 0 as infitiy
    
    T0s <- big.matrix(n.perm,length(pow),type = "double")
    
    Ts.abs <- abs(Ts)
    
    Res = calcT0Wsim3(as.matrix(CovSsqrt), as.matrix(weight),as.matrix(pow),as.matrix(Ts.abs),T0s@address,n.perm)
    
    minp0 = Res$minp0
    pPerm0 = Res$pPerm0
    
    Paspu <- (sum(minp0 <= min(pPerm0)) + 1) / (n.perm+1)
    pvs <- c(pPerm0, Paspu)
    
    Ts <- c(Ts, min(pPerm0))
    pow[pow==0] <- Inf
    names(Ts) <- c(paste("SPU(", pow,")", sep=""), "aSPU")
    names(pvs) <- names(Ts)
    
    return(list(Ts = Ts, pvs = pvs))
}

