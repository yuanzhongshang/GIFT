#conditional model (call)
library(MASS)

JointSum=function(B,S,N,XX=diag(1,nrow=1)){ #x first
  B = as.matrix(B)
  S = as.matrix(S)
  N = as.matrix(N)
  
  n=max(N)
  
  nrx = dim(B)[1]
  yy1 = NULL
  
  for(rx in 1:nrx){

    Dj = XX[rx,rx]
    Sj = S[rx]
    Bj = B[rx]
    
    yy1 = c(yy1, N[rx,1] / n * (n-1) * Dj * Sj^2 + Dj * Bj^2)
    
  }
  
  yy = median(yy1)

  xxd = diag(XX)
  XY1 = xxd * B      #vector

  B = c(XY1)
  
  nx = ncol(XX)
  sA = ginv(XX)

  beta= sA %*% B
  
  sigma2=(yy - t(beta) %*% B ) / (n- nx)
  
  se = sigma2[1,1] * sA   #cov

  pvalue = 2 * pnorm(abs(beta/sqrt(diag(se))), lower.tail=FALSE)
  
  return(list(beta=beta,cov=se,pvalue=pvalue,sigma2=sigma2))
}



JointRidge = function(B,S,N,XX=diag(1,nrow=1),lambda= 0.1){ #x first
    
    B = as.matrix(B)
    S = as.matrix(S)
    N = as.matrix(N)
    
    n=max(N)
    
    nrx = dim(B)[1]
    yy1 = NULL
    
    for(rx in 1:nrx){
        
        Dj = XX[rx,rx]
        Sj = S[rx]
        Bj = B[rx]
        
        yy1 = c(yy1, N[rx,1] / n * (n-1) * Dj * Sj^2 + Dj * Bj^2)
        
    }
    
    yy = median(yy1)
    
    xxd = diag(XX)
    XY1 = xxd * B      #vector
    
    B = c(XY1)
    
    
    sA = ginv( (XX + diag(lambda,dim(XX)[1],dim(XX)[2])) )
    
    beta = sA %*% B
    
    eign = eigen(XX)$values
    df = sum(eign/(eign + lambda))
    
    sigma2 = (yy - t(beta) %*% B ) / (n- df )
    
    se = sigma2[1,1] * sA %*% XX %*% sA   #cov
    
    pvalue = 2 * pnorm(abs(beta/sqrt(diag(se))) ,lower.tail=FALSE)
    
    return(list(beta=beta,cov=se,pvalue=pvalue,sigma2=sigma2,yy = yy,xy = XY1,df = df))
}


AIC <- function(B,S,N,XX, lambda.set = seq(0.01,0.5,length.out= 100) ) {
    
    res = matrix(NA,length(lambda.set),3)
    
    n = max(N)
    for(i in 1:length(lambda.set)) {
        tmp = JointRidge(B,S,N,XX,lambda.set[i])
        yy = tmp$yy
        xy = tmp$xy
        beta = tmp$beta
        df = tmp$df
        res[i,1] = lambda.set[i]
        RSS = yy - 2 * t(beta) %*% xy + t(beta) %*% XX %*% beta
        res[i,2] = n * log(RSS) + 2 * df
        res[i,3] = n * log(RSS) + df * log(n)
    }
    colnames(res) = c("lambda","AIC","BIC")
    return(res)
}


conRidge = function(B,S,N,XX,lambda.set= seq(0.01,0.1,length.out= 100) ){ #x first
    AIC.res = AIC(B,S,N,XX,lambda.set)
    AIC.res = AIC.res[AIC.res[,2] == min(AIC.res[,2]),]
    
    AIC.res = as.matrix(AIC.res)
    lambda = AIC.res[1,1]
    
    res = JointRidge(B,S,N,XX,lambda)
    return(res)
}
