#SMI
SMI = function(BM,SM,N,Bref,lambda = 0.1,mult=30){
  q=length(BM)
  grp=1:q
  n=N

  if(mult>1){
  bdie=matrix(nrow=mult,ncol=q)
  bdie2=matrix(nrow=mult,ncol=q)
  withindie=matrix(0,nrow=q,ncol=q)
  withindie2=matrix(0,nrow=q,ncol=q)
  Bref2=Bref

  for(ii in 1:mult){
    sasa2=sample(1:nrow(Bref2),nrow(Bref2),replace=TRUE)
    if(ii==1){
      sasa2=1:nrow(Bref2)
    }
    Bref2s=Bref2[sasa2,]
    kk2=JointRidge(BM,SM,N=n,XX=cov(Bref2s),lambda)
    bdie2[ii,]=kk2$beta
    if(ii==1){
      withindie2=kk2$cov
    }
  }

  betweendie2=cov(bdie2[-1,])
  vdie2=betweendie2*(1+1/mult)+withindie2
  bb30=as.matrix(bdie2[1,])
  cc30=vdie2
  bb3=as.matrix(bb30[grp])
  cc3=cc30[grp,grp]
  }else{
    kk2=JointRidge(BM,SM,N=n,XX=cov(Bref),lambda)
    bb3=kk2$beta
    cc3=kk2$cov
  }
  ll=t(bb3)%*%ginv(cc3)%*%bb3
  dff=qr(cc3)$rank
  pp=pchisq(ll,df=dff,lower.tail=FALSE)
  return(list(beta=bb3,cov=cc3,chisq=ll,df=dff,pvalue=pp))
}



SMI2 = function(BM,SM,N,Bref,mult=30){ # regularized version, a older version; we find a better way to regularized later
    q=length(BM)
    grp=1:q
    n=N
    
    if(mult>1){
        bdie=matrix(nrow=mult,ncol=q)
        bdie2=matrix(nrow=mult,ncol=q)
        withindie=matrix(0,nrow=q,ncol=q)
        withindie2=matrix(0,nrow=q,ncol=q)
        Bref2=Bref
        
        for(ii in 1:mult){
            sasa2=sample(1:nrow(Bref2),nrow(Bref2),replace=TRUE)
            if(ii==1){
                sasa2=1:nrow(Bref2)
            }
            Bref2s=Bref2[sasa2,]
            XX.tmp = cov(Bref2s)
            XX.tmp = XX.tmp +  diag(0.1,dim(XX.tmp)[1],dim(XX.tmp)[2])
            kk2=JointSum(BM,SM,N=n,XX=XX.tmp)
            bdie2[ii,]=kk2$beta
            if(ii==1){
                withindie2=kk2$cov
            }
        }
        
        betweendie2=cov(bdie2[-1,])
        vdie2=betweendie2*(1+1/mult)+withindie2
        bb30=as.matrix(bdie2[1,])
        cc30=vdie2
        bb3=as.matrix(bb30[grp])
        cc3=cc30[grp,grp]
    }else{
        XX.tmp = cov(Bref)
        XX.tmp = XX.tmp + diag(0.1,dim(XX.tmp)[1],dim(XX.tmp)[2])
        kk2=JointSum(BM,SM,N=n,XX=XX.tmp)
        bb3=kk2$beta
        cc3=kk2$cov
    }
    
    ll=t(bb3)%*%ginv(cc3)%*%bb3
    dff=qr(cc3)$rank
    pp=pchisq(ll,df=dff,lower.tail=FALSE)
    return(list(beta=bb3,cov=cc3,chisq=ll,df=dff,pvalue=pp))
}

