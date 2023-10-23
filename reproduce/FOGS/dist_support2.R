### SSU #################

allele.qc = function(a1,a2,ref1,ref2) {
    ref = ref1
    flip = ref
    flip[ref == "A"] = "T"
    flip[ref == "T"] = "A"
    flip[ref == "G"] = "C"
    flip[ref == "C"] = "G"
    flip1 = flip
    
    ref = ref2
    flip = ref
    flip[ref == "A"] = "T"
    flip[ref == "T"] = "A"
    flip[ref == "G"] = "C"
    flip[ref == "C"] = "G"
    flip2 = flip;
    
    snp = list()
    snp[["keep"]] = !((a1=="A" & a2=="T") | (a1=="T" & a2=="A") | (a1=="C" & a2=="G") | (a1=="G" & a2=="C")) & (a1=="A" | a1=="T" | a1=="G" | a1=="C") & (a2=="A" | a2=="T" | a2=="G" | a2=="C") & (ref1=="A" | ref1=="T" | ref1=="G" | ref1=="C") & (ref2=="A" | ref2=="T" | ref2=="G" | ref2=="C")
    snp[["flip"]] = (a1 == ref2 & a2 == ref1) | (a1 == flip2 & a2 == flip1)
    
    return(snp)
}


SumSqU<-function(U, CovS){
    if (is.null(dim(CovS))) {# only one-dim:
        Tscore<- sum(U^2 /CovS)
        if (is.na(Tscore) || is.infinite(Tscore) || is.nan(Tscore)) Tscore<-0
        pTg1<-as.numeric(1-pchisq(Tscore, 1))
    }
    else {
        #it's possible U=0 and Cov(U)=0:
        if (all(abs(U)<1e-20)) pTg1<-1 else{
            Tg1<- t(U) %*% U
            ##distr of Tg1 is sum of cr Chisq_1:
            cr<-eigen(CovS, only.values=TRUE)$values
            ##approximate the distri by alpha Chisq_d + beta:
            alpha1<-sum(cr*cr*cr)/sum(cr*cr)
            beta1<-sum(cr) - (sum(cr*cr)^2)/(sum(cr*cr*cr))
            d1<-(sum(cr*cr)^3)/(sum(cr*cr*cr)^2)
            alpha1<-as.double(alpha1)
            beta1<-as.double(beta1)
            d1<-as.double(d1)
            pTg1<-as.numeric(1-pchisq((Tg1-beta1)/alpha1, d1))
        }
    }
    return(pTg1)
}

##########SumTest########################
Sum<-function(U, CovS){
    #it's possible U=0 and Cov(U)=0:
    if (all(abs(sum(U))<1e-20)) pTsum<-1 else{
        a<-rep(1, length(U))
        Tsum<- sum(U)/(sqrt(as.numeric(t(a) %*% CovS %*% (a))))
        pTsum <- as.numeric( 1-pchisq(Tsum^2, 1) )
    }
    out = list(pval = pTsum, T = Tsum)
    
    return(out)
}

##########UminP Test########################
UminPd<-function(U, CovS){
    
    if (is.null(dim(CovS))) {# only one-dim:
        Tu<- sum(U^2 /CovS)
        if (is.na(Tu) || is.infinite(Tu) || is.nan(Tu)) Tu<-0
        pTu<-as.numeric(1-pchisq(Tu, 1))
    }
    else{
        ####it's POSSIBLR Ui=0 and CovS[i,i]=0!
        Tu<-as.vector(abs(U)/(sqrt(diag(CovS)) + 1e-20) )
        k<-length(U)
        V <- matrix(0,nrow=k, ncol=k)
        for(i in 1:k){
            for(j in 1:k){
                if (abs(CovS[i,j])>1e-20)
                V[i,j] <- CovS[i,j]/sqrt(CovS[i,i]*CovS[j,j])
                else   V[i,j] <- 1e-20
            }
        }
        pTu <- as.numeric(PowerUniv(Tu,V))
    }
    
    pTu
}


PowerUniv <- function(U,V){
    n <- dim(V)[1]
    
    x <- as.numeric(max(abs(U)))
    TER <- as.numeric(1-pmvnorm(lower=c(rep(-x,n)),upper=c(rep(x,n)),mean=c(rep(0,n)),sigma=V))
    
    return(TER)
}

calc.pvalue <- function(cur.Z) {
    # Compute LD matrix
    if ( length(cur.Z) == 0 ) {
        cat( "WARNING : " , unlist(wgtlist0[w,]) , " had no overlapping SNPs\n")
        cur.FAIL = TRUE
    } else {
        cur.LD = t(cur.genos) %*% cur.genos / (nrow(cur.genos)-1)
        
        cur.miss = is.na(cur.Z)
        # Impute missing Z-scores
        if ( sum(cur.miss) != 0 ) {
            if ( sum(!cur.miss) == 0 ) {
                cat( "WARNING : " , unlist(wgtlist0[w,]) , " had no overlapping GWAS Z-scores\n")
                cur.FAIL = TRUE
            } else {
                cur.wgt =  cur.LD[cur.miss,!cur.miss] %*% solve( cur.LD[!cur.miss,!cur.miss] + 0.1 * diag(sum(!cur.miss)) )
                cur.impz = cur.wgt %*% cur.Z[!cur.miss]
                cur.r2pred = diag( cur.wgt %*% cur.LD[!cur.miss,!cur.miss] %*% t(cur.wgt) )
                cur.Z[cur.miss] = cur.impz / sqrt(cur.r2pred)
            }
        }
        
        if ( !cur.FAIL ) {
            # Compute TWAS Z-score
            
            cur.twasz = wgt.matrix[,mod.best] %*% cur.Z
            cur.twasr2pred = wgt.matrix[,mod.best] %*% cur.LD %*% wgt.matrix[,mod.best]
            
            non_zero_ind <- (wgt.matrix[,mod.best] !=0)
            
            ### standardize weights （not useful, just keep it)
            diag_element <- as.vector(wgt.matrix[non_zero_ind,mod.best])
            diag_sd <- diag_element/sum(abs(diag_element))
            weight_diag <- diag(diag_sd,nrow = length(diag_sd))
            
            Zstat.w <- weight_diag %*% cur.Z[non_zero_ind]
            corSNP.w <- weight_diag %*% cur.LD[non_zero_ind, non_zero_ind] %*% t(weight_diag)
            
            U = cur.Z
            U = as.matrix(U)
            
            V = cur.LD
            weight = wgt.matrix[,mod.best]
            
            name = rownames(V)
            rownames(U) = name
            
            PaSPU <- aSPU(U, V,weight,pow = c(1:6, Inf), n.perm = 1e3)$pvs
            
            if(min(PaSPU)< 5e-4){
                PaSPU <- aSPU(U, V,weight,pow = c(1:6, Inf), n.perm = 1e5)$pvs
            }
            
            if(min(PaSPU) < 5e-5){
                PaSPU <- aSPU(U, V,weight,pow = c(1:6, Inf), n.perm = 1e6)$pvs
            }
            
            #For simulation, comment out this time connsuimg part
            if(min(PaSPU) < 5e-6){
                PaSPU <- aSPU(U, V,weight,pow = c(1:6, Inf), n.perm = 1e7)$pvs
            }
            
            SumRes = Sum(U=Zstat.w, CovS=corSNP.w)
            pSum = SumRes$pval
            pSSU = SumSqU(U=Zstat.w, CovS=corSNP.w)
            pUminP = 1 - (1 - min(pSum,pSSU))^2
            Pout = as.numeric(c(PaSPU[length(PaSPU)], pSum, pSSU, pUminP, PaSPU,SumRes$T))
        }
    }
    Pout
}
