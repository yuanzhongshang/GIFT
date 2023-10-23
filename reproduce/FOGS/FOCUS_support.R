LOG_INFO <- "INFO"
LOG_WARNING <- "WARNING"
LOG_ERROR <- "ERROR"

LOG <- function(x, level=LOG_INFO) {
    str <- paste0("[%D %H:%M:%S - ", level, "]")
    cat(format(Sys.time(), str), x, "\n")
}

log_sum_exp <- function(a, b) {
    x <- c(a, b)
    # Computes log(sum(exp(x))
    # Uses offset trick to avoid numeric overflow: http://jblevins.org/notes/log-sum-exp
    if ( max(abs(x)) > max(x) ) {
        offset <- min(x)
    } else {
        offset <- max(x)
    }
    log(sum(exp(x - offset))) + offset
}

annotate_cred_set <- function(df, prb=0.90) {
    # the call to abs  function may seem weird, but there is numerical issues with 1 - cumsum for strange reason.
    # w/o it sometimes small negative numbers appear and throw off CS computation
    df %>% group_by(BLOCK) %>%
    mutate(NPIP=PIP / sum(PIP)) %>%
    arrange(BLOCK, NPIP) %>%
    group_by(BLOCK) %>%
    mutate(IN.CRED.SET=abs((1 - cumsum(NPIP))) <= prb) %>%
    select(-NPIP)
}

get_independent <- function(CHR, ID, P0, P1, regions) {
    # Sometimes labels are the same bc of bugs...
    P1 = P0 + 1
    
    # compute overlaps
    g_ranges <- GRanges(seqnames = CHR, ranges= IRanges::IRanges(start = P0, end = P1, names = ID))
    r_ranges <- GRanges(seqnames = regions$CHR, ranges= IRanges::IRanges(start = regions$START, end = regions$STOP))
    overlaps <- findOverlaps(g_ranges, r_ranges, select="arbitrary", maxgap=1e4)
    
    # some annotations don't overlap exactly
    # just get nearest for now
    if (sum(is.na(overlaps)) > 0) {
        subset <- g_ranges[is.na(overlaps)]
        miss <- nearest(subset, r_ranges, select="arbitrary")
        overlaps[is.na(overlaps)] <- miss
    }
    
    # prettify regions and output them as 'blocks'
    pranges <- paste0(regions$CHR, ":", regions$START, "..", regions$STOP)
    pranges[overlaps]
}

get_local_params <- function(wgt.mat, cur.ID, cur.MODEL, genos) {
    # load weights into matrix after QCing...
    Wlist <- lapply(1:length(cur.ID), function(i) {
        wgt.matrix = wgt.mat[wgt.mat[,2]==cur.ID[i],]
        
        # Match up the SNPs and weights
        m = match(  wgt.matrix[,1] , genos$bim[,2] )
        m.keep = !is.na(m)
        wgt.matrix = wgt.matrix[m.keep,,drop=F]
        
        cur.genos = genos$bed[,m[m.keep]]
        cur.bim = genos$bim[m[m.keep],]
        
        # Flip WEIGHTS for mismatching alleles
        qc = allele.qc( wgt.matrix[,5] , wgt.matrix[,4] , cur.bim[,5] , cur.bim[,6] )
        wgt.matrix[qc$flip,"weight"] = -1 * wgt.matrix[qc$flip,"weight"]
        
        # Predict into reference
        mod = cur.MODEL[i]
        
        wgt <- tibble(SNP=wgt.matrix[,1], WGT=as.double(wgt.matrix[, mod]))
        colnames(wgt)[2] <- cur.ID[i]
        wgt
    })
    
    W <- purrr::reduce(Wlist, full_join, by="SNP")
    snps <- W$SNP
    W$SNP <- NULL
    W[is.na(W)] <- 0
    W <- as.matrix(W)
    
    # check if we have single gene
    if (is.null(dim(W)) && length(W) > 0) {
        W <- t(t(W))
    }
    rownames(W) <- snps
    
    # scale weights and compute LD
    m <- match(rownames(W), genos$bim[, 2])
    X <- genos$bed[,m]
    S <- apply(X %*% W, 2, sd)
    if (length(S) == 1) {
        flags <- S != 0
        if (S != 0) {
            S <- t(t(1 / S))
        }
    } else {
        # drop genes with 0 genetic covariance
        flags <- S != 0
        S <- S[flags]
        W <- W[, flags]
        S <- diag(1 / S)
    }
    SW <- W %*% S
    LD <- cor(X)
    
    return (list(SW=SW, LD=LD, FLAGS=flags))
}

fine_map <- function(cur.Z, cur.ID, cur.W, cur.LD, prb, prior_chisq, intercept=F, posterior_check=0, tol=2.220446e-14, verbose=F) {
    
    # check length
    m <- length(cur.Z)
    if (m > 1) {
        S <- t(cur.W) %*% cur.LD %*% cur.W
    } else {
        S <- t(1)
    }
    rownames(S) <- cur.ID
    colnames(S) <- cur.ID
    
    if (m > 1) {
        # this isn't the optimal solution, but _is_ unbiased.
        # optimal would require computing GLS during REML solution
        # and potentially induce large overhead
        if (intercept) {
            tmp_v <- rowSums(t(cur.W) %*% cur.LD)
            inter.est <- (t(tmp_v) %*% cur.Z) / sum(tmp_v ** 2)
            inter.s2 <- sum((cur.Z - tmp_v * inter.est)^2) / (m - 1)
            inter.se <- sqrt(inter.s2 / sum(tmp_v ** 2))
            inter.z <- inter.est / inter.se
            inter.p <- 2 * pnorm(abs(inter.z), lower.tail=F)
            mu <- tmp_v * inter.est
            Zp <- (cur.Z - mu)
        } else {
            mu <- rep(0, m)
            inter.z <- 0
            inter.p <- 1
            Zp <- cur.Z
        }
        
        pip <- rep(0, m)
        null.pip <- m * log(1 - prb)
        log.marginal <- null.pip # initialize with NULL prior-weighted BF
        
        k <- min(5, m) # this needs to be a parameter at some point
        pset <- unlist(sapply(1:k, function(x) combn(m, x, simplify=F)), recursive=F)
        for (idx_set in pset) {
            # only need genes in the causal configuration using FINEMAP BF trick
            cur.S <- S[idx_set, idx_set]
            cur.Zp <- Zp[idx_set]
            
            # compute SVD for robust estimation
            # if rank deficient, drop corresponding eigenvectors/values
            res <- svd(cur.S)
            keep <- res$d > tol
            nc <- sum(keep)
            cur.eig <- res$d[keep]
            cur.U <- res$u[,keep]
            cur.chi2 <- prior_chisq / nc
            cur.scaled.Zp <- (t(cur.Zp) %*% cur.U)^2
            
            # log BF + log prior
            cur.log.BF <- 0.5 * -sum(log(1 + cur.chi2 * cur.eig)) +
            0.5 *  sum((cur.chi2 / (1 + cur.chi2 * cur.eig)) * cur.scaled.Zp) +
            nc * log(prb) + (m - nc) * log(1 - prb)
            
            # keep track for marginal likelihood
            log.marginal <- log_sum_exp(log.marginal, cur.log.BF)
            
            # marginalize the posterior for marginal-posterior-inclusion probability (pip) on causals
            for (idx in idx_set) {
                if (pip[idx] == 0) {
                    pip[idx] <- cur.log.BF
                } else {
                    pip[idx] <- log_sum_exp(pip[idx], cur.log.BF)
                }
            }
        }
        
        # convert logpips to pips
        pip <- exp(pip - log.marginal)
        null.pip <- exp(null.pip - log.marginal)
        
        if (posterior_check > 0) {
            res <- svd(S)
            eigs <- res$d
            U <- res$u
            sim.mu <- mu
            simulation <- bind_rows(lapply(1:posterior_check,
            function(x) {
                sim.caus <- rbinom(n=m, size=1, prob=pip)
                sim.chi2 <- ifelse(sum(sim.caus) > 0, prior_chisq / sum(sim.caus), 0)
                sim.V <- S %*% diag(prior_chisq * sim.caus) %*% S + S
                sim.Zs <- rmvn(n=1, mu=sim.mu, sigma=sim.V)
                #sim.Zs <- t(rnorm(n=m, mean=sim.mu, sd=sqrt(eigs * prior_chisq * sim.caus + 1))) %*% U
                data.frame(ID=cur.ID,
                DATA.Z=cur.Z,
                CAUS=t(t(sim.caus)),
                SIM.Z=t(sim.Zs))
            }))
        } else {
            simulation <- NA
        }
        
    } else {
        marginal <- 0
        inter.z <- NA
        inter.p <- NA
        BF <- dnorm(cur.Z, mean=0, sd=sqrt(prior_chisq + 1)) / dnorm(cur.Z, mean=0, sd=1) * prb
        marginal <-  (1 - prb) + BF
        pip <- BF / marginal
        null.pip <- (1 - prb) / marginal
        log.marginal <- log(marginal)
        Zp <- cur.Z
        if (posterior_check > 0) {
            sim.caus <- rbinom(posterior_check, size=1, prob=pip)
            sim.Z <- rnorm(posterior_check, mean=0, sd=sqrt(prior_chisq * sim.caus + 1))
            simulation <- data.frame(ID=rep(cur.ID, posterior_check),
            CAUS=t(t(sim.caus)),
            DATA.Z=rep(cur.Z, posterior_check),
            SIM.Z=t(t(sim.Z)))
        } else {
            simulation <- NA
        }
    }
    
    return (list(RESID.Z=Zp, INTER.Z=inter.z, INTER.P=inter.p, PIP=pip, NULL.PIP=null.pip,
    MARG.LOG.LIKE=log.marginal, SIM=simulation))
}



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
    snp[["keep"]] = !((a1=="A" & a2=="T") | (a1=="T" & a2=="A") | (a1=="C" & a2=="G") | (a1=="G" & a2=="C"))
    snp[["flip"]] = (a1 == ref2 & a2 == ref1) | (a1 == flip2 & a2 == flip1)
    
    return(snp)
}

FOCUS <- function(TWAS.P,TWAS.Z,wgt.mat) {
    
    TWAS.P = out.res[, "SUM"]
    TWAS.Z = out.res[, "ZSum"]
    wgtlist$TWAS.P = TWAS.P
    wgtlist$TWAS.Z = TWAS.Z
    genos = genos2 #read_plink(paste0(outd,"/LDref_tmp_",job),impute="avg")
    genos$bed = scale(genos$bed)
    
    
    # Revise here to get the results..
    
    wgt.file = paste(opt$weights_dir,"/",wgtlist$WGT,sep='')
    cur.FILE <- wgt.file
    
    
    cur.MODEL <- rep(3,dim(wgtlist)[1])
    cur.ID <- wgtlist$ID
    cur.CHR <- wgtlist$CHR
    cur.P0 <- wgtlist$P0
    cur.P1 <- wgtlist$P1
    cur.Z <- wgtlist$TWAS.Z
    cur.P <- wgtlist$TWAS.P
    
    
    
    params <- get_local_params(wgt.mat, cur.ID, cur.MODEL, genos)
    cur.LD <- params$LD
    cur.SW <- params$SW
    cur.keep <- params$FLAGS
    
    # ffs just make a data.frame already...
    #cur.FILE <- cur.FILE[cur.keep] #R.data
    cur.MODEL <- cur.MODEL[cur.keep]
    cur.ID <- cur.ID[cur.keep]
    cur.CHR <- cur.CHR[cur.keep]
    cur.P0 <- cur.P0[cur.keep]
    cur.P1 <- cur.P1[cur.keep]
    cur.Z <- cur.Z[cur.keep]
    cur.P <- cur.P[cur.keep]
    
    
    posterior_check = 0
    use_intercept = F
    res <- fine_map(cur.Z,cur.ID,cur.SW,cur.LD,prb,prior_chisq,intercept=use_intercept,posterior_check=0,verbose=T,tol=tol)
    
    
    cur.INTERCEPT.Z <- res$INTER.Z
    cur.INTERCEPT.P <- res$INTER.P
    cur.RESID.Z <- res$RESID.Z
    cur.PIP <- res$PIP
    cur.NULL.PIP <- res$NULL.PIP
    
    
    tbl <-  data.frame (
                        ID=c(paste0("NULL.", i), paste0("INTERCEPT.", i), as.character(cur.ID)),
                        CHR=chr,
                        P0=c(0, 0, cur.P0),
                        P1=c(0, 0, cur.P1),
                        Z=c(0, cur.INTERCEPT.Z, cur.Z),
                        P=c(1, cur.INTERCEPT.P, cur.P),
                        RESID.Z=c(0, 0, cur.RESID.Z),
                        PIP=c(cur.NULL.PIP, 0, cur.PIP))
    return(tbl)
}
