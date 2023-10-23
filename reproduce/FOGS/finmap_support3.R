
PatchUp <- function(M) {
  for (p in 1:ncol(M)) {
    # Get the current column vector
    vec <- M[, p]
    
    # Skip if no NA is found
    if (sum(is.na(vec)) == 0) {
      next()
    }

    # Impute by averaging
    M[is.na(vec), p] <- mean(vec, na.rm = TRUE)
  }
  
  return(M)
}


pre.process <- function(keep.ambigous=TRUE) {
    # Match summary data to input, record NA where summary data is missing
    m = match( genos$bim[,2] , sumstat$SNP )
    sum.missing = is.na(m)
    sumstat = sumstat[m,]
    sumstat$SNP = genos$bim[,2]
    sumstat$A1[ sum.missing ] = genos$bim[sum.missing,5]
    sumstat$A2[ sum.missing ] = genos$bim[sum.missing,6]
    
    # QC / allele-flip the input and output
    qc = allele.qc( sumstat$A1 , sumstat$A2 , genos$bim[,5] , genos$bim[,6] )
    
    # Flip Z-scores for mismatching alleles
    sumstat$Z[ qc$flip ] = -1 * sumstat$Z[ qc$flip ]
    sumstat$beta[ qc$flip ] = -1 * sumstat$beta[ qc$flip ]
    
    sumstat$A1[ qc$flip ] = genos$bim[qc$flip,5]
    sumstat$A2[ qc$flip ] = genos$bim[qc$flip,6]
    
    # Remove strand ambiguous SNPs (if any)
    if ( sum(!qc$keep) > 0 & keep.ambigous) {
        genos$bim = genos$bim[qc$keep,]
        genos$bed = genos$bed[,qc$keep]
        sumstat = sumstat[qc$keep,]
    }
    
    ## get the pruned SNP set
    prune.cutoff = 0.9
    
    prune.snp = SNP
    
    used.snp = wgtlist[wgtlist[,"gene"] ==wgtlist0[w,"gene"],1]
    
    prune.snp = prune.snp[!prune.snp%in% used.snp]
    prune.snp = prune.snp[prune.snp %in% sumstat[,"SNP"]]
    
    prune.snp.length.preproc = length(prune.snp)

    snp.list.name = paste(outd,"/",wgtlist0$gene[w],"_",job,"_SNPlist.txt",sep="")
    write.table(prune.snp,snp.list.name,quote = FALSE,row.names = FALSE,col.names = FALSE)
    
    system(paste0("plink --bfile ", outd,"/LDref_tmp_",job," --extract ", snp.list.name," --make-bed --out ",outd,"/",wgtlist0$gene[w],"_",job))
    system(paste0("plink --bfile ", outd,"/",wgtlist0$gene[w],"_",job," --indep-pairwise 1000 1 ",prune.cutoff," --out ",outd,"/",wgtlist0$gene[w]))
    
    prune.snp = read.table(paste(outd,"/",wgtlist0$gene[w],".prune.in",sep=""),stringsAsFactors =FALSE)
    
    prune.snp = prune.snp[,1]
    system(paste0("rm ",outd,"/",wgtlist0$gene[w],"*"))
    
    
    twas.weight.snp.set = used.snp
    # remove the SNPs with NAN, Inf, etc...
    qual.sumstat =  sumstat[sumstat[,"SNP"] %in% prune.snp, ]
    qual.sumstat = qual.sumstat[!duplicated(qual.sumstat[,"SNP"]),]
    rownames(qual.sumstat) = qual.sumstat[,"SNP"]
    qual.sumstat = qual.sumstat[prune.snp,]
    prune.snp = prune.snp[rowSums(is.na(qual.sumstat))==0]
    
    qual.sumstat =  sumstat[sumstat[,"SNP"] %in%twas.weight.snp.set, ]
    qual.sumstat = qual.sumstat[!duplicated(qual.sumstat[,"SNP"]),]
    rownames(qual.sumstat) = qual.sumstat[,"SNP"]
    qual.sumstat = qual.sumstat[twas.weight.snp.set,]
    
    twas.weight.snp.set = twas.weight.snp.set[rowSums(is.na(qual.sumstat))==0]
    
    # remove pruned SNPs set that are highly correlated with the twas.weight.snp.set
    out = list(prune.snp.length.preproc = prune.snp.length.preproc,prune.snp =prune.snp,twas.weight.snp.set = twas.weight.snp.set)
    return(out)
}
