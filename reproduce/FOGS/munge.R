# check if the required packages have been installed
#list.of.packages <- c("data.table","optparse","Rcpp","RcppArmadillo","mvtnorm","BEDMatrix","bigmemory","dplyr","mvnfast")
#new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
#if(length(new.packages)) install.packages(new.packages)

# To run FOCUS; make sure install GenomicRanges by BiocManager::install("GenomicRanges")

library(data.table)


munge_sumstat <- function(sumstats,sampleN=NULL) {
    
    #sumstats = "/gpfs/research/chongwu/shared/summary_statistics/COVID19/release5/processed/ANA_B2_eur_V5.txt"

    
    if (grepl(".txt",sumstats)) {
        raw = fread(sumstats)
        raw = as.data.frame(raw)
    } else {
        raw = fread(sumstats)
        raw = as.data.frame(raw)
        warning("We try to read GWAS summary data other txt format. The data may not read properly")
        
        cat("Following are the first few lines of the dataset:", sep = "\n")
        print(head(raw))
    }

    header.inner = colnames(raw)
    # Initialize the header
    header.inner <- tolower(header.inner)

    # SNP
    try.snp <- c("snp", "markername", "snpid", "rs", "rsid", "rs_number", "snps")
    header.inner[header.inner %in% try.snp] <- "SNP"

    # A1
    try.a1 <- c("a1", "allele1", "allele_1", "effect_allele", "reference_allele", "inc_allele", "ea", "ref", "a1lele1", "al1ele1")
    header.inner[header.inner %in% try.a1] <- "A1"

    # A2
    try.a2 <- c("a2", "allele2", "allele_2", "other_allele", "non_effect_allele", "dec_allele", "nea", "alt", "a0")
    header.inner[header.inner %in% try.a2] <- "A2"

    # Z-score
    try.z <- c("zscore", "z-score", "gc_zscore", "z")
    header.inner[header.inner %in% try.z] <- "Z"


    try.chromosome <- c("chrom", "ch", "chr", "chromosome","#chr")
    header.inner[header.inner %in% try.chromosome] <- "CHR"

    # P
    try.p <- c("pvalue", "p_value", "pval", "p_val", "gc_pvalue", "p")
    header.inner[header.inner %in% try.p] <- "P"


    # Beta
    try.beta <- c("b", "beta", "effects", "effect","all_inv_var_meta_beta")
    header.inner[header.inner %in% try.beta] <- "BETA"

    # Odds ratio
    try.or <- c("or")
    header.inner[header.inner %in% try.or] <- "ODDS_RATIO"

    # Log odds
    try.logodds <- c("log_odds", "logor", "log_or")
    header.inner[header.inner %in% try.logodds] <- "LOG_ODDS"

    # Standard error
    try.se <- c("se", "sebeta", "beta_se","all_inv_var_meta_sebeta")
    header.inner[header.inner %in% try.se] <- "SE"

    # sample size
    try.samplesize <- c("samplesize", "N","all_meta_sample_N")
    
    if(sum(header.inner %in% try.samplesize)==1) {
        header.inner[header.inner %in% try.samplesize] <- "N"
    } else {
        raw$N = sampleN
    }

    colnames(raw) <- header.inner


    list.coerce <- c("Z", "BETA", "ODDS_RATIO", "LOG_ODDS", "SE","N")

    for (i in 1:length(header.inner)) {
        if (header.inner[i] %in% list.coerce) {
            if (class(raw[, header.inner[i]]) != "numeric") {
                class(raw[, header.inner[i]]) <- "numeric"
                cat(paste0("Column ", header.inner[i], " has wrong class and has been coerced to numeric."), sep = "\n")
                cat("=============================================================================================================", sep = "\n")
            }
        }
    }

    # Missing z-score?
    calculate.z <- FALSE
    if (!("Z" %in% header.inner)) {
        warning("No Z score column, we calculate one based on the information we have")
        
        if ("BETA" %in% header.inner & "SE" %in% header.inner) {
            raw['Z'] <- raw$BETA / raw$SE
            calculate.z <- TRUE
        } else if ("ODDS_RATIO" %in% header.inner & "SE" %in% header.inner) {
            raw['Z'] <- log(raw$ODDS_RATIO) / raw$SE
            calculate.z <- TRUE
        } else if ("LOG_ODDS" %in% header.inner & "SE" %in% header.inner) {
            raw['Z'] <- raw$LOG_ODDS / raw$SE
            calculate.z <- TRUE
        } else if ("BETA" %in% header.inner & "P" %in% header.inner) {
            raw['Z'] <- sign(raw$BETA) * abs(qnorm(raw$P / 2))
            calculate.z <- TRUE
        } else if ("ODDS_RATIO" %in% header.inner & "P" %in% header.inner) {
            raw['Z'] <- sign(log(raw$ODDS_RATIO)) * abs(qnorm(raw$P / 2))
            calculate.z <- TRUE
        } else if ("LOG_ODDS" %in% header.inner & "P" %in% header.inner) {
            raw['Z'] <- sign(raw$ODDS_RATIO) * abs(qnorm(raw$P / 2))
            calculate.z <- TRUE
        } else {
            stop("I can't calculate z-score based on the information I have. SAD FACE EMOJI.", sep = "\n")
        }
        
        if (sum(is.na(raw$Z)) != 0) {
            n.start <- nrow(raw)
            
            raw <- raw[!is.na(raw$Z),]
            n.end <- nrow(raw)
            cat(paste0(n.start - n.end, " rows removed for having invalid z-score!"), sep = "\n")
        }
        cat("=============================================================================================================", sep = "\n")
    }

    if(sum(header.inner %in% c("SNP","A1","A2","CHR")) !=4) {
        stop("We tried our best to match the colnames of the summary data. Please revise the colnames to the following format. You can also report in GitHub to make the lists more comphrensive:\n SNP: snp, markername, snpid, rs, rsid, rs_number, snps\n A1: a1, allele1, allele_1, effect_allele, reference_allele, inc_allele, ea, ref, a1lele1, al1ele1\n A2: a2, allele2, allele_2, other_allele, non_effect_allele, dec_allele, nea, alt, a0\n Z: zscore, z-score, gc_zscore, z\n CHR: chrom, ch, chr, chromosome")
    }
    
    col.name = colnames(raw)
    col.name[col.name=="pos"] = "POS"
    col.name[col.name=="BETA"] = "beta"
    
    col.name[col.name=="SE"] = "se"
    
    colnames(raw) = col.name
    
    raw = raw[,c("CHR","POS","SNP","A1","A2","beta","se","N","Z")]
    
    # remove the deletion or insertion
    indx = nchar(raw$A1)==1 & nchar(raw$A2) ==1
    
    raw = raw[indx,]
    
    # remove duplicated SNPs
    raw = raw[!duplicated(raw$SNP),]
    
    # return the processed data
    return(raw)
}


# use Zichen's code to prepare this;
# also remove the ambinguous SNPs and nchar !=1

####################################
# Example codes                 ####
####################################
#sumstats = "/gpfs/research/chongwu/shared/summary_statistics/COVID19/release5/processed/ANA_B2_eur_V5.txt"
#data = munge_sumstat(sumstats)
#write.table(data,"processed_data.txt",col.names=TRUE,row.names=FALSE,quote=FALSE)
####################################
####################################

#3. sumstat.orgin[1:3,]
#  CHR    POS         SNP A1 A2     beta      se      N          Z
#1   1 528642  rs76388980  G  A  0.28847 0.72139 900687  0.3998808
#3   1 565987 rs564223368  C  T -1.40180 1.73360 900687 -0.8086064

