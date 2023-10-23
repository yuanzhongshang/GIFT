# check if the required packages have been installed
#list.of.packages <- c("data.table","optparse","Rcpp","RcppArmadillo","mvtnorm","BEDMatrix","bigmemory","dplyr","mvnfast")
#new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
#if(length(new.packages)) install.packages(new.packages)

# To run FOCUS; make sure install GenomicRanges by BiocManager::install("GenomicRanges")

#setwd("/gpfs/research/chongwu/Chong/MWAS/Finemapping")

library(Rcpp)
library(RcppArmadillo)
library(bigmemory)

library(mvtnorm)
sourceCpp("GMaSPU_support.cpp")
source("aSPU.R")
source("dist_support2.R")
source("JointRidge.R")
source("SMI.R")
source("FOCUS_support.R")
source("finmap_support3.R") # cutoff 0.9

suppressMessages(library("BEDMatrix"))
suppressMessages(library("optparse"))
library(data.table)

suppressMessages(library("dplyr"))
suppressMessages(library("GenomicRanges"))
suppressMessages(library("mvnfast"))


conditional.method <- "ridgeCOJO"
fine.model <- "reduced"

#########################################
## Inputs                             ###
#########################################
#loci.indx = 14
#chr.id = 21

option_list = list(
make_option("--refld", action="store", default=NA, type='character',
help="LD reference (e.g., 1000 Genomes) path [required]"),
make_option("--outd", action="store", default=NA, type='character',
help="The temporary directory for intermediate results [required]"),
make_option("--loci", action="store", default=NA, type='character',
help="The directory for the LDetect defined loci [required]"),
make_option("--weights", action="store", default=NA, type='character',
help="The eQTL-derived weights file  [required]"),
make_option("--genelist", action="store", default=NA, type='character',
help="The corresponding gene list  [required]"),
make_option("--sumstat", action="store", default=NA, type='character',
help="summary statistics (txt format) [required]"),
make_option("--saveprefix", action="store", default=NA, type='character',
help="The prefix for the output [required]"),
make_option("--chr_id", action="store", default=-1, type='character',
help="The chromosome ID. We recommend parallel the computations by chromosomes and locus id [required]"),
make_option("--locus_id", action="store", default=NA, type='character',
help="The locus ID [required]")
)

opt = parse_args(OptionParser(option_list=option_list))

#--refld /gpfs/research/chongwu/shared/1000Genomes/1000G.EUR.ALLSNP.QC.CHR --outd /gpfs/research/chongwu/Chong/MWAS/Finemapping/Blood --loci /gpfs/research/chongwu/shared/LDetect_LD_regions/EUR/ --weights /gpfs/research/chongwu/Chong/Application/COVID19/Finemapping/Whole_Blood_weights.txt --genelist /gpfs/research/chongwu/Chong/Application/COVID19/Finemapping/Whole_Blood_UT_gene_list.txt --sumstat /gpfs/research/chongwu/Chong/Application/COVID19/Finemapping/processed_data.txt --saveprefix Blood_UT --chr_id 3 --locus_id 32

#opt = list(
#refld = "/gpfs/research/chongwu/shared/1000Genomes/1000G.EUR.ALLSNP.QC.CHR",
#outd = "/gpfs/research/chongwu/Chong/MWAS/Finemapping/Blood",
#loci = "/gpfs/research/chongwu/shared/LDetect_LD_regions/EUR/",
#weights = "/gpfs/research/chongwu/Chong/Application/COVID19/Finemapping/Blood_JTI_weights.txt",
#genelist = "/gpfs/research/chongwu/Chong/Application/COVID19/Finemapping/Blood_JTI_gene_list.txt",
#sumstat = "/gpfs/research/chongwu/Chong/Application/COVID19/Finemapping/processed_data.txt",
#saveprefix = "test",
#chr_id = 3,
#locus_id = 32
#)

loci.indx = as.numeric(opt$locus_id)
chr.id = as.numeric(opt$chr_id)
job = paste("CHR",chr.id,"Loci",loci.indx,sep="")



system(paste("mkdir -p ",opt$outd,sep=""))

loci.dir = opt$loci
outd = opt$outd


wgtlist = as.data.frame(fread(opt$weights))
#wgtlist
#rsid            gene      weight ref_allele eff_allele
#1 rs141364387 ENSG00000169583 -0.19646764          C          T


wgtlist0 = as.data.frame(fread(opt$genelist))
#wgtlist0
#       left     right  chr            gene
#1 169818772 169863408 chr1 ENSG00000000457
#2 169631245 169823221 chr1 ENSG00000000460

wgtlist0[,3] = as.numeric(gsub("chr","",wgtlist0[,3]))
wgtlist0 = wgtlist0[!is.na(wgtlist0[,1]),]


# load summary statistics
sumstat.orgin = as.data.frame(fread(opt$sumstat))
# sumstat.orgin[1:3,]
#  CHR    POS         SNP A1 A2     beta      se      N          Z
#1   1 528642  rs76388980  G  A  0.28847 0.72139 900687  0.3998808
#3   1 565987 rs564223368  C  T -1.40180 1.73360 900687 -0.8086064


cat("Load the data\n")

#34602206 34637980
##############################################################
# Step 1: find the prunned SNP for each gene in a locus    ###
##############################################################
tmp.loci = fread(paste(loci.dir, "fourier_ls-chr",chr.id,".bed",sep=""))
tmp.loci = as.data.frame(tmp.loci)

tmp.loci[, 1] = as.numeric(gsub("chr", "", tmp.loci[, 1]))

risk.region.used = tmp.loci[loci.indx,,drop=F]
wgtlist0 = wgtlist0[wgtlist0[, 3] == risk.region.used[1, 1],]
wgtlist0 = wgtlist0[(wgtlist0[,1] < risk.region.used[1, 3] & wgtlist0[,1] > risk.region.used[1, 2]),]

colnames(wgtlist0) = c("P0","P1","CHR","gene")

wgtlist = wgtlist[wgtlist[,"gene"] %in% wgtlist0[,4],]
wgt.mat = wgtlist

SNP = unique(wgtlist[,1])


start.p0 = min(wgtlist0$P0) - 500 * 1000

if (start.p0 < 0) {
    start.p0 = 0
}
start.p1 = max(wgtlist0$P1) + 500 * 1000
start.chr = wgtlist0$CHR[1]

chr.id = risk.region.used[1, 1]
system(paste0("plink --bfile ", opt$refld, " --chr ", chr.id, " --from-kb ", start.p0 / 1e3, " --to-kb ", start.p1 / 1e3, " --maf 0.01 --make-bed --out ", opt$outd, "/LDref_tmp_", job))


# Load in reference data
# make the genos format as the one by plink2R
#genos = read_plink(paste0(opt$outd, "/LDref_tmp_", job), impute = "avg")
geno = BEDMatrix(paste0(opt$outd, "/LDref_tmp_", job))
geno <- as.matrix(geno)

geno = PatchUp(geno)

snp.name = gsub("_.*","",colnames(geno))
colnames(geno) = snp.name

bim <- fread(paste0(opt$outd, "/LDref_tmp_", job,".bim"))
bim <- as.data.frame(bim)

fam <- fread(paste0(opt$outd, "/LDref_tmp_", job,".fam"))
fam <- as.data.frame(fam)

genos = list(bed = geno,fam = fam, bim = bim)

tmp.bim = genos$bim
snp.inf = tmp.bim[, 2]

sumstat.orgin = sumstat.orgin[sumstat.orgin[,1]==chr.id,]

m = match(snp.inf, sumstat.orgin$SNP)
## For each wgt file:
sumstat.orgin = sumstat.orgin[m,]


genos2 = genos

out.res <- as.data.frame(matrix(NA, nrow(wgtlist0), 38))

if (dim(wgtlist0)[1] < 2) {
    out.res[1,] <- -1
    out.file <- paste(outd, "/out_", loci.indx, ".rds", sep = "")
    saveRDS(out.res, out.file)
}

outd = opt$outd


for (w in 1:nrow(wgtlist0)) {
    #out.fun <- function(w) {
    tryCatch({
        sumstat = sumstat.orgin
        genos = genos2 #read_plink(paste0(outd,"/LDref_tmp_",job),impute="avg")
        
        save.name = paste(outd, "/snp_", loci.indx, "_", w, ".rds", sep = "")
        out = pre.process(keep.ambigous=FALSE)
        
        prune.snp.length.preproc = out$prune.snp.length.preproc
        prune.snp = out$prune.snp
        twas.weight.snp.set = out$twas.weight.snp.set
        
        #if (length(prune.snp) >= 2) {
            saveRDS(out, save.name)
        #}
        cat(length(prune.snp),"\n")
    }, error = function(e) { cat("ERROR :", conditionMessage(e), "\n") })
}


##########################################
## Step 2: Run FOGS                    ###
##########################################
keep.ambigous = FALSE
causal.gene.id <- 99 # dummy variable to make sure the code run smoothly
colnames(out.res) <- c("CHR", "ID", "P0", "P1", "n.SNP", "n.condSNP", "casual.id", "afSPU", "fSUM", "fSSU", "afSSU", "fSPU(1)", "fSPU(2)", "fSPU(3)", "fSPU(4)", "fSPU(5)", "fSPU(6)", "fSPU(Inf)", "afSPU", "ZfSum", "aSPU", "SUM", "SSU", "aSSU", "SPU(1)", "SPU(2)", "SPU(3)", "SPU(4)", "SPU(5)", "SPU(6)", "SPU(Inf)", "aSPU", "ZSum", "Focus-Sum", "Focus-SSU", "Focus-aSPU", "Focus-aSSU", "runing_time")

for (w in 1:nrow(wgtlist0)) {
    # out.fun <- function(w) {
    tryCatch({
        sumstat <- sumstat.orgin
        start.time = proc.time()[3]
        genos <- genos2 # read_plink(paste0(outd,"/LDref_tmp_",job),impute="avg")
        
        save.name <- paste(outd, "/snp_", loci.indx, "_", w, ".rds", sep = "")
        
        out <- readRDS(save.name)
        
        prune.snp.length.preproc <- out$prune.snp.length.preproc
        prune.snp <- out$prune.snp
        twas.weight.snp.set <- out$twas.weight.snp.set
        
        
        # quality control
        m <- match(genos$bim[, 2], sumstat$SNP)
        sum.missing <- is.na(m)
        sumstat <- sumstat[m,]
        sumstat$SNP <- genos$bim[, 2]
        sumstat$A1[sum.missing] <- genos$bim[sum.missing, 5]
        sumstat$A2[sum.missing] <- genos$bim[sum.missing, 6]
        
        # QC / allele-flip the input and output
        qc <- allele.qc(sumstat$A2, sumstat$A1, genos$bim[, 5], genos$bim[, 6])
        
        # Flip Z-scores for mismatching alleles
        sumstat$Z[qc$flip] <- -1 * sumstat$Z[qc$flip]
        sumstat$beta[qc$flip] <- -1 * sumstat$beta[qc$flip]
        
        sumstat$A1[qc$flip] <- genos$bim[qc$flip, 5]
        sumstat$A2[qc$flip] <- genos$bim[qc$flip, 6]
        
        # Remove strand ambiguous SNPs (if any)
        if (sum(!qc$keep) > 0 & keep.ambigous) {
            genos$bim <- genos$bim[qc$keep,]
            genos$bed <- genos$bed[, qc$keep]
            sumstat <- sumstat[qc$keep,]
        }
        
        # make sure prune SNPs is in the summary data and corresponding reference
        cur.genos <- genos$bed
        prune.snp <- prune.snp[prune.snp %in% colnames(cur.genos)]
        twas.weight.snp.set <- twas.weight.snp.set[twas.weight.snp.set %in% colnames(cur.genos)]
        
        prune.snp <- prune.snp[prune.snp %in% sumstat$SNP]
        twas.weight.snp.set <- twas.weight.snp.set[twas.weight.snp.set %in% sumstat$SNP]
        
        # make sure prune SNPs is in the summary data and corresponding reference
        cur.genos <- genos$bed
        prune.snp <- prune.snp[prune.snp %in% colnames(cur.genos)]
        twas.weight.snp.set <- twas.weight.snp.set[twas.weight.snp.set %in% colnames(cur.genos)]
        
        
        # calcualte the conditional score for each SNPs
        sumstat.final <- NULL
        for (i in 1:length(twas.weight.snp.set)) {
            snp.used <- c(twas.weight.snp.set[i], prune.snp)
            cur.genos <- genos$bed
            cur.genos <- cur.genos[, snp.used]
            
            sumstat.tmp2 <- sumstat[sumstat[, "SNP"] %in% snp.used,]
            
            rownames(sumstat.tmp2) <- sumstat.tmp2[, "SNP"]
            sumstat.tmp2 <- sumstat.tmp2[snp.used,]
            
            B <- sumstat.tmp2[, "beta"]
            S <- sumstat.tmp2[, "se"]
            N <- sumstat.tmp2[, "N"]
            
            cur.genos <- cur.genos - matrix(rep(colMeans(cur.genos), each = dim(cur.genos)[1]), dim(cur.genos)[1], dim(cur.genos)[2])
            
            XX <- cov(cur.genos)
            
            if (conditional.method == "regCOJO") {
                XX <- XX + diag(0.1, dim(XX)[1], dim(XX)[2])
                tmp <- JointSum(B, S, N, XX)
            }
            
            if (conditional.method == "ridgeCOJO") {
                lambda <- 0.1
                tmp <- JointRidge(B, S, N, XX, lambda)
            }
            
            beta.tmp <- tmp$beta
            se.tmp <- sqrt(diag(tmp$cov))
            z.tmp <- beta.tmp / se.tmp
            
            sumstat.tmp2$cond.z <- z.tmp
            sumstat.final <- rbind(sumstat.final, sumstat.tmp2[1,])
        }
        
        # sumstat.final2 = sumstat.final
        sumstat.tmp <- sumstat.final
        
        tmp.index <- genos$bim[, 2] %in% sumstat.tmp[, "SNP"]
        genos$bim <- genos$bim[tmp.index,,drop=F]
        genos$bed <- genos$bed[, tmp.index,drop=F]
        
        wgt.matrix = wgtlist[wgtlist[,"gene"] ==wgtlist0[w,"gene"],]
        
        # Match up the SNPs and weights
        m <- match(wgt.matrix[, 1], genos$bim[, 2])
        m.keep <- !is.na(m)
        wgt.matrix <- wgt.matrix[m.keep,,drop=F]
        
        cur.genos <- scale(genos$bed[, m[m.keep]])
        cur.bim <- genos$bim[m[m.keep],]
        # Flip WEIGHTS for mismatching alleles
        qc <- allele.qc(wgt.matrix[, 5], wgt.matrix[, 4], cur.bim[, 5], cur.bim[, 6])
        wgt.matrix[qc$flip,"weight"] <- -1 * wgt.matrix[qc$flip,"weight"]
        
        
        cur.FAIL <- FALSE
        
        # Match up the SNPs and the summary stats
        m <- match(cur.bim[, 2], sumstat.tmp$SNP)
        
        cur.Z <- sumstat.tmp$cond.z[m]
        Z <- sumstat.tmp$Z[m]
        
        mod.best = 3
        if(length(cur.Z)==1) {
            tmp.res <- c(wgtlist0$CHR[w], as.character(wgtlist0$gene[w]), wgtlist0$P0[w], wgtlist0$P1[w], length(cur.Z), length(prune.snp), causal.gene.id)
            out.res[w, 1:length(tmp.res)] <- tmp.res
            out.res[w, c("afSPU","ZfSum")] <- c(pnorm(abs(cur.Z),lower.tail=F)*2, cur.Z)
            out.res[w, c("SUM","ZSum")] <- c(pnorm(abs(Z),lower.tail=F)*2, Z)
                        
        } else {
            fine.p <- calc.pvalue(cur.Z)
            org.p <- calc.pvalue(Z)
             
            tmp.res <- c(wgtlist0$CHR[w], as.character(wgtlist0$gene[w]), wgtlist0$P0[w], wgtlist0$P1[w], length(cur.Z), length(prune.snp), causal.gene.id)
            out.res[w, 1:length(tmp.res)] <- tmp.res
            out.res[w, (length(tmp.res) + 1):(length(tmp.res) + length(fine.p))] <- fine.p
            out.res[w, (length(tmp.res) + length(fine.p) + 1):(length(tmp.res) + length(fine.p) + length(org.p))] <- org.p
            out.res[w, 38] = proc.time()[3] - start.time
             
        }
 
        rm(genos)
        rm(cur.genos)
        cat(w)
        # return(out.res)
    }, error = function(e) {
        cat("ERROR :", conditionMessage(e), "\n")
    })
}


# Run FOCUS with default settings
max_iter <- 10

prior_chisq <- 40
prb <- 1e-3
tol <- 2.220446e-14

wgt.matrix = wgtlist

wgtlist <- wgtlist0
tmp.na <- is.na(out.res[, 33])
causal.gene.id <- causal.gene.id - sum(tmp.na[1:(causal.gene.id - 1)])

out.res <- out.res[!is.na(out.res[, 33]),]
out.res[, "casual.id"] <- causal.gene.id
# cascual id should be channged

chr <- unique(out.res[, 1])
chr <- as.numeric(chr[1])

rownames(out.res) <- paste(out.res[, 2], out.res[, 3], out.res[, 4], sep = "-")
rownames(wgtlist) <- paste(wgtlist[, 4], wgtlist[, 1], wgtlist[, 2], sep = "-")

colnames(wgtlist) = c("P0","P1","CHR","ID")

wgtlist <- wgtlist[rownames(wgtlist) %in% rownames(out.res),]

rownames(wgtlist) = 1:dim(wgtlist)[1]
tryCatch({
    tmp <- FOCUS(out.res[, "SUM"], out.res[, "ZSum"],wgt.mat)
    rownames(tmp) <- paste(tmp[, "ID"], tmp[, "P0"], tmp[, "P1"], sep = "-")
    tmp <- tmp[rownames(out.res),]
    out.res[, "Focus-Sum"] <- tmp[, "PIP"]
    
}, error = function(e) {
    cat("ERROR :", conditionMessage(e), "\n")
})

out.res = out.res[,c("CHR", "ID", "P0", "P1", "n.SNP", "n.condSNP",  "afSPU",  "SUM",  "Focus-Sum", "runing_time")]
colnames(out.res) = c("CHR", "ID", "P0", "P1", "n.SNP", "n.condSNP",  "FOGS-aSPU",   "TWAS",  "Focus", "Runtime(s)")


out.file <- paste(opt$saveprefix,"CHR_",chr.id,"_Locus", loci.indx, ".txt", sep = "")
write.table(out.res, out.file,quote=F,row.names=F,col.names=T)
