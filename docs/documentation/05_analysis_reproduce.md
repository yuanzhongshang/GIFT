---
layout: page
title: Analysis Reproduce
description: ~
---

## Generate the simulation data
We randomly selected 50 regions from 1,533 regions. The 1,533 region information is [here](https://github.com/yuanzhongshang/GIFT/blob/main/reproduce/LDetectregion1533.txt). The full information about the region can be found in [here](https://github.com/yuanzhongshang/GIFT/blob/main/reproduce/LDetect). For each region in turn, we performed 20 simulation replicates, resulting in a total of 1,000 simulation replicates per setting.
For each region, we conducted the simulations based on the realistic genotypes from GEUVADIS (n1=465) and UK Biobank (n2=5,000). Take a region on chr 4 for example. This region includes seven genes: C4orf29, C4orf33, LARP1B, PGRMC2, PHF17, RP11-420A23.1, and SCLT1. We set C4orf29 as the causal gene with the effect size being sqrt(0.1).

### Generate the individual-level data
```r
library("BEDMatrix")
library(mvtnorm)
library(GIFT)

##load the genotypes from GEUVADIS

###load the directory containing the files to be processed only (e.g., plink binary format)
dir <- getwd()
filelocation <- paste0(dir,"/reproduce/simulation_data_generate/GEUVADIS")

###load the directory of plink exe file
plinkexe <- "plink"

###pre-process the file to be a list including gene names vector, cis-genotype matrix and pindex
convert <- pre_process_individual(filelocation, plinkexe)
gene <- convert$gene
Zx <- convert$Z
pindex <- convert$pindex

###please see https://yuanzhongshang.github.io/GIFT/documentation/04_GIFT_Example.html for the details of pre-poccess of genotype data.  

##load the genotypes from UK Biobank
filelocation <- paste0(dir,"/reproduce/simulation_data_generate/UKB")

###pre-process the file to be a list including gene names vector, cis-genotype matrix and pindex
convert <- pre_process_individual(filelocation, plinkexe)
Zy <- convert$Z

##set the correlation of the residual of gene expressions

###We can set the correlation based on the real data or some specific correlation structure 
###Option 1: the correlation below is realistic from GEUVADIS
R <- as.matrix(read.table(paste0(dir,"/reproduce/simulation_data_generate/R.txt")))

###Option 2: generate the correlation from an exponential covariance structure,e.g.ρ=0.9
AR = function(rho, p) {
  outer(1:p, 1:p, FUN = function(x, y) rho^(abs(x - y))) 
}
R <- AR(0.9, length(gene))

##generate the gene expression data 

###set the heritability of gene (PVEzx),here we fixed heritability to be 0.1 for each gene.
PVEzx <- rep(0.1, length(gene))

###generate the effect size from the normal distribution
b=list()
for(i in 1:length(gene)){
  set.seed(i)
  b[[i]] <- matrix(rnorm(pindex[i], 0, sd = sqrt(PVEzx[i]/pindex[i])), pindex[i], 1)
}

###generate the residuals from the multi-normal distribution
set.seed(123)
E <- rmvnorm(dim(Zx)[1], mean = rep(0, length(gene)), sigma = (sqrt(1-PVEzx) %*% t(sqrt(1-PVEzx)))*R)

###sum
pindexsum <- c(0,cumsum(pindex))
X=matrix(0,dim(Zx)[1],length(gene))
for(i in 1:length(gene)){
  X[,i] <- Zx[,(pindexsum[i]+1):pindexsum[i+1]] %*% b[[i]] + E[,i]
}

###set the heritability (PVEzy),here we set RASA1 as the causal gene with the effect size being sqrt(0.1).
PVEzy <- c(0.01,0,0,0,0,0,0)
casual_effect <- sqrt(PVEzy/PVEzx)

###generate the phenotype
set.seed(456)
if (sum(PVEzy) == 0){
  Y <- rnorm(dim(Zy)[2], 0, sqrt(1)) 
} else {
  lp_z=0
  for(i in 1:length(gene)){
    lp_z <- lp_z + casual_effect[i]* as.matrix(Zy[,(pindexsum[i]+1):pindexsum[i+1]]) %*% b[[i]]
  }
  Y <- lp_z + rnorm(dim(Zy)[1], 0, sqrt(1 - sum(PVEzy)))
}

X=scale(X)
Y=scale(Y)
Zx=scale(Zx)
Zy=scale(Zy)

###save these variables
setwd(dir)
save(X, Y, Zx, Zy, gene, pindex, file = "./reproduce/simulation_data_generate/data_generate_individual.RData")
```

### Convert the individual-level data into the summary statistics
```r
###calculate the z-score from GEUVADIS
n1=dim(Zx)[1]
Zscore1 =NULL
for(i in 1:length(gene)){
  Zscore1 <- cbind(Zscore1,(1/sqrt(n1-1))*t(Zx) %*% X[,i])
}

###calculate the LD matrix from GEUVADIS
LDmatrix1 <- (1/(n1-1))*t(Zx) %*% Zx

###calculate the z-score from UK Biobank
n2=dim(Zy)[1]

###calculate the LD matrix from UK Biobank
LDmatrix2 <- (1/(n2-1))*t(Zy) %*% Zy
Zscore2 <- (1/sqrt(n2-1))* (t(Zy) %*% Y)

###save these variables
save(Zscore1, LDmatrix1, Zscore2, LDmatrix2, R, n1, n2, gene, pindex, file = "./reproduce/simulation_data_generate/data_generate_summary.RData")
```

## Run GIFT
```r
library(GIFT)
```
### Using individual-level data as input 
```r
##load the simulation data or the real data with the specific format from individual-level data.
###Here we used the simulation data above.
load("./reproduce/simulation_data_generate/data_generate_individual.RData")

##run GIFT
result <- GIFT_individual(X, Y, Zx, Zy, gene, pindex, maxiter=1000, tol=1e-4, pleio=0, ncores=1, filter=F)

result
         gene   causal_effect            p
1     C4orf29    0.3940422901 0.0001300531
2     C4orf33   -0.0274405121 0.7029942822
3      LARP1B   -0.1511359367 0.2229299816
4      PGRMC2    0.1500225583 0.0690925279
5       PHF17   -0.0574943210 0.1795472908
6 RP11-420A23.1  0.0756822259 0.3894757754
7       SCLT1   -0.0004761227 0.9933530763
```

### Using summary statistics as input
```r
##load the simulation data or the real data with the specific format from summary statistics.
###Here we used the simulation data above.
load("./reproduce/simulation_data_generate/data_generate_summary.RData")

##run GIFT
result <- GIFT_summary(Zscore1, Zscore2, LDmatrix1, LDmatrix2, n1, n2, gene, pindex, R=R, maxiter=1000, tol=1e-4, pleio=0, ncores=1, in_sample_LD=T, filter=F)

result
         gene    causal_effect           p
1     C4orf29     0.411634628 9.111960e-06
2     C4orf33     0.007711403 7.982367e-01
3      LARP1B    -0.107352053 2.677252e-01
4      PGRMC2     0.135373699 4.580617e-02
5       PHF17    -0.046618238 1.993550e-01
6 RP11-420A23.1   0.041653980 4.335418e-01
7       SCLT1    -0.013595634 7.343946e-01
```

## Run FOCUS
FOCUS takes GWAS summary statistics, reference LD, eQTL weight database as inputs. Here we used the same data as above. If you do not have the in-sample LD, you may use the LD from 1000 Genomes reference: [https://alkesgroup.broadinstitute.org/FUSION/LDREF.tar.bz2](https://alkesgroup.broadinstitute.org/FUSION/LDREF.tar.bz2).

```r
library("data.table")
library("gtools")
library(RSQLite)

##load the simulation data or the real data with the specific format from individual-level data.
###Here we used the simulation data.
dir=getwd()
load("./reproduce/simulation_data_generate/data_generate_individual.RData")

##If you need compute own weights, the following additional steps are required.

###integrate the gene expression into the .fam file 
setwd("./reproduce/FOCUS")
system("mkdir weight")
setwd(gsub("/reproduce/FOCUS", "", getwd()))
###FUSION.weights.R is from TWAS/Fusion website:http://gusevlab.org/projects/fusion/
for(i in 1:length(gene)){
  fam <- read.table(paste0("./reproduce/simulation_data_generate/GEUVADIS/",gene[i],".fam"))
  fam <- cbind(fam[,1:5],X[,i])
  write.table(fam,paste0("./reproduce/simulation_data_generate/GEUVADIS/",gene[i],".fam"), quote=F,sep=" ",col.names=F,row.names=F)
  system(paste0("Rscript FUSION.weights.R --bfile ./reproduce/simulation_data_generate/GEUVADIS/",gene[i]," --tmp ./reproduce/FOCUS/tmp --out ./reproduce/FOCUS/weight/",gene[i]," --verbose 0 --models bslmm,lasso,top1,enet,blup"))
}

##the function outputs above and the downloaded existing weight often in the .wgt.RDat format (for example, [http://gusevlab.org/projects/fusion/](http://gusevlab.org/projects/fusion/))

###conduct .db

####model.txt
setwd("./reproduce/FOCUS/weight")
ot <- mixedsort(list.files(".",full.names=F))
n <- length(ot)
model <- cbind.data.frame(c(1:n), rep("bslmm",n), rep(1,n), c(1:n))
colnames(model) <- c("id", "inference", "ref_id", "mol_id")
setwd(gsub("/weight", "", getwd()))
system("mkdir db")
setwd("./db")
write.table(model, "model.txt", row.names = F,quote = F, col.names = T, sep = "\t")

####weight.txt
wt = att = list()
setwd(gsub("/db", "/weight", getwd()))
for(i in 1:n){
    onfi <- ot[i]
    load(onfi)
    genefile <- gsub(".wgt.RDat","",onfi)
    bim <- as.data.frame(fread(paste0(dir,"/reproduce/simulation_data_generate/GEUVADIS/",genefile,".bim"),header=F))
    p <- nrow(bim)
    ##Here we choose the BSLMM weight.
    all <- cbind.data.frame(bim[,2],bim[,1],bim[,4:6],wgt.matrix[,colnames(wgt.matrix)=="bslmm"],rep(NA,p),rep(i,p))
    colnames(all) <- c("snp","chrom","pos","effect_allele","alt_allele","effect","std_error","model_id")
    wt[[i]] <- all
    cv <- cbind.data.frame(c("cv.R2","cv.R2.pval"),cv.performance[,1],rep(i,2))
    colnames(cv) <- c("attr_name","value","model_id")
    att[[i]] <- cv
  }
rb <- Reduce(rbind,wt)
N <- nrow(rb)
All <- cbind.data.frame(c(1:N),rb)
colnames(All)[1] <- "id"
setwd(gsub("/weight", "/db", getwd()))
write.table(All, "weight.txt", row.names = F,quote = F, col.names = T, sep = "\t")

####modelattribute.txt
modelatt <- Reduce(rbind, att)
N <- nrow(modelatt)
att.txt <- cbind.data.frame(c(1:N), modelatt)
colnames(att.txt)[1] <- "id"
write.table(att.txt, "modelattribute.txt", row.names = F, quote = F, col.names = T,sep = "\t")

####molecularfeature.txt
o1 <- gsub(".wgt.RDat","",ot)
infos=as.data.frame(fread(paste0(dir,"/reproduce/simulation_data_generate/geneinfo.txt"),header=T))

###We focused on protein coding genes and lincRNAs that are annotated in GENCODE (release 12).
###please see https://github.com/yuanzhongshang/GIFT/blob/main/reproduce/gencodev12.tsv for the information of 15,577 genes included in this research.
info <- infos[infos$genetype2 %in% o1,]
sub <- info[order(match(info$genetype2,o1)),]
mole <- cbind.data.frame(c(1:n), sub[,5], rep("<NA>",n), sub[,4], sub[,3], sub[,6], sub[,1:2])
colnames(mole) <- c("id", "ens_gene_id", "ens_tx_id", "mol_name", "type", "chrom", "tx_start", "tx_stop")
write.table(mole, "molecularfeature.txt",row.names = F, quote = F, col.names = T, sep = "\t")

####refpanel.txt
myref <- cbind.data.frame(1, "GEUVADIS", "blood", "rnaseq")
colnames(myref) <- c("id", "ref_name", "tissue", "assay")
write.table(myref, "refpanel.txt", row.names = F, quote = F, col.names = T,sep = "\t")
  
####construct FOCUS.db file
model <- as.data.frame(fread("model.txt",header = T))
modelattribute <- as.data.frame(fread("modelattribute.txt", header = T))
molecularfeature <- as.data.frame(fread("molecularfeature.txt", header = T))
refpanel <- as.data.frame(fread("refpanel.txt",header = T))
weight <- as.data.frame(fread("weight.txt",header = T))

setwd(gsub("/db", "", getwd()))
mydb <- dbConnect(RSQLite::SQLite(), "weight.db")
dbWriteTable(mydb, "model", model)
dbWriteTable(mydb, "modelattribute", modelattribute)
dbWriteTable(mydb, "molecularfeature",molecularfeature)
dbWriteTable(mydb, "refpanel", refpanel)
dbWriteTable(mydb, "weight", weight)

##perform GWAS using PLINK or process the obtained GWAS summary statistics
##If you perform GWAS using other software, please provide a .txt file contains columns: "CHR", "SNP", "BP", "A1", "A2", "Z", "P", "N".
##conduct the pheotype file used in plink
FID <- read.table(paste0(dir,"/reproduce/simulation_data_generate/Zy.fam"))$V1
IID <- FID
data <- cbind(FID,IID,Y)
colnames(data)[3] <- "V1"
write.table(data,"gwasY.txt",quote=F,row.names=F)

setwd(dir)
system("plink --bfile ./reproduce/simulation_data_generate/Zy --pheno  ./reproduce/FOCUS/gwasY.txt --pheno-name V1  --allow-no-sex --assoc --out ./reproduce/FOCUS/gwas")
data <- fread("./reproduce/FOCUS/gwas.qassoc")
bim <- fread("./reproduce/simulation_data_generate/Zy.bim")
summary <- cbind(data[,1:3], bim[,5:6], data[,8], data[,9], data[,4])
colnames(summary) <- c("CHR", "SNP", "BP", "A1", "A2", "Z", "P", "N")
write.table(summary,"./reproduce/FOCUS/gwas.txt",quote=F,row.names=F)

###clean GWAS summary data
system("focus munge ./reproduce/FOCUS/gwas.txt --output ./reproduce/FOCUS/GWAS.cleaned")

##run FOCUS
system("focus finemap ./reproduce/FOCUS/GWAS.cleaned.sumstats.gz ./reproduce/simulation_data_generate/Zy ./reproduce/FOCUS/weight.db --locations 37:EUR --chr 4 --start 128996665 --stop 130591885 --prior-prob ./reproduce/gencodev12.tsv --p-threshold 1 --out ./reproduce/FOCUS/result")

result <- read.table("./reproduce/FOCUS/result.focus.tsv", header = T)

result
block	ens_gene_id	ens_tx_id	mol_name	tissue	ref_name	type	chrom	tx_start	tx_stop	block_genes	trait	inference_pop1	inter_z_pop1	cv.R2_pop1	cv.R2.pval_pop1	ldregion_pop1	twas_z_pop1	pips_pop1	in_cred_set_pop1
4:128996665-4:130591885	ENSG00000164074.8	<NA>	C4orf29	blood	GEUVADIS	protein_coding	4	128786435	129060866	7	trait	bslmm	NA	0.0318458223286902	6.44829004061019e-05	4:129000314-4:130133779	4.77	0.998	1
4:128996665-4:130591885	ENSG00000077684.10	<NA>	PHF17	blood	GEUVADIS	protein_coding	4	129630779	129896379	7	trait	bslmm	NA	0.191771884804571	1.99043122525914e-23	4:129000314-4:130133779	-1.66	0.123	1
4:128996665-4:130591885	ENSG00000138709.13	<NA>	LARP1B	blood	GEUVADIS	protein_coding	4	128882423	129244086	7	trait	bslmm	NA	0.0551955272133811	1.78083420363514e-07	4:129000314-4:130133779	-3.1	0.115	1
4:128996665-4:130591885	ENSG00000151466.7	<NA>	SCLT1	blood	GEUVADIS	protein_coding	4	129686076	130114764	7	trait	bslmm	NA	0.0959445604427946	5.09993685259739e-12	4:129000314-4:130133779	-1.11	0.0683	1
4:128996665-4:130591885	ENSG00000164040.10	<NA>	PGRMC2	blood	GEUVADIS	protein_coding	4	129090392	129309984	7	trait	bslmm	NA	0.0742207822860967	1.40211140963927e-09	4:129000314-4:130133779	-1.96	0.051	0
4:128996665-4:130591885	ENSG00000251432.2	<NA>	RP11-420A23.1	blood	GEUVADIS	lincRNA	4	129113906	129540549	7	trait	bslmm	NA	-0.00166001190601528	0.630987463676111	4:129000314-4:130133779	-1.81	0.0486	0
4:128996665-4:130591885	ENSG00000151470.7	<NA>	C4orf33	blood	GEUVADIS	protein_coding	4	129914472	130134487	7	trait	bslmm	NA	0.0123982415772871	0.00928127102852506	4:129000314-4:130133779	0.106	0.0408	0
4:128996665-4:130591885	NULL.MODEL	NA	NULL	NA	NA	NULL	4	NA	NA	7	trait	NA	NA	NA	NA	4:129000314-4:130133779	0	0.000349	0
```

## Run FOGS
FOGS takes eQTL-derived weights, reference LD, GWAS summary statistics and gene list as inputs.
```r
library("data.table")
library(RSQLite)
library(gtools)

##load the pindex from the simulation data.
##You can also load from the weight file directly.
dir=getwd()
load("./reproduce/simulation_data_generate/data_generate_individual.RData")

##prepare the eQTL-derived weights
weights <- "./reproduce/FOCUS/weight.db"
weightsave <- paste0(dir,"/reproduce/FOGS/weight.txt")
sqlite.driver <- dbDriver("SQLite")
db <- dbConnect(sqlite.driver,dbname = weights)
dbListTables(db)
weights <- dbReadTable(db, "weight")

setwd("./reproduce/FOCUS/weight")
ot <- mixedsort(list.files(".",full.names=F))
dn <- NULL
for(i in 1:length(gene)){
  onfi <- ot[i]
  gen1 <- gsub(".wgt.RDat","",onfi)
  dn <- c(dn, rep(gen1, pindex[which(gene %in% gen1)]))
}
weightsn <- cbind(weights, dn)
weightsn <- weightsn[,c(2,10,7,6,5)]
colnames(weightsn) <- c("rsid", "gene", "weight", "ref_allele", "eff_allele")
write.table(weightsn, weightsave, col.names = T, row.names = F, quote = F)

weights <- paste0(dir,"/reproduce/FOGS/weight.txt")

##prepare the summary statistics
##If you perform GWAS using other software, please provide a .txt file contains columns: "CHR", "SNP", "POS", "A1", "A2", "beta", "se", "Z", "P", "samplesize".
setwd(paste0(dir, "/reproduce/FOGS"))
source("munge.R")

data <- fread(paste0(dir,"/reproduce/FOCUS/gwas.qassoc"))
bim <- fread(paste0(dir,"/reproduce/simulation_data_generate/Zy.bim"))
summary <- cbind(data[,1:3], bim[,5:6], data[,5:6], data[,8], data[,9],data[,4])
colnames(summary) <- c("CHR", "SNP", "POS", "A1", "A2", "beta", "se", "Z", "P", "samplesize")
write.table(summary, "gwas.txt", quote = F, row.names = F)

sumstats <- paste0(dir,"/reproduce/FOGS/gwas.txt")
data_processed <- munge_sumstat(sumstats)
write.table(data_processed, "gwas_processed_data.txt", col.names = T, row.names = F, quote = F)

sumstat <- paste0(dir,"/reproduce/FOGS/gwas_processed_data.txt")

##load the directory of LD reference
##Here we used the same data as above.
refld <- paste0(dir,"/reproduce/simulation_data_generate/Zy")

##load the directory of gene list
genelist <- paste0(dir,"/reproduce/FOGS/genelist15577.txt")

##match the locus in LDetect region used in the FOGS
ldetect <- read.table(paste0(dir,"/reproduce/FOGS/LDetectregion.txt"),header=T)
chr_id <- 4
locus_id <- which(ldetect[ldetect$chr=="chr4",]$start==128996665)

##other inputs
outd <- paste0(dir,"/reproduce/FOGS/result")
saveprefix <- paste0(dir,"/reproduce/FOGS/result")
loci <-paste0(dir,"/reproduce/LDetect/EUR/")

##run FOGS
system(paste0("Rscript FOGS.R --refld ",refld," --outd ",outd," --loci ",loci," --weights ",weights," --genelist ",genelist," --sumstat ",sumstat," --saveprefix ",saveprefix," --chr_id ",chr_id," --locus_id ",locus_id))

result <- read.table(paste0(dir,"/reproduce/FOGS/resultCHR_4_Locus",locus_id,".txt"),header = T)

result
CHR ID P0 P1 n.SNP n.condSNP FOGS-aSPU TWAS Focus Runtime(s)
4 C4orf33 129914472 130134487 239 112 0.948051948051948 0.912378129020417 0.00015761489003644 24.856
4 RP11-420A23.1 129113906 129540549 252 91 0.68031968031968 0.0666979609501761 0.000805957175765842 18.123
4 SCLT1 129686076 130114764 324 101 0.405594405594406 0.261608296740732 0.000289547991162607 27.113
4 PHF17 129630779 129896379 163 120 0.175824175824176 0.0951010742555474 0.00060903461396202 17.343
4 PGRMC2 129090392 129309984 103 131 0.010989010989011 0.0493022814419665 0.00102967099297041 12.53
##Of note, FOGS fails to perform one gene when the cis-SNPs of one gene are fully contained in the cis-SNPs of another gene. 
##Here cis-SNPs of C4orf29 and LARP1B are included in that of C4orf33. Thus, the result does not contain these two genes.
```

## Run MV-IWAS
MV-IWAS takes GWAS summary statistics, reference LD, eQTL weight database as inputs. Here we used the same data as above.
```r
library(MVIWAS)

##load the simulation data or the real data with the specific format from the summary statistics
###Here we used the simulation data.
dir=getwd()
load("./reproduce/simulation_data_generate/data_generate_summary.RData")

##prepare the eQTL-derived weights
p=sum(pindex)

##Here we used the BSLMM weight for the fair comparison.
betax <- NULL
setwd("./reproduce/FOCUS/weight")
for(i in 1:length(gene)){
  if(i == length(gene)){
    load(paste0(gene[i],".wgt.RDat"))
    betax1 <- wgt.matrix[,colnames(wgt.matrix)=="bslmm"]
    betax <- c(betax,betax1)
  }else{
   load(paste0(gene[i],".wgt.RDat"))
   betax1 <- wgt.matrix[,colnames(wgt.matrix)=="bslmm"]
   betax <- c(betax, betax1, rep(0,p))
  }
}
betax <- matrix(betax, p, length(gene))

betay <- Zscore2/sqrt(n2-1)
se_betaZY <- matrix(1/sqrt(n2-1),p,1)

##run MV-IWAS
result <- mv_iwas_summ(betay, se_betaZY, betax, LDmatrix2, n2, "Continuous", n_case = NULL, n_control = NULL)
result$TERM <- gene
write.table(result,paste0(dir,"/reproduce/MVIWAS/result.txt"), quote = F, row.names = F)

result
    MODEL          TERM         BETA         SE          Z            P NSNP
1 MV-IWAS       C4orf29  0.221597362 0.05474690  4.0476699 5.173003e-05 1231
2 MV-IWAS       C4orf33 -0.008441355 0.06985365 -0.1208434 9.038150e-01 1231
3 MV-IWAS        LARP1B -0.105409347 0.06662337 -1.5821678 1.136113e-01 1231
4 MV-IWAS        PGRMC2  0.055406960 0.05005174  1.1069938 2.682966e-01 1231
5 MV-IWAS         PHF17 -0.024894469 0.02305678 -1.0797030 2.802745e-01 1231
6 MV-IWAS RP11-420A23.1  0.062661326 0.12442069  0.5036246 6.145252e-01 1231
7 MV-IWAS         SCLT1 -0.005552580 0.03470806 -0.1599796 8.728972e-01 1231
  DiseaseType
1  Continuous
2  Continuous
3  Continuous
4  Continuous
5  Continuous
6  Continuous
7  Continuous

##Note that, matrix objects now also inherit from class "array" since R 4.0.0. mv_iwas_summ() contains some class checks which should modify to avoid the bug.
##If you used the R verion over 4.0.0, source("./reproduce/MVIWAS/mv_iwas_summ.R")
```
## Run Marginal GWAS
We used PLINK-1.9 to perform the GWAS analysis. 
```r
library("data.table")

##conduct the pheotype file used in plink
dir <- getwd()
FID <- read.table(paste0(dir,"/reproduce/simulation_data_generate/Zy.fam"))$V1
IID <- FID
data <- cbind(FID,IID,Y)
colnames(data)[3] <- "V1"
write.table(data,"gwasY.txt",quote=F,row.names=F)

setwd(dir)
system("plink --bfile ./reproduce/simulation_data_generate/Zy --pheno  ./reproduce/FOCUS/gwasY.txt --pheno-name V1  --allow-no-sex --assoc --out ./reproduce/FOCUS/gwas")
data <- fread("./reproduce/FOCUS/gwas.qassoc")

head(data)
   CHR         SNP        BP NMISS     BETA      SE        R2      T         P
1:   4  rs34350456 129000314  5000 -0.07705 0.02115 0.0026490 -3.644 0.0002717
2:   4   rs1993722 129002172  5000 -0.07705 0.02115 0.0026490 -3.644 0.0002717
3:   4   rs7684955 129007262  5000 -0.07705 0.02115 0.0026490 -3.644 0.0002717
4:   4 rs116213612 129011061  5000  0.06143 0.03414 0.0006473  1.799 0.0720400
5:   4   rs2306054 129012181  5000 -0.07705 0.02115 0.0026490 -3.644 0.0002717
6:   4   rs1064205 129012638  5000 -0.07705 0.02115 0.0026490 -3.644 0.0002717
```

## Run Marginal TWAS
### Using individual-level data as input 
```r
library(gtools)
library("BEDMatrix")

##Here we also used the BSLMM weight above.
dir=getwd()
setwd("./reproduce/FOCUS/weight")

ot <- mixedsort(list.files(".",full.names=F))
gene <- NULL
for(i in 1:length(ot)){
  onfi <- ot[i]
  gene[i] <- gsub(".wgt.RDat","",onfi)
}

##impute function
##We impute the genotypes from the UK Biobank.
impu <- function(x) {
  x[is.na(x)] <- mean(x, na.rm = TRUE)
  x
}

##load the phenotype
setwd(dir)
load("./reproduce/simulation_data_generate/data_generate_individual.RData")

Z <- NULL
P <- NULL
setwd("./reproduce/FOCUS/weight")
for(i in 1:length(gene)){
  load(ot[i])
  beta <- wgt.matrix[,colnames(wgt.matrix)=="bslmm"]
  Zy <- BEDMatrix(paste0(dir,"/reproduce/simulation_data_generate/UKB/",gene[i]))
  Zy <- apply(Zy ,2, impu)
  p <- dim(Zy)[2]
  tmp <- cbind(Zy,Y)
  tmp <- scale(na.omit(tmp))
  Zy <- tmp[,1:p]
  Y <- tmp[,p+1]
  pred <- Zy %*% beta
  fit <- lm(Y ~ pred)
  Z[i] <- coefficients(summary(fit))[2,3]
  P[i] <- coefficients(summary(fit))[2,4]
}
result <- data.frame(gene, Z, P)
write.table(result,paste0(dir, "/reproduce/TWAS/result_individual.txt"), quote = F, row.names = F)

result
           gene          Z            P
1       C4orf29  4.9493385 7.691249e-07
2       C4orf33  0.1629571 8.705588e-01
3        LARP1B -3.1384610 1.708262e-03
4        PGRMC2 -1.9652850 4.943655e-02
5         PHF17 -1.6692476 9.513096e-02
6 RP11-420A23.1 -1.8330907 6.684853e-02
7         SCLT1 -1.1222759 2.617990e-01
```

### Using summary statistics as input
```r
##load the simulation data or the real data with the specific format from the summary statistics
###Here we used the simulation data.
dir=getwd()
load("./reproduce/simulation_data_generate/data_generate_summary.RData")

Z <- NULL
P <- NULL
pindexsum <- c(0,cumsum(pindex))
setwd("./reproduce/FOCUS/weight")
for(i in 1:length(gene)){
  load(ot[i])
  beta <- wgt.matrix[,colnames(wgt.matrix)=="bslmm"]
  Z[i] <- sum(beta * Zscore2[(pindexsum[i]+1):pindexsum[i+1]])/sqrt(sum(beta*(LDmatrix2[(pindexsum[i]+1):pindexsum[i+1],(pindexsum[i]+1):pindexsum[i+1]] %*% beta)))
  P[i] <- 2*pnorm(-abs(Z[i]), 0, 1)
}
result <- data.frame(gene, Z, P)
write.table(result,paste0(dir, "/reproduce/TWAS/result_summary.txt"), quote = F, row.names = F)

result
           gene         Z            P
1       C4orf29  4.937748 7.902983e-07
2       C4orf33  0.162973 8.705397e-01
3        LARP1B -3.135687 1.714523e-03
4        PGRMC2 -1.964723 4.944635e-02
5         PHF17 -1.668949 9.512741e-02
6 RP11-420A23.1 -1.832658 6.685343e-02
7         SCLT1 -1.122247 2.617575e-01
```
