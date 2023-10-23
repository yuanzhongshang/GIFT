---
layout: page
title: Analysis Reproduce
description: ~
---

## Generate the simulation data
We randomly selected 50 regions from 1,533 regions. The region information is [here](https://github.com/yuanzhongshang/GIFT/blob/main/reproduce/LDetectregion1533.txt). The full information about the region can be found in [https://bitbucket.org/nygcresearch/ldetect-data/src/master/](https://bitbucket.org/nygcresearch/ldetect-data/src/master/). For each region in turn, we performed 20 simulation replicates, resulting in a total of 1,000 simulation replicates per setting.
For each region, we conducted the simulations based on the realistic genotypes from GEUVADIS (n1=465) and UK Biobank (n2=5,000). Take a region on chr 4 for example. This region includes seven genes: C4orf29, C4orf33, LARP1B, PGRMC2, PHF17, RP11-420A23.1, and SCLT1. We set C4orf29 as the causal gene with the effect size being sqrt(0.1).

### Generate the individual level data
```r
library("data.table")
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
PVEzy <- c(0,0,0.01,0)
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
save(X, Y, Zx, Zy, gene, pindex, file = "./reproduce/simulation_data_generate/data_generate_individual.Rdata")
```

### Convert the individual level data into the summary statistics
```r
###calculate the z-score from GEUVADIS
n1=dim(Zx)[1]
Zscore1 =NULL
for(i in 1:length(gene)){
  Zscore1 <- rbind(Zscore1,(1/sqrt(n1-1))*t(Zx[,(pindexsum[i]+1):pindexsum[i+1]]) %*% X[,i])
}

###calculate the LD matrix from GEUVADIS
LDmatrix1 <- (1/(n1-1))*t(Zx) %*% Zx

###calculate the z-score from UK Biobank
n2=dim(Zy)[1]

###calculate the LD matrix from UK Biobank
LDmatrix2 <- (1/(n2-1))*t(Zy) %*% Zy
Zscore2 <- (1/sqrt(n2-1))* (t(Zy) %*% Y)

###save these variables
save(Zscore1, LDmatrix1, Zscore2, LDmatrix2, R, n1, n2, gene, pindex, file = "./reproduce/simulation_data_generate/data_generate_summary.Rdata")
```

## Run GIFT
```r
library(GIFT)
```
### Using individual-level data as input 
```r
##load the simulation data or the real data with the specific format from individual level data.
###Here we used the simulation data above.
load("./reproduce/simulation_data_generate/data_generate_individual.Rdata")

##run GIFT
result <- GIFT_individual(X, Y, Zx, Zy, gene, pindex, maxiter=1000, tol=1e-4, ncores=1)

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
load("./reproduce/simulation_data_generate/data_generate_summary.Rdata")

##run GIFT
result <- GIFT_summary(Zscore1, Zscore2, LDmatrix1, LDmatrix2, R, n1, n2, gene, pindex, maxiter=1000, tol=1e-4, ncores=1)

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
FOCUS takes GWAS summary statistics, reference LD, eQTL weight database as input. Here we used the same data as above. If you do not have the in-sample LD, you may use the LD from 1000 Genomes reference: [https://alkesgroup.broadinstitute.org/FUSION/LDREF.tar.bz2](https://alkesgroup.broadinstitute.org/FUSION/LDREF.tar.bz2).

```r
library("data.table")
library("gtools")

##load the simulation data or the real data with the specific format from individual level data.
###Here we used the simulation data.
dir=getwd()
load("./reproduce/simulation_data_generate/data_generate_individual.Rdata")

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

##the function outputs above and the downloaded existing weight often in the .wgt.RDat format

###conduct .db
library(RSQLite)
library(gtools)

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
    genefile < -gsub(".wgt.RDat","",onfi)
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
colnames(mole)= <- ("id", "ens_gene_id", "ens_tx_id", "mol_name", "type", "chrom", "tx_start", "tx_stop")
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
summar <- cbind(data[,1:3], bim[,5:6], data[,8], data[,9], data[,4])
colnames(summary) <- c("CHR", "SNP", "BP", "A1", "A2", "Z", "P", "N")
write.table(summary,"./reproduce/FOCUS/gwas.txt",quote=F,row.names=F)

###clean GWAS summary data
system("focus munge ./reproduce/FOCUS/gwas.txt --output ./reproduce/FOCUS/GWAS.cleaned")

##Run FOCUS
system("focus finemap ./reproduce/FOCUS/GWAS.cleaned.sumstats.gz ./reproduce/simulation_data_generate/Zy ./reproduce/FOCUS/weight.db --locations 37:EUR --chr 4 --start 128996665 --stop 130591885 --prior-prob ./reproduce/gencodev12.tsv --p-threshold 1 --out ./reproduce/FOCUS/result")

result <- read.table("./reproduce/FOCUS/result", header = T)

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
