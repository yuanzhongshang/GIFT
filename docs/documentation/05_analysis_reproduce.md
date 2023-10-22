---
layout: page
title: Analysis Reproduce
description: ~
---

## Generate the simulation data
We randomly selected 50 regions from 1,533 regions. The region information is [here](https://github.com/yuanzhongshang/GIFT/blob/main/reproduce/LDetectregion1533.txt). For each region in turn, we performed 20 simulation replicates, resulting in a total of 1,000 simulation replicates per setting.
For each region, we conducted the simulations based on the realistic genotypes from GEUVADIS (n1=465) and UK Biobank (n2=5,000). Take a region on chr 5 for example. This region includes four genes: RASA1, COX7C, CCNH and TMEM161B. We set RASA1 as the causal gene with the effect size being sqrt(0.1).

### Generate the individual level data
```r
##load the genotypes from GEUVADIS

###load the directory containing the files to be processed only (e.g., plink binary format)
filelocation <- "./example/simulation/individual/pre_process/plink_binary"

###load the directory of plink exe file
plinkexe <- "plink"

###pre-process the file to be a list including gene names vector, cis-genotype matrix and pindex
convert <- pre_process_individual(filelocation, plinkexe)
gene <- convert$gene
Zx <- convert$Z
pindex <- convert$pindex

###please see https://yuanzhongshang.github.io/GIFT/documentation/04_GIFT_Example.html for the details of pre-poccess of genotype data.  

##load the genotypes from UK Biobank
load("./reproduce/simulation_data_generate/Zy.Rdata")

##set the correlation of the residual of gene expressions

###We can set the correlation based on the real data or some specific correlation structure 
###Option 1: the correlation below is realistic from GEUVADIS
R <- as.matrix(read.table("./example/simulation/summary/R.txt"))

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
  b[[i]]=matrix(rnorm(pindex[i], 0, sd = sqrt(PVEzx[i]/pindex[i])), pindex[i], 1)
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
    lp_z =lp_z + casual_effect[i]* as.matrix(Zy[,(pindexsum[i]+1):pindexsum[i+1]]) %*% b[[i]]
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
