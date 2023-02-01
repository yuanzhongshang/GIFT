---
layout: page
title: Example Analysis
description: ~
---
This tutorial is the example analysis with GIFT for the individual-level data and summary statistics, respectively. Before runing the tutorial, make sure that the GIFT package is installed. Installation instructions see the [link](https://github.com/yuanzhongshang/GIFT/documentation/02_installation.html)

## For the individual-level data
The example data for runing the tutorial can be downloaded in this [page](https://github.com/yuanzhongshang/GIFT/documentation/03_data.html)
Here are the details about the required data input illustrated. 
### 1. standardized cis-genotype matrix in eQTL data, e.g.,
```r
#### load the simulated scaled genenotype matrix in eQTL data,
Zx<-read.table("Zx.txt")
Zx<-as.matrix(Zx)
```

### 2. standardized cis-genotype matrix in GWAS data,  e.g.,
```r
#### load the simulated scaled genenotype matrix in GWAS data
Zy<-read.table("Zy.txt")
Zy<-as.matrix(Zy)
```

### 3. complex gene expression matrix,  e.g.,
```r
#load the simulated exposure or gene expression vector
X<-read.table("X.txt")
X<-as.matrix(X)
```

### 4. standarized trait vector,  e.g.,
```r
#load the simulated phenotype vector
Y<-read.table("Y.txt")$V1
```

## Conditional fine-mapping for TWAS analysis
```r
library(GIFT)
``` 
The function `GIFT_individual` is for conditional fine-mapping for in TWAS with individual-level data. The essential inputs are:
- X: The complex gene expression matrix, each column is the standardized specific expression. 
- Y: The standarized trait vector.
- Zx: The standardized cis-genotype matrix in eQTL data.
- Zy: The standardized cis-genotype matrix in GWAS data.
- pindex: A vector with each element represents the number of cis-SNPs for each gene.
- max_iterin: The maximum iteration, which can be determined by users. Default is 1000. 
- epsin: The convergence tolerance of the absolute value of the difference between the nth and (n+1)th log likelihood, which can be determined by users. Default is 1e-4. 
- Cores: The number of cores used in analysis. If the number of cores is greater than 1, analysis will perform with fast parallel computing. The function mclapply() depends on another R package "parallel" in Linux. Default is 1.

```r
result<-GIFT_individual(X, Y, Zx, Zy, pindex, max_iterin =1000,epsin=1e-4,Cores=1)

```
The result is a list of estimated parameters including the causal effects and p values for each gene in a focal region. 
