---
layout: page
title: Example Analysis
description: ~
---
This tutorial is the example analysis with GIFT for the individual-level data and summary statistics, respectively. Before runing the tutorial, make sure that the GIFT package is installed. Installation instructions see the [link](https://yuanzhongshang.github.io/GIFT/documentation/02_installation.html)

## For the individual-level data
The example data for runing the tutorial can be downloaded in this [page](https://yuanzhongshang.github.io/GIFT/documentation/03_data.html)
Here are the details about the required data input illustrated. 
### 1. Standardized cis-genotype matrix in eQTL data, e.g.,
```r
#### load the simulated scaled genenotype matrix in eQTL data,
Zx<-read.table("Zx.txt")
Zx<-as.matrix(Zx)
```

### 2. Standardized cis-genotype matrix in GWAS data,  e.g.,
```r
#### load the simulated scaled genenotype matrix in GWAS data
Zy<-read.table("Zy.txt")
Zy<-as.matrix(Zy)
```

### 3. Complex gene expression matrix,  e.g.,
```r
#### load the simulated exposure or gene expression vector
X<-read.table("X.txt")
X<-as.matrix(X)
```

### 4. Standarized trait vector,  e.g.,
```r
#### load the simulated phenotype vector
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
pindex=c(24,33)
result<-GIFT_individual(X, Y, Zx, Zy, pindex, max_iterin =1000,epsin=1e-4,Cores=1)
```
The result is a list of estimated parameters including the causal effects and p values for each gene in a focal region. 

## For the summary statistics
The example data for runing the tutorial can be downloaded in this [page](https://yuanzhongshang.github.io/GIFT/documentation/03_data.html)
Here are the details about the required data input illustrated. 
### 1. Zscore matrix of the cis-SNP effect size matrix, e.g.,
```r
#### load the Zscore matrix for the cis-SNP effect size from the eQTL data
Zscore1<-read.table("Zscore1.txt")
Zscore1<-as.matrix(Zscore1)
```

### 2. Zscore vector of the cis-SNP effect size vector for one specific trait in GWAS data,  e.g.,
```r
#### load the Zscore vector for the cis-SNP effect size from the GWAS data
Zscore2<-read.table("Zscore2.txt")$V1
```

### 3. LD matrix in eQTL data,  e.g.,
```r
#### load the LD matrix for the cis-SNPs in the eQTL data
LDmatrix1<-read.table("LDmatrix1.txt")
LDmatrix1<-as.matrix(LDmatrix1)
```

### 4. LD matrix in GWAS data,  e.g.,
```r
#### load the LD matrix for the cis-SNPs in the GWAS data
LDmatrix2<-read.table("LDmatrix2.txt")
LDmatrix2<-as.matrix(LDmatrix2)
```
### 5. Estimated correlated matrix of gene expressions,  e.g.,
```r
#### load the estimated correlated matrix of gene expressions
R<-read.table("R.txt")
R<-as.matrix(R)
``` 

## Conditional fine-mapping for TWAS analysis
```r
library(GIFT)
``` 
The function `GIFT_summary` is for conditional fine-mapping for in TWAS with summary statistics. The essential inputs are:
- Zscore_1: The Zscore matrix of the cis-SNP effect size matrix, each column for one specific gene in eQTL data.
- Zscore_2: The Zscore vector of the cis-SNP effect size vector for one specific trait in GWAS data.
- Sigma1: The LD matrix in eQTL data.
- Sigma2: The LD matrix in GWAS data, both Sigma1 and Sigma2 are often the same from the reference panel.
- R: The estimated correlation matrix of gene expressions.
- n1: The sample size of eQTL data.
- n2: The sample size of GWAS data.
- pindex: A vector with each element represents the number of cis-SNPs for each gene.
- max_iterin: The maximum iteration, which can be determined by users. Default is 1000. 
- epsin: The convergence tolerance of the absolute value of the difference between the nth and (n+1)th log likelihood, which can be determined by users. Default is 1e-4. 
- Cores: The number of cores used in analysis. If the number of cores is greater than 1, analysis will perform with fast parallel computing. The function mclapply() depends on another R package "parallel" in Linux. Default is 1.

```r
pindex=c(24,33)
n1=465
n2=5000
result<-GIFT_summary(Zscore1, Zscore2, LDmatrix1, LDmatrix2, R, n1, n2, pindex, max_iterin =1000,epsin=1e-4, Cores=1)

```
The result is a list of estimated parameters including the causal effects and p values for each gene in a focal region. 
