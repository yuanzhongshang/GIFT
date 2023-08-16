---
layout: page
title: Example Analysis
description: ~
---
This tutorial is the example analysis with GIFT for the individual-level data and summary statistics, respectively. Before runing the tutorial, make sure that the GIFT package is installed. Installation instructions see the [link](https://yuanzhongshang.github.io/GIFT/documentation/02_installation.html).

## For the individual-level data
The example data for runing the tutorial can be downloaded in this [page](https://yuanzhongshang.github.io/GIFT/documentation/03_data.html).
Here are the details about the required data input illustrated. 
### 1. Standardized cis-genotype matrix in eQTL data, e.g.,
```r
#### load the simulated scaled genenotype matrix in eQTL data
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

### Conditional fine-mapping for TWAS analysis
```r
library(GIFT)
``` 
The function `GIFT_individual` is for conditional fine-mapping for in TWAS with individual-level data. The essential inputs are:
- X: The complex gene expression matrix, each column is the standardized specific expression. 
- Y: The standarized trait vector.
- Zx: The standardized cis-genotype matrix in eQTL data.
- Zy: The standardized cis-genotype matrix in GWAS data.
- gene: The gene name vector.
- pindex: A vector with each element represents the number of cis-SNPs for each gene.
- max_iterin: The maximum iteration, which can be determined by users. Default is 1000. 
- epsin: The convergence tolerance of the absolute value of the difference between the nth and (n+1)th log likelihood, which can be determined by users. Default is 1e-4. 
- Cores: The number of cores used in analysis. If the number of cores is greater than 1, analysis will perform with fast parallel computing. The function mclapply() depends on another R package "parallel" in Linux. Default is 1.

```r
gene=c("RASA1", "COX7C", "CCNH", "TMEM161B")
pindex=c(63, 23, 41, 96)
result<-GIFT_individual(X, Y, Zx, Zy, gene, pindex, max_iterin =1000, epsin=1e-4, Cores=1)
```
The result is a data frame including the causal effect estimates and p values for each gene in a focal region. 
```r
result
      gene causal_effect            p
1    RASA1    0.35441121 7.132522e-06
2    COX7C    0.02106124 8.028767e-01
3     CCNH    0.01116380 8.598955e-01
4 TMEM161B   -0.03302958 3.257469e-01
```

## For the summary statistics
The example data for runing the tutorial can be downloaded in this [page](https://yuanzhongshang.github.io/GIFT/documentation/03_data.html).
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

### Conditional fine-mapping for TWAS analysis
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
- gene: The gene name vector.
- pindex: A vector with each element represents the number of cis-SNPs for each gene.
- max_iterin: The maximum iteration, which can be determined by users. Default is 1000. 
- epsin: The convergence tolerance of the absolute value of the difference between the nth and (n+1)th log likelihood, which can be determined by users. Default is 1e-4. 
- Cores: The number of cores used in analysis. If the number of cores is greater than 1, analysis will perform with fast parallel computing. The function mclapply() depends on another R package "parallel" in Linux. Default is 1.

```r
gene=c("RASA1", "COX7C", "CCNH", "TMEM161B")
pindex=c(63, 23, 41, 96)
n1=465
n2=5000
result<-GIFT_summary(Zscore1, Zscore2, LDmatrix1, LDmatrix2, R, n1, n2, gene, pindex, max_iterin =1000, epsin=1e-4, Cores=1)

```
The result is a data frame including the causal effect estimates and p values for each gene in a focal region. 
```r
result
      gene causal_effect            p
1    RASA1    0.35439571 7.088416e-06
2    COX7C    0.02104398 8.029700e-01
3     CCNH    0.01114944 8.609294e-01
4 TMEM161B   -0.03302900 3.257526e-01
```

## For the two-stage version
The example data for runing the tutorial can be downloaded in this [page](https://yuanzhongshang.github.io/GIFT/documentation/03_data.html).
Here are the details about the required data input illustrated. 
### 1. Weight matrix for the cis-SNP effect size, e.g.,
```r
#### load the weight matrix for the cis-SNP effect size from the eQTL data
betax<-read.table("betax.txt")
betax<-as.matrix(betax)
```

### 2. Beta vector and the corresponding se vector for the cis-SNP effect size vector for one specific trait in GWAS data,  e.g.,
```r
#### load the beta vector and the corresponding se vector for the cis-SNP effect size from the GWAS data
betay<-read.table("betay.txt")
betay<-as.matrix(betay)

se_betay<-read.table("se_betay.txt")
se_betay<-as.matrix(se_betay)
```

### 3. LD matrix in GWAS data,  e.g.,
```r
#### load the LD matrix for the cis-SNPs in the GWAS data
Sigma<-read.table("Sigma.txt")
Sigma<-as.matrix(Sigma)
```
### 4. Sample size n from GWAS data,  e.g.,
```r
n=5000
``` 

### Conditional fine-mapping for TWAS analysis
```r
library(GIFT)
``` 
The function `GIFT_two_stage_summ` is for conditional fine-mapping for in TWAS with summary statistics. The essential inputs are:
- betax: The weight matrix for the cis-SNP effect size from the eQTL data.
- betay: The beta vector of the cis-SNP effect size vector for one specific trait in GWAS data.
- se_betay: The se vector of the cis-SNP effect size vector for one specific trait in GWAS data.
- Sigma: The LD matrix in GWAS data, Sigma is often from the reference panel.
- n: The sample size of GWAS data.
- gene: The vector of gene names.

```r
gene=c("RASA1", "COX7C", "CCNH", "TMEM161B")
result<-GIFT_two_stage_summ(betax, betay, se_betay, Sigma, n, gene)
```
The result is a data frame including the z-scores and p values for each gene in a focal region. 
```r
result
      gene         Z         P NSNP
1    RASA1  1.473215 0.1406930  223
2    COX7C -1.216562 0.2237709  223
3     CCNH  1.251621 0.2107080  223
4 TMEM161B -0.758564 0.4481134  223
```

## Visualization for the GIFT result
We visualize the GIFT result incorporating with the marginal GWAS and TWAS results. Here, we load the analyzed results of GWAS and TWAS directly using the example data. The example data for runing the tutorial can be downloaded in this [page](https://yuanzhongshang.github.io/GIFT/documentation/03_data.html).
```r
####load the GWAS results
GWASresult=read.table("GWASresult.txt",header=T)
GWASresult=GWASresult[,c(1,3,5)]
GWASresult$index="GWAS"
colnames(GWASresult)=c("X","BP","P","index")

####load the TWAS results
TWASresult=read.table("TWASresult.txt",header=T)
TWASresult$BP=apply(TWASresult[,c(3,4)],1,mean)
TWASresult=TWASresult[,c(1,6,5)]
TWASresult$index="TWAS"
colnames(TWASresult)=c("X","BP","P","index")

####load the GIFT results
GIFTresult=result
GIFTresult$BP=TWASresult$BP
GIFTresult=GIFTresult[,c(1,4,3)]
GIFTresult$index="GIFT"
colnames(GIFTresult)=c("X","BP","P","index")

####visualize the result by Manhattan plot
data=rbind(GWASresult,TWASresult)
data=rbind(data,GIFTresult)
data$BP=data$BP/1000000
data$index=factor(data$index,levels=c("GWAS","TWAS","GIFT"))

library(ggplot2)
ggplot(data) +
  labs(x="Position on Chr5(Mb)", y=expression(paste(-log[10]," (p-value)"))) +
  geom_point(aes(x=BP, y=-log10(P),color=index,shape=index,size=index)) +
  theme_bw() +
  theme( 
    legend.title = element_blank(),
    legend.position="top",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(), 
    axis.line.y = element_line(color = "black", linetype ="solid"), 
    axis.line.x = element_line (color = "black",linetype = "solid") 
  )+geom_hline(yintercept =-log10(0.05),lty="dashed")+
  scale_discrete_manual(values=c("grey","#377EB8","#F23557"), aesthetics = 'colour')+
  scale_shape_manual(values=c(19,15,18))+  scale_size_manual(values=c(1,1.5,2))+theme(panel.grid=element_blank())
``` 
Here is an example output:
![GIFT\_pipeline](visualization.png)
