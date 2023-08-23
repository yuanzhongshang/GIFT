---
layout: page
title: Example Analysis
description: ~
---
GIFT examines one genomic region at a time, explicitly models the gene expression correlation and cis-SNP LD across different genes in the region and accounts for the uncertainty in the constructed GReX, and carries out TWAS conditional analysis in a maximum likelihood framework. This tutorial is the example analysis with GIFT. Before runing the tutorial, make sure that the GIFT package is installed. Installation instructions see the [link](https://yuanzhongshang.github.io/GIFT/documentation/02_installation.html). The example data for runing the tutorial can be downloaded in this [page](https://yuanzhongshang.github.io/GIFT/documentation/03_data.html).

## Simulation
We generated the simulation data from the realistic genotypes from GEUVADIS (n1=465) and UK Biobank (n2=5,000) in a region of chr 5. This region includes four genes: RASA1, COX7C, CCNH and TMEM161B. We set RASA1 to be causal with the effect size being sqrt(0.1).

### GIFT: Using individual-level data as input
The function `GIFT_individual` is for conditional fine-mapping for in TWAS with individual-level data. The essential inputs are:
- X: The complex gene expression matrix containing all gene expressions in a specific region, each column is the standardized specific expression. 
- Y: The standarized trait vector.
- Zx: The standardized cis-genotype matrix for all genes in a specific region from eQTL data.
- Zy: The standardized cis-genotype matrix for all genes in a specific region from GWAS data.
- gene: The gene name vector represents the genes in a specific region.
- pindex: A vector with each element represents the number of cis-SNPs for each gene.
The optional inputs are:
- maxiter: The maximum iteration, which can be determined by users. Default is 1000. 
- tol: The convergence tolerance of the absolute value of the difference between the nth and (n+1)th log likelihood, which can be determined by users. Default is 1e-4. 
- ncores: The number of cores used in analysis. If the number of cores is greater than 1, analysis will perform with fast parallel computing. The function mclapply() depends on another R package "parallel" in Linux. Default is 1.

#### Step 1: pre-process the genotype data in various formats.
The function `pre_process_individual` can convert common genotype formats to match GIFT input. Specifically, this function is fexible to handle plink binary format (.bim/.fam./.bed), vcf, ped/map format, csv, and tsv file. Here, we provide various genotype formats from GEUVADIS data in [page](https://yuanzhongshang.github.io/GIFT/documentation/03_data.html) for example. Note that, in this step, cis-genotype matrix has been standardized to have a mean of zero and standard derivation of one. 
```r
#### load the directory contains the files to be processed only (e.g., plink binary format)
filelocation <- "./simulation/individual/pre_process/plink_binary"
#### load the directory of plink exe file
plinkexe <- "plink"
#### pre-process the file to a list including gene names vector, cis-genotype matrix and pindex
convert <- pre_process_individual(filelocation, plinkexe)
gene <- convert$gene
Zx <- convert$Z
pindex <- convert$pindex
```

#### Step 2: Read other required data.
The required example data can be downloaded in this [page](https://yuanzhongshang.github.io/GIFT/documentation/03_data.html). 
```r
#### load the Rdata file containing X, Zy and Y
load("individual_data.Rdata")
```

#### Step 3: Perform conditional fine-mapping for TWAS analysis.
```r
library(GIFT)
result <- GIFT_individual(X, Y, Zx, Zy, gene, pindex, maxiter=1000, tol=1e-4, ncores=1)
```
The result is a data frame including the causal effect estimates and p values for each gene in a focal region. 
```r
result
      gene causal_effect            p
1     CCNH    0.01209639 8.594982e-01
2    COX7C    0.02111053 8.028104e-01
3    RASA1    0.35309347 7.131827e-06
4 TMEM161B   -0.03306309 3.257628e-01
```

### GIFT: Using summary statistics as input
The function `GIFT_summary` is for conditional fine-mapping for in TWAS with summary statistics. The essential inputs are:
- Zscore_1: Zscore matrix of the cis-SNP effect size matrix, each column from eQTL data.
- Zscore_2: Zscore vector of the cis-SNP effect size vector from GWAS data.
- Sigma1: LD matrix from eQTL data.
- Sigma2: LD matrix from GWAS data, both Sigma1 and Sigma2 are often the same from the reference panel.
- R: Estimated correlation matrix of gene expressions.
- n1: Sample size of eQTL data.
- n2: Sample size of GWAS data.
- gene: The gene name vector.
- pindex: A vector with each element represents the number of cis-SNPs for each gene.
The optional inputs are:
- maxiter: The maximum iteration, which can be determined by users. Default is 1000. 
- tol: The convergence tolerance of the absolute value of the difference between the nth and (n+1)th log likelihood, which can be determined by users. Default is 1e-4. 
- ncores: The number of cores used in analysis. If the number of cores is greater than 1, analysis will perform with fast parallel computing. The function mclapply() depends on another R package "parallel" in Linux. Default is 1.

#### Step 1: pre-process the summary statistics in various formats.
The function `pre_process_summary` can convert common summary statistics and LD matrix format to match GIFT input. Specifically, this function is fexible to handle association test output from plink (.qassoc), GEMMA (.assoc.txt) and SAIGE (.txt). While, this function is also fexible to handle LD matrix either from matrix or a long format such as h5 format. Here, we provide various formats from example data in [page](https://yuanzhongshang.github.io/GIFT/documentation/03_data.html). Note that, the summary statistics version of GIFT often requires the in-sample LD matrix. If we cannot have the in-sample LD matrix, it can be also calculated from the reference panel data (e.g., 1,000 Genomes project), it would be better to ensure the ethnicity of the reference panel is consistent with that of the analyzed data. 
```r
#### load the directory contains files of summary statistics from eQTL data only (e.g., the SAIGE output)
eQTLfilelocation <- "./simulation/summary/pre_process/saige/eQTL"
#### load summary statistics from GWAS data (e.g., the SAIGE output)
GWASfile <- "./summary/pre_process/saige/GWAS.txt"
#### load LD matrix from eQTL data and GWAS data (e.g., a long format: h5 format)
eQTLLDfile <- "./summary/pre_process/LDmatrix1.h5"
GWASLDfile <- "./summary/pre_process/LDmatrix2.h5"
#### load the SNP list and cis-SNP number for each gene in a region
snplist <- read.table("./simulation/summary/pre_process/snplist.txt")$V1
pindex <- c(41, 23, 63, 96)
#### pre-process the file to a list including gene names vector, z-score matrix and LD matrix of eQTL data and GWAS data
convert <- pre_process_summary(eQTLfilelocation, eQTLLDfile, GWASfile, GWASLDfile, snplist, pindex)
gene <- convert$gene
Zscore1 <- convert$Zscore1
Zscore2 <- convert$Zscore2
LDmatrix1 <- convert$LDmatrix1
LDmatrix2 <- convert$LDmatrix2
```

#### Step 2: Read other required data.
```r
### input the sample sizes of eQTL data and GWAS data
n1 <- 465
n2 <- 5000
### load the estimated correlated matrix of gene expressions
R <- as.matrix(read.table("./simulation/summary/R.txt"))
```

#### Step 3: Perform conditional fine-mapping for TWAS analysis.
```r
library(GIFT)
result <- GIFT_summary(Zscore1, Zscore2, LDmatrix1, LDmatrix2, R, n1, n2, gene, pindex, maxiter=1000, tol=1e-4, ncores=1)
```
The result is a data frame including the causal effect estimates and p values for each gene in a focal region. 
```r
result
      gene causal_effect            p
1     CCNH    0.01160943 8.651676e-01
2    COX7C    0.02205361 7.943974e-01
3    RASA1    0.35339282 8.671951e-06
4 TMEM161B   -0.03373786 3.197526e-01 
```

### Two-stage version of GIFT: Using pre-trained weights and summary statistics as input
The function `GIFT_two_stage_summ` is for conditional fine-mapping for in TWAS using pre-trained weights and summary statistics. The two-stage version of GIFT is not only computationally efficient but also allows us to make use of existing gene expression prediction models for more convenient TWAS fine-mapping. The example data for runing the tutorial can be downloaded in this [page](https://yuanzhongshang.github.io/GIFT/documentation/03_data.html). The essential inputs are:
- betax: Weight matrix for the cis-SNP effect size from eQTL data.
- betay: Beta vector of the cis-SNP effect size vector from GWAS data.
- se_betay: Se vector of the cis-SNP effect size vector from GWAS data.
- Sigma: LD matrix from GWAS data.
- n: Sample size of GWAS data.
- gene: The vector of gene names.
  
#### Step 1: Read the eQTL weight.
Gene expression prediction is a critical component of two-stage TWAS methods. Examples of prediction models include lasso and elastic net (enet) as implemented in prediXcan, Best Linear Unbiased Prediction (BLUP), the top SNPs (top1) and Bayesian sparse linear mixed model (BSLMM) as implemented in TWAS/FUSION, latent Dirichlet process regression (DPR) as implemented in both DPR and TIGAR. For a specific region, the weights from all genes consist a block diagonal matrix. Each block is a vector of cis-SNP effect using the same gene expression prediction model. The function `weightconvert` can convert a list including weights for multiple genes into a required block diagonal matrix.
```r
#### load the weight matrix from the eQTL data (e.g., BLUP)
setwd("./simulation/two_stage/weights")
CCNH <- as.matrix(read.table("CCNHweight.txt"))
COX7C <- as.matrix(read.table("COX7Cweight.txt"))
RASA1 <- as.matrix(read.table("RASA1weight.txt"))
TMEM161B <- as.matrix(read.table("TMEM161Bweight.txt"))
weightlist <- list(CCNH = CCNH, COX7C = COX7C, RASA1 = RASA1, TMEM161B = TMEM161B)
betax <- weightconvert(weightlist)
```

#### Step 2: Read the beta vector, corresponding se vector and LD matrix from GWAS data.
The function `pre_process_twostage` can convert common summary statistics and LD matrix format to match GIFT input. Specifically, this function is fexible to handle association test output from plink (.qassoc), GEMMA (.assoc.txt) and SAIGE (.txt). While, this function is also fexible to handle LD matrix either from a matrix or a long format such as h5 format. The example data is same as the section from GWAS data above in [page](https://yuanzhongshang.github.io/GIFT/documentation/03_data.html). Note that, the two-stage version of GIFT often requires the in-sample LD matrix. If we cannot have the in-sample LD matrix, it can be also calculated from the reference panel data (e.g., 1,000 Genomes project), it would be better to ensure the ethnicity of the reference panel is consistent with that of the analyzed data. 
```r
#### load summary statistics from GWAS data (e.g., the plink output)
GWASfile="./summary/pre_process/plink/GWAS.qassoc"
#### load LD matrix from eQTL data and GWAS data (e.g., a matrix)
GWASLDfile="./summary/pre_process/LDmatrix2.txt"
#### load the SNP list and cis-SNP number for each gene in a region
snplist=read.table("./simulation/summary/pre_process/snplist.txt")$V1
#### pre-process the file to a list including the beta vector, corresponding se vector and LD matrix from GWAS data
convert <- pre_process_twostage(GWASfile, GWASLDfile, snplist)
betay <- as.matrix(convert$beta)
se_betay <- as.matrix(convert$se)
Sigma <- convert$LDmatrix
```

#### Step 3: Read other required data.
```r
#### load the sample size from GWAS data
n=5000
#### load the gene name vector.
gene=c("CCNH", "COX7C", "RASA1", "TMEM161B")
``` 

#### Step 4: Conditional fine-mapping for TWAS analysis.
```r
library(GIFT)
result<-GIFT_two_stage_summ(betax, betay, se_betay, Sigma, n, gene)
```
The result is a data frame including the z-scores and p values for each gene in a focal region. 
```r
      gene          z         p
1     CCNH  0.9994516 0.3175760
2    COX7C -1.1635643 0.2446006
3    RASA1  1.5391589 0.1237655
4 TMEM161B -0.7940135 0.4271876
```

### Visualization for the GIFT result
We visualize the GIFT result incorporating with the marginal GWAS and TWAS results in a Manhattan plot. Here, we load the analyzed results of GWAS and TWAS directly using the example data. The example data for runing the tutorial can be downloaded in this [page](https://yuanzhongshang.github.io/GIFT/documentation/03_data.html).
```r
#### load the GWAS results (e.g., the GEMMA output)
GWASresult=read.table("./summary/pre_process/gemma/GWAS.assoc.txt",header=T)
GWASresult=GWASresult[,c(2,3,11)]
GWASresult$index="GWAS"
colnames(GWASresult)=c("X","BP","P","index")
#### load the TWAS results (e.g., using the BSLMM weight)
TWASresult=read.table("./simulation/visualization/TWASresult.txt",header=T)
TWASresult$BP=apply(TWASresult[,c(3,4)],1,mean)
TWASresult=TWASresult[,c(1,6,5)]
TWASresult$index="TWAS"
colnames(TWASresult)=c("X","BP","P","index")
#### load the GIFT results
GIFTresult=result
GIFTresult$BP=TWASresult$BP
GIFTresult=GIFTresult[,c(1,4,3)]
GIFTresult$index="GIFT"
colnames(GIFTresult)=c("X","BP","P","index")
#### visualize the result by Manhattan plot
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

## Real data application
We used GIFT to perform the condition TWAS fine-mapping in a region on chr 9 (107,581,749-109,298,754) for HDL. This region includes four genes: NIPSNAP3A, NIPSNAP3B, ABCA1, SLC44A1, FSD1L, FKTN, TAL2, and TMEM38B. The data for runing the tutorial can be downloaded in this [page](https://yuanzhongshang.github.io/GIFT/documentation/03_data.html).
```r
#### load the required data
load(/realdata/.realdataRdata)
#### perform conditional fine-mapping for TWAS analysis
library(GIFT)
library(parallel)
result <- GIFT_individual(X, Y, Zx, Zy, gene, pindex, maxiter=1000, tol=1e-4, ncores=8)
result
       gene causal_effect             p
1     ABCA1 -1.857260e+00 2.938669e-179
2 NIPSNAP3A -1.500607e+01  1.971649e-01
3      TAL2 -3.555218e-02  4.157664e-01
4     FSD1L  7.122696e-04  1.000000e+00
5   SLC44A1 -6.702161e-02  1.171462e-01
6   TMEM38B -2.780012e-03  7.819523e-01
7 NIPSNAP3B -2.223423e+00  8.610865e-01
8      FKTN  3.041027e-02  9.781747e-01
#### visualize the result by Manhattan plot
GIFTresult=result
GIFTresult$BP=TWASresult$BP
GIFTresult=GIFTresult[,c(1,4,3)]
GIFTresult$index="GIFT"
colnames(GIFTresult)=c("X","BP","P","index")
data=rbind(GWASresult,TWASresult)
data=rbind(data,GIFTresult)
data$BP=data$BP/1000000
data$index=factor(data$index,levels=c("GWAS","TWAS","GIFT"))

library(ggplot2)
ggplot(data) +
  labs(x="Position on Chr9(Mb)", y=expression(paste(-log[10]," (p-value)"))) +
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
Here is a real data output:
![GIFT\_pipeline](visualization_real.png)

After performing the conditional fine-mapping analysis for HDL across all GWAS risk regions with TWAS significant genes, GIFT produces calibrated p-values for the conditional TWAS tests, see the quantile–quantile plot of -log10 p-values below. In addition, we also provided a Manhattan plot to show different levels of evidence: a gene is “Known” (red) if its association with the trait has been previously reported and well documented; a gene is a significant “TWAS” gene (blue) if its marginal TWAS p-value is below the Bonferroni corrected transcriptome-wide threshold; a genes is significant “GWAS” gene (purple) if its marginal GWAS p-value is below the usual genome-wide threshold 5×10-8 or previously reported; otherwise, a gene is denoted as “NA” (brown). Details are in the Supplementary Tables 4.
![GIFT\_pipeline](realdata.png)
