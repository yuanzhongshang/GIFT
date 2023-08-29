---
layout: page
title: Example Analysis
description: ~
---
GIFT examines one genomic region at a time, explicitly models the gene expression correlation and cis-SNP LD across different genes in the focal region and accounts for the uncertainty in the constructed GReX, and carries out TWAS conditional analysis in a maximum likelihood framework. This tutorial is the example analysis with GIFT. Before running the tutorial, make sure that the GIFT package is successfully installed. Please see the [link](https://yuanzhongshang.github.io/GIFT/documentation/02_installation.html) for the installation instructions. The example data for the tutorial can be downloaded in this [page](https://yuanzhongshang.github.io/GIFT/documentation/03_data.html).

## A simulated data example
We conducted the simulations based on the realistic genotypes from GEUVADIS (n1=465) and UK Biobank (n2=5,000) in a region on chr 5. This region includes four genes: RASA1, COX7C, CCNH and TMEM161B. We set RASA1 as the causal gene with the effect size being sqrt(0.1).

### GIFT: Using individual-level data as input
The function `GIFT_individual` is the main function for GIFT with individual-level data. The essential inputs are:
- X: The standardized gene expression matrix for all genes in a specific region from eQTL data. 
- Y: The standardized trait vector from GWAS data.
- Zx: The standardized cis-genotype matrix for all genes in a specific region from eQTL data.
- Zy: The standardized cis-genotype matrix for all genes in a specific region from GWAS data.
- gene: The gene name vector represents the genes in a specific region, the order of the gene name should be consistent with that in X.
- pindex: A vector with each element representing the number of cis-SNPs for each gene.

The optional inputs are:
- maxiter: The user-defined maximum iteration, with the default to be 1000.
- tol: The user-defined convergence tolerance of the absolute value of the difference between the nth and (n+1)th log likelihood, with the default value as 1e-4. 
- ncores: The number of cores used in analysis, with the default to be 1. The analysis will be performed with parallel computing once the number of cores is greater than 1. Of note, the incorporated function mclapply() depends on another R package "parallel" in Linux. 

#### Step 1: Pre-process the genotype data with different formats.
The function `pre_process_individual` is able to convert different genotype data formats to GIFT inputs. In particular, this function is flexible to handle plink binary format (.bim/.fam./.bed), vcf, ped/map format, csv, and tsv file. Here, we take various genotype data formats from GEUVADIS data in [page](https://yuanzhongshang.github.io/GIFT/documentation/03_data.html) for example. Of note, in this step, cis-genotype matrix has been standardized to have a mean of zero and standard derivation of one. 
```r
library(GIFT)
#### load the directory containing the files to be processed only (e.g., plink binary format)
filelocation <- "./simulation/individual/pre_process/plink_binary"
#### load the directory of plink exe file
plinkexe <- "plink"
#### pre-process the file to be a list including gene names vector, cis-genotype matrix and pindex
convert <- pre_process_individual(filelocation, plinkexe)
gene <- convert$gene
Zx <- convert$Z
pindex <- convert$pindex
```

#### Step 2: Read other datasets.
The required example data can be downloaded in this [page](https://yuanzhongshang.github.io/GIFT/documentation/03_data.html). 
```r
#### load the Rdata file containing X, Zy and Y
load("./simulation/individual/individual_data.Rdata")
```

#### Step 3: Perform conditional fine-mapping for TWAS analysis.
```r
result <- GIFT_individual(X, Y, Zx, Zy, gene, pindex, maxiter=1000, tol=1e-4, ncores=1)
```
The result is a data frame including the causal effect estimates and p values for each gene within a focal region. 
```r
result
      gene causal_effect            p
1     CCNH    0.01209639 8.594982e-01
2    COX7C    0.02111053 8.028104e-01
3    RASA1    0.35309347 7.131827e-06
4 TMEM161B   -0.03306309 3.257628e-01
```

### GIFT: Using summary statistics as input
The function `GIFT_summary` is the main function for GIFT with summary statistics. The essential inputs are:
- Zscore_1: Zscore matrix of the cis-SNP effect size for all genes in a specific region from eQTL data.
- Zscore_2: Zscore vector of the cis-SNP effect size for all genes in a specific region from GWAS data.
- Sigma1: LD matrix from eQTL data.
- Sigma2: LD matrix from GWAS data, both Sigma1 and Sigma2 are often from the same reference panel.
- R: Estimated correlation matrix of gene expressions.
- n1: Sample size of eQTL data.
- n2: Sample size of GWAS data.
- gene: The gene name vector, the order of the gene name should be consistent with that in Zscore_1.
- pindex: A vector with each element represents the number of cis-SNPs for each gene.

The optional inputs are:
- maxiter: The user-defined maximum iteration, with the default to be 1000.
- tol: The user-defined convergence tolerance of the absolute value of the difference between the nth and (n+1)th log likelihood, with the default value as 1e-4. 
- ncores: The number of cores used in analysis, with the default to be 1. The analysis will be performed with parallel computing once the number of cores is greater than 1. Of note, the incorporated function mclapply() depends on another R package "parallel" in Linux. 

#### Step 1: Pre-process the summary statistics with different formats.
The function `pre_process_summary` is able to convert different summary statistics and LD matrix data formats to GIFT inputs. In particular, this function is flexible to handle association test output from plink (.qassoc), GEMMA (.assoc.txt) and SAIGE (.txt). While, this function is also flexible to handle LD matrix either from matrix or a long format such as h5 format. Here, we take various data formats from example data in [page](https://yuanzhongshang.github.io/GIFT/documentation/03_data.html). Note that, the summary statistics version of GIFT often requires the in-sample LD matrix. If the in-sample LD matrix is not available, it can be also calculated from the reference panel data (e.g., 1,000 Genomes project). It would be better to ensure the ethnicity of the reference panel is consistent with that of the analyzed data. 
```r
#### load the directory containing files of summary statistics from eQTL data only (e.g., the SAIGE output)
eQTLfilelocation <- paste0(getwd(), "/simulation/summary/pre_process/saige/eQTL")
#### load the directory of summary statistics from GWAS data (e.g., the SAIGE output)
GWASfile <- "./simulation/summary/pre_process/saige/GWAS.txt"
#### load the directory of LD matrix from eQTL data and GWAS data (e.g., a long format: h5 format)
eQTLLDfile <- "./simulation/summary/pre_process/LDmatrix1.h5"
GWASLDfile <- "./simulation/summary/pre_process/LDmatrix2.h5"
#### load the SNP list and cis-SNP number for each gene in a region
snplist <- read.table("./simulation/summary/pre_process/snplist.txt")$V1
pindex <- c(41, 23, 63, 96)
#### pre-process the file to be a list including gene names vector, z-score matrix and LD matrix of eQTL data and GWAS data
convert <- pre_process_summary(eQTLfilelocation, eQTLLDfile, GWASfile, GWASLDfile, snplist, pindex)
gene <- convert$gene
Zscore1 <- convert$Zscore1
Zscore2 <- convert$Zscore2
LDmatrix1 <- convert$LDmatrix1
LDmatrix2 <- convert$LDmatrix2
```

#### Step 2: Read other datasets.
```r
### input the sample sizes of eQTL data and GWAS data
n1 <- 465
n2 <- 5000
### load the estimated correlated matrix of gene expressions
setwd(gsub("/simulation/summary/pre_process/saige/eQTL", "", getwd()))
R <- as.matrix(read.table("./simulation/summary/R.txt"))
```

#### Step 3: Perform conditional fine-mapping for TWAS analysis.
```r
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
The function `GIFT_two_stage_summ` is developed for conditional fine-mapping in TWAS with pre-trained weights and summary statistics. The two-stage version of GIFT is not only computationally efficient but also allows us to make use of existing gene expression prediction models for more convenient TWAS fine-mapping. The example data for runing the tutorial can be downloaded in this [page](https://yuanzhongshang.github.io/GIFT/documentation/03_data.html). The essential inputs are:
- betax: Weight matrix for the cis-SNP effect size from eQTL data.
- betay: Beta vector of the cis-SNP effect size vector from GWAS data.
- se_betay: Se vector of the cis-SNP effect size vector from GWAS data.
- Sigma: LD matrix from GWAS data.
- n: Sample size of GWAS data.
- gene: The vector of gene names.
  
#### Step 1: Read the eQTL weight.
Gene expression prediction is the key for two-stage TWAS methods. The commonly used prediction models include lasso and elastic net (enet) as implemented in prediXcan, Best Linear Unbiased Prediction (BLUP), the top SNPs (top1) and Bayesian sparse linear mixed model (BSLMM) as implemented in TWAS/FUSION, latent Dirichlet process regression (DPR) as implemented in both DPR and TIGAR. For a specific region, the weights from all genes can be represented to be a block diagonal matrix. The function `weightconvert` is able to convert a list including weights for multiple genes into a required block diagonal matrix for GIFT.
```r
#### load the weight matrix from the eQTL data (e.g., BLUP)
setwd("./simulation/two_stage/weights")
CCNH <- as.matrix(read.table("CCNHweight.txt"))
COX7C <- as.matrix(read.table("COX7Cweight.txt"))
RASA1 <- as.matrix(read.table("RASA1weight.txt"))
TMEM161B <- as.matrix(read.table("TMEM161Bweight.txt"))
weightlist <- list(CCNH = CCNH, COX7C = COX7C, RASA1 = RASA1, TMEM161B = TMEM161B)
#### convert the weights from all genes into a block diagonal matrix
betax <- weightconvert(weightlist)
```

#### Step 2: Read the beta vector, corresponding se vector and LD matrix from GWAS data.
The function `pre_process_twostage` is able to convert common summary statistics and LD matrix data formats to GIFT inputs. Specifically, this function is flexible to handle output from plink (.qassoc), GEMMA (.assoc.txt) and SAIGE (.txt). Meanwhile, this function is also flexible to handle LD matrix either from a matrix or a long format such as h5 format. The example data is the same as above in [page](https://yuanzhongshang.github.io/GIFT/documentation/03_data.html). Note that, the two-stage version of GIFT often requires the in-sample LD matrix. If the in-sample LD matrix is not available, it can be also calculated from the reference panel data (e.g., 1,000 Genomes project). It would be better to ensure the ethnicity of the reference panel is consistent with that of the analyzed data. 
```r
#### load the directory of summary statistics from GWAS data (e.g., the plink output)
GWASfile <- "./simulation/summary/pre_process/plink/GWAS.qassoc"
#### load LD matrix from eQTL data and GWAS data (e.g., a matrix)
GWASLDfile <- "./simulation/summary/pre_process/LDmatrix2.txt"
#### load the SNP list and cis-SNP number for each gene in a region
snplist <- read.table("./simulation/summary/pre_process/snplist.txt")$V1
#### pre-process the file to be a list including the beta vector, corresponding se vector and LD matrix from GWAS data
convert <- pre_process_twostage(GWASfile, GWASLDfile, snplist)
betay <- as.matrix(convert$beta)
se_betay <- as.matrix(convert$se)
Sigma <- convert$LDmatrix
```

#### Step 3: Read other datasets.
```r
#### load the sample size from GWAS data
n=5000
#### load the gene name vector.
gene=c("CCNH", "COX7C", "RASA1", "TMEM161B")
``` 

#### Step 4: Conditional fine-mapping for TWAS analysis.
```r
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
The GIFT result can be visualized with the marginal GWAS and TWAS results in a Manhattan plot. Here, we load the results from GWAS and TWAS directly using the example data. The example data for runing the tutorial can be downloaded in this [page](https://yuanzhongshang.github.io/GIFT/documentation/03_data.html).
```r
#### load the GWAS results (e.g., the GEMMA output)
GWASresult=read.table("./simulation/summary/pre_process/gemma/GWAS.assoc.txt",header=T)
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

## A real data example
We used GIFT to perform the condition TWAS fine-mapping in a region on chr 9 (107,581,749-109,298,754) for HDL. This region includes eight genes: NIPSNAP3A, NIPSNAP3B, ABCA1, SLC44A1, FSD1L, FKTN, TAL2, and TMEM38B. The data for runing the tutorial can be downloaded in this [page](https://yuanzhongshang.github.io/GIFT/documentation/03_data.html).
```r
#### load the required data
load("./realdata/realdata.Rdata")
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

After performing the conditional fine-mapping analysis for HDL across all GWAS risk regions with TWAS significant genes, GIFT produces calibrated p-values for the conditional TWAS tests, please see the quantile–quantile plot of -log10 p-values below. A Manhattan plot can be also provided to show different levels of evidence: a gene is “Known” (red) if its association with the trait has been previously reported and well documented; a gene is a significant “TWAS” gene (blue) if its marginal TWAS p-value is below the Bonferroni corrected transcriptome-wide threshold; a genes is significant “GWAS” gene (purple) if its marginal GWAS p-value is below the usual genome-wide threshold 5×10-8 or previously reported; otherwise, a gene is denoted as “NA” (brown). Details are in the Supplementary Tables 4.
![GIFT\_pipeline](realdata.png)
