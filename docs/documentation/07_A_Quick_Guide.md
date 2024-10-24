---
layout: page
title: A Quick Guide
---

GIFT is a Gene-based Integrative Fine-mapping for performing conditional TWAS analysis. GIFT examines one genomic region at a time, jointly models the GReX of all genes residing in the focal region, and carries out TWAS conditional analysis in a maximum likelihood framework. We used LD blocks defined by LDetect: [hg19] (https://github.com/yuanzhongshang/GIFT/tree/main/reproduce/LDetect), [hg38](https://github.com/mancusolab/ma-focus/tree/master/pyfocus/data/ld_blocks). For each region, we focus on protein-coding genes and long intergenic noncoding RNAs (lincRNAs) that are annotated in GENCODE (release 12). The latest release can also be downloaded [here](https://www.gencodegenes.org/human/). Other definitions can also be used to determine LD blocks. 

Data Availability
-------------------
The commonly used eQTL summary statistics includes: genome-wide eQTL summary statistics from GEUVADIS data in [dropbox](https://www.dropbox.com/scl/fo/4nqcmkblerspfmva5stwf/ANHZU_kX2AlveEEbx9DKbZU?rlkey=qjcxprlk83t7pw8ka2ne2v4w9&dl=0), cis-eQTL and trans-eQTL summary statistics from [eQTLGen Consortium](https://www.eqtlgen.org/phase1.html), and cis-eQTL mapping summary statistics for African American and European American from [GENOA](https://xiangzhou.github.io/resources/).

The commonly used GWAS summary statistics can be downloaded in [GWAS Catalog](https://www.ebi.ac.uk/gwas/), [GWAS ATLAS](https://atlas.ctglab.nl/), [IEU Open GWAS Project](https://gwas.mrcieu.ac.uk/), [FinnGen](https://www.finngen.fi/en/access_results), [Biobank Japan](https://pheweb.jp/) and [Neale Lab UK Biobank GWAS results](https://www.nealelab.is/uk-biobank).

The commonly used reference panel are [1000 Genomes project data](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502), and [UK Biobank LD reference](https://uchicago.app.box.com/s/jqocacd2fulskmhoqnasrknbt59x3xkn).

Data Preparation
-------------------
Since summary statistics may be available in different formats, we provide a guideline for data preparation. 
  * First, ensure the same genome builds for three datasets including eQTL, GWAS, reference panel ([code for lifting over genome builds](https://github.com/yuanzhongshang/GIFT/issues/12)).
  * Second, perform quality control on each dataset, including removing strand ambiguous variants, multi-allelic variants, and SNPs with MAF below the threshold.
  * Third, obtain the overlapping SNP list with alleles in the target region. If summary statistics are only available for each gene in the cis-region, combine the set of cis-SNPs for each gene and overlap them with the SNP list from the GWAS dataset and reference panel.
  * Fourth, extract SNPs from the SNP list, align the reference alleles, and switch the sign of beta or z score based on the allele. Note that, missing SNPs effect on gene in the eQTL dataset will automatically imputed with zeros in the `pre_process_summary` function. 
Following these steps, you can obtain data in the format like [here](https://yuanzhongshang.github.io/GIFT/documentation/03_data.html#gift-using-summary-statistics-as-input), and proceed to pre-process the data into GIFT input format.

Running GIFT with Summary Statistics
-------------------
We recommend performing the GIFT analysis using the summary statistic version. If you have individual-level data, please [convert it into into the corresponding summary statistics](https://yuanzhongshang.github.io/GIFT/documentation/05_analysis_reproduce.html#convert-the-individual-level-data-into-the-summary-statistics). First, ensure GIFT is successfully installed. For more details, refer to the (installation guide)[https://yuanzhongshang.github.io/GIFT/documentation/02_installation.html].
```r
#devtools::install_github('yuanzhongshang/GIFT')
library(GIFT)
```
#### Step 1: Pre-process the summary statistics with different formats.
In this guild, we use the simulations based on the realistic genotypes from GEUVADIS (n1=465) and UK Biobank (n2=5,000) in a region on chr 5 as an example. This region includes four genes: RASA1, COX7C, CCNH and TMEM161B. We set RASA1 as the causal gene with the effect size being sqrt(0.1).Here, we convert different summary statistics (plink (.qassoc), GEMMA (.assoc.txt) and SAIGE (.txt)) and LD matrix data formats (matrix (.txt) or h5 format (.h5)) to GIFT inputs. 
```r
library(GIFT)
dir <- getwd()

#### load the directory containing files of summary statistics from eQTL data only (e.g., the SAIGE output)
eQTLfilelocation <- paste0(dir, "/example/simulation/summary/pre_process/saige/eQTL")

#### load the directory of summary statistics from GWAS data (e.g., the SAIGE output)
GWASfile <- paste0(dir,"/example/simulation/summary/pre_process/saige/GWAS.txt")

#### load the directory of LD matrix from eQTL data and GWAS data (e.g., a long format: h5 format)
eQTLLDfile <- paste0(dir, "/example/simulation/summary/pre_process/LDmatrix1.h5")
GWASLDfile <- paste0(dir, "/example/simulation/summary/pre_process/LDmatrix2.h5")

#### load the SNP list and cis-SNP number for each gene in a region
snplist <- read.table(paste0(dir, "/example/simulation/summary/pre_process/snplist.txt"))$V1
pindex <- c(41, 23, 63, 96)

#### pre-process the file to be a list including gene names vector, z-score matrix and LD matrix of eQTL data and GWAS data
convert <- pre_process_summary(eQTLfilelocation, eQTLLDfile, GWASfile, GWASLDfile, snplist, pindex)
gene <- convert$gene
Zscore1 <- convert$Zscore1
Zscore2 <- convert$Zscore2
Sigma1 <- convert$LDmatrix1
Sigma2 <- convert$LDmatrix2
n1 <- convert$n1
n2 <- convert$n2
```

### Step 2: Perform conditional fine-mapping for TWAS analysis.
```r
result <- GIFT_summary(Zscore1, Zscore2, Sigma1, Sigma2, n1, n2, gene, pindex, R=NULL, maxiter=100, tol=1e-3, pleio=0, ncores=1, in_sample_LD=T, filter=T, split=5)
```
The result is a data frame including the causal effect estimates and p values for each gene in a focal region. 
```r
result
      gene causal_effect            p
1     CCNH   0.002940875 9.712126e-01
2    COX7C   0.028515942 7.435854e-01
3    RASA1   0.364477685 1.337356e-05
4 TMEM161B  -0.034292965 3.134303e-01
```
Note that, the summary statistics version of GIFT often requires the in-sample LD matrix. If the in-sample LD matrix is not available, it can be also calculated from the reference panel data. It would be better to ensure the ethnicity of the reference panel is consistent with that of the analyzed data, details in [here](https://yuanzhongshang.github.io/GIFT/documentation/06_Summary_statistic_issues.html). If in-sample LD was not used, the LD matrix is regularized to be (1-s1)\*Sigma1+s1\*E and (1-s2)\*Sigma2+s2\*E where both s1 and s2 are estimated by [estimate_s_rss](https://stephenslab.github.io/susieR/reference/estimate_s_rss.html) in susieR. A grid search algorithm is performed over the range from 0.1 to 1 once the estimation from susieR does not work well. 
```r
### load the LD matrix from 1,000 Genomes project
LD <- as.matrix(read.table(paste0(dir, "/example/simulation/summary/LDmatrix10000G.txt")))

result <- GIFT_summary(Zscore1, Zscore2, LD, LD, n1, n2, gene, pindex, R=NULL, maxiter=100, tol=1e-3, pleio=0, ncores=1, in_sample_LD=F, filter=T, split=5)
result
      gene causal_effect            p
1     CCNH   0.017590501 8.296347e-01
2    COX7C   0.007763401 8.633912e-01
3    RASA1   0.319433644 7.881041e-05
4 TMEM161B  -0.064599635 3.086057e-02
```