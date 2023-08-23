---
layout: page
title: Data Input
description: ~
---
The following are the links for the example dataset used in the GIFT. 

## Simulation
### GIFT: Using individual-level data as input:
  Different genotype data formats:
  * [plink binary (.bim/.fam./.bed)](https://github.com/yuanzhongshang/GIFT/tree/main/example/simulation/individual/pre_process/plink_binary)
  * [vcf](https://github.com/yuanzhongshang/GIFT/tree/main/example/simulation/individual/pre_process/vcf)
  * [ped/map](https://github.com/yuanzhongshang/GIFT/tree/main/example/simulation/individual/pre_process/pedmap)
  * [csv](https://github.com/yuanzhongshang/GIFT/tree/main/example/simulation/individual/pre_process/csv)
  * [tsv](https://github.com/yuanzhongshang/GIFT/tree/main/example/simulation/individual/pre_process/tsv)
 
 Other dataset:
  * [individual_data.Rdata](https://github.com/yuanzhongshang/GIFT/blob/main/example/simulation/individual/individual_data.Rdata)
  
### GIFT: Using summary statistics as input:
  Different summary statistics data formats:
  * [plink binary (.bim/.fam./.bed)](https://github.com/yuanzhongshang/GIFT/tree/main/example/simulation/summary/pre_process/plink)
  * [GEMMA](https://github.com/yuanzhongshang/GIFT/tree/main/example/simulation/summary/pre_process/gemma)
  * [SAIGE](https://github.com/yuanzhongshang/GIFT/tree/main/example/simulation/summary/pre_process/saige)

  Different LD matrix data formats:
  * [matrix format from eQTL data](https://github.com/yuanzhongshang/GIFT/blob/main/example/simulation/summary/pre_process/LDmatrix1.txt)
  * [matrix format from GWAS data](https://github.com/yuanzhongshang/GIFT/blob/main/example/simulation/summary/pre_process/LDmatrix2.txt)
  * [long format from eQTL data (.h5)](https://github.com/yuanzhongshang/GIFT/blob/main/example/simulation/summary/pre_process/LDmatrix1.h5)
  * [long format from GWAS data (.h5)](https://github.com/yuanzhongshang/GIFT/blob/main/example/simulation/summary/pre_process/LDmatrix2.h5)

  Other dataset:
  * [SNP list](https://github.com/yuanzhongshang/GIFT/blob/main/example/simulation/summary/pre_process/snplist.txt)
  * [estimated correlated matrix of gene expressions](https://github.com/yuanzhongshang/GIFT/blob/main/example/simulation/summary/R.txt)
  
### Two-stage version of GIFT: Using pre-trained weights and summary statistics as input:
  eQTL weight:
  * [the BLUP weight](https://github.com/yuanzhongshang/GIFT/tree/main/example/simulation/two_stage/weights)

  Different summary statistics data formats:
  * [plink binary (.bim/.fam./.bed)](https://github.com/yuanzhongshang/GIFT/blob/main/example/simulation/summary/pre_process/plink/GWAS.qassoc)
  * [GEMMA](https://github.com/yuanzhongshang/GIFT/blob/main/example/simulation/summary/pre_process/gemma/GWAS.assoc.txt)
  * [SAIGE](https://github.com/yuanzhongshang/GIFT/blob/main/example/simulation/summary/pre_process/saige/GWAS.txt)  
  
### For the visualization:
  * [TWAS result](https://github.com/yuanzhongshang/GIFT/blob/main/example/simulation/visualization/TWASresult.txt)

## Real data application
  * [realdata.RData](https://github.com/yuanzhongshang/GIFT/blob/main/example/realdata/realdata.RData)
