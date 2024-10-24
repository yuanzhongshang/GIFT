#' @title The pro-process function converts common genotype formats to match GIFT input using individual-level data.
#' @description pre_process_individual is flexible to handle plink binary format (.bim/.fam./.bed), vcf, ped/map format, csv, and tsv file.
#' @param filelocation A file location only contains the genotype files to be processed.
#' @param plinkexe A path for the executable plink software.
#' @return A list including the gene names, z-scores and cis-SNP numbers for each gene. 

pre_process_individual <- function(filelocation, plinkexe="plink"){
  setwd(filelocation)
  file <- list.files(filelocation)
  gene <- unique(sub("\\.[^.]*$", "", file))
  
  if(length(grep(".ped",file)) > 0){
    system("mkdir ./plink_binary")
    setwd("./plink_binary") 
    for(i in 1:length(gene)){
      system(paste0(plinkexe," --file ",filelocation,"/",gene[i]," --make-bed --out ",gene[i]))
    }
  }
  file <- list.files()
  if(length(grep(".bim",file)) > 0){
    k <- file[grep(".bim",file)]
    k <- gsub(".bim","",k)
    if (!requireNamespace("BEDMatrix", quietly = TRUE)) {
      install.packages("BEDMatrix")
    }
    library("BEDMatrix")
    Z=NULL
    p=NULL
    for(i in 1:length(k)){
      Z1<-BEDMatrix(k[i])
      Z<-cbind(Z,as.matrix(Z1))
      p[i] <- dim(Z1)[2]
    }
  }
  
  
  if(length(grep("vcf",file)) > 0){
    if (!requireNamespace("VariantAnnotation", quietly = TRUE)) {
      if (!require("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
      BiocManager::install("VariantAnnotation")
    }
    library(VariantAnnotation)
    Z=NULL
    p=NULL
    for(i in 1:length(file)){
      Z1 <-readVcf(file[i])
      Z1 <- geno(Z1)
      Z1 <- t(as.matrix(Z1@listData$GT))
      Z1[Z1 == "0/0"] <- 0
      Z1[Z1 == "0/1"] <- 1
      Z1[Z1 == "1/0"] <- 1
      Z1[Z1 == "1/1"] <- 2
      Z1[Z1 == "./."] <- NA
      Z1 <- matrix(as.numeric(Z1),dim(Z1)[1],dim(Z1)[2])
      Z<-cbind(Z,as.matrix(Z1))
      p[i] <- dim(Z1)[2]
    }
  }
  
  if(length(grep(".csv",file)) > 0){
    Z=NULL
    p=NULL
    for(i in 1:length(file)){
      Z1 <- read.csv(file[i], header = TRUE, row.names = 1)
      Z <- cbind(Z,as.matrix(Z1))
      p[i] <- dim(Z1)[2]
    }
  }
  
  if(length(grep(".tsv",file)) > 0){
    Z=NULL
    p=NULL
    for(i in 1:length(file)){
      Z1 <- read.table(file[i], header = TRUE, sep = "\t", row.names = 1)
      Z <- cbind(Z,as.matrix(Z1))
      p[i] <- dim(Z1)[2]
    }
  }
  Z <- scale(na.omit(Z))
  r <- list(gene = gene, Z = Z, pindex = p)
  return(r)
}


#' @title The pro-process function converts common summary statistics and LD matrix format to match GIFT input using summary statistics.
#' @description pre_process_summary flexible to handle association test output from plink (.qassoc), GEMMA (.assoc.txt) and SAIGE (.txt); and LD matrix either from matrix or a long format such as h5 format. 
#' @param eQTLfilelocation A file location only contains the summary statistics files from eQTL data.
#' @param eQTLLDfile A file path for the LD matrix from eQTL data.
#' @param GWASfile A file path for the summary statistics from GWAS data.
#' @param GWASLDfile A file path for the LD matrix from GWAS data.
#' @param snplist A vector represents the cis-SNP list for all genes in a region.
#' @param pindex A vector with each element represents the number of cis-SNPs for each gene.
#' @return A list including the gene names, cis-SNP numbers for each gene, z-scores from eQTL data and GWAS data, LD matrix from eQTL data and GWAS data. 

pre_process_summary <- function(eQTLfilelocation, eQTLLDfile, GWASfile, GWASLDfile, snplist, pindex){
  
  setwd(eQTLfilelocation)
  eQTLfile <- list.files(eQTLfilelocation)
  eQTLgene <- gsub("eQTL","",eQTLfile)
  split_character <- "."
  gene <- unique(sub("\\.[^.]*$", "", eQTLgene))
  pindexcum <- c(0,cumsum(pindex))
  n1 <- NULL
  
  cat("Processing the eQTL data:\n")
  eQTLz <- NULL
  for(i in 1:length(gene)){
    cat(gene[i],":\n")
    eQTLdata=read.table(eQTLfile[i],header = T)
    eQTLz1 <- rep(0,length(snplist))
    ##plink format
    if(sum(colnames(eQTLdata) %in% "T")!=0){
      cat("Reading the plink format summary statistics: ")
      eQTLz1 <- eQTLdata[match(snplist,eQTLdata$SNP),]$T
	  n1[i] <- median(eQTLdata[match(snplist,eQTLdata$SNP),]$NMISS)
      if(sum(is.na(eQTLz1))>0){
		cat("You input only the cis-SNP for each gene, setting the z-scores of SNPs other than the cis-SNP for each gene to zero. ")
		tmp=eQTLz1[(pindexcum[i]+1):pindexcum[i+1]]
		eQTLz1 <- rep(0,length(snplist))
		eQTLz1[(pindexcum[i]+1):pindexcum[i+1]]=tmp
      }
      cat("Success!\n")
    }
    ##gemma format
    if(sum(colnames(eQTLdata) %in% "beta")!=0){
      cat("Reading the gemma format summary statistics: ")
      eQTLz1 <- eQTLdata[match(snplist,eQTLdata$rs),]$beta/eQTLdata[match(snplist,eQTLdata$rs),]$se
	  n1[i] <- median(eQTLdata[match(snplist,eQTLdata$rs),]$n_obs)
      if(sum(is.na(eQTLz1))>0){
  	cat("You input only the cis-SNP for each gene, setting the z-scores of SNPs other than the cis-SNP for each gene to zero. ")
    	tmp=eQTLz1[(pindexcum[i]+1):pindexcum[i+1]]
    	eQTLz1 <- rep(0,length(snplist))
    	eQTLz1[(pindexcum[i]+1):pindexcum[i+1]]=tmp
      }
      cat("Success!\n")
    }
    ##SAIGE format
    if(sum(colnames(eQTLdata) %in% "Tstat")!=0){
      cat("Reading the SAIGE format summary statistics: ")
      eQTLz1 <- eQTLdata[match(snplist,eQTLdata$MarkerID),]$BETA/eQTLdata[match(snplist,eQTLdata$MarkerID),]$SE
	  n1[i] <- median(eQTLdata[match(snplist,eQTLdata$MarkerID),]$N)
      if(sum(is.na(eQTLz1))>0){
    	 cat("You input only the cis-SNP for each gene, setting the z-scores of SNPs other than the cis-SNP for each gene to zero. ")
    	 tmp=eQTLz1[(pindexcum[i]+1):pindexcum[i+1]]
    	 eQTLz1 <- rep(0,length(snplist))
    	 eQTLz1[(pindexcum[i]+1):pindexcum[i+1]]=tmp
       }
       cat("Success!\n")
    }
    eQTLz <- cbind(eQTLz,eQTLz1)
  }
  n1 <- median(n1)
  if(sum(is.na(eQTLz))>0){
    cat("Unmatched SNP list. Please check your SNP list. \n")
  }
  
  cat("Processing the GWAS data:\n")
  GWASdata <- read.table(GWASfile,header=T)
  GWASz <- NULL
  n2 <- NULL
  
  ##plink format
  if(sum(colnames(GWASdata) %in% "T")!=0){
    cat("Reading the plink format summary statistics: ")
    GWASz <- GWASdata[match(snplist,GWASdata$SNP),]$T
	n2 <- median(GWASdata[match(snplist,GWASdata$SNP),]$NMISS)
    cat("Success!\n")
  }  
  ##gemma format
  if(sum(colnames(GWASdata) %in% "beta")!=0){
    cat("Reading the gemma format summary statistics: ")
    GWASz <- GWASdata[match(snplist,GWASdata$rs),]$beta/GWASdata[match(snplist,GWASdata$rs),]$se
	n2 <- median(GWASdata[match(snplist,GWASdata$rs),]$n_obs)
    cat("Success!\n")
  }
  ##saige format
  if(sum(colnames(GWASdata) %in% "Tstat")!=0){
    cat("Reading the SAIGE format summary statistics: ")
    GWASz <- GWASdata[match(snplist,GWASdata$MarkerID),]$BETA/GWASdata[match(snplist,GWASdata$MarkerID),]$SE
	n2 <- median(GWASdata[match(snplist,GWASdata$MarkerID),]$N)
    cat("Success!\n")
  }
  if(sum(is.na(GWASz))>0){
    cat("Unmatched SNP list. Please check your SNP list. \n")
  }
  
  cat("Processing the LD matrix:\n")
  if(sum(grep(".txt",eQTLLDfile))!=0){ 
    cat("Reading the txt format LD matrix from eQTL data: ")
    LDmatrix1 <- read.table(eQTLLDfile)
  }
  if(sum(grep(".h5",eQTLLDfile))!=0){     
    if (!requireNamespace("rhdf5", quietly = TRUE)) {
      if (!require("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
        BiocManager::install("rhdf5")
    }
    library(rhdf5)
    cat("Reading the h5 format LD matrix from eQTL data: ")
    LDmatrix1 <- h5read(eQTLLDfile, "matrix")
  }
  LDmatrix1 <- as.matrix(LDmatrix1)
  if(dim(LDmatrix1)[1]!=length(snplist)){
    cat("Converting the LD matrix of the unique SNP list to the LD matrix of the input SNP list.")
    LDmatrix1 <- LDmatrix1[match(snplist,colnames(LDmatrix1)),match(snplist,colnames(LDmatrix1))]
    cat("Success!\n")
  }else{
    cat("Success!\n")
  }
  if(sum(is.na(LDmatrix1))>0){
    cat("LD matrix cannot contain missing values or missing column names. Please check your LD matrix. \n")
  }
  
  if(sum(grep(".txt",GWASLDfile))!=0){  
    cat("Reading the txt format LD matrix from GWAS data: ")  
    LDmatrix2 <- read.table(GWASLDfile)
  }
  if(sum(grep(".h5",GWASLDfile))!=0){     
    if (!requireNamespace("rhdf5", quietly = TRUE)) {
      if (!require("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
        BiocManager::install("rhdf5")
    }
    library(rhdf5)
    cat("Reading the h5 format LD matrix from GWAS data: ")  
    LDmatrix2 <- h5read(GWASLDfile, "matrix")
  }  
  LDmatrix2 <- as.matrix(LDmatrix2)
  if(dim(LDmatrix2)[1]!=length(snplist)){
    cat("Converting the LD matrix of the unique SNP list to the LD matrix of the input SNP list.")
    LDmatrix2 <- LDmatrix2[match(snplist,colnames(LDmatrix2)),match(snplist,colnames(LDmatrix2))]
    cat("Success!\n")
  }else{
    cat("Success!\n")
  }
  if(sum(is.na(LDmatrix2))>0){
    cat("LD matrix cannot contain missing values or missing column names. Please check your LD matrix. \n")
  }
  
  r <- list(gene = gene, pindex = pindex, Zscore1 = eQTLz, Zscore2 = GWASz, LDmatrix1 = LDmatrix1, LDmatrix2 = LDmatrix2, n1 = n1, n2 = n2)
  return(r)
}


#' @title The pro-process function converts the weights from all genes into a required input of two-stage version of GIFT.
#' @description weightconvert converts a list including the weights for multiple genes into a required block diagonal matrix. 
#' @param weightlist A list contains the weights for all genes in a region.
#' @return A weight matrix to match the two-stage version of GIFT input. 
weightconvert <- function(weightlist){
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    install.packages("Matrix")
  }
  library(Matrix)  
  weightinput <- as.matrix(bdiag(weightlist))
  return(weightinput)
}


#' @title The pro-process function converts common summary statistics and LD matrix format from GWAS data to match a required input of two-stage version of GIFT.
#' @description pre_process_twostage flexible to handle association test output from plink (.qassoc), GEMMA (.assoc.txt) and SAIGE (.txt); and LD matrix either from matrix or a long format such as h5 format. 
#' @param GWASfile A file path for the summary statistics from GWAS data.
#' @param GWASLDfile A file path for the LD matrix from GWAS data.
#' @param snplist A vector represents the cis-SNP list for all genes in a region.
#' @return A list including the beta vector, corresponding se vector and LD matrix from GWAS data.

pre_process_twostage <- function(GWASfile, GWASLDfile, snplist){

  cat("Processing the GWAS data:\n")
  GWASdata <- read.table(GWASfile,header=T)
  beta <- NULL
  se <- NULL
  ##plink format
  if(sum(colnames(GWASdata) %in% "T")!=0){
    cat("Reading the plink format summary statistics: ")
    beta <- GWASdata[match(snplist,GWASdata$SNP),]$BETA
    se <- GWASdata[match(snplist,GWASdata$SNP),]$SE
    cat("Success!\n")
  }
  ##gemma format
  if(sum(colnames(GWASdata) %in% "beta")!=0){
    cat("Reading the gemma format summary statistics: ")
    beta <- GWASdata[match(snplist,GWASdata$rs),]$beta
    se <- GWASdata[match(snplist,GWASdata$rs),]$se
    cat("Success!\n")
  }
  ##saige format
  if(sum(colnames(GWASdata) %in% "Tstat")!=0){
    cat("Reading the SAIGE format summary statistics: ")
    beta <- GWASdata[match(snplist,GWASdata$MarkerID),]$BETA
    se <- GWASdata[match(snplist,GWASdata$MarkerID),]$SE
    cat("Success!\n")
  }

  cat("Processing the LD matrix:\n")
  if(sum(grep(".txt",GWASLDfile))!=0){
    cat("Reading the txt format LD matrix from eQTL data: ")
    LDmatrix <- read.table(GWASLDfile)
  }
  
  if(sum(grep(".h5",GWASLDfile))!=0){     
    if (!requireNamespace("rhdf5", quietly = TRUE)) {
      if (!require("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
        BiocManager::install("rhdf5")
    }
    library(rhdf5)
    cat("Reading the h5 format LD matrix from eQTL data: ")
    LDmatrix <- h5read(GWASLDfile, "matrix")
  }  
  LDmatrix <- as.matrix(LDmatrix)
  if(dim(LDmatrix)[1]!=length(snplist)){
    cat("Converting the LD matrix of the unique SNP list to the LD matrix of the input SNP list.")
    LDmatrix <- LDmatrix[match(snplist,colnames(LDmatrix)),match(snplist,colnames(LDmatrix))]
    cat("Success!\n")
  }else{
    cat("Success!\n")
  }
  if(sum(is.na(LDmatrix))>0){
    cat("LD matrix cannot contain missing values or missing column names. Please check your LD matrix. \n")
  }
  
  r <- list(beta = beta, se = se, LDmatrix = LDmatrix)
  return(r)
}
