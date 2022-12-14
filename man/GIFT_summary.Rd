% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GIFT_summary.R
\name{GIFT_summary}
\alias{GIFT_summary}
\title{The main function for conditional fine-mapping for in transcriptome-wide association studies with summary-level data}
\usage{
GIFT_summary(
  Zscore1,
  Zscore2,
  Sigma1,
  Sigma2,
  R,
  n1,
  n2,
  pindex,
  max_iterin = 1000,
  epsin = 1e-04,
  Cores = 1
)
}
\arguments{
\item{Sigma1}{the LD matrix in eQTL data}

\item{Sigma2}{the LD matrix in GWAS data,both Sigma2sin and sigma1sin are often the same from the reference panel}

\item{R}{the estimated correlation matrix of gene expressions}

\item{n1}{the sample size of eQTL data}

\item{n2}{the sample size of GWAS data}

\item{pindex}{a vector with each element represents the number of cis-SNPs for each gene.}

\item{max_iterin}{The maximum iteration, which can be determined by users.}

\item{epsin}{The convergence tolerance of the absolute value of the difference  between the nth and (n+1)th log likelihood, which can be determined by users.}

\item{Cores}{The number of cores used in analysis. If the number of cores is greater than 1, analysis will perform with fast parallel computing. The function mclapply() depends on another R package "parallel" in Linux.}

\item{Zscore_1}{the Zscore matrix of the cis-SNP effect size matrix, each column for one specific gene in eQTL data}

\item{Zscore_2}{the Zscore vector of the cis-SNP effect size vector for one specific trait in GWAS data}
}
\value{
A list of estimated parameters including the p values for the gene-based test. 
\item{causal_effect}{The estimates of causal effect for each gene in a specific region}
\item{gene_based_test_pvalue}{The p values for each gene by the gene-based test}
}
\description{
GIFT_individual applies a likelihood-based approach, accounting for the correlated cis-SNPs and genes in a region
}
