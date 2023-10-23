#define ARMA_NO_DEBUG

#include <armadillo>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo, BH, bigmemory)]]

#include <Rcpp.h>
#include <vector>

using namespace Rcpp;
using namespace arma;

#include <bigmemory/BigMatrix.h>
// [[Rcpp::plugins(cpp11)]]


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export()]]
Rcpp::List calcT0Wsim1 (arma::mat& CvSqrt,arma::mat& weight, arma::vec powV, int nperm) {
    
    //const int n = X1.n_rows;
    const int npow = powV.n_rows;
    const int k = CvSqrt.n_rows;
    // containers
    arma::mat T0s1(nperm,npow);
    T0s1.fill(0);
    
    for (int i = 0; i < nperm; i++) {
        arma::mat U00 = arma::randn(k,1);
        arma::mat U0 = CvSqrt * U00;
        arma::mat U0tmp = weight % U0;
        
        for (int j = 0; j < npow; j++) {
            if (powV[j] == 0) {
                arma::mat tmpU01 = abs(U0tmp);
                T0s1(i,j) = tmpU01.max();
            } else {
                T0s1(i,j) = accu(pow(U0tmp,powV[j]));
            }
        }
    }
    Rcpp::List res;
    res["T0"] =T0s1;
    
    return(res);
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export()]]
arma::vec avg_rank(arma::vec x) {
    
    arma::uvec w  = arma::stable_sort_index(x,"descend");
    R_xlen_t sz = x.size();
    arma::vec r(sz);
    
    for (R_xlen_t n, i = 0; i < sz; i += n) {
        n = 1;
        while (i + n < sz && x[w[i]] == x[w[i + n]]) ++n;
        for (R_xlen_t k = 0; k < n; k++) {
            r[w[i + k]] = i + (n + 1) / 2.;
        }
    }
    
    return r;
    
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export()]]
Rcpp::List calcT0Wsim3 (arma::mat& CvSqrt,arma::mat& weight, arma::vec powV, arma::vec Tsabs, SEXP pBigMat, int nperm) {
    
    //const int n = X1.n_rows;
    const int npow = powV.n_rows;
    const int k = CvSqrt.n_rows;
    
    // containers
    XPtr<BigMatrix> xpMat(pBigMat);
    arma::mat T0s = arma::Mat<double> ( (double *)xpMat->matrix(), xpMat->nrow(), xpMat->ncol(), false);
    
    for (int i = 0; i < nperm; i++) {
        arma::mat U00 = arma::randn(k,1);
        arma::mat U0 = CvSqrt * U00;
        arma::mat U0tmp = weight % U0;
        
        for (int j = 0; j < npow; j++) {
            if (powV[j] == 0) {
                arma::mat tmpU01 = abs(U0tmp);
                T0s(i,j) = tmpU01.max();
            } else {
                T0s(i,j) = accu(pow(U0tmp,powV[j]));
            }
        }
    }
    
    T0s = arma::abs(T0s);
    
    arma::vec minp0(nperm);
    arma::vec P0s(nperm);
    arma::vec pPerm0(npow);
    
    for (int i = 0; i < npow; i++) {
        
        // Calculate p-value
        arma::vec T0stmp = T0s.col(i);
        int tmp3 = 0;
        arma::vec rankarma  = avg_rank(T0stmp);
        
        for( int tt=0 ; tt < nperm ; tt++) {
            if( Tsabs(i)  <= T0stmp(tt) ) {
                tmp3++;
            }
            
            P0s(tt) = (double) rankarma(tt) / (double) nperm;
        }
        
        
        if(i == 0) {
            minp0 = P0s;
        } else {
            for( int ii = 0; ii < nperm; ii++) {
                if( minp0(ii) > P0s(ii) ) {
                    minp0(ii) = P0s(ii);
                }
            }
        }
        pPerm0(i) = (double) tmp3/ (double) nperm;
        
    }
    
    return Rcpp::List::create(Rcpp::Named("minp0") = minp0,
                              Rcpp::Named("pPerm0") = pPerm0 );
    
}






// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export()]]
Rcpp::List calcT0Wsim4 (arma::mat& CvSqrt,arma::mat& weight, arma::vec Tsabs, arma::vec powV, SEXP pBigMat, SEXP pBigMat3, int nperm) {
    
    //const int n = X1.n_rows;
    const int npow = powV.n_rows;
    const int nweight = weight.n_cols;
    const int k = CvSqrt.n_rows;
    // containers
    XPtr<BigMatrix> xpMat(pBigMat);
    arma::mat T0s = arma::Mat<double> ( (double *)xpMat->matrix(), xpMat->nrow(), xpMat->ncol(), false);
    
    XPtr<BigMatrix> xpMat3(pBigMat3);
    arma::mat minp0_sign = arma::Mat<double> ( (double *)xpMat3->matrix(), xpMat3->nrow(), xpMat3->ncol(), false);
    
    for (int i = 0; i < nperm; i++) {
        arma::mat U00 = arma::randn(k,1);
        arma::mat U0 = CvSqrt * U00;
        
        for (int k = 0; k < nweight; k++) {
            arma::mat weight_tmp = weight.col(k);
            arma::mat U0tmp = weight_tmp % U0;
            
            for (int j = 0; j < npow; j++) {
                if (powV[j] == 0) {
                    arma::mat tmpU01 = abs(U0tmp);
                    T0s(i,j * nweight + k) = tmpU01.max();
                } else {
                    T0s(i,j * nweight + k) = accu(pow(U0tmp,powV[j]));
                }
            }
        }
    }
    
    for (int k =0; k<nweight; k++) {
        minp0_sign.col(k) = T0s.col(k);
    }
    
    int npw = npow * nweight;
    arma::mat cov_res(npw,npw);
    cov_res.fill(0);
    cov_res = cov(T0s);
    
    T0s = arma::abs(T0s);
    
    arma::vec minp0(nperm);
    arma::vec P0s(nperm);
    arma::mat pPerm0(npow + 1,nweight);
    pPerm0.fill(-1);
    
    for(int k = 0; k < nweight; k++) {
        
        double minP_tmp = 1.0;
        
        for (int j = 0; j < npow; j++) {
            
            // Calculate p-value
            arma::vec T0stmp = T0s.col(j * nweight + k);
            int tmp3 = 0;
            arma::vec rankarma  = avg_rank(T0stmp);
            
            for( int tt=0 ; tt < nperm ; tt++) {
                if( Tsabs(j * nweight + k)  <= T0stmp(tt) ) {
                    tmp3++;
                }
                
                P0s(tt) = (double) rankarma(tt) / (double) nperm;
            }
            
            if(j == 0) {
                minp0 = P0s;
            } else {
                for( int ii = 0; ii < nperm; ii++) {
                    if( minp0(ii) > P0s(ii) ) {
                        minp0(ii) = P0s(ii);
                    }
                }
            }
            double tmp_pvalue = (double) tmp3 / (double) nperm;
            pPerm0(j,k) = tmp_pvalue;
            
            if (tmp_pvalue <= minP_tmp) {
                minP_tmp = tmp_pvalue;
            }
        }
        
        int count_pvalue = 0;
        for(int j=0; j < nperm; j++) {
            if( minP_tmp > minp0(j) ) {
                count_pvalue = count_pvalue + 1;
            }
        }
        
        pPerm0(npow,k) = (double) (count_pvalue + 1) / (double) (nperm + 1);
        
    }
    
    return Rcpp::List::create(Rcpp::Named("pPerm0") = pPerm0,
                              Rcpp::Named("cov") = cov_res );
    
}





// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export()]]
Rcpp::List calcT0WsimV3 (arma::mat& CvSqrt,arma::mat& weight, arma::vec powV, SEXP pBigMat, SEXP pBigMat3,int nperm) {
    
    //const int n = X1.n_rows;
    const int npow = powV.n_rows;
    const int nweight = weight.n_cols;
    const int k = CvSqrt.n_rows;
    // containers
    XPtr<BigMatrix> xpMat(pBigMat);
    arma::mat T0s = arma::Mat<double> ( (double *)xpMat->matrix(), xpMat->nrow(), xpMat->ncol(), false);
    
    XPtr<BigMatrix> xpMat3(pBigMat3);
    arma::mat minp0_sign = arma::Mat<double> ( (double *)xpMat3->matrix(), xpMat3->nrow(), xpMat3->ncol(), false);
    
    for (int i = 0; i < nperm; i++) {
        arma::mat U00 = arma::randn(k,1);
        arma::mat U0 = CvSqrt * U00;
        
        for (int k = 0; k < nweight; k++) {
            arma::mat weight_tmp = weight.col(k);
            arma::mat U0tmp = weight_tmp % U0;
            
            for (int j = 0; j < npow; j++) {
                if (powV[j] == 0) {
                    arma::mat tmpU01 = abs(U0tmp);
                    T0s(i,j * nweight + k) = tmpU01.max();
                } else {
                    T0s(i,j * nweight + k) = accu(pow(U0tmp,powV[j]));
                }
            }
        }
    }
    
    for (int k =0; k<nweight; k++) {
        minp0_sign.col(k) = T0s.col(k);
    }
    
    int npw = npow * nweight;
    arma::mat cov_res(npw,npw);
    cov_res.fill(0);
    cov_res = arma::cov(T0s);
    
    arma::mat mean_res(npw,1);
    mean_res.fill(0);
    mean_res = arma::mean(T0s,0);
    
    return Rcpp::List::create(Rcpp::Named("cov") = cov_res,
                              Rcpp::Named("mean") = mean_res );
    
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export()]]
void calc_test_ch(arma::mat& cov, arma::mat& weight, arma::vec powV, SEXP pBigMat, SEXP pBigMat2,int nperm) {
    
    //const int n = X1.n_rows;
    const int npow = powV.n_rows;
    const int nweight = weight.n_cols;
    
    // containers
    XPtr<BigMatrix> xpMat(pBigMat);
    arma::mat T0s = arma::Mat<double> ( (double *)xpMat->matrix(), xpMat->nrow(), xpMat->ncol(), false);
    
    XPtr<BigMatrix> xpMat2(pBigMat2);
    arma::mat T1s = arma::Mat<double> ( (double *)xpMat2->matrix(), xpMat2->nrow(), xpMat2->ncol(), false);
    
    arma::vec pPerm0(npow);
    
    for (int i = 0; i < npow; i++) {
        
        arma::mat subcov(nweight,nweight);
        subcov.fill(0);
        
        subcov = cov.submat(i*nweight,i*nweight,((i+1)*nweight -1),((i+1)*nweight -1));
        subcov = subcov.i();
        
        for (int j =0; j < nperm; j++) {
            arma::vec tmp_ts = T0s.row(j);
            arma::mat tmp_ts2(1,nweight);
            tmp_ts2 = tmp_ts.subvec(i*nweight,((i+1) * nweight -1));
            arma::mat tmp = tmp_ts2.t() * subcov * tmp_ts2;
            T1s(j,i) = tmp(0,0);
        }
    }
}




// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export()]]
void calc_test_ch_V3(arma::mat& cov,arma::vec & mean, arma::mat& weight, arma::vec powV, SEXP pBigMat, SEXP pBigMat2,int nperm) {
    
    //const int n = X1.n_rows;
    const int npow = powV.n_rows;
    const int nweight = weight.n_cols;
    
    // containers
    XPtr<BigMatrix> xpMat(pBigMat);
    arma::mat T0s = arma::Mat<double> ( (double *)xpMat->matrix(), xpMat->nrow(), xpMat->ncol(), false);
    
    XPtr<BigMatrix> xpMat2(pBigMat2);
    arma::mat T1s = arma::Mat<double> ( (double *)xpMat2->matrix(), xpMat2->nrow(), xpMat2->ncol(), false);
    
    arma::vec pPerm0(npow);
    
    for (int i = 0; i < npow; i++) {
        
        arma::mat subcov(nweight,nweight);
        subcov.fill(0);
        
        subcov = cov.submat(i*nweight,i*nweight,((i+1)*nweight -1),((i+1)*nweight -1));
        subcov = subcov.i();
        
        arma::mat tmp_mean(1,nweight);
        tmp_mean = mean.subvec(i*nweight,((i+1) * nweight -1));
        for (int j =0; j < nperm; j++) {
            arma::vec tmp_ts = T0s.row(j);
            arma::mat tmp_ts2(1,nweight);
            
            tmp_ts2 = tmp_ts.subvec(i*nweight,((i+1) * nweight -1));
            tmp_ts2 = tmp_ts2 - tmp_mean;
            arma::mat tmp = tmp_ts2.t() * subcov * tmp_ts2;
            T1s(j,i) = tmp(0,0);
        }
    }
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export()]]
void calc_test_pan(arma::mat& cov, arma::mat& weight, arma::vec powV, SEXP pBigMat,SEXP pBigMat2, int nperm) {
    
    //const int n = X1.n_rows;
    const int npow = powV.n_rows;
    const int nweight = weight.n_cols;
    
    // containers
    XPtr<BigMatrix> xpMat(pBigMat);
    arma::mat T1s = arma::Mat<double> ( (double *)xpMat->matrix(), xpMat->nrow(), xpMat->ncol(), false);
    
    XPtr<BigMatrix> xpMat2(pBigMat2);
    arma::mat T2s = arma::Mat<double> ( (double *)xpMat2->matrix(), xpMat2->nrow(), xpMat2->ncol(), false);
    
    arma::vec pPerm0(npow);
    
    for (int i = 0; i < npow; i++) {
        
        arma::mat subcov(nweight,nweight);
        subcov.fill(0);
        
        subcov = cov.submat(i*nweight,i*nweight,((i+1)*nweight -1),((i+1)*nweight -1));
        
        for (int j =0; j < nperm; j++) {
            arma::vec tmp_ts = T2s.row(j);
            arma::mat tmp_ts2(1,nweight);
            tmp_ts2 = tmp_ts.subvec(i*nweight,((i+1) * nweight -1));
            arma::mat tmp = tmp_ts2.t() * subcov * tmp_ts2;
            T1s(j,i) = tmp(0,0);
        }
    }
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export()]]
arma::mat calc_test_ts(arma::mat& cov, arma::mat& weight, arma::vec powV, arma::vec Ts) {
    
    //const int n = X1.n_rows;
    const int npow = powV.n_rows;
    const int nweight = weight.n_cols;
    
    arma::mat Ts2(npow,1);
    Ts2.fill(0);
    
    for (int i = 0; i < npow; i++) {
        arma::mat subcov(nweight,nweight);
        subcov.fill(0);
        
        subcov = cov.submat(i*nweight,i*nweight,((i+1)*nweight -1),((i+1)*nweight -1));
        
        arma::mat tmp_ts2(1,nweight);
        tmp_ts2 = Ts.subvec(i*nweight,((i+1)*nweight -1));
        arma::mat tmp =  tmp_ts2.t() * subcov * tmp_ts2;
        Ts2(i,0) = tmp(0,0);
    }
    
    return Ts2;
}



// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export()]]
arma::mat calc_test_ts_V3(arma::mat& cov, arma::vec & mean, arma::mat& weight, arma::vec powV, arma::vec Ts) {
    
    //const int n = X1.n_rows;
    const int npow = powV.n_rows;
    const int nweight = weight.n_cols;
    
    arma::mat Ts2(npow,1);
    Ts2.fill(0);
    
    for (int i = 0; i < npow; i++) {
        arma::mat subcov(nweight,nweight);
        subcov.fill(0);
        
        subcov = cov.submat(i*nweight,i*nweight,((i+1)*nweight -1),((i+1)*nweight -1));
        subcov = subcov.i();
        
        arma::mat tmp_mean(1,nweight);
        tmp_mean = mean.subvec(i*nweight,((i+1) * nweight -1));
        
        arma::mat tmp_ts2(1,nweight);
        tmp_ts2 = Ts.subvec(i*nweight,((i+1)*nweight -1));
        tmp_ts2 = tmp_ts2 - tmp_mean;
        
        arma::mat tmp =  tmp_ts2.t() * subcov * tmp_ts2;
        Ts2(i,0) = tmp(0,0);
    }
    
    return Ts2;
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export()]]
Rcpp::List calc_p_ch (arma::vec powV, arma::vec Tsabs, SEXP pBigMat2, int nperm) {
    
    //const int n = X1.n_rows;
    const int npow = powV.n_rows;
    
    // containers
    XPtr<BigMatrix> xpMat2(pBigMat2);
    arma::mat T1s = arma::Mat<double> ( (double *)xpMat2->matrix(), xpMat2->nrow(), xpMat2->ncol(), false);
    
    arma::vec minp0(nperm);
    arma::vec P0s(nperm);
    arma::vec pPerm0(npow);
    
    for (int i = 0; i < npow; i++) {
        
        // Calculate p-value
        arma::vec T1stmp = T1s.col(i);
        int tmp3 = 0;
        
        arma::vec rankarma  = avg_rank(T1stmp);
        
        for( int tt=0 ; tt < nperm ; tt++) {
            if( Tsabs(i)  <= T1stmp(tt) ) {
                tmp3++;
            }
            P0s(tt) = (double) rankarma(tt) / (double) nperm;
        }
        
        
        if(i == 0) {
            minp0 = P0s;
        } else {
            for( int ii = 0; ii < nperm; ii++) {
                if( minp0(ii) > P0s(ii) ) {
                    minp0(ii) = P0s(ii);
                }
            }
        }
        pPerm0(i) = (double) tmp3 / (double) nperm;
        
    }
    
    return Rcpp::List::create(Rcpp::Named("minp0") = minp0,
                              Rcpp::Named("pPerm0") = pPerm0 );
    
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export()]]
Rcpp::List calc_p_ch_pan (arma::vec powV, arma::vec Tsabs, SEXP pBigMat, int nperm) {
    
    //const int n = X1.n_rows;
    const int npow = powV.n_rows;
    
    // containers
    XPtr<BigMatrix> xpMat(pBigMat);
    arma::mat T1s = arma::Mat<double> ( (double *)xpMat->matrix(), xpMat->nrow(), xpMat->ncol(), false);
    
    
    arma::vec minp0(nperm);
    arma::vec P0s(nperm);
    arma::vec pPerm0(npow);
    
    for (int i = 0; i < npow; i++) {
        
        // Calculate p-value
        arma::vec T1stmp = T1s.col(i);
        int tmp3 = 0;
        
        arma::vec rankarma  = avg_rank(T1stmp);
        
        for( int tt=0 ; tt < nperm ; tt++) {
            if( Tsabs(i)  <= T1stmp(tt) ) {
                tmp3++;
            }
            P0s(tt) = (double) rankarma(tt) / (double) nperm;
        }
        
        
        if(i == 0) {
            minp0 = P0s;
        } else {
            for( int ii = 0; ii < nperm; ii++) {
                if( minp0(ii) > P0s(ii) ) {
                    minp0(ii) = P0s(ii);
                }
            }
        }
        pPerm0(i) = (double) tmp3 / (double) nperm;
        
    }
    
    return Rcpp::List::create(Rcpp::Named("minp0") = minp0,
                              Rcpp::Named("pPerm0") = pPerm0 );
    
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export()]]
Rcpp::List calc_p_pan (arma::vec Tsabs,SEXP pBigMat,SEXP pBigMat2,SEXP pBigMat3, int nTs,int nweight, int nperm) {
    
    // containers
    XPtr<BigMatrix> xpMat(pBigMat);
    arma::mat T0s = arma::Mat<double> ( (double *)xpMat->matrix(), xpMat->nrow(), xpMat->ncol(), false);
    
    XPtr<BigMatrix> xpMat2(pBigMat2);
    arma::mat T2s = arma::Mat<double> ( (double *)xpMat2->matrix(), xpMat2->nrow(), xpMat2->ncol(), false);
    
    
    XPtr<BigMatrix> xpMat3(pBigMat3);
    arma::mat minp0_sign = arma::Mat<double> ( (double *)xpMat3->matrix(), xpMat3->nrow(), xpMat3->ncol(), false);
    
    
    // T0s = arma::abs(T0s);
    
    arma::vec minp0(nperm);
    arma::vec P0s(nperm);
    arma::vec pPerm0(nTs);
    
    
    for (int i = 0; i < nTs; i++) {
        
        // Calculate p-value
        arma::vec T1stmp = T0s.col(i);
        int tmp3 = 0;
        
        arma::vec rankarma  = avg_rank(T1stmp);
        
        int minp0_index = i - floor(i/nweight) * nweight;
        
        for( int tt=0 ; tt < nperm ; tt++) {
            if( Tsabs(i)  <= T1stmp(tt) ) {
                tmp3++;
            }
            
            NumericVector tmp_res;
            tmp_res = (double) rankarma(tt) / (double) nperm;
            
            int p0_sign = -1;
            if (minp0_sign(tt,minp0_index) >0) {
                p0_sign = 1;
            }
            
            tmp_res = qnorm(1 - tmp_res/2) * p0_sign;
            T2s(tt,i) = tmp_res[0];
        }
        
        pPerm0(i) = (double) tmp3 / (double) nperm;
    }
    
    arma::mat cov_res(nTs,nTs);
    cov_res.fill(0);
    cov_res = arma::cov(T2s);
    
    return Rcpp::List::create(Rcpp::Named("pPerm0") = pPerm0,
                              Rcpp::Named("cov") = cov_res );
    
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export()]]
arma::vec calc_aSPU (arma::mat p_cov, SEXP pBigMat2, int nperm) {
    
    // containers
    XPtr<BigMatrix> xpMat2(pBigMat2);
    arma::mat minp0 = arma::Mat<double> ( (double *)xpMat2->matrix(), xpMat2->nrow(), xpMat2->ncol(), false);
    
    arma::vec res(nperm);
    
    for (int i = 0; i < nperm; i++) {
        arma::mat tmp =  minp0.row(i) * p_cov * minp0.row(i).t();
        res(i) = tmp(0,0);
    }
    
    return res;
}



// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export()]]
Rcpp::List daSPU_calcT0Wsim4 (arma::mat& CvSqrt,arma::mat& weight, arma::vec powV, arma::vec Tsabs, SEXP pBigMat, SEXP pBigMat2,SEXP pBigMat3,int nperm) {
    
    //const int n = X1.n_rows;
    const int npow = powV.n_rows;
    const int nweight = weight.n_cols;
    const int k = CvSqrt.n_rows;
    // containers
    XPtr<BigMatrix> xpMat(pBigMat);
    arma::mat T0s = arma::Mat<double> ( (double *)xpMat->matrix(), xpMat->nrow(), xpMat->ncol(), false);
    
    // containers
    XPtr<BigMatrix> xpMat2(pBigMat2);
    arma::mat minp0 = arma::Mat<double> ( (double *)xpMat2->matrix(), xpMat2->nrow(), xpMat2->ncol(), false);
    
    // containers
    XPtr<BigMatrix> xpMat3(pBigMat3);
    arma::mat minp0_sign = arma::Mat<double> ( (double *)xpMat3->matrix(), xpMat3->nrow(), xpMat3->ncol(), false);
    
    for (int i = 0; i < nperm; i++) {
        arma::mat U00 = arma::randn(k,1);
        arma::mat U0 = CvSqrt * U00;
        
        for (int k = 0; k < nweight; k++) {
            arma::mat weight_tmp = weight.col(k);
            arma::mat U0tmp = weight_tmp % U0;
            
            for (int j = 0; j < npow; j++) {
                if (powV[j] == 0) {
                    arma::mat tmpU01 = abs(U0tmp);
                    T0s(i,k * npow + j) = tmpU01.max();
                } else {
                    T0s(i,k * npow + j) = accu(pow(U0tmp,powV[j]));
                }
            }
        }
    }
    
    for (int k =0; k<nweight; k++) {
        minp0_sign.col(k) = T0s.col(k * npow);
    }
    
    T0s = arma::abs(T0s);
    
    arma::vec minp0_tmp(nperm);
    
    arma::vec P0s(nperm);
    arma::vec pPerm0(npow* nweight);
    
    for (int k = 0; k < nweight; k++) {
        
        for (int i = (npow*k); i < (npow* (k + 1)); i++) {
            
            // Calculate p-value
            arma::vec T0stmp = T0s.col(i);
            int tmp3 = 0;
            arma::vec rankarma  = avg_rank(T0stmp);
            
            for( int tt=0 ; tt < nperm ; tt++) {
                if( Tsabs(i)  <= T0stmp(tt) ) {
                    tmp3++;
                }
                
                P0s(tt) = (double) rankarma(tt) / (double) nperm;
            }
            
            
            if(i == 0) {
                minp0_tmp = P0s;
            } else {
                for( int ii = 0; ii < nperm; ii++) {
                    if( minp0_tmp(ii) > P0s(ii) ) {
                        minp0_tmp(ii) = P0s(ii);
                    }
                }
            }
            pPerm0(i) = (double) tmp3/ (double) nperm;
        }
        
        minp0.col(k) = minp0_tmp;
    }
    
    for (int k = 0; k < nweight; k++) {
        for (int i = 0; i <nperm; i++) {
            NumericVector tmp_res;
            
            tmp_res = minp0(i,k);
            int p0_sign = -1;
            if (minp0_sign(i,k) >0) {
                p0_sign = 1;
            }
            tmp_res = qnorm(1 - tmp_res/2) * p0_sign;
            minp0(i,k) = tmp_res[0];
        }
    }
    
    arma::mat cov_res(nweight,nweight);
    cov_res.fill(0);
    cov_res = arma::cov(minp0);
    
    return Rcpp::List::create(Rcpp::Named("pPerm0") = pPerm0,
                              Rcpp::Named("cov") = cov_res );
    
    
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export()]]
arma::vec daSPU_calc_aSPU (arma::mat p_cov, SEXP pBigMat2, int nperm) {
    
    // containers
    XPtr<BigMatrix> xpMat2(pBigMat2);
    arma::mat minp0 = arma::Mat<double> ( (double *)xpMat2->matrix(), xpMat2->nrow(), xpMat2->ncol(), false);
    
    arma::vec res(nperm);
    
    for (int i = 0; i < nperm; i++) {
        arma::mat tmp =  minp0.row(i) * p_cov * minp0.row(i).t();
        res(i) = tmp(0,0);
    }
    
    return res;
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export()]]
Rcpp::List calcT0Wsim2 (arma::mat& CvSqrt,arma::mat& weight, arma::vec pow1, arma::vec pow2, int nperm) {
    
    //const int n = X1.n_rows;
    const int n_pow1 = pow1.size();
    const int n_pow2 = pow2.size();
    const int n_weight = weight.n_cols;
    
    const int k = CvSqrt.n_rows;
    // containers
    arma::mat T0s1(n_pow1,n_weight);
    T0s1.fill(0);
    arma::mat T0s(nperm,n_pow1 * n_pow2);
    T0s.fill(0);
    
    for (int i = 0; i < nperm; i++) {
        arma::mat U00 = arma::randn(k,1);
        arma::mat U0 = CvSqrt * U00;
        
        for (int i2 = 0; i2 < n_weight; i2 ++) {
            arma::mat U0tmp = weight.col(i2) % U0;
            for (int j1 = 0; j1 < n_pow1; j1++) {
                if( pow1[j1] > 0) {
                    arma::vec tmp1 = pow(U0tmp,pow1[j1]);
                    double tmp2 = sum(tmp1);
                    
                    if( tmp2 > 0 ) {
                        T0s1(j1,i2) = pow(std::abs(tmp2), 1/pow1[j1]);
                    } else {
                        T0s1(j1,i2) = -pow(std::abs(tmp2), 1/pow1[j1]);
                    }
                    
                } else {
                    arma::vec tmp1 = abs(U0tmp);
                    T0s1(j1,i2) = max(tmp1);
                }
            }
        }
        
        for (int j1 = 0; j1 < n_pow1; j1++) {
            for (int j2 = 0; j2 < n_pow2; j2++) {
                if (pow2[j2] > 0) {
                    arma::mat tmp3 = pow(T0s1.row(j1),pow2[j2]);
                    T0s(i, j2 * n_pow1 + j1) = accu(tmp3);
                } else {
                    T0s(i, j2 * n_pow1 + j1) = max(abs(T0s1.row(j1)));
                }
            }
        }
    }
    
    
    Rcpp::List res;
    res["T0s"] =T0s;
    
    return(res);
}





// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export()]]
arma::vec daSPU_calcT0 (arma::mat& CvSqrt,arma::mat& weight, arma::vec pow1, arma::vec pow2,  arma::vec Tsabs,SEXP pBigMat, int nperm) {
    
    //const int n = X1.n_rows;
    const int n_pow1 = pow1.size();
    const int n_pow2 = pow2.size();
    const int n_pow = n_pow1 * n_pow2;
    const int n_weight = weight.n_cols;
    
    const int k = CvSqrt.n_rows;
    // containers
    arma::mat T0s1(n_pow1,n_weight);
    T0s1.fill(0);
    
    XPtr<BigMatrix> xpMat(pBigMat);
    arma::mat T0s = arma::Mat<double> ( (double *)xpMat->matrix(), xpMat->nrow(), xpMat->ncol(), false);
    
    for (int i = 0; i < nperm; i++) {
        arma::mat U00 = arma::randn(k,1);
        arma::mat U0 = CvSqrt * U00;
        
        for (int i2 = 0; i2 < n_weight; i2 ++) {
            arma::mat U0tmp = weight.col(i2) % U0;
            for (int j1 = 0; j1 < n_pow1; j1++) {
                if( pow1[j1] > 0) {
                    arma::vec tmp1 = pow(U0tmp,pow1[j1]);
                    double tmp2 = sum(tmp1);
                    
                    if( tmp2 > 0 ) {
                        T0s1(j1,i2) = pow(std::abs(tmp2), 1/pow1[j1]);
                    } else {
                        T0s1(j1,i2) = -pow(std::abs(tmp2), 1/pow1[j1]);
                    }
                    
                } else {
                    arma::vec tmp1 = abs(U0tmp);
                    T0s1(j1,i2) = max(tmp1);
                }
            }
        }
        
        for (int j1 = 0; j1 < n_pow1; j1++) {
            for (int j2 = 0; j2 < n_pow2; j2++) {
                if (pow2[j2] > 0) {
                    arma::mat tmp3 = pow(T0s1.row(j1),pow2[j2]);
                    T0s(i, j2 * n_pow1 + j1) = accu(tmp3);
                } else {
                    T0s(i, j2 * n_pow1 + j1) = max(abs(T0s1.row(j1)));
                }
            }
        }
    }
    
    //
    T0s = arma::abs(T0s);
    
    arma::vec minp0(nperm);
    arma::vec P0s(nperm);
    arma::vec pPerm0(n_pow + 1);
    double minP_tmp = 1.0;

    for(int k = 0; k < n_pow; k++) {
        
        
        // Calculate p-value
        arma::vec T0stmp = T0s.col(k);
        int tmp3 = 0;
        arma::vec rankarma  = avg_rank(T0stmp);
        
        for( int tt=0 ; tt < nperm ; tt++) {
            if( Tsabs(k)  <= T0stmp(tt) ) {
                tmp3++;
            }
            
            P0s(tt) = (double) rankarma(tt) / (double) nperm;
        }
        
        if(k == 0) {
            minp0 = P0s;
        } else {
            for( int ii = 0; ii < nperm; ii++) {
                if( minp0(ii) > P0s(ii) ) {
                    minp0(ii) = P0s(ii);
                }
            }
        }
        double tmp_pvalue = (double) tmp3 / (double) nperm;
        pPerm0(k) = tmp_pvalue;
        
        if (tmp_pvalue <= minP_tmp) {
            minP_tmp = tmp_pvalue;
        }
        
    }
    
    int count_pvalue = 0;
    for(int j=0; j < nperm; j++) {
        if( minP_tmp >= minp0(j) ) {
            count_pvalue = count_pvalue + 1;
        }
    }
    
    pPerm0(n_pow) = (double) (count_pvalue + 1) / (double) (nperm + 1);

    
    return pPerm0;
}






