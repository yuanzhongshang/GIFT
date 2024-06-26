// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// GIFT_individualcpp
List GIFT_individualcpp(arma::vec X, arma::vec y, arma::mat Zx, arma::mat Zy, arma::mat Omega, arma::vec constralpha, arma::vec constrp, int maxIter, double epsStopLogLik);
RcppExport SEXP _GIFT_GIFT_individualcpp(SEXP XSEXP, SEXP ySEXP, SEXP ZxSEXP, SEXP ZySEXP, SEXP OmegaSEXP, SEXP constralphaSEXP, SEXP constrpSEXP, SEXP maxIterSEXP, SEXP epsStopLogLikSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Zx(ZxSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Zy(ZySEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Omega(OmegaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type constralpha(constralphaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type constrp(constrpSEXP);
    Rcpp::traits::input_parameter< int >::type maxIter(maxIterSEXP);
    Rcpp::traits::input_parameter< double >::type epsStopLogLik(epsStopLogLikSEXP);
    rcpp_result_gen = Rcpp::wrap(GIFT_individualcpp(X, y, Zx, Zy, Omega, constralpha, constrp, maxIter, epsStopLogLik));
    return rcpp_result_gen;
END_RCPP
}
// GIFT_individualcpppleio
List GIFT_individualcpppleio(arma::vec X, arma::vec y, arma::mat Zx, arma::mat Zy, arma::mat Omega, arma::vec constralpha, arma::vec constrp, int maxIter, double epsStopLogLik, arma::vec pleioindex);
RcppExport SEXP _GIFT_GIFT_individualcpppleio(SEXP XSEXP, SEXP ySEXP, SEXP ZxSEXP, SEXP ZySEXP, SEXP OmegaSEXP, SEXP constralphaSEXP, SEXP constrpSEXP, SEXP maxIterSEXP, SEXP epsStopLogLikSEXP, SEXP pleioindexSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Zx(ZxSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Zy(ZySEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Omega(OmegaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type constralpha(constralphaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type constrp(constrpSEXP);
    Rcpp::traits::input_parameter< int >::type maxIter(maxIterSEXP);
    Rcpp::traits::input_parameter< double >::type epsStopLogLik(epsStopLogLikSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type pleioindex(pleioindexSEXP);
    rcpp_result_gen = Rcpp::wrap(GIFT_individualcpppleio(X, y, Zx, Zy, Omega, constralpha, constrp, maxIter, epsStopLogLik, pleioindex));
    return rcpp_result_gen;
END_RCPP
}
// GIFT_summarycpp
List GIFT_summarycpp(arma::mat betax, arma::vec betay, arma::mat Sigma1s, arma::mat Sigma2s, arma::mat Omega, arma::vec constralpha, arma::vec constrp, int k, int n1, int n2, int maxIter, double epsStopLogLik);
RcppExport SEXP _GIFT_GIFT_summarycpp(SEXP betaxSEXP, SEXP betaySEXP, SEXP Sigma1sSEXP, SEXP Sigma2sSEXP, SEXP OmegaSEXP, SEXP constralphaSEXP, SEXP constrpSEXP, SEXP kSEXP, SEXP n1SEXP, SEXP n2SEXP, SEXP maxIterSEXP, SEXP epsStopLogLikSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type betax(betaxSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type betay(betaySEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Sigma1s(Sigma1sSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Sigma2s(Sigma2sSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Omega(OmegaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type constralpha(constralphaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type constrp(constrpSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< int >::type n1(n1SEXP);
    Rcpp::traits::input_parameter< int >::type n2(n2SEXP);
    Rcpp::traits::input_parameter< int >::type maxIter(maxIterSEXP);
    Rcpp::traits::input_parameter< double >::type epsStopLogLik(epsStopLogLikSEXP);
    rcpp_result_gen = Rcpp::wrap(GIFT_summarycpp(betax, betay, Sigma1s, Sigma2s, Omega, constralpha, constrp, k, n1, n2, maxIter, epsStopLogLik));
    return rcpp_result_gen;
END_RCPP
}
// GIFT_summarycpppleio
List GIFT_summarycpppleio(arma::mat betax, arma::vec betay, arma::mat Sigma1s, arma::mat Sigma2s, arma::mat Omega, arma::vec constralpha, arma::vec constrp, int k, int n1, int n2, int maxIter, double epsStopLogLik, arma::vec pleioindex);
RcppExport SEXP _GIFT_GIFT_summarycpppleio(SEXP betaxSEXP, SEXP betaySEXP, SEXP Sigma1sSEXP, SEXP Sigma2sSEXP, SEXP OmegaSEXP, SEXP constralphaSEXP, SEXP constrpSEXP, SEXP kSEXP, SEXP n1SEXP, SEXP n2SEXP, SEXP maxIterSEXP, SEXP epsStopLogLikSEXP, SEXP pleioindexSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type betax(betaxSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type betay(betaySEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Sigma1s(Sigma1sSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Sigma2s(Sigma2sSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Omega(OmegaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type constralpha(constralphaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type constrp(constrpSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< int >::type n1(n1SEXP);
    Rcpp::traits::input_parameter< int >::type n2(n2SEXP);
    Rcpp::traits::input_parameter< int >::type maxIter(maxIterSEXP);
    Rcpp::traits::input_parameter< double >::type epsStopLogLik(epsStopLogLikSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type pleioindex(pleioindexSEXP);
    rcpp_result_gen = Rcpp::wrap(GIFT_summarycpppleio(betax, betay, Sigma1s, Sigma2s, Omega, constralpha, constrp, k, n1, n2, maxIter, epsStopLogLik, pleioindex));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_GIFT_GIFT_individualcpp", (DL_FUNC) &_GIFT_GIFT_individualcpp, 9},
    {"_GIFT_GIFT_individualcpppleio", (DL_FUNC) &_GIFT_GIFT_individualcpppleio, 10},
    {"_GIFT_GIFT_summarycpp", (DL_FUNC) &_GIFT_GIFT_summarycpp, 12},
    {"_GIFT_GIFT_summarycpppleio", (DL_FUNC) &_GIFT_GIFT_summarycpppleio, 13},
    {NULL, NULL, 0}
};

RcppExport void R_init_GIFT(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
