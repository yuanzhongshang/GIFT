# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

GIFT_individualcpp <- function(X, y, Zx, Zy, Omega, constralpha, constrp, maxIter, epsStopLogLik) {
    .Call(`_GIFT_GIFT_individualcpp`, X, y, Zx, Zy, Omega, constralpha, constrp, maxIter, epsStopLogLik)
}

GIFT_individualcpppleio <- function(X, y, Zx, Zy, Omega, constralpha, constrp, maxIter, epsStopLogLik, pleioindex) {
    .Call(`_GIFT_GIFT_individualcpppleio`, X, y, Zx, Zy, Omega, constralpha, constrp, maxIter, epsStopLogLik, pleioindex)
}

GIFT_summarycpp <- function(betax, betay, Sigma1s, Sigma2s, Omega, constralpha, constrp, k, n1, n2, maxIter, epsStopLogLik) {
    .Call(`_GIFT_GIFT_summarycpp`, betax, betay, Sigma1s, Sigma2s, Omega, constralpha, constrp, k, n1, n2, maxIter, epsStopLogLik)
}

GIFT_summarycpppleio <- function(betax, betay, Sigma1s, Sigma2s, Omega, constralpha, constrp, k, n1, n2, maxIter, epsStopLogLik, pleioindex) {
    .Call(`_GIFT_GIFT_summarycpppleio`, betax, betay, Sigma1s, Sigma2s, Omega, constralpha, constrp, k, n1, n2, maxIter, epsStopLogLik, pleioindex)
}

