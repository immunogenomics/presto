#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <iostream>
#include <vector>
#include <algorithm>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(Rcpp]]


using namespace std;


// [[Rcpp::export]]
arma::cube collapse_mats(vector<arma::mat> mat_list, unsigned n) {
//     arma::cube result = arma::zeros<arma::cube>(10, 10, n);
    arma::cube result = arma::zeros<arma::cube>(mat_list[0].n_rows, mat_list[0].n_cols, n);
    for (int i = 0; i < n; i++) {
        result.slice(i) = mat_list[i];
    }
    return result;
}


// [[Rcpp::export]]
arma::mat collapse_vecs(vector<arma::vec> vec_list, unsigned n) {
//     arma::cube result = arma::zeros<arma::cube>(10, 10, n);
    arma::mat result = arma::zeros<arma::mat>(vec_list[0].n_elem, n);
    for (int i = 0; i < n; i++) {
        result.col(i) = vec_list[i];
    }
    return result;
}

