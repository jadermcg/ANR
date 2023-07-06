#pragma once
#include <RcppArmadillo.h>
#include <strings.h>
#include <vector>
#include <omp.h>

//[[Rcpp::depends(RcppArmadillo)]]

Rcpp::List anr(const std::vector<std::string> &fasta, arma::mat alpha, const arma::mat &beta, const double cutoff, int niter, double w);
Rcpp::List loganr(const std::vector<std::string> &fasta, arma::mat alpha, const arma::mat &beta, const double cutoff, int niter, double w);
