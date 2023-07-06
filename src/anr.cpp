#include "utils.h"
#include "prob_utils.h"
#include "anr.h"

//' Runs Expectation Maximization ANR and reestimates model parameters.
//'@name anr
//'@param fasta Dataset of sequences.
//'@param alpha PWM model.
//'@param beta 0-order Markov Chain.
//'@param cutoff Cutoff for EM convergence.
//'@param niter Maximum number of iterations.
//'@param w Priori probability for motif belongs to a positive model.
//'@return Updated PWM model.
//[[Rcpp::export]]
Rcpp::List anr(const std::vector<std::string> &fasta, arma::mat alpha, const arma::mat &beta, const double cutoff, int niter, double w = 0.5) {
  /**
   * Parameters
   */
  int n = fasta.size();
  int t = fasta[0].size();
  int k = alpha.n_cols;
  int m = t - k + 1;
  
  /**
   * Convergence control
   */
  std::vector<double> convergence;
  std::vector<double> changes;
  convergence.push_back(-std::numeric_limits<double>::infinity());
  changes.push_back(0);
  
  /**
   * Posteriori
   */
  arma::vec z(m);
  
  /**
   * Model to reestimate
   */
  arma::mat new_alpha(4, k);
  double n1 = 0.0;
  
  while (true) {
    new_alpha.fill(1e-100);
    n1 = 0.0;
    for (const auto &seq : fasta) { // Foreach sequence
      z.fill(0.0);
      /**
       * E-STEP
       */
      
      for (int j = 0; j < m; ++j) {
        const auto &kmer = seq.substr(j, k);
        double a = w * probSeqGivenAlpha(kmer, alpha);
        double b = (1-w) * probSeqGivenBeta(kmer, beta);
        z[j] = a / (a + b);
        n1 += z[j];
      }
      
      /**
       * M-STEP
       */
      
      soft_update(new_alpha, seq, z.t());
    }
    
    w = n1 / (n*m);
    alpha = new_alpha / n1;
    
    /**
     * Convergence control
     */
    if (hasConverged(cutoff, niter, alpha, convergence, changes)) break;
    
    /**
     * Next iteration
     */
    --niter;
    
  }
  
  return Rcpp::List::create(alpha, convergence, changes, w);
}

//' Runs Expectation Maximization ANR and reestimates model parameters.
//'@name loganr
//'@param fasta Dataset of sequences.
//'@param alpha PWM model.
//'@param beta 0-order Markov Chain.
//'@param cutoff Cutoff for EM convergence.
//'@param niter Maximum number of iterations.
//'@param w Priori probability for motif belongs to a positive model.
//'@return Updated PWM model.
//[[Rcpp::export]]
Rcpp::List loganr(const std::vector<std::string> &fasta, arma::mat alpha, const arma::mat &beta, const double cutoff, int niter, double w = 0.5) {
  /**
   * Parameters
   */
  int n = fasta.size();
  int t = fasta[0].size();
  int k = alpha.n_cols;
  int m = t - k + 1;
  
  /**
   * Convergence control
   */
  std::vector<double> convergence;
  std::vector<double> changes;
  convergence.push_back(-std::numeric_limits<double>::infinity());
  changes.push_back(0);
  
  arma::mat new_alpha(4, k);
  arma::vec z(m);
  double n1 = 0.0;
  
  while (true) {
    new_alpha.fill(1e-100);
    n1 = 0.0;
    for (const auto &seq : fasta) { // Foreach sequence
      z.fill(0.0);
      /**
       * E-STEP
       */
      
      for (int j = 0; j < m; ++j) {
        const auto &kmer = seq.substr(j, k);
        double a = std::log(w) + probSeqGivenAlphaLog(kmer, alpha);
        double b = std::log(1-w) + probSeqGivenBetaLog(kmer, beta);
        z[j] = std::exp(a - std::log(std::exp(a) + std::exp(b)));
        n1 += z[j];
      }
      
      /**
       * M-STEP
       */
      
      soft_update(new_alpha, seq, z.t());
    }
    
    w = n1 / (n*m);
    alpha = new_alpha / n1;
    
    /**
     * Convergence control
     */
    if (hasConverged(cutoff, niter, alpha, convergence, changes)) break;
    
    /**
     * Next iteration
     */
    --niter;
    
  }
  
  return Rcpp::List::create(alpha, convergence, changes, w);
}

