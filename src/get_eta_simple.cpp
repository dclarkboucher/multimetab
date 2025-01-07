#include <RcppArmadillo.h>
#include <RcppDist.h>
using namespace Rcpp;

// [[Rcpp::export]]
List get_prior_stat(double omega, arma::vec etas, arma::mat R,
                    int draws, int burnin, int thin, double quantile){

  int p = R.n_rows;
  int n_eta = etas.n_elem;
  int total_kept = draws / thin;
  int p1 = 0;
  int p0 = 0;

  // Initial value for gamma
  arma::uvec gamma(p);
  arma::uvec gamma_initial(p);
  int i;
  int j;
  for (int j = 0; j < p; j++){
    gamma_initial(j) = int(R::rbinom(1.0, 0.5));
  }
  int gamma_sum_initial = int(sum(gamma_initial));

  // Set up
  int r_sub;
  int move;

  double eta;
  arma::uvec which_p1;
  arma::uvec which_p0;
  double u;
  double prior_ratio;
  double prop_ratio;
  arma::uvec model_size(total_kept);
  arma::uvec model_size_order(total_kept);
  arma::vec size_quantile(n_eta);
  arma::vec size_q50(n_eta);
  arma::vec size_mean(n_eta);

  // Rprintf("n_eta: %u \n", n_eta);
  for (int k = 0; k < n_eta; k++){

    eta = etas(k);
    gamma = gamma_initial;
    p1 = gamma_sum_initial;

    for (int r = 0; r < draws + burnin; r++){

      // Determine move type
      which_p1 = find(gamma == 1);
      which_p0 = find(gamma == 0);
      p0 = p - p1;
      u = R::runif(0.0,1.0);

      if (p1 == 0){
        move = 0; // Add operation

      } else if (p1 == p){
        move = 1; // Remove operation

      } else{
        move = int(R::runif(0,3)); // Random operation

      }

      if (move == 0){ // Add

        j = which_p0(int(R::runif(0, p0)));
        prior_ratio =
          exp(omega + 2 * eta * sum(R.col(j) % gamma));
        prop_ratio =
          ((p0) * 1.0) / ((p1 + 1) * 1.0);

        if (p1 + 1 == p) prop_ratio = prop_ratio * 3.0;
        if (p1 == 0) prop_ratio = prop_ratio / 3.0;
        if( u < prior_ratio * prop_ratio){
          gamma(j) = 1;
          p1 = p1 + 1;
        }

      } else if (move == 1){

        j = which_p1(int(R::runif(0, p1)));
        prior_ratio =
          exp(-omega - 2 * eta * sum(R.col(j) % gamma));
        prop_ratio =
          (p1 * 1.0) / (1.0 * (p0 + 1)) ;

        if (p1 == 1) prop_ratio = prop_ratio * 3.0;
        if (p1 == p) prop_ratio = prop_ratio / 3.0;
        if( u < prior_ratio * prop_ratio){
          gamma(j) = 0;
          p1 = p1 - 1;
        }

      } else {
        j = which_p0(int(R::runif(0, p0)));
        i = which_p1(int(R::runif(0, p1)));
        prior_ratio =
          exp(2 * eta * (sum((R.col(j) - R.col(i)) % gamma) - R(i,j)));

        if( u < prior_ratio){
          gamma(j) = 1;
          gamma(i) = 0;
        }

      }

      r_sub = r - burnin + 1;
      if ((r_sub >= 1) & (r_sub % thin == 0)){
        r_sub = r_sub / thin - 1;
        model_size(r_sub) = p1;

      }
    }

    size_mean(k) = mean(model_size);
    size_q50(k) = median(model_size);
    model_size_order = sort(model_size);
    size_quantile(k) = model_size_order(int(total_kept * quantile) - 1);
  }

  return(
    List::create(
      Named("eta") = etas,
      Named("mean") = size_mean,
      Named("median") = size_q50,
      Named("quantile") = size_quantile
    )
  );

}
