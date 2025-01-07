#include <RcppArmadillo.h>
#include <RcppDist.h>
using namespace Rcpp;

arma::vec arma_ifelse(arma::uvec cond, const arma::vec& a, const arma::vec& b){

  int n = cond.n_elem;
  int n1 = a.n_elem;
  int n2 = b.n_elem;

  if ((n != n1) || (n != n2))
    Rcpp::stop("Ifelse failed. Incompatible dimensions.");

  arma::vec out(n, arma::fill::zeros);
  for (int i = 0; i < n; i++){
    if (cond(i)){
      out(i) = a(i);
    } else{
      out(i) = b(i);
    }
  }
  return out;
}


void update_beta0(
    double& beta0,
    const arma::vec& y_star,
    const arma::vec& u,
    const arma::vec& X_beta,
    const double& sigma2,
    const double& delta,
    const arma::vec& z,
    const double& nu2_0
){

  double A = 1/nu2_0 + sum(u) / sigma2;
  double a = sum(u % (y_star - delta * z - X_beta)) / sigma2;
  beta0 = R::rnorm(a / A, pow(A, -0.5));
}

void update_u(
    arma::vec& u,
    const arma::uvec& which_missing,
    const double& psi,
    const arma::vec& y_star,
    const double& rho
){

  int n_miss = which_missing.n_elem;
  for (int i = 0; i < n_miss; i++){
    u(which_missing(i)) = 0;
    if (y_star(which_missing(i)) < psi){
      u(which_missing(i)) = R::rbinom(1, rho);
    }

  }
}

void update_rho(
    double& rho,
    const arma::vec& u,
    const double& rho_0,
    const double& rho_1
){

  int n = u.n_elem;
  int n1 = sum(u);
  rho = R::rbeta(rho_0 + n1, rho_1 + n - n1);

}

void update_sigma2(
    double& sigma2,
    const arma::vec& y_star,
    const arma::vec& u,
    const arma::vec& current_mean,
    const double& xi_0,
    const double& sigma2_0
){

  double shape_new = (xi_0 + sum(u))/2.0;
  double SSR = sum(u % pow(y_star - current_mean, 2.0));
  double scale_new = (xi_0 * sigma2_0 + SSR) / 2.0;
  sigma2 = scale_new / R::rgamma(shape_new, 1.0);

}

void update_z(
    arma::vec& z,
    const arma::vec& y_star,
    const arma::vec& u,
    const double& beta0,
    const arma::vec& X_beta,
    const double& sigma2,
    const double& delta
){

  int n = X_beta.n_elem;
  int ny = y_star.n_elem;

  if (n != ny)
    Rcpp::stop("Incompatible dimensions");

  arma::vec a = delta * (1/sigma2) * (y_star - X_beta - beta0);
  double A_inv = 1 / (1 + pow(delta,2.0) / sigma2 );
  for (int i = 0; i < n; i++){

    z(i) = r_truncnorm(A_inv * a(i), pow(A_inv, 0.5), 0, pow(10.0, 10.0));

  }
}

void update_delta(
    double& delta,
    const arma::vec& y_star,
    const arma::vec& u,
    const double& beta0,
    const arma::vec& X_beta,
    const double& sigma2,
    const arma::vec& z,
    const double& nu2_d){

  double B = 1/nu2_d + sum(pow(z,2.0) % u) / sigma2;
  double b = sum(u % z % (y_star - X_beta - beta0)) / sigma2;
  delta = R::rnorm(b / B, pow(B, -0.5));

}

void update_y_star(
    arma::vec& y_star,
    const arma::uvec& which_missing,
    const double& psi,
    const arma::vec& u,
    const arma::vec& current_mean,
    const double& sigma2){

  int n_update = which_missing.n_elem;
  int update_index;
  double max;
  for (int i = 0; i < n_update; i++){

    update_index = which_missing(i);
    max = pow(10,10);

    if (u(update_index) == 1){
      max = psi;
    }

    y_star(update_index) =
      r_truncnorm(current_mean(update_index), pow(sigma2, 0.5),
                  -pow(10.0,10.0), max);
  }
}

void refine_beta(
    arma::vec& beta,
    const arma::vec& y_star,
    const double& beta0,
    const arma::mat& X,
    const arma::mat& XtX,
    const arma::vec& gamma,
    const double& sigma2,
    const double& delta,
    const arma::vec& z,
    const arma::mat& nu2_inv_mat
){
  int p1 = int(sum(gamma));
  if (p1 > 0){
    // might be able to simply this code to reduce copying
    arma::uvec which_beta = find(gamma == 1);
    arma::vec y_use = y_star - beta0 - delta * z;
    arma::mat L_0 = nu2_inv_mat.submat(which_beta, which_beta) +
      XtX.submat(which_beta, which_beta) / sigma2;
    L_0 = L_0.i();
    beta.elem(which_beta) =
      mvnrnd(L_0 * (X.cols(which_beta).t() * y_use / sigma2), L_0);
  }
}


double expit(double x){

  return(1 / (1 + exp(-x)));

}

void update_gamma_beta(
    arma::vec& beta,
    arma::vec& gamma,
    int& decision,
    int& move,
    const arma::vec& y_star,
    const arma::vec& u,
    const arma::mat& X,
    const arma::vec& current_mean,
    const arma::vec& sigma_vec,
    const arma::vec& nu2,
    const arma::mat& R,
    const double& omega,
    const double& eta,
    const arma::vec& proposal_mean,
    const double& theta2
){

  // Determine move
  int p = gamma.n_elem;
  int p1 = int(sum(gamma));
  if (p1 == 0){
    move = 0; // Add operation

  } else if (p1 == p){
    move = 1; // Remove operation

  } else{
    move = int(R::runif(0,3)); // Random operation

  }

  // Preliminary things
  arma::uvec which_gamma_1 = find(gamma == 1);
  arma::uvec which_gamma_0 = find(gamma == 0);
  double beta_j_proposal;
  double log_lik_ratio;
  double prior_ratio;
  double proposal_ratio;
  double alpha;
  double decider = R::runif(0,1);
  decision = 0;
  int j;
  int k;
  List out;

  if (move == 0){ // Add

    j = which_gamma_0(int(R::runif(0, p - p1)));
    beta_j_proposal = R::rnorm(proposal_mean(j), pow(theta2, 0.5));
    log_lik_ratio =
      sum(u % (arma::log_normpdf(y_star, current_mean + beta_j_proposal * X.col(j), sigma_vec) -
      arma::log_normpdf(y_star, current_mean , sigma_vec)));

    prior_ratio =
      exp(omega + 2 * eta * sum(R.col(j) % gamma)) *
      arma::normpdf(beta_j_proposal, 0.0, pow(nu2(j), 0.5));
    proposal_ratio =
      ((p - p1) * 1.0) / ((p1 + 1) * 1.0) / arma::normpdf(beta_j_proposal, proposal_mean(j), pow(theta2, 0.5));

    if (p1 + 1 == p) proposal_ratio = proposal_ratio * 3.0;
    if (p1 == 0) proposal_ratio = proposal_ratio / 3.0;

    alpha = prior_ratio * proposal_ratio * exp(log_lik_ratio);
    if (decider < alpha){
      beta(j) = beta_j_proposal;
      gamma(j) = 1;
      decision = 1;
    }

  } else if (move == 1){ // Remove

    j = which_gamma_1(int(R::runif(0, p1)));

    log_lik_ratio =
      sum(u % (arma::log_normpdf(y_star, current_mean - beta(j) * X.col(j), sigma_vec) -
      arma::log_normpdf(y_star, current_mean , sigma_vec)));

    prior_ratio =
      exp(-omega - 2 * eta * sum(R.col(j) % gamma)) /
        arma::normpdf(beta(j), 0.0, pow(nu2(j), 0.5));

    proposal_ratio =
      (p1 * 1.0) / (1.0 * (p - p1 + 1)) * arma::normpdf(beta(j), proposal_mean(j), pow(theta2, 0.5));

    if (p1 == 1) proposal_ratio = proposal_ratio * 3.0;
    if (p1 == p) proposal_ratio = proposal_ratio / 3.0;

    alpha = prior_ratio * proposal_ratio * exp(log_lik_ratio);
    if (decider < alpha){
      beta(j) = 0;
      gamma(j) = 0;
      decision = 1;
    }

  } else { // Swap

    j = which_gamma_0(int(R::runif(0, p - p1)));
    k = which_gamma_1(int(R::runif(0, p1)));

    beta_j_proposal = R::rnorm(proposal_mean(j), pow(theta2, 0.5));

    log_lik_ratio =
      sum(u % (
          arma::log_normpdf(
            y_star,
            current_mean + beta_j_proposal * X.col(j) - beta(k) * X.col(k),
            sigma_vec
          ) -
            arma::log_normpdf(y_star, current_mean , sigma_vec))
      );

    prior_ratio =
      exp(2 * eta * (sum((R.col(j) - R.col(k)) % gamma) - R(k,j))) *
      arma::normpdf(beta_j_proposal, 0.0, pow(nu2(j), 0.5)) /
        arma::normpdf(beta(k), 0.0, pow(nu2(k), 0.5));

    proposal_ratio =
      arma::normpdf(beta(k), proposal_mean(k), pow(theta2, 0.5)) /
        arma::normpdf(beta_j_proposal, proposal_mean(j), pow(theta2, 0.5));

    alpha = prior_ratio * proposal_ratio * exp(log_lik_ratio);
    if (decider < alpha){
      beta(j) = beta_j_proposal;
      gamma(j) = 1;
      beta(k) = 0;
      gamma(k) = 0;
      decision = 1;
    }

  }

}

// [[Rcpp::export]]
void bvs_mcmc(
    arma::vec y, arma::mat X, double psi, List hyper_params, arma::mat R, double theta2,
    int reps, int burnin, int thinning, bool infer_delta, bool refine_betas, bool adaptive,
    bool vs, double adapt_prop, arma::vec& beta0_samples, arma::mat& beta_samples,
    arma::mat& gamma_samples, arma::vec& sigma2_samples,
    arma::vec& delta_samples, arma::vec& rho_samples,
    arma::vec& acceptance, arma::vec& moves, int& pt_check, int pmax, int pmax_draws
){

  // Set up
  int n = y.n_elem;
  int p = X.n_cols;
  arma::uvec which_missing = find(y < psi);
  arma::uvec which_present = find(y >= psi);
  int n_zero = which_missing.n_elem;
  arma::mat XtX = X.t() * X;
  int total_draws = burnin + reps;
  int r_sub;
  if (!vs){
    adaptive = vs;
    refine_betas = !vs;
  }
  int decision = 1;
  int move = 5;

  // Hyper-parameters
  double rho_0 = hyper_params["rho_0"];
  double rho_1 = hyper_params["rho_1"];
  double nu2_0 = hyper_params["nu2_0"];
  double nu2_d = hyper_params["nu2_d"];
  arma::vec nu2 = hyper_params["nu2"];
  double xi_0 = hyper_params["xi_0"];
  double sigma2_0 = hyper_params["sigma2_0"];
  double omega = hyper_params["omega"];
  double eta = hyper_params["eta"];
  arma::mat nu2_inv_mat(p, p, arma::fill::zeros);
  nu2_inv_mat.diag() =  1/nu2;

  // Initial values
  arma::vec y_star = y;
  arma::vec proposal_mean(p, arma::fill::zeros);
  double rho = rho_0 / (rho_0 + rho_1);
  arma::vec u(n, arma::fill::ones);
  for (int i = 0; i < n_zero; i++) u(which_missing(i)) = R::rbinom(1, rho);
  double beta0 = R::rnorm(mean(y.elem(which_present)), 0.5);
  double delta = 0;
  arma::vec z(n, arma::fill::randn); z = abs(z);
  if (infer_delta){
    delta = mean(y.elem(which_present)) - beta0;
  }
  arma::vec gamma(p, arma::fill::ones);
  arma::vec beta(p, arma::fill::zeros);
  for (int j = 0; j < p; j++){
    if (vs){
      gamma(j) = R::rbinom(1, expit(omega));
    }

    if (gamma(j) == 1) beta(j) = R::rnorm(0.0, 0.5);
  }

  arma::vec X_beta = X * beta;
  double sigma2 = var(y.elem(which_present));
  arma::vec current_mean = X_beta + beta0 + delta * z;
  arma::vec sigma_vec(n, arma::fill::value(pow(sigma2, 0.5)));

  // Adaptation
  // (adapt_prop: proportion of burn-in to use in adaptation)
  int adapt_point = int(burnin * (1 - adapt_prop));
  arma::mat adapt_beta_storage(p, burnin - adapt_point);
  arma::vec adapt_beta_mean(p, arma::fill::zeros);
  arma::vec beta_positive_recent = beta;
  arma::vec beta_vec_temp(burnin - adapt_point);
  arma::uvec which_betas(burnin - adapt_point);

  // Things for phase transition check
  int size_counter;

  // Main sampler
  for (int r = 0; r < total_draws; r++){

    // Update y_star (HAS TO COME FIRST DUE TO INITIALIZATION)
    update_y_star(y_star, which_missing, psi, u, current_mean, sigma2);

    // Update u
    update_u(u, which_missing, psi, y_star, rho);

    // Update rho
    update_rho(rho, u, rho_0, rho_1);

    // Update delta
    if (infer_delta){
      update_z(z, y_star, u, beta0, X_beta, sigma2, delta);

      update_delta(delta, y_star, u, beta0, X_beta, sigma2,
                     z, nu2_d);
    }

    // Update beta0
    update_beta0(beta0, y_star, u, X_beta, sigma2, delta, z, nu2_0);

    // Update sigma2
    current_mean = X_beta + beta0 + delta * z;
    update_sigma2(sigma2, y_star, u, current_mean, xi_0, sigma2_0);

    // Adaptation check
    // Throughout burnin: most recent non-zero beta is proposal mean
    // After burnin: proposal mean is the average nonzero beta from the last
    // 25% of the burnin draws
    // Code for this is wrong - maybe an indexing error
    if (adaptive){

      // Storage of betas
      if ((r > adapt_point) & (r <= burnin)){

        adapt_beta_storage.row(r - adapt_point - 1) = beta.t();
      }

      // Iterative update of proposal mean
      if (r < burnin){
        beta_positive_recent = arma_ifelse(beta == 0, beta_positive_recent, beta);
        proposal_mean = beta_positive_recent;
      }

      // Final choice of proposal mean
      if (r == burnin){

        for (int j = 0; j < p; j++){

          proposal_mean(j) = 0;
          beta_vec_temp = adapt_beta_storage.row(j);
          which_betas = find(beta_vec_temp != 0);
          if (which_betas.n_elem > 0){
            proposal_mean(j) = mean(beta_vec_temp.elem(which_betas));
          }


        }
      }
    }

    // Update beta and gamma
    if (vs){
      sigma_vec = arma::ones<arma::vec>(n) * pow(sigma2,0.5);
      update_gamma_beta(
        beta,gamma,decision,move,y_star,u,X,current_mean,sigma_vec,nu2,R,
        omega,eta,proposal_mean,theta2
      );
    }

    // Check for possible phase transition
    if (int(sum(gamma)) > pmax){
      size_counter++;
    } else {
      size_counter = 0;
    }

    if (size_counter > pmax_draws){
      pt_check = 1;
      break;
    }
    // Refine betas
    if (refine_betas){
      refine_beta(beta,y_star,beta0,X,XtX,gamma,sigma2,
                    delta,z,nu2_inv_mat);
    }

    X_beta = X * beta;
    current_mean = X_beta + beta0 + delta * z;

    // Save results
    r_sub = r - burnin + 1;
    if ((r_sub >= 1) & (r_sub % thinning == 0)){
      r_sub = r_sub / thinning - 1;
      beta0_samples(r_sub) = beta0;
      delta_samples(r_sub) = delta;
      sigma2_samples(r_sub) = sigma2;
      rho_samples(r_sub) = rho;
      acceptance(r_sub) = decision;
      moves(r_sub) = move;
      for (int j = 0; j < p; j++){
        beta_samples(r_sub, j) = beta(j);
        gamma_samples(r_sub, j) = gamma(j);

      }

    }

  }
}







