// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace std;

//' @useDynLib mvregrp
//' @importFrom Rcpp evalCpp

static double const log2pi = std::log(2.0 * M_PI);

//' @export
// [[Rcpp::export]]
arma::mat generateSymmetricMatrix(const arma::vec& offDiagonal, int n) {
 arma::mat A = arma::eye(n, n);
 arma::uvec lower_indices = arma::trimatl_ind( size(A), -1 );
 A(lower_indices) = offDiagonal;
 arma::mat symmetricMatrix = arma::symmatl(A); // create a full symmetric matrix from the sparse matrix
 return symmetricMatrix;
}
//' @export
// [[Rcpp::export]]
double my_dmvnorm(arma::vec x, arma::mat Sigma){
   int xdim = x.n_elem;
   arma::mat rooti = arma::trans(arma::inv(trimatu(arma::chol(Sigma))));
   double rootisum = arma::sum(log(rooti.diag()));
   arma::vec z = rooti * x ;   
   double loglik = -(double)xdim/2.0 * log2pi - 0.5 * arma::accu(z%z) + rootisum;
   return loglik;
}
//' @export
// [[Rcpp::export]]
double calc_loglik(const arma::vec &sh_len_g1_indices, 
                   const arma::vec &sh_h_mapper, const arma::vec &z_sh_ind, 
              const arma::vec &Zbeta, const arma::vec &VH_all, const double sigma2_v){
  // loglik <- 0
  // for(sh_loc in which(sh_len!=1)){
  //   h_loc <- which(sh_h_mapper == sh_loc)
  //   z_loc <- which(z_covs[,1] == sh_loc)
  //   rho <- Zbeta[z_loc]
  //   vs <- VH_all[h_loc]
  //   R <- diag(length(h_loc))
  //   R[lower.tri(R)] <- rho
  //   R <- R + t(R) - diag(diag(R))
  //   Sigma_vs <- sigma2_v * R
  //   loglik <- loglik + dmvnorm(vs, sigma = Sigma_vs, log = T)
  // }
  double loglik = 0.0;
  for(auto sh_loc : sh_len_g1_indices ){
    arma::uvec h_loc = arma::find(sh_h_mapper == sh_loc);
    arma::uvec z_loc = arma::find(z_sh_ind == sh_loc);
    arma::vec rho = Zbeta(z_loc);
    arma::vec vs = VH_all(h_loc);
    arma::mat Sigma_vs = sigma2_v * generateSymmetricMatrix(rho, h_loc.n_elem);
    if(!Sigma_vs.is_sympd())
      return 1;
    loglik += my_dmvnorm(vs, Sigma_vs);
  }
  return loglik;
}
//' @export
// [[Rcpp::export]]
void update_VH_multi(const arma::vec &sh_len_g1_indices, const arma::vec &hit_len,
                          const arma::vec &sh_h_mapper, const arma::vec &z_sh_ind,
                          const arma::vec &Zbeta, arma::vec &VH_all, const arma::vec &temp,
                          const double sigma2_e_inv, const double sigma2_v_inv){
//   for(sh_loc in which(sh_len!=1)){
//     h_loc <- which(sh_h_mapper == sh_loc)
//     Sigma_s_inv <- diag(hit_len[h_loc] * sigma2_e_inv)
// # Sigma_vs_inv <- sigma2_v_inv * diag(length(h_loc))
//     z_loc <- which(z_sh_ind == sh_loc)
//     rho <- Zbeta[z_loc]
//     R <- diag(length(h_loc))
//     R[lower.tri(R)] <- rho
//     R <- R + t(R) - diag(diag(R))
//     Sigma_vs_inv <- sigma2_v_inv * solve(R)
//     Sigma_vs_cond <- solve(Sigma_s_inv + Sigma_vs_inv)
//     mu_vs_cond <- Sigma_vs_cond %*% temp[h_loc,] * sigma2_e_inv
//     VH_all[h_loc,] <- mu_vs_cond + t(rmvnorm(n=1, sigma = Sigma_vs_cond))
//   }
  for(auto sh_loc : sh_len_g1_indices ){
    arma::uvec h_loc = arma::find(sh_h_mapper == sh_loc);
    arma::uvec z_loc = arma::find(z_sh_ind == sh_loc);
    arma::vec rho = Zbeta(z_loc);
    arma::vec vs = VH_all(h_loc);
    arma::mat Sigma_s_inv = sigma2_e_inv * arma::diagmat(hit_len(h_loc));
    arma::mat Sigma_vs_inv = sigma2_v_inv * arma::inv_sympd(generateSymmetricMatrix(rho, h_loc.n_elem));
    arma::mat Sigma_vs_cond = arma::inv_sympd(Sigma_s_inv + Sigma_vs_inv);
    arma::vec mu_vs_cond = Sigma_vs_cond * temp(h_loc) * sigma2_e_inv;
    VH_all(h_loc) = arma::mvnrnd(mu_vs_cond, Sigma_vs_cond).as_col();
  }
  // return VH_all;
}




