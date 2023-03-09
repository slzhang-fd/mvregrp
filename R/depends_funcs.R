#' @importFrom truncnorm rtruncnorm
#' @importFrom mvtnorm rmvnorm
#' @noRd
#' general function for sampling random effects
sample_lv_ge <- function(X_all, sigma2_e_inv, sigma2_x_inv, x_len, temp){
  # temp <- rowsum(Y_star - Xbeta - U_all[i_ind,] - V_all[j_ind,], jt_ind, reorder = T)
  for(pp in unique(x_len)){
    sigma2_X_cond <- 1.0 / (pp * sigma2_e_inv + sigma2_x_inv)
    x_loc <- which(x_len==pp)
    ## temp is already sum
    mu_X_cond <- temp[x_loc,] * sigma2_e_inv * sigma2_X_cond
    X_all[x_loc,] <- mu_X_cond + sqrt(sigma2_X_cond) * rnorm(length(x_loc))
  }
  X_all
}
sample_lv_grp <- function(VH_all, sigma2_e_inv, sigma2_v_inv,
                          sh_len, sh_h_mapper, hit_len, temp){
  # temp <- rowsum(Y_star - Xbeta - U_all[i_ind], hit_ind, reorder = T)
  
  ## calculate posterior covariance matrices for v_s of superhousehold groups
  ## for n_s = 1 (superhousehold contains only 1 house)
  sh_singles_loc <- which(sh_len==1)
  sh_single_h_loc <- match(sh_singles_loc, sh_h_mapper)
  sigma2_v_cond_singles <- 1.0 / (hit_len[sh_single_h_loc] * sigma2_e_inv
                                  + sigma2_v_inv)
  mu_v_cond_singles <- temp[sh_single_h_loc] * sigma2_e_inv * sigma2_v_cond_singles
  VH_all[sh_single_h_loc] <- mu_v_cond_singles +
    sqrt(sigma2_v_cond_singles) * rnorm(length(sh_singles_loc))
  
  ## for n_s > 1 (superhousehold contains more than 1 house)
  ## loop each superhousehold loc
  for(sh_loc in which(sh_len!=1)){
    h_loc <- which(sh_h_mapper == sh_loc)
    Sigma_s_inv <- diag(hit_len[h_loc]) * sigma2_e_inv
    Sigma_vs_inv <- sigma2_v_inv * diag(length(h_loc))
    Sigma_vs_cond <- solve(Sigma_s_inv + Sigma_vs_inv)
    mu_vs_cond <- temp[h_loc] %*% Sigma_s_inv %*% Sigma_vs_cond
    
    VH_all[h_loc] <- mu_vs_cond + rmvnorm(length(h_loc), sigma = Sigma_vs_cond)
  }
  VH_all
}
#' @noRd
sample_params_he <- function(x_covs, XtX, Y_star, U_all, VH_all, coeffs,
                             sigma2_e, sigma2_u, sigma2_v,
                             i_ind, hit_ind, sh_ind){
  ## sample coeffs, sigma2_e, sigma2_u, sigma2_v
  p <- nrow(coeffs)
  ## sample coeffs
  coeffs <- sample_coeffs(x_covs, XtX, Y_star - U_all[i_ind,] - VH_all[hit_ind,],
                          matrix(1.0 / sigma2_e))
  
  ## sample sigma2_e
  sigma2_e <- sample_sigma2(Y_star - x_covs%*% coeffs - U_all[i_ind,] - VH_all[hit_ind,])
  
  ## sample sigma2_u
  sigma2_u <- sample_sigma2(U_all)
  
  ## sample sigma2_v
  sigma2_v <- sample_sigma2(VH_all)
  
  ## sample R_coeffs (R_vs)
  
  list('coeffs' = coeffs,
       'sigma2_e' = sigma2_e,
       'sigma2_u' = sigma2_u,
       'sigma2_v' = sigma2_v)
}