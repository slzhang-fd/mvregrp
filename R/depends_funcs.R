#' @importFrom truncnorm rtruncnorm
#' @importFrom mvtnorm rmvnorm
#' @noRd
sample_lv_ge <- function(X_all, sigma2_e_inv, sigma2_x_inv, x_len, temp){
  # temp <- rowsum(Y_star - Xbeta - U_all[i_ind,] - V_all[j_ind,], jt_ind, reorder = T)
  for(pp in unique(x_len)){
    sigma2_X_cond <- 1.0 / (pp * sigma2_e_inv + sigma2_x_inv)
    x_loc <- which(x_len==pp)
    ## temp is already sum
    mu_X_cond <- temp[x_loc,] * sigma2_e_inv * sigma2_X_cond
    X_all[x_loc,] <- mu_X_cond + rnorm(n=length(x_loc), sd=sqrt(sigma2_X_cond))
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
  mu_v_cond_singles <- temp[sh_single_h_loc,] * sigma2_e_inv * sigma2_v_cond_singles
  VH_all[sh_single_h_loc,] <- mu_v_cond_singles +
    rnorm(n=length(sh_singles_loc), sd=sqrt(sigma2_v_cond_singles))
  
  ## for n_s > 1 (superhousehold contains more than 1 house)
  ## loop each superhousehold loc
  ## currently set correlation matrix as identity
  for(sh_loc in which(sh_len!=1)){
    h_loc <- which(sh_h_mapper == sh_loc)
    Sigma_s_inv <- diag(hit_len[h_loc] * sigma2_e_inv)
    Sigma_vs_inv <- sigma2_v_inv * diag(length(h_loc))
    Sigma_vs_cond <- solve(Sigma_s_inv + Sigma_vs_inv)
    mu_vs_cond <- Sigma_vs_cond %*% temp[h_loc,] * sigma2_e_inv
    # cat(length(h_loc), dim(mu_vs_cond), "\n")
    VH_all[h_loc,] <- mu_vs_cond + t(rmvnorm(n=1, sigma = Sigma_vs_cond))
  }
  VH_all
}
#' @noRd
sample_params_he <- function(x_covs, XtX, Y_star, U_all, VH_all, coeffs,
                             sigma2_e, sigma2_u, sigma2_v,
                             i_ind, hit_ind, cor_step_size){
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
#' @noRd
my_seconds_to_period = function(x) {
  days = round(x %/% (60 * 60 * 24))
  hours = round((x - days*60*60*24) %/% (60 * 60))
  minutes = round((x - days*60*60*24 - hours*60*60) %/% 60)
  seconds = round(x - days*60*60*24 - hours*60*60 - minutes*60)
  days_str = ifelse(days == 0, "", paste0(days, "d "))
  hours_str = ifelse((hours == 0 & days == 0), "", paste0(hours, "hour "))
  minutes_str = ifelse((minutes == 0 & days == 0 & hours == 0), "", paste0(minutes, "min "))
  seconds_str = paste0(seconds, "s")
  final_str = paste0(days_str, hours_str, minutes_str, seconds_str)
  return(final_str)
}
#' @noRd
sample_coeffs <- function(X, XtX, Z, Sigma_inv){
  K <- ncol(Z)
  N <- nrow(Z)
  # cat("\n dim(Z)=", dim(Z),"\n")
  sig2_beta0 <- 100
  V_beta <- kronecker(Sigma_inv, XtX)
  diag(V_beta) <- diag(V_beta) + 1.0/sig2_beta0
  V_beta <- solve(V_beta)
  mu_beta <- V_beta %*% kronecker(Sigma_inv, t(X)) %*% as.vector(Z)
  params <- mu_beta + t(rmvnorm(1, sigma = V_beta))
  return(matrix(params, ncol = K))
}
#' @importFrom MCMCpack rinvgamma
#' @noRd
sample_sigma2 <- function(X){
  ## inverse Gamma IG(alpha, beta) distribution for sigma2
  alpha <- 0.01
  beta <- 0.01
  N <- nrow(X)
  S <- sum(t(X) %*% X)
  return( rinvgamma(n=1, alpha + 0.5 * N, beta + 0.5 * S) )
}