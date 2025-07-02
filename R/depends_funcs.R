#' @importFrom truncnorm rtruncnorm
#' @importFrom mvtnorm rmvnorm
#' @noRd
sample_lv_ge <- function(X_all, sigma2_e_inv, sigma2_x_inv, x_len, temp) {
  # temp <- rowsum(Y_star - Xbeta - U_all[i_ind,] - V_all[j_ind,], jt_ind, reorder = T)
  for (pp in unique(x_len)) {
    sigma2_X_cond <- 1.0 / (pp * sigma2_e_inv + sigma2_x_inv)
    x_loc <- which(x_len == pp)
    ## temp is already sum
    mu_X_cond <- temp[x_loc, ] * sigma2_e_inv * sigma2_X_cond
    X_all[x_loc, ] <- mu_X_cond + rnorm(n = length(x_loc), sd = sqrt(sigma2_X_cond))
  }
  X_all
}
#' @noRd
calcu_Deviance <- function(Y_star, x_covs, z_covs, i_ind, hit_ind, sh_len, sh_h_mapper,
                           mean_coeffs, corr_coeffs,
                           sigma2_u, sigma2_v, sigma2_e,
                           U_all, VH_all) {
  loglik <- sum(dnorm(Y_star,
    mean = x_covs %*% mean_coeffs + U_all[i_ind, ] + VH_all[hit_ind, ],
    sd = sqrt(sigma2_e)
  ), log = TRUE)
  # loglik <- loglik + sum(dnorm(U_all, mean = 0, sd = sqrt(sigma2_u), log = TRUE))

  # Zbeta <- z_covs[, -(1:3), drop = FALSE] %*% corr_coeffs
  # sh_singles_loc <- which(sh_len == 1)
  # sh_single_h_loc <- match(sh_singles_loc, sh_h_mapper)
  # loglik <- loglik + calc_loglik(
  #   which(sh_len != 1), sh_h_mapper,
  #   z_covs[, 1], Zbeta, VH_all, sigma2_v
  # ) +
  #   sum(dnorm(VH_all[sh_single_h_loc], mean = 0, sd = sqrt(sigma2_v), log = TRUE))
  return(-2 * loglik)
}
#' @noRd
calcu_Deviance_area <- function(Y_star, x_covs, z_covs, i_ind, area_ind, hit_ind, sh_len, sh_h_mapper,
                                mean_coeffs, corr_coeffs,
                                sigma2_u, sigma2_w, sigma2_v, sigma2_e,
                                U_all, W_all, VH_all) {
  loglik <- sum(dnorm(Y_star,
    mean = x_covs %*% mean_coeffs + U_all[i_ind, ] + VH_all[hit_ind, ] + W_all[area_ind, ],
    sd = sqrt(sigma2_e)
  ), log = TRUE)
  # loglik <- loglik + sum(dnorm(U_all, mean = 0, sd = sqrt(sigma2_u), log = T))
  # loglik <- loglik + sum(dnorm(W_all, mean = 0, sd = sqrt(sigma2_w), log = T))

  # Zbeta <- z_covs[, -(1:3), drop = FALSE] %*% corr_coeffs
  # sh_singles_loc <- which(sh_len == 1)
  # sh_single_h_loc <- match(sh_singles_loc, sh_h_mapper)
  # loglik <- loglik + calc_loglik(
  #   which(sh_len != 1), sh_h_mapper,
  #   z_covs[, 1], Zbeta, VH_all, sigma2_v
  # ) +
  #   sum(dnorm(VH_all[sh_single_h_loc], mean = 0, sd = sqrt(sigma2_v), log = T))
  return(-2 * loglik)
}
#' @noRd
calcu_Deviance_area_nohe <- function(Y_star, x_covs, i_ind, area_ind,
                                mean_coeffs,
                                sigma2_u, sigma2_w, sigma2_e,
                                U_all, W_all) {
  loglik <- sum(dnorm(Y_star,
                      mean = x_covs %*% mean_coeffs + U_all[i_ind, ] + W_all[area_ind, ],
                      sd = sqrt(sigma2_e)
  ), log = TRUE)
  # loglik <- loglik + sum(dnorm(U_all, mean = 0, sd = sqrt(sigma2_u), log = T))
  # loglik <- loglik + sum(dnorm(W_all, mean = 0, sd = sqrt(sigma2_w), log = T))
  
  # Zbeta <- z_covs[, -(1:3), drop = FALSE] %*% corr_coeffs
  # sh_singles_loc <- which(sh_len == 1)
  # sh_single_h_loc <- match(sh_singles_loc, sh_h_mapper)
  # loglik <- loglik + calc_loglik(
  #   which(sh_len != 1), sh_h_mapper,
  #   z_covs[, 1], Zbeta, VH_all, sigma2_v
  # ) +
  #   sum(dnorm(VH_all[sh_single_h_loc], mean = 0, sd = sqrt(sigma2_v), log = T))
  return(-2 * loglik)
}
sample_lv_grp <- function(VH_all, sigma2_e_inv, sigma2_v_inv,
                          sh_len, sh_h_mapper, hit_len, Zbeta,
                          z_sh_ind, temp) {
  # temp <- rowsum(Y_star - Xbeta - U_all[i_ind], hit_ind, reorder = T)
  ## calculate posterior covariance matrices for v_s of superhousehold groups
  ## for n_s = 1 (superhousehold contains only 1 house)
  sh_singles_loc <- which(sh_len == 1)
  sh_single_h_loc <- match(sh_singles_loc, sh_h_mapper)
  sigma2_v_cond_singles <- 1.0 / (hit_len[sh_single_h_loc] * sigma2_e_inv
    + sigma2_v_inv)
  mu_v_cond_singles <- temp[sh_single_h_loc, ] * sigma2_e_inv * sigma2_v_cond_singles
  VH_all[sh_single_h_loc, ] <- mu_v_cond_singles +
    rnorm(n = length(sh_singles_loc), sd = sqrt(sigma2_v_cond_singles))

  ## for n_s > 1 (superhousehold contains more than 1 house)
  ## loop each superhousehold loc
  ## currently set correlation matrix as identity
  # for(sh_loc in which(sh_len!=1)){
  #   h_loc <- which(sh_h_mapper == sh_loc)
  #   Sigma_s_inv <- diag(hit_len[h_loc] * sigma2_e_inv)
  #   # Sigma_vs_inv <- sigma2_v_inv * diag(length(h_loc))
  #   z_loc <- which(z_sh_ind == sh_loc)
  #   rho <- Zbeta[z_loc]
  #   R <- diag(length(h_loc))
  #   R[lower.tri(R)] <- rho
  #   R <- R + t(R) - diag(diag(R))
  #   Sigma_vs_inv <- sigma2_v_inv * solve(R)
  #   Sigma_vs_cond <- solve(Sigma_s_inv + Sigma_vs_inv)
  #   mu_vs_cond <- Sigma_vs_cond %*% temp[h_loc,] * sigma2_e_inv
  #   VH_all[h_loc,] <- mu_vs_cond + t(rmvnorm(n=1, sigma = Sigma_vs_cond))
  # }
  update_VH_multi(
    which(sh_len != 1), hit_len, sh_h_mapper, z_sh_ind, Zbeta, VH_all,
    temp, sigma2_e_inv, sigma2_v_inv
  )
  VH_all
}
#' @noRd
sample_params_he <- function(x_covs, z_covs, XtX, Y_star, U_all, VH_all,
                             mean_coeffs, corr_coeffs,
                             sigma2_e, sigma2_u, sigma2_v,
                             i_ind, hit_ind, sh_len, sh_h_mapper,
                             cor_step_size, corr_vs_diag) {
  ## sample coeffs, sigma2_e, sigma2_u, sigma2_v
  p <- nrow(mean_coeffs)
  ## sample coeffs
  mean_coeffs <- sample_coeffs(
    x_covs, XtX, Y_star - U_all[i_ind, ] - VH_all[hit_ind, ],
    matrix(1.0 / sigma2_e)
  )

  ## sample sigma2_e
  sigma2_e <- sample_sigma2(Y_star - x_covs %*% mean_coeffs - U_all[i_ind, ] - VH_all[hit_ind, ])

  ## sample sigma2_u
  sigma2_u <- sample_sigma2(U_all)

  ## sample sigma2_v
  # sigma2_v <- sample_sigma2(VH_all)
  sigma2_v <- sample_sigma2_v_MH(
    corr_coeffs, z_covs, sigma2_v, VH_all,
    sh_len, sh_h_mapper, cor_step_size
  )

  ## sample R_coeffs (R_vs)
  if (!corr_vs_diag) {
    corr_coeffs <- sample_corr_coeffs_MH(
      corr_coeffs, z_covs, sigma2_v, VH_all,
      sh_len, sh_h_mapper, cor_step_size
    )
  }

  list(
    "mean_coeffs" = mean_coeffs,
    "corr_coeffs" = corr_coeffs,
    "sigma2_e" = sigma2_e,
    "sigma2_u" = sigma2_u,
    "sigma2_v" = sigma2_v
  )
}
#' @noRd
sample_params_he_area <- function(x_covs, z_covs, XtX, Y_star, U_all, VH_all, W_all,
                                  mean_coeffs, corr_coeffs,
                                  sigma2_e, sigma2_u, sigma2_v,
                                  i_ind, area_ind, hit_ind, sh_len, sh_h_mapper,
                                  cor_step_size, corr_vs_diag) {
  ## sample coeffs, sigma2_e, sigma2_u, sigma2_v
  p <- nrow(mean_coeffs)
  ## sample coeffs
  mean_coeffs <- sample_coeffs(
    x_covs, XtX, Y_star - U_all[i_ind, ] - VH_all[hit_ind, ] - W_all[area_ind, ],
    matrix(1.0 / sigma2_e)
  )

  ## sample sigma2_e
  sigma2_e <- sample_sigma2(Y_star - x_covs %*% mean_coeffs - U_all[i_ind, ] - VH_all[hit_ind, ] - W_all[area_ind, ])

  ## sample sigma2_u
  sigma2_u <- sample_sigma2(U_all)

  ## sample sigma2_w
  sigma2_w <- sample_sigma2(W_all)

  ## sample sigma2_v
  # sigma2_v <- sample_sigma2(VH_all)
  sigma2_v <- sample_sigma2_v_MH(
    corr_coeffs, z_covs, sigma2_v, VH_all,
    sh_len, sh_h_mapper, cor_step_size
  )

  ## sample R_coeffs (R_vs)
  if (!corr_vs_diag) {
    corr_coeffs <- sample_corr_coeffs_MH(
      corr_coeffs, z_covs, sigma2_v, VH_all,
      sh_len, sh_h_mapper, cor_step_size
    )
  }

  list(
    "mean_coeffs" = mean_coeffs,
    "corr_coeffs" = corr_coeffs,
    "sigma2_e" = sigma2_e,
    "sigma2_u" = sigma2_u,
    "sigma2_v" = sigma2_v,
    "sigma2_w" = sigma2_w
  )
}
#' @noRd
sample_params_he_area_nohe <- function(x_covs, XtX, Y_star, U_all, W_all,
                                  mean_coeffs,
                                  sigma2_e, sigma2_u,
                                  i_ind, area_ind) {
  ## sample coeffs, sigma2_e, sigma2_u, sigma2_v
  p <- nrow(mean_coeffs)
  ## sample coeffs
  mean_coeffs <- sample_coeffs(
    x_covs, XtX, Y_star - U_all[i_ind, ] - W_all[area_ind, ],
    matrix(1.0 / sigma2_e)
  )
  
  ## sample sigma2_e
  sigma2_e <- sample_sigma2(Y_star - x_covs %*% mean_coeffs - U_all[i_ind, ] - W_all[area_ind, ])
  
  ## sample sigma2_u
  sigma2_u <- sample_sigma2(U_all)
  
  ## sample sigma2_w
  sigma2_w <- sample_sigma2(W_all)
  
  ## sample sigma2_v
  # sigma2_v <- sample_sigma2(VH_all)
  # sigma2_v <- sample_sigma2_v_MH(
  #   corr_coeffs, z_covs, sigma2_v, VH_all,
  #   sh_len, sh_h_mapper, cor_step_size
  # )
  
  ## sample R_coeffs (R_vs)
  # if (!corr_vs_diag) {
  #   corr_coeffs <- sample_corr_coeffs_MH(
  #     corr_coeffs, z_covs, sigma2_v, VH_all,
  #     sh_len, sh_h_mapper, cor_step_size
  #   )
  # }
  
  list(
    "mean_coeffs" = mean_coeffs,
    # "corr_coeffs" = corr_coeffs,
    "sigma2_e" = sigma2_e,
    "sigma2_u" = sigma2_u,
    "sigma2_w" = sigma2_w
  )
}
#' @importFrom matrixcalc is.positive.definite
#' @noRd
sample_corr_coeffs_MH <- function(corr_coeffs, z_covs, sigma2_v, VH_all,
                                  sh_len, sh_h_mapper, cor_step_size) {
  Zbeta <- z_covs[, -(1:3), drop = FALSE] %*% corr_coeffs
  loglik <- calc_loglik(
    which(sh_len != 1), sh_h_mapper,
    z_covs[, 1], Zbeta, VH_all, sigma2_v
  )
  # loglik <- 0
  # for(sh_loc in which(sh_len!=1)){
  #   h_loc <- which(sh_h_mapper == sh_loc)
  #   z_loc <- which(z_covs[,1] == sh_loc)
  #   rho <- Zbeta[z_loc]
  #   vs <- VH_all[h_loc]
  #   R <- diag(length(h_loc))
  #   R[lower.tri(R)] <- rho
  #   R <- R + t(R) - diag(diag(R))
  #   Sigma_vs <- sigma2_v * R
  #   loglik <- loglik + dmvnorm(vs, sigma = Sigma_vs, log = T)
  # }
  # cat("cpp_loglik:", loglik1, "R_loglik:", loglik)
  for (l in 1:length(corr_coeffs)) {
    # nonspd_flag <- FALSE
    corr_coeffs_new <- corr_coeffs
    corr_coeffs_new[l] <- corr_coeffs_new[l] + cor_step_size[l] * rnorm(1)
    Zbeta_new <- z_covs[, -(1:3), drop = FALSE] %*% corr_coeffs_new
    # loglik_new <- 0
    # for(sh_loc in which(sh_len!=1)){
    #   h_loc <- which(sh_h_mapper == sh_loc)
    #   z_loc <- which(z_covs[,1] == sh_loc)
    #   rho_new <- Zbeta_new[z_loc]
    #   vs <- VH_all[h_loc]
    #   R_new <- diag(length(h_loc))
    #   R_new[lower.tri(R_new)] <- rho_new
    #   R_new <- R_new + t(R_new) - diag(diag(R_new))
    #   if(is.positive.definite(R_new)){
    #     Sigma_vs_new <- sigma2_v * R_new
    #     loglik_new <- loglik_new + dmvnorm(vs, sigma = Sigma_vs_new, log = T)
    #   }
    #   else
    #     nonspd_flag <- TRUE
    # }
    loglik_new <- calc_loglik(
      which(sh_len != 1), sh_h_mapper,
      z_covs[, 1], Zbeta_new, VH_all, sigma2_v
    )
    if (loglik_new < 0) {
      ## accept probability
      alpha <- min(1, exp(loglik_new - loglik))
      if (runif(1) < alpha) {
        corr_coeffs <- corr_coeffs_new
        loglik <- loglik_new
      }
    }
  }
  return(corr_coeffs)
}
#' @noRd
sample_sigma2_v_MH <- function(corr_coeffs, z_covs, sigma2_v, VH_all,
                               sh_len, sh_h_mapper, cor_step_size) {
  Zbeta <- z_covs[, -(1:3), drop = FALSE] %*% corr_coeffs
  sh_singles_loc <- which(sh_len == 1)
  sh_single_h_loc <- match(sh_singles_loc, sh_h_mapper)
  loglik <- calc_loglik(
    which(sh_len != 1), sh_h_mapper,
    z_covs[, 1], Zbeta, VH_all, sigma2_v
  ) +
    sum(dnorm(VH_all[sh_single_h_loc], mean = 0, sd = sqrt(sigma2_v), log = T))

  corr_cov_num <- length(corr_coeffs)
  sigma2_v_new <- sigma2_v
  sigma2_v_new <- sigma2_v_new + cor_step_size[corr_cov_num + 1] * rnorm(1)
  if (sigma2_v_new < 0) {
    return(sigma2_v)
  }
  loglik_new <- calc_loglik(
    which(sh_len != 1), sh_h_mapper,
    z_covs[, 1], Zbeta, VH_all, sigma2_v_new
  ) +
    sum(dnorm(VH_all[sh_single_h_loc], mean = 0, sd = sqrt(sigma2_v_new), log = T))
  if (loglik_new < 0) {
    ## accept probability
    alpha <- min(1, exp(loglik_new - loglik))
    if (runif(1) < alpha) {
      sigma2_v <- sigma2_v_new
      loglik <- loglik_new
    }
  }
  return(sigma2_v)
}
#' @noRd
my_seconds_to_period <- function(x) {
  days <- round(x %/% (60 * 60 * 24))
  hours <- round((x - days * 60 * 60 * 24) %/% (60 * 60))
  minutes <- round((x - days * 60 * 60 * 24 - hours * 60 * 60) %/% 60)
  seconds <- round(x - days * 60 * 60 * 24 - hours * 60 * 60 - minutes * 60)
  days_str <- ifelse(days == 0, "", paste0(days, "d "))
  hours_str <- ifelse((hours == 0 & days == 0), "", paste0(hours, "hour "))
  minutes_str <- ifelse((minutes == 0 & days == 0 & hours == 0), "", paste0(minutes, "min "))
  seconds_str <- paste0(seconds, "s")
  final_str <- paste0(days_str, hours_str, minutes_str, seconds_str)
  return(final_str)
}
#' @noRd
sample_coeffs <- function(X, XtX, Z, Sigma_inv) {
  K <- ncol(Z)
  N <- nrow(Z)
  # cat("\n dim(Z)=", dim(Z),"\n")
  sig2_beta0 <- 100
  V_beta <- kronecker(Sigma_inv, XtX)
  diag(V_beta) <- diag(V_beta) + 1.0 / sig2_beta0
  V_beta <- solve(V_beta)
  mu_beta <- V_beta %*% kronecker(Sigma_inv, t(X)) %*% as.vector(Z)
  params <- mu_beta + t(rmvnorm(1, sigma = V_beta))
  return(matrix(params, ncol = K))
}
#' @importFrom MCMCpack rinvgamma
#' @noRd
sample_sigma2 <- function(X) {
  ## inverse Gamma IG(alpha, beta) distribution for sigma2
  alpha <- 0.01
  beta <- 0.01
  N <- nrow(X)
  S <- sum(t(X) %*% X)
  return(rinvgamma(n = 1, alpha + 0.5 * N, beta + 0.5 * S))
}
#' @export
res_compare_plot <- function(refit_res, true_val, burn_in) {
  refit_res <- refit_res[(burn_in + 1):nrow(refit_res), ]
  yl <- apply(refit_res, 2, quantile, 0.025)
  yu <- apply(refit_res, 2, quantile, 0.975)
  ## to remove warnings of arrows commmand
  yu <- pmax(yl + mean(abs(yl)) * 0.002, yu)
  x <- seq_len(ncol(refit_res))
  plot(x, true_val,
    ylim = range(c(yl, yu)), pch = 19, cex = .3, xaxt = "n",
    xlab = "parameters", ylab = "95% credible interval and true value"
  )
  axis(1, at = x)
  arrows(x, yl, x, yu, angle = 90, code = 3, length = .04)
}
