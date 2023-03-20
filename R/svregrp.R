#' Gibbs sampling of the grouped household effect model
#'
#' This function takes response, multilevel indices as input and
#' returns MCMC sample chains of model parameters and latent variables.
#'
#' @param Y_star A NN * 1 matrix containing continuous responses.
#' @param x_covs A NN * p matrix containing covariates for the mean model.
#' @param z_covs A NN * q matrix containing covariates for the correlation model.
#' @param i_ind A vector containing indices for individual.
#' @param sh_ind A vector containing indices for super-household.
#' @param hit_ind A vector containing indices for unique household.
#' @param max_steps Length of MCMC chains to be drawn.
#' @param cor_step_size A vector length q specifying step sizes for q coefficients 
#'        of correlation model.
#' @param corr_vs_diag Whether set the grouped household effects to be independent,
#'        default value is FALSE.
#' @return A dataframe with columns of MCMC chains of sampled parameters and latent variables.
#' @export
svregrp_Gibbs <- function(Y_star, x_covs, z_covs,
                          i_ind, sh_ind, hit_ind, 
                          max_steps, cor_step_size,
                          corr_vs_diag = FALSE){
  K <- ncol(Y_star)
  ## record number
  rcd_num <- nrow(Y_star)
  ## individual number
  ind_num <- max(i_ind)
  ## number of super-household
  sh_num <- max(sh_ind)
  ## number of households
  hit_num <- max(hit_ind)
  
  mean_cov_num <- ncol(x_covs)
  corr_cov_num <- ncol(z_covs) - 3
  
  ## table for unique superhousehold id -> household id
  tb_sh_hit <- unique(data.frame(sh_ind, hit_ind))
  sh_h_mapper <- tb_sh_hit$sh_ind
  
  ## length for each superhouseholds
  sh_len <- aggregate(tb_sh_hit$sh_ind, by=list(tb_sh_hit$sh_ind), length)$x
  ## Functions to map ind to # of records
  u_len <- aggregate(i_ind, by=list(i_ind), length)$x
  hit_len <- aggregate(hit_ind, by=list(hit_ind), length)$x
  
  ## initialize random variables
  # Y_star <- matrix(0, rcd_num, K)
  U_all <- matrix(rnorm(ind_num*K),ind_num,K)
  VH_all <- matrix(rnorm(hit_num*K),hit_num,K)
  
  ## initialize parameters
  mean_coeffs <- matrix(0, mean_cov_num, K)
  corr_coeffs <- matrix(0, corr_cov_num, K)
  # corr_coeffs <- matrix(c(0.1, 0.01, 0.1, 0.1, 0.1, -0.01, 0.02), ncol = 1)
  sigma2_e <- sigma2_v <- sigma2_u <- 1
  
  mean_coeffs_all <- matrix(0, max_steps, mean_cov_num*K)
  corr_coeffs_all <- matrix(0, max_steps, corr_cov_num*K)
  sigma2_e_all <- sigma2_v_all <- sigma2_u_all <- rep(0, max_steps)
  XtX <- t(x_covs)%*%x_covs
  
  ## for plotting progress bar
  init <- numeric(max_steps)
  end <- numeric(max_steps)
  extra <- 6
  width <- 30
  time <- remainining <- 0
  rejection_rate <- rep(0, 7)
  for(iter in 1:max_steps){
    init[iter] <- Sys.time()
    step <- round(iter / max_steps * (width - extra))
    text <- sprintf('\r|%s%s|% 3s%% | Execution time:%s | Estimated time remaining:%s | rejection rate:%s       ',
                    strrep('=', step), strrep(' ', width - step - extra),
                    round(iter / max_steps * 100), my_seconds_to_period(time),
                    my_seconds_to_period(remainining), paste(round(rejection_rate,3), collapse=","))
    cat(text)
    sigma2_e_inv <- 1.0 / sigma2_e
    sigma2_u_inv <- 1.0 / sigma2_u
    sigma2_v_inv <- 1.0 / sigma2_v
    
    Xbeta <- x_covs %*% mean_coeffs
    Zbeta <- z_covs[,4:(3+corr_cov_num)] %*% corr_coeffs
    ## stochastic E step
    # sample U
    U_all <- sample_lv_ge(U_all, sigma2_e_inv, sigma2_u_inv, u_len,
                          rowsum(Y_star - Xbeta - VH_all[hit_ind,],
                                 i_ind, reorder = T))
    # # sample V
    VH_all <- sample_lv_grp(VH_all, sigma2_e_inv, sigma2_v_inv,
                            sh_len, sh_h_mapper, hit_len, Zbeta, z_covs[,1],
                            rowsum(Y_star - Xbeta - U_all[i_ind,],
                                   hit_ind, reorder = T))
    # update parameters
    params <- sample_params_he(x_covs, z_covs, XtX, Y_star, U_all, VH_all, 
                               mean_coeffs, corr_coeffs,
                               sigma2_e, sigma2_u, sigma2_v, 
                               i_ind, hit_ind, sh_len, sh_h_mapper,
                               cor_step_size, corr_vs_diag)
    mean_coeffs <- params$mean_coeffs
    corr_coeffs <- params$corr_coeffs
    sigma2_u <- params$sigma2_u
    sigma2_v <- params$sigma2_v
    sigma2_e <- params$sigma2_e
    
    # store results
    mean_coeffs_all[iter,] <- c(mean_coeffs)
    corr_coeffs_all[iter,] <- c(corr_coeffs)
    sigma2_e_all[iter] <- sigma2_e
    sigma2_u_all[iter] <- sigma2_u
    sigma2_v_all[iter] <- sigma2_v
    
    ## progress bar
    end[iter] <- Sys.time()
    time <- round(sum(end - init), 0)
    est <- max_steps * (mean(end[end != 0] - init[init != 0])) - time
    remainining <- round(est, 0)
    cat(if (iter == max_steps) '\n' else '\r')
    if(!corr_vs_diag){
      if(iter %% 100 == 0)
        rejection_rate = colMeans(diff(corr_coeffs_all[(iter-99):iter,])==0)
    }
  }
  return(list('mean_coeffs_all'=mean_coeffs_all,
              'corr_coeffs_all'=corr_coeffs_all,
              'sigma2_e_all'=sigma2_e_all,
              'sigma2_u_all'=sigma2_u_all,
              'sigma2_v_all'=sigma2_v_all,
              'U_all'=U_all,
              'VH_all'=VH_all))
}
