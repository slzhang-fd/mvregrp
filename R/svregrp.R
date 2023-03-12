#' @export
svregrp_Gibbs <- function(Y_star, x_covs, i_ind, sh_ind,
                          hit_ind, 
                          max_steps, cor_step_size = 0.01){
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
  # corr_cov_num <- ncol(z_covs)
  
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
  coeffs <- matrix(0, mean_cov_num, K)
  sigma2_e <- sigma2_v <- sigma2_u <- 1
  
  coeffs_all <- matrix(0, max_steps, mean_cov_num*K)
  sigma2_e_all <- sigma2_v_all <- sigma2_u_all <- rep(0, max_steps)
  # Y_pos_ind <- Y==1
  # Y_zero_ind <- Y==0
  XtX <- t(x_covs)%*%x_covs
  
  ## for plotting progress bar
  init <- numeric(max_steps)
  end <- numeric(max_steps)
  extra <- 6
  width <- 30
  time <- remainining <- 0
  rejection_rate <- 0
  for(iter in 1:max_steps){
    init[iter] <- Sys.time()
    step <- round(iter / max_steps * (width - extra))
    text <- sprintf('|%s%s|% 3s%% | Execution time:%s | Estimated time remaining:%s | rejection rate:%.3f       ',
                    strrep('=', step), strrep(' ', width - step - extra),
                    round(iter / max_steps * 100), my_seconds_to_period(time),
                    my_seconds_to_period(remainining), rejection_rate)
    cat(text)
    sigma2_e_inv <- 1.0 / sigma2_e
    sigma2_u_inv <- 1.0 / sigma2_u
    sigma2_v_inv <- 1.0 / sigma2_v
    
    Xbeta <- x_covs %*% coeffs
    ## stochastic E step
    # sample U
    U_all <- sample_lv_ge(U_all, sigma2_e_inv, sigma2_u_inv, u_len,
                          rowsum(Y_star - Xbeta - VH_all[hit_ind,],
                                 i_ind, reorder = T))
    # # sample V
    VH_all <- sample_lv_grp(VH_all, sigma2_e_inv, sigma2_v_inv,
                            sh_len, sh_h_mapper, hit_len,
                            rowsum(Y_star - Xbeta - U_all[i_ind,],
                                   hit_ind, reorder = T))
    # update parameters
    params <- sample_params_he(x_covs, XtX, Y_star, U_all, VH_all, coeffs,
                               sigma2_e, sigma2_u, sigma2_v, 
                               i_ind, hit_ind, cor_step_size)
    coeffs <- params$coeffs
    sigma2_u <- params$sigma2_u
    sigma2_v <- params$sigma2_v
    sigma2_e <- params$sigma2_e
    
    # store results
    coeffs_all[iter,] <- c(coeffs)
    sigma2_e_all[iter] <- sigma2_e
    sigma2_u_all[iter] <- sigma2_u
    sigma2_v_all[iter] <- sigma2_v
    
    ## progress bar
    end[iter] <- Sys.time()
    time <- round(sum(end - init), 0)
    est <- max_steps * (mean(end[end != 0] - init[init != 0])) - time
    remainining <- round(est, 0)
    cat(if (iter == max_steps) '\n' else '\r')
    # if(iter %% 100 == 0 && K > 1)
    #   rejection_rate = mean(diff(Sigma_e_all[(iter-99):iter,2])==0)
    
  }
  return(list('coeffs_all'=coeffs_all,
              'sigma2_e_all'=sigma2_e_all,
              'sigma2_u_all'=sigma2_u_all,
              'sigma2_v_all'=sigma2_v_all,
              'U_all'=U_all,
              'VH_all'=VH_all))
}
