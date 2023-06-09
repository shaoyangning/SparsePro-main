## AnnotateSusie function
## S. Ning
## 20230602
#' Run susie+annotate in batches
#' @param list_dat_gwas A list of susie-formatted datasets in batches. 
#'  Each batch dataset is a list of SNP's `z_score` from GWAS and corresponding LD matrix `ld_corr`.
#' @param list_snp_gwas A list of SNP names in batches, corresponding to `list_dat_gwas`. 
#' @param dat_annot_susie A dataframe/matrix of SNP annotations, with SNPs matching the order of `list_dat_gwas`. 
#' @param seed A number. Random seed. Default 1001
#' @param n_iter A number. Max number of iteration allowed. Default 20.
#' @param ind_prev_fit. Logic for whether using previous iteration to initialize Susie. Default `TRUE`.
#' @param prior_hyper. A hyperparameters for prior for all SNPs. If `NULL`, all equal to 1.
#' 
#' @returns A list with components \tabular{ll}{
#'    \code{susie_est} \tab A list of final susie fits in batches \cr
#'    \tab \cr
#'    \code{pip_est} \tab A vector of final pip estimates. \cr
#'    \tab \cr
#'    \code{prob_annot_est} \tab A vector of final annotation prior estimate. \cr
#'    \tab \cr
#'    \code{list_susie} \tab A list of susie fit in each iteration. \cr
#'    \tab \cr
#'    \code{list_pois} \tab A list of prior GLM fit in each iteration.  \cr
#'    \tab \cr
#'    \code{list_elbo} \tab A list of ELBO in each iteration.  \cr
#'    \tab \cr
#'    \code{diff_alpha} \tab A matrix difference in alphas by each iteration \cr
#'    \tab \cr
#'    \code{diff_alpha_sum} \tab A matrix difference in sums of alphas by each iteration \cr
#' }


annotate_susie_beta <- function(list_dat_gwas, list_snp_gwas, dat_annot_susie, 
                                seed = 1001, n_iter=20,
                                ind_prev_fit = T, prior_hyper = NULL){
  num_snp_cum <- c(0, cumsum(lapply(list_snp_gwas, length)))
  if(length(prior_hyper) == 0){
    prior_cur_all <- rep(1, times = nrow(dat_annot_susie))
  }

  list_susie <- list()  # list with length n_iter, each is a list of n_batch models
  list_pois <- list() # list with length n_iter
  
  niter_susie <- matrix(NA, n_iter, n_batch)
  diff_alpha <- matrix(NA, n_iter, n_batch)
  diff_alpha_sum <- matrix(NA, n_iter, n_batch)
  list_elbo <- matrix(NA, n_iter, n_batch)
  
  time_log_iter <- numeric(n_iter)
  time_log_glm <- numeric(n_iter)
  time_log_susie <- matrix(NA, n_iter, n_batch)
  
  
  set.seed(seed)
  list_susie_cur <- vector(mode='list', length=n_batch)
  list_susie_prev <- vector(mode='list', length=n_batch)
  alpha_cur <- Inf
  time_cur <- Sys.time()
  
  options(warn=-1)
  
  for(i in 1:n_iter){
    for(j in 1:n_batch){
      time_susie <- Sys.time()
      susie_init <- list_susie_cur[[j]]
      prior_cur <- prior_cur_all[(num_snp_cum[j]+1):num_snp_cum[j+1]]
      if(length(susie_init)!=0){
        susie_init$pi <- prior_cur[[j]]
      }
      susie_cur <- susie_rss(z = list_dat_gwas[[j]]$z.score, 
                             R = list_dat_gwas[[j]]$ld.cor, 
                             prior_weights = prior_cur,
                             s_init = list_susie_prev[[j]],
                             refine=TRUE)
      diff_alpha[i,j] <- sum((list_susie_cur[[j]]$alpha - susie_cur$alpha)^2)
      if(i > 1){
        diff_alpha_sum[i,j] <- sum((apply(list_susie_cur[[j]]$alpha,2,sum) - apply(susie_cur$alpha,2,sum))^2)
      }
      list_susie_cur[[j]] <- susie_cur
      list_elbo[i,j] <- susie_cur$elbo[length(susie_cur$elbo)]
      niter_susie[i,j] <- susie_cur$niter
      time_log_susie[i,j] <- difftime(Sys.time(), time_susie, units = "min")
      gc()
      #if(j%%10 == 0){
      #  print(paste("iter:", i , "batch:", j))
      #}
    }
    list_susie[[i]] <- list_susie_cur
    alpha_cur <- do.call(cbind,lapply(list_susie_cur, function(x){x$alpha}))
    
    alpha_sum_cur <- apply(alpha_cur, 2, sum)
    
    dat_curr <- data.frame(y = alpha_sum_cur, dat_annot_susie)
    
    ## pois regression for prior update
    time_glm <- Sys.time()
    glm_curr <- glm(y~., data = dat_curr, family = "poisson") 
    time_log_glm[i] <- difftime(Sys.time(), time_glm, units = "min")
    
    prior_cur_all <- predict(glm_curr, type = "response")
    ## better way to put it into prior of susie
    list_pois[[i]] <- glm_curr
    
    print(paste0("iter:", i, "done"))
    if(i > 1){
      print(paste0("elbo diff:", sum((list_elbo[i,] - list_elbo[i-1,])^2)))
    }
    
    time_prev <- time_cur
    time_cur <- Sys.time()
    time_log_iter[i] <- difftime(time_cur, time_prev, units = "min")
    gc()
  }
  susie_est <- list_susie_cur
  pip_est <- do.call(c,lapply(susie_est, function(x){x$pip}))
  prob_annot_est <- prior_cur_all
  gc()
  return(list(susie_est, pip_est, prob_annot_est, 
              list_elbo = list_elbo, list_susie = list_susie, list_pois = list_pois, time_log_iter = time_log_iter,
              diff_alpha = diff_alpha, diff_alpha_sum=diff_alpha_sum))
}


