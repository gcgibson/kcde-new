#####################'
#' Leave one season out cross validation to estimate kernel parameters
#' Inputs:
#' ts_obj = time series object containing a time series under 
#'              ts_obj$val
#'          and a season indicator under
#'               ts_obj$season
#'          and a calendar date under
#'               ts_obj$date
#' k = number of neighbors to use
#' h = horizon 
#' num_lags = 4
#' sigma = bandwith
#' Outputs:
####################'

possible_sigmas <- c(1e-20,1e-19,1e-18,1e-17,1e-16,1e-15,1e-14,1e-13,1e-12,1e-11,1e-10,1e-9,1e-8,1e-7,1e-6,1e-5,1e-4,1e-3,1e-2,1e-1,1e-0,1e1,1e2,1e3,1e4)
possible_etas  <- c(1e-30,1e-20,1e-10,1e0,1e10,1e20,1e30)


source('R/predict.R')

#' @rdname leave_one_season_out_cv
#' @export leave_one_season_out_cv
leave_one_season_out_cv <- function(ts_obj,k,h,num_lags,sd){
  
  sigma_eta_results <- rep(NA,length(possible_sigmas)*length(possible_etas))
  sigma_eta_idx <-1
  seasons <- unique(ts_obj$season)
  sigma_eta_combinations <-matrix(NA,nrow=length(possible_sigmas)*length(possible_etas),ncol=2)
  for (sigma in possible_sigmas){
    for (eta in possible_etas){
        ls <- rep(NA,nrow(ts_obj))
        ls_idx<-1
        for (season in seasons){
          test_data_from_other_season <- ts_obj[ts_obj$season != season,]$val
          test_data_from_this_season <- ts_obj[ts_obj$season == season,]$val
          if (length(test_data_from_this_season) >=h){
          for (obs_idx in 1:(length(test_data_from_this_season)-h)){
            pred_i <- predict_seasonal_kcde(c(test_data_from_other_season,test_data_from_this_season[1:obs_idx]),k,h,num_lags,sigma=sigma,eta = eta,sd=sd)
            ls[ls_idx] <- max(log(sum(round(pred_i,1) <= round(test_data_from_this_season[obs_idx+h],1) + .5 &round(pred_i,1) >= round(test_data_from_this_season[obs_idx+h],1) - .5 )/length(pred_i)),-10000)
            ls_idx <- ls_idx + 1
          }
          }
        }
        sigma_eta_combinations[sigma_eta_idx,] <- c(sigma,eta)
        sigma_eta_results[sigma_eta_idx] <- mean(ls,na.rm=T)
        sigma_eta_idx <- sigma_eta_idx +1

    }
  }
  ret_list <- list()
  ret_list[[1]] <- sigma_eta_results
  ret_list[[2]] <- sigma_eta_combinations
  return (ret_list)
}


#' @title get_k_step_ahead_bw: The my happy function
#' @param x numeric number
#' @param ... other arguments
#'
#' @rdname get_k_step_ahead_bw
#' @export get_k_step_ahead_bw
get_k_step_ahead_bw <- function(data,nneighbor,week,h,num_lags){
  best_sigmas <- rep(NA,1)
  max_sigma_idx <-  which.max(exp(leave_one_season_out_cv(data,nneighbor,h,num_lags)))
  best_sigmas[h] <-possible_sigmas[max_sigma_idx]
  return (best_sigmas)
}
