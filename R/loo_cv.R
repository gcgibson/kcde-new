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

possible_sigmas <- c(1e-1,1e0,1e1)
possible_etas  <-c(1e-1,1e0,1e1)
possible_sds <- c(.01,.05,.1,.15)

source('R/utils.R')
source('R/predict.R')

#' @rdname leave_one_season_out_cv
#' @export leave_one_season_out_cv
leave_one_season_out_cv <- function(par,ts_obj,k,h,num_lags){
  library(quantmod)
  seasons <- unique(ts_obj$season)
  lambda <- 1
  ls <- rep(NA, nrow(ts_obj))
  ls_idx <- 1
  for (season in seasons){
          test_data_from_other_season <- ts_obj[ts_obj$season != season,]$val
          test_data_from_this_season <- ts_obj[ts_obj$season == season,]$val
          if (length(test_data_from_this_season) >=h){
          for (obs_idx in 1:(length(test_data_from_this_season)-h)){
            pred_i <- predict_seasonal_kcde(c(test_data_from_other_season,test_data_from_this_season[1:obs_idx]),k,h,num_lags,sigma=par[1],eta = par[2],sd=par[3],lambda=lambda,doBC=F)
            ls[ls_idx] <- max(log(sum(round(pred_i,1) == round(test_data_from_this_season[obs_idx+h],1)  )/length(pred_i)),-10000)
            ls_idx <- ls_idx + 1
          }
        }
  }
  
  return (-mean(ls,na.rm = T))
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
