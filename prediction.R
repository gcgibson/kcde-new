library(tsfknn)
library(quantmod)
# Utility function for RBF similarity
rbf <- function(x, y, sigma = 1)
{
  exp(- sigma * (x - y) ^ 2)
}



# prediction function using KCDE 
predict_kcde <- function(ts,k,h,num_lags,sigma,epiweek){
  design_matrix <- matrix(NA,nrow=length(ts),ncol=num_lags)
  
  for (i in 0:(num_lags-1)){
    design_matrix[,(4-i)] <- (Lag(ts,i))
  }
  
  design_matrix_complete <- design_matrix[complete.cases(design_matrix),]
  

  test_x <- tail(ts,num_lags)
  
  
  sim_mat <- apply(design_matrix_complete,1,function(x){rbf(x,test_x,sigma = sigma)})
  similarities <- colSums(sim_mat)/sum(colSums(sim_mat))
  names(similarities) <- 1:length(similarities)
  
  similarities <- similarities[1:(length(similarities)-h)]
  
  
  
  top_k_similar <- tail(sort(similarities),k)
  top_k_similar_ys <-as.numeric(names(top_k_similar))
  
  pred_density <- sample(ts[top_k_similar_ys+num_lags-1+h],10000,prob=top_k_similar,replace=T)

  return(pred_density)
}

loo_cv <- function(ts,k,h,num_lags,epiweek){
  possible_sigmas <- c(1e-10,.00001,.001,.01,.1,1,10,100,1000)
  sigma_results <- rep(NA,length(possible_sigmas))
  sigma_idx <-1
  for (sigma in possible_sigmas){
    mse <- rep(NA,length(ts))
    for (i in 1:length(ts)){
      ts_minus_i <- ts[-i]
      pred_i <- predict_kcde(ts_minus_i,k,h,num_lags,sigma,epiweek)
      mse[i] <- max(log(sum(round(pred_i,1) == round(ts[i],1))/length(pred_i)),-10)
    }
    sigma_results[sigma_idx] <- mean(mse)
    sigma_idx <- sigma_idx +1
  }
  return (sigma_results)
}

leave_one_season_out_cv <- function(ts_obj,k,h,num_lags,epiweek){
  possible_sigmas <- c(1e-10,.00001,.001,.01,.1,1,10,100,1000)
  sigma_results <- rep(NA,length(possible_sigmas))
  sigma_idx <-1
  seasons <- unique(ts_obj$season)
  for (sigma in possible_sigmas){
    ls <- rep(NA,nrow(ts_obj))
    ls_idx<-1
    for (season in seasons){
      test_data_from_other_season <- ts_obj[ts_obj$season != season,]$unweighted_ili
      test_data_from_this_season <- ts_obj[ts_obj$season == season,]$unweighted_ili
      
      for (obs_idx in 1:(length(test_data_from_this_season)-h)){
        pred_i <- predict_kcde(c(test_data_from_other_season,test_data_from_other_season[1:obs_idx]),k,h,num_lags,sigma,epiweek)
       # p <- ggplot(data=data.frame(x=pred_i),aes(x=x))+ geom_histogram()
        #ggsave(paste0("pred_density-",substr(season,1,4),"-",obs_idx),plot=p,device = "png")
        ls[ls_idx] <- max(log(sum(round(pred_i,1) <= round(test_data_from_this_season[obs_idx+h],1) + .5 &round(pred_i,1) >= round(test_data_from_this_season[obs_idx+h],1) - .5 )/length(pred_i)),-10)
        ls_idx <- ls_idx + 1
      }
    }
    sigma_results[sigma_idx] <- mean(ls,na.rm=T)
    sigma_idx <- sigma_idx +1
  }
  return (sigma_results)
}


get_k_step_ahead_bw <- function(data,nneighbor,week){
  possible_sigmas <- c(1e-10,.00001,.001,.01,.1,1,10,100,1000)
  
  best_sigmas <- rep(NA,33)
  for (h in 1:33){
    max_sigma_idx <-  which.max(exp(leave_one_season_out_cv(flu_data_ma,300,h,4,flu_nat_week)))
    best_sigmas[h] <-possible_sigmas[max_sigma_idx]
  }
  return (best_sigmas)
}
#flu <-read.csv("/Users/Graham/Desktop/cdcfluutils/data-raw/flu_data.csv")
#flu_nat <-flu[flu$region=="National",]$weighted_ili
#flu_nat_week <-flu[flu$region=="National",]$week
library(cdcfluutils)
flu_data <- download_and_preprocess_state_flu_data()
flu_data_ma <- flu_data[flu_data$region == "Massachusetts",]

#debug(predict_kcde)
pred_dist <- predict_kcde(rep(1:4,10),10,6,4,10,flu_nat_week)
hist(pred_dist)
mean(pred_dist)
#var(pred_dist)

#
#debug(leave_one_season_out_cv)
#print(exp(leave_one_season_out_cv(flu_data_ma,300,,4,flu_nat_week)))

#undebug(get_k_step_ahead_bw)


#debug(leave_one_season_out_cv)
#debug(predict_kcde)

get_k_step_ahead_bw(flu_data_ma,300,flu_nat_week)
