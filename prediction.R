library(tsfknn)
library(quantmod)

possible_sigmas <- c(1e-20,1e-15,1e-10,1e-5,1e-2)

# Utility function for RBF similarity
rbf <- function(x, y, sigma = 1)
{
  exp(- sigma * (x - y) ^ 2)
}


periodic <- function(x, y,h=1, rho=pi/52,eta=1000)
{
  exp(- (sin((rho^h)*(x-y))^2)/(2*(eta^h)))
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
  print (ts)
 # print(ts[top_k_similar_ys+num_lags-1+h])
  print("--------")
  pred_density <- sample(ts[top_k_similar_ys+num_lags-1+h],10000,prob=top_k_similar,replace=T)

  return(pred_density)
}



predict_seasonal_kcde <- function(ts,k,h,num_lags,sigma,epiweek){
  design_matrix <- matrix(NA,nrow=length(ts),ncol=num_lags)
  
  for (i in 0:(num_lags-1)){
    design_matrix[,(4-i)] <- (Lag(ts,i))
  }
  
  design_matrix_complete <- design_matrix[complete.cases(design_matrix),]
  
  
  test_x <- tail(ts,num_lags)
  
  #covariate similarities
  sim_mat <- apply(design_matrix_complete,1,function(x){rbf(x,test_x,sigma = sigma)})
  similarities <- colSums(sim_mat)/sum(colSums(sim_mat))
  names(similarities) <- 1:length(similarities)
  
  similarities <- similarities[1:(length(similarities)-h)]
  periodic_similarities<- unlist(lapply(1:length(ts),function(x){periodic(x,length(ts)+1)}))
  periodic_similarities <- periodic_similarities[num_lags:(length(periodic_similarities)-h)]
  
  similarities <- similarities*periodic_similarities
  
  top_k_similar <- tail(sort(similarities),k)
  top_k_similar_ys <-as.numeric(names(top_k_similar))
 
  
  
  pred_density <- sample(ts[top_k_similar_ys+num_lags-1+h],10000,prob=top_k_similar,replace=T)
  
  return(pred_density)
}


loo_cv <- function(ts,k,h,num_lags,epiweek){
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
        pred_i <- predict_seasonal_kcde(c(test_data_from_other_season,test_data_from_this_season[1:obs_idx]),k,h,num_lags,sigma,epiweek)
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

  best_sigmas <- rep(NA,1)
  for (h in 1:1){
    max_sigma_idx <-  which.max(exp(leave_one_season_out_cv(data,nneighbor,h,4,flu_nat_week)))
    best_sigmas[h] <-possible_sigmas[max_sigma_idx]
  }
  return (best_sigmas)
}




library(HIDDA.forecasting)

chili_data <- HIDDA.forecasting::CHILI
data <- data.frame(unweighted_ili=chili_data)
data$date <- rownames(as.data.frame(chili_data))
data$season <-c(rep(1,35),rep(2:17,each=52),rep(18,20))
data_train <- data[1:(nrow(data)-213),]
data_test <- data[(nrow(data)-213):nrow(data),]



#debug(predict_seasonal_kcde)
periodic_dist <- predict_seasonal_kcde(data_train$unweighted_ili,k = 10,h=1,num_lags = 4,sigma = 1e-5)

### 1  week ahead 
#debug(leave_one_season_out_cv)
bw_optim <-get_k_step_ahead_bw(data_train,50,flu_nat_week)
truth <- rep(NA,nrow(data_test))
pred_dists <- matrix(NA,ncol=10000,nrow=nrow(data_test))

for (test_idx in 1:nrow(data_test)){
  pred_dists[test_idx,] <- predict_seasonal_kcde(c(data_train$unweighted_ili,data_test$unweighted_ili[1:test_idx]),k = 50,h = 1,num_lags = 4,sigma = 5e-10)
  truth[test_idx] <- data_test$unweighted_ili[test_idx]
}
#truth <- data_test$unweighted_ili[1]
scores <- HIDDA.forecasting::scores_sample(truth,pred_dists)

colMeans(scores)

## long term stuff
truth <- rep(NA,30*4)
pred_dists <- matrix(NA,ncol=10000,nrow=30*4)

truth_idx <- 1


for (test_idx in 1:30){
  pred_dists[truth_idx,] <- predict_seasonal_kcde(data_train$unweighted_ili,k = 200,h = test_idx,num_lags = 4,sigma = bw_optim)
  truth[truth_idx] <- data_test$unweighted_ili[test_idx]
  truth_idx <- truth_idx+1
}


for (test_idx in 1:30){
  pred_dists[truth_idx,] <- predict_seasonal_kcde(c(data_train$unweighted_ili,data_test$unweighted_ili[1:53]),k = 200,h = test_idx,num_lags = 4,sigma = bw_optim)
  truth[truth_idx] <- data_test$unweighted_ili[test_idx + 53]
  truth_idx <- truth_idx+1
}

for (test_idx in 1:30){
  pred_dists[truth_idx,] <- predict_seasonal_kcde(c(data_train$unweighted_ili,data_test$unweighted_ili[1:105]),k = 200,h = test_idx,num_lags = 4,sigma = bw_optim)
  truth[truth_idx] <- data_test$unweighted_ili[test_idx + 105]
  truth_idx <- truth_idx+1
}

for (test_idx in 1:30){
  pred_dists[truth_idx,] <- predict_seasonal_kcde(c(data_train$unweighted_ili,data_test$unweighted_ili[1:157]),k = 200,h = test_idx,num_lags = 4,sigma = bw_optim)
  truth[truth_idx] <- data_test$unweighted_ili[test_idx + 157]
  truth_idx <- truth_idx+1
}





scores <- HIDDA.forecasting::scores_sample(truth,pred_dists)

colMeans(scores)
