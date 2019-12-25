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
  
  similarities <- similarities[1:(length(similarities)-num_lags)]
  
  
  
  top_k_similar <- tail(sort(similarities),k)
  top_k_similar_ys <-as.numeric(names(top_k_similar))
  
  pred_density <- sample(ts[top_k_similar_ys+num_lags-1+h],10000,prob=top_k_similar,replace=T)

  return(pred_density)
}

loo_cv <- function(ts,k,h,num_lags,epiweek){
  possible_sigmas <- c(.001,.01,.1,1,10,100,1000)
  sigma_results <- rep(NA,length(possible_sigmas))
  sigma_idx <-1
  for (sigma in possible_sigmas){
    mse <- rep(NA,length(ts))
    for (i in 1:length(ts)){
      ts_minus_i <- ts[-i]
      sigma <-1
      pred_i <- predict_kcde(ts_minus_i,k,h,num_lags,sigma,epiweek)
      mse[i] <- (mean(pred_i)-ts[i])^2 
    }
    sigma_results[sigma_idx] <- mean(mse)
    sigma_idx <- sigma_idx +1
  }
  return (sigma_results)
}

flu <-read.csv("/Users/gcgibson/cdcfluforecasts/data-raw/flu_data.csv")
flu_nat <-flu[flu$region=="National",]$weighted_ili
flu_nat_week <-flu[flu$region=="National",]$week

#debug(predict_kcde)
#pred_dist <- predict_kcde(rep(1:4,10),10,1,4,10,flu_nat_week)
#hist(pred_dist)
#mean(pred_dist)
#var(pred_dist)


debug(loo_cv)
print(loo_cv(head(flu_nat,100),10,1,4,flu_nat_week))

