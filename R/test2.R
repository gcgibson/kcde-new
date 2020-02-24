# The Stan model statement:
cat(
  '
data {
  int<lower=1> N1;
  real x1[N1];
  vector[N1] y1;
  int<lower=1> N2;
  real x2[N2];
  }
  transformed data {
  real delta = 1e-9;
  int<lower=1> N = N1 + N2;
  real x[N];
  for (n1 in 1:N1) x[n1] = x1[n1];
  for (n2 in 1:N2) x[N1 + n2] = x2[n2];
  }
  parameters {
  real<lower=0> rho;
  real<lower=0> alpha;
  real<lower=0> sigma;
  vector[N] eta;
  }
  transformed parameters {
  vector[N] f;
  {
  matrix[N, N] L_K;
  matrix[N, N] K = cov_exp_quad(x, alpha, rho);
  
  // diagonal elements
  for (n in 1:N)
  K[n, n] = K[n, n] + delta;
  
  L_K = cholesky_decompose(K);
  f = L_K * eta;
  }
  }
  model {
  rho ~ inv_gamma(5, 5);
  alpha ~ normal(0, 1);
  sigma ~ normal(0, .1);
  eta ~ normal(0, 1);
  
  y1 ~ normal(f[1:N1], sigma);
  }
  generated quantities {
    vector[N2] y2;
    for (n2 in 1:N2)
    y2[n2] = normal_rng(f[N1 + n2], sigma);
  }', 
  file = "fit.stan", sep="", fill=T)


fit_and_predict_gp <- function(train_data){
  # FITTING
  # For stan model we need the following variables:
  train_data$diff <- c(0,diff(train_data$wili))
  train_data$diff2 <- c(0,0,diff(train_data$wili,lag = 2))
  
  train_data$week <- unlist(lapply(train_data$epiweek,function(x){return (substr(x,5,7))}))
  train_data$season <- unlist(lapply(train_data$epiweek,function(x){
    if (as.numeric(substr(x,5,7)) <=20){
      return (as.numeric(substr(x,1,4))-1)
    } else{
      return (as.numeric(substr(x,1,4)))
    }
    }))

  #we need to upvote this seasons data
  test_season <- tail(unique(train_data$season),1)
  subset_to_repeat <- train_data[train_data$season == test_season,]
  
  train_data_full <- rbind(train_data,rbindlist(replicate(n = 100, expr = subset_to_repeat, simplify = FALSE)
))
  
  library(reticulate)
  use_python("/Users/gcgibson/anaconda/bin/python2.7")
  source_python("/Users/gcgibson/kcde-new/R/gp.py")
  preds <- fit_and_predict(as.numeric(as.character(train_data_full$week)),train_data_full$diff,1:52)
  preds2 <- fit_and_predict(as.numeric(as.character(train_data_full$week)),train_data_full$diff2,1:52)
  ret_list <- list()
  ret_list[[1]] <-preds
  ret_list[[2]] <- preds2
  ggplot(train_data,aes(x=as.numeric(as.character(week)),y=diff,col=as.factor(season))) + geom_line() +
    geom_line(data=data.frame(x=1:52,y=preds[[1]],season=10),aes(x=x,y=y,col=season))
  return (ret_list)
}

