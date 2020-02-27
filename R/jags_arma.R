fit_and_predict_jags <- function(data,epiweeks,params_ar){
  epiweeks <- as.numeric(epiweeks)
model.loc = "arma.txt"
model_code <- cat('
                  model {
                  # Set up residuals
                  for(t in 1:max(p,q)) {
                  eps[t] <- d[t] - alpha
                  } 
                  
                  
                  X[1] <- 0 

                  for (t in 2:T){
                    X[t] ~ dnorm(X[t-1],10);
                  }

                  S[1] <- 0
                  
                  for (t in 2:53){
                      S[t] ~ dnorm(S[t-1],100);
                  }
                  
                  for (t in (max(p,q)+1):T) {
                  d[t] ~ dnorm(alpha + ar_mean[t] + ma_mean[t]  + X[t] + S[week[t]] , sigma)
                  ma_mean[t] <- inprod(theta, eps[(t-q):(t-1)])
                  ar_mean[t] <- inprod(phi, d[(t-p):(t-1)])
                  eps[t] <- d[t] - alpha - ar_mean[t] - ma_mean[t]
                  }
                  # Likelihood
                  # Priors
                  alpha ~ dnorm(0.0,0.01)
                  alpha_s ~ dnorm(0.0,0.01)
                  sigma ~ dgamma(1,1);

                  gamma_1 ~ dunif(0, 10)
                  gamma_2 ~ dunif(0, 1)
                    for (i in 1:q) {
                        theta[i] ~ dnorm(0.0,0.01)
                    }
                    for(i in 1:p) {
                        phi[i] ~ dnorm(0.0,0.01)
                    }
                  }',file=model.loc)
library(R2jags)

#perform box-cox transformation
bc_gamma  <- .5
bc_params <- car::powerTransform(data + bc_gamma, family = "bcPower")
bc_params <- list(
  lambda = bc_params$lambda,
  gamma = bc_gamma)

transformed_data <- car::bcPower(
  U = data + bc_params$gamma,
  lambda = bc_params$lambda)

#need a utility function to get all weeks ahead
weeks_to_model <- c(44:52,1:17)
current_week <- tail(epiweeks,1)
current_week_idx <- which(weeks_to_model==current_week)
weeks_to_simulate <- weeks_to_model[5:length(weeks_to_model)]

jags.data = list(d= c(transformed_data,rep(NA,length(weeks_to_simulate))),p=params_ar[1],q=params_ar[2],T=length(transformed_data)+length(weeks_to_simulate),week=c(epiweeks,weeks_to_simulate))
jags.params = c("d")
mod_ar1_intercept = jags(jags.data, parameters.to.save = jags.params, 
    model.file = model.loc, n.chains = 3, n.burnin = 5000, n.thin = 1, 
    n.iter = 10000, DIC = TRUE)
posterior <-mod_ar1_intercept$BUGSoutput$sims.matrix[,paste0("d[",length(data):(length(data)+length(weeks_to_simulate)),"]")][1:10000,]

return (invert_bc_transform(posterior,bc_params$lambda,bc_params$gamma))
}


invert_bc_transform <- function(b, lambda, gamma) {
  ## Two steps: 1) undo box-cox 2) subtract offset gamma
  
  ## 1) undo box-cox
  if(abs(lambda) <= 1e-10) {
    z <- exp(b)
  } else {
    z <- (lambda * b + 1)^(1 / lambda)
  }
  
  ## 2) Subtract gamma
  return(z - gamma)
}