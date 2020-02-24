# The Stan model statement:
cat(
  '
  functions {
  
  // This largely follows the deSolve package, but also includes the x_r and x_i variables.
  // These variables are used in the background.
  
  real[] SI(real t,
  real[] y,
  real[] params,
  real[] x_r,
  int[] x_i) {
  
  real dydt[3];
  
  dydt[1] = - params[1] * y[1] * y[2];
  dydt[2] = params[1] * y[1] * y[2] - params[2] * y[2];
  dydt[3] = params[2] * y[2];
  
  return dydt;
  }
  
  }
  
  data {
  int<lower = 1> n_obs; // Number of days sampled
  int<lower = 1> n_params; // Number of model parameters
  int<lower = 1> n_difeq; // Number of differential equations in the system
  int<lower = 1> n_sample; // Number of hosts sampled at each time point.
  int<lower = 1> n_fake; // This is to generate "predicted"/"unsampled" data
  int<lower = 1> num_seasons;
  int n_obs_total;
  real y_obs[n_obs_total]; // The binomially distributed data
  real t0; // Initial time point (zero)
  real ts[n_obs]; // Time points that were sampled
  int y_obs_idx[n_obs_total];
  int y_season_idx[n_obs_total];
  int length_of_seasons[num_seasons];
  real fake_ts[n_fake]; // Time points for "predicted"/"unsampled" data
  }
  
  transformed data {
  real x_r[0];
  int x_i[0];
  }
  
  parameters {
  real<lower = 0> params[n_params]; // Model parameters
  real<lower = 0, upper = 1> S0; // Initial fraction of hosts susceptible
  real sqrtQ;
  real sqrtQ2;

  vector[52] x;
  matrix[52,6] season_offset;
  }
  
  transformed parameters{
  real y_hat[n_obs, n_difeq]; // Output from the ODE solver
  real y0[n_difeq]; // Initial conditions for both S and I
  
  y0[1] = S0;
  y0[2] = 1 - S0;
  y0[3] = 0;
  
  y_hat = integrate_ode_rk45(SI, y0, t0, ts, params, x_r, x_i);
  
  }
  
  model {
  params ~ normal(0, 2); //constrained to be positive
  S0 ~ normal(0.5, 0.5); //constrained to be 0-1.
  

 

  for (i in 1:n_obs_total){
      y_obs[i] ~ normal(to_vector(y_hat[, 2])[y_obs_idx[i]],.00001); //y_hat[,2] are the fractions infected from the ODE solver
    }
  }
  
  generated quantities {
  // Generate predicted data over the whole time series:
  real fake_I[n_fake, n_difeq];

  fake_I = integrate_ode_rk45(SI, y0, t0, fake_ts, params, x_r, x_i);
  
  }
  
  ', 
  file = "SI_fit.stan", sep="", fill=T)

# FITTING
# For stan model we need the following variables:
t_max <- 52
stan_d = list(n_obs = 52,
              n_params = 2,
              n_difeq = 3,
              n_sample = 1000,
              n_fake = length(1:t_max),
              season_idx = rep(1:52),
              y_obs =train_data,
              y_obs_idx = rep(1:52,length.out=length(train_data)),
              y_season_idx = rep(1:6,each=52)[1:length(train_data)],
              length_of_seasons = c(52,52,52,52,52,2),
              num_seasons =6,
              n_obs_total = length(train_data),
              t0 = 0,
              ts = 1:52/10,
              fake_ts = c(1:52)/10)

# Which parameters to monitor in the model:
params_monitor = c("y_hat", "y0", "params", "fake_I","x","season_offset")

# Test / debug the model:
test = stan("SI_fit.stan",
            data = stan_d,
            pars = params_monitor,
            chains = 1, iter = 10)

# Fit and sample from the posterior
mod = stan(fit = test,
           data = stan_d,
           pars = params_monitor,
           chains = 1,
           warmup = 500,
           iter = 1500)

# You should do some MCMC diagnostics, including:
#traceplot(mod, pars="lp__")
#traceplot(mod, pars=c("params", "y0"))
#summary(mod)$summary[,"Rhat"]

# These all check out for my model, so I'll move on.

# Extract the posterior samples to a structured list:
posts <- extract(mod)
plot(colMeans(posts$fake_I[,,2]) + colMeans(posts$x)+ colMeans(posts$season_offset[,,5]),ylim=c(0,.13),type='p')
lines(augmented_data[,5],col='red')
