import pymc3 as pm
import numpy as np
import pandas as pd
from theano import shared, tensor as tt

def norm_cdf(z):
    return 0.5 * (1 + tt.erf(z / np.sqrt(2)))

def stick_breaking(v):
    return v * tt.concatenate([tt.ones_like(v[:, :1]),
                               tt.extra_ops.cumprod(1 - v, axis=1)[:, :-1]],
                              axis=1)
                              
def run_dpp2(x_lidar,y_lidar,x_new):                      
  N = len(y_lidar)
  K = 20

  x_lidar = shared(np.array(x_lidar).reshape((-1,1)), broadcastable=(False, True))
  with pm.Model() as model:
      alpha = pm.Normal('alpha', 0., 5., shape=K)
      beta = pm.Normal('beta', 0., 5., shape=K)
      v = norm_cdf(alpha + beta * x_lidar)
      w = pm.Deterministic('w', stick_breaking(v))
      
      gamma = pm.Normal('gamma', 0., 10., shape=K)
      delta = pm.Normal('delta', 0., 10., shape=K)
      mu = pm.Deterministic('mu', gamma + delta * x_lidar)
      
  with model:
      tau = pm.Gamma('tau', 1., 1., shape=K)
      tau_det = pm.Deterministic('tau_det', gamma + delta * x_lidar)

      obs = pm.NormalMixture('obs', w, mu, tau=tau, observed=y_lidar)
      
  SAMPLES = 20000
  BURN = 10000
  
  with model:
      step = pm.Metropolis()
      trace = pm.sample(SAMPLES, step, chains=1, tune=BURN)
      
  PP_SAMPLES = 10000
  
  lidar_pp_x = np.array(x_new).reshape(-1)
  x_lidar.set_value(lidar_pp_x[:, np.newaxis])
  
  with model:
      pp_trace = pm.sample_ppc(trace, PP_SAMPLES)
      
  return pp_trace
  
