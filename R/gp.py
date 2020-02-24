import pymc3 as pm

import theano
import theano.tensor as tt
import numpy as np
import GPy 
Xu_init = 10*np.random.rand(20)

def fit_and_predict(X,y,X_new):
  X= np.array(X).reshape((-1,1))
  y= np.array(y).reshape((-1,1))
  X_new = np.array(X_new).reshape((-1,1))
 
  kernel = GPy.kern.RBF(input_dim=1, variance=1., lengthscale=1.)
  m = GPy.models.GPRegression(X,y,kernel)
  m.optimize(messages=True)
  
  return (m.predict(X_new))

fit_and_predict(np.arange(1,52),np.arange(1,52),np.arange(1,52))
