import GPy
import numpy as np

def fit_gp(X,y):
  X = np.array(X).reshape((-1,1))
  y = np.array(y).reshape((-1,1))
  
  kernel = GPy.kern.RBF(input_dim=1, variance=1., lengthscale=1.)
  m = GPy.models.GPRegression(X,y,kernel)
  #m.het_Gauss.variance = .05
  #m.het_Gauss.variance.fix()
  m.optimize()
  return (m)
def predict_gp(m,X):
  X = np.array(X).reshape((-1,1))
  
  return (m._raw_predict(X))
  
