import numpy as np
from keras.models import Sequential
from keras.layers.core import Dense, Activation
from keras.optimizers import RMSprop
from rbflayer import RBFLayer, InitCentersRandom, InitCentersKMeans
import matplotlib.pyplot as plt
from keras.layers import LSTM
import sys

def load_data():

    data = np.loadtxt("data/data.txt")
    X = data[:, :-1]
    y = data[:, -1:]
    return X, y


def fit(X,y,X_test):
  
    model = Sequential()
   # rbflayer = RBFLayer(X.shape[0],
     #                   initializer=InitCentersRandom(X),
     #                   betas=1,
     #                   input_shape=(X.shape[1],))
   # model.add(rbflayer)
  #  original_shape= X.shape
   # X = X.reshape((-1,1,original_shape[1]))
  #  model.add(LSTM(4, input_shape=(1, X.shape[2])))
    #model.add(Dense(1))
   # model.add(Activation('Linear'))
    
    model.add(Dense(1,  input_shape=(X.shape[1],)))
    #model.add(Dropout(0.2))
    #model.add(Dense(1))

    model.compile(loss='mean_squared_error',
                  optimizer=RMSprop())

    model.fit(X, y,
              batch_size=X.shape[0],
              epochs=2000,
              verbose=0)
    
    return (model)

def predict_rbf(X_test,model):
    return (model.predict(np.array(X_test).reshape((1,-1))))
    
