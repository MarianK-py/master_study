import sys
import numpy as np
import scipy
import matplotlib.pyplot as plt

# Autor: Marian Kravec

def load(filename):
  f = open(filename)
  X = []
  Y = []
  for l in f:
    its = l.strip().split()
    x = float(its[0])
    y = float(its[1])
    X.append(x)
    Y.append(y)

  return np.array(X), np.array(Y)

# Tuto funkciu treba doprogramovat
# X - 1rozmerne pole vstupov
# Y - 1rozmerne pole ocakavanych vystupov
#
# Fill out this function
# X - 1 dimensional array of inputs
# Y - 1 dimensional array of expected outputs
def fit(X, Y):
  # we have problem y = A*e^(-x/B) first lets calculate natural logarithm of this expression
  # we get:
  # ln(y) = ln(A*e^(-x/B))
  # ln(y) = ln(A)+ln(e^(-x/B))
  # ln(y) = ln(A)+(-x/B)
  # ln(y) = ln(A)+(-1/B)*x
  # now let's substitute few things (y_log = ln(y), theta0 = ln(A), theta1 = -1/B
  # y_log = theta0 + theta1*x -> standard linear regression

  # to make solution more interesting (and because we use quite small datasets)
  # we do cross-validation, we split dataset into 4 parts (number of parts can be changed with parameter k)
  k = 4
  n = Y.shape[0]
  testSize = int(n/k)
  bestErr = float("inf")
  for i in range(k):
      # for each of k runs we prepare out training and testing dataset
      testStart = i*testSize
      testStop = (i+1)*testSize
      Xtrain = np.hstack((X[:testStart], X[testStop:]))
      Ytrain = np.hstack((Y[:testStart], Y[testStop:]))
      Xtest = X[testStart:testStop]
      Ytest = Y[testStart:testStop]
      # we transform out data based on substitution we
      y_log = np.log(Ytrain)
      X_mod = np.c_[np.ones(Xtrain.shape[0]), Xtrain]
      # because we have small dataset we compute theta using closed-form solution to linear regression
      # we hope it will work fine
      theta = np.dot(np.dot(np.linalg.inv(np.dot(X_mod.T, X_mod)), X_mod.T), y_log)
      # after we compute theta by solving linear equation we can compute
      # parameter A and B based on our substitutions we used
      Atrain = np.e**theta[0]
      Btrain = -1/theta[1]

      # even thought we train our model on modified logaritmized version of problem
      # (which means that which means that we optimoze difference of logarithms not values themselves)
      # I decided that I would choose best model using standard MSE of model as originaly described in assignment
      err = np.sum(((predict(Xtest, Atrain, Btrain))-Ytest)**2)/testSize
      if err < bestErr:
          bestErr = err
          A = Atrain
          B = Btrain

  return A, B

# Pomocna funkcia (pre dane, X, A, B vrati predikcie)
# Helper function, for given A, B returns predictions
def predict(X, A, B):
  return A * np.exp(-X/B)

#Avlado, Bvlado = 7, 7
#Avlado, Bvlado = 10, 10
#Avlado, Bvlado = 23, 47
#Avlado, Bvlado = 47, 23
#X, Y = load("test-"+str(Avlado)+"-"+str(Bvlado)+".txt")

X, Y = load(sys.argv[1])

A, B = fit(X, Y)

"""
Abest, Bbest = 1,1
err = np.sum((predict(X, 1, 1)-Y)**2)
for Atry in range(1,100):
    for Btry in range(1,100):
        errNew = np.sum((predict(X, Atry, Btry)-Y)**2)
        if errNew < err:
            Abest, Bbest = (Atry, Btry)
            err = errNew

fig, ax = plt.subplots(1,2)

ax[0].scatter(X,Y)
ax[0].plot(np.sort(X),predict(np.sort(X),A,B), label="Regression")
ax[0].plot(np.sort(X),predict(np.sort(X),Abest,Bbest), label="Grid Search (just very basic)")
ax[0].plot(np.sort(X),predict(np.sort(X),Avlado, Bvlado), label="Teoretical optimum (from Vlado)")
ax[0].legend(loc="upper right")

ax[1].scatter(X,Y)
ax[1].set_yscale('log')
ax[1].plot(np.sort(X),predict(np.sort(X),A,B), label="Regression")
ax[1].plot(np.sort(X),predict(np.sort(X),Abest,Bbest), label="Grid Search (just very basic)")
ax[1].plot(np.sort(X),predict(np.sort(X),Avlado, Bvlado), label="Teoretical optimum (from Vlado)")

plt.show()
"""

print("A: %.3f, B: %.3f" % (A, B))
#print("A: %.3f, B: %.3f" % (Abest, Bbest))
