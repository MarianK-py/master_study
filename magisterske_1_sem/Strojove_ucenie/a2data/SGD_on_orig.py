import sys
import numpy as np
import scipy
import matplotlib.pyplot as plt

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

def error(Y, X, A, B, size):
    return np.sum(((predict(X, A, B))-Y)**2)/size
def logError(Y, X, A, B, size):
    return np.sum((np.log(predict(X, A, B))-np.log(Y))**2)/size

# Tuto funkciu treba doprogramovat
# X - 1rozmerne pole vstupov
# Y - 1rozmerne pole ocakavanych vystupov
#
# Fill out this function
# X - 1 dimensional array of inputs
# Y - 1 dimensional array of expected outputs
def fit(X, Y):
    # we have problem y = A*e^(-x/B)
    # we get:
    # y = A*e^(-1/B)*e^(x)
    # y = e^(log(A))*e^(-1/B)*e^(x)
    # y = e^(-log(A)/B)*e^(x)
    # now let's substitute few things (C = e^(-log(A)/B), x_mod = e^x, theta1 = -1/B
    # y_log = theta0 + theta1*x -> standard linear regression
    n = Y.shape[0]
    train  = int(0.9*n)
    Xtrain = X[:train]
    Ytrain = Y[:train]
    Xtest = X[train:]
    Ytest = Y[train:]

    #y_log = np.log(Ytrain)
    X_mod = np.exp(Xtrain)
    theta0 = np.ones(1)
    alfa = 1
    err = np.sum((np.dot(X_mod, theta0) - Ytrain)**2)

    for i in range(10000):
        grad = np.dot(X_mod.T, (np.dot(X_mod, theta0) - Ytrain))
        theta1 = theta0 - alfa*grad
        errNew = np.sum((np.dot(X_mod, theta1) - Ytrain)**2)
        if errNew < err:
            theta0 = theta1
            err = errNew
        else:
            alfa /= 2
    # after we compute theta by solving linear equation we can compute
    # parameter A and B based on our substitutions we used
    logTheta = np.log(theta0[0])
    A = np.e**theta0[0]
    B = -1/theta0[1]

    #err = np.sum(((predict(Xtest, Atrain, Btrain))-Ytest)**2)/testSize
    #err = error(Ytest, Xtest, Atrain, Btrain, testSize)
    #err = logError(Ytest, Xtest, Atrain, Btrain, train)
    """
    y_log = np.log(Y)
    X_mod = np.c_[np.ones(X.shape[0]), X]
    theta = np.dot(np.dot(np.linalg.inv(np.dot(X_mod.T, X_mod)), X_mod.T), y_log)
    # after we compute theta by solving linear equation we can compute
    # parameter A and B based on our substitutions we used
    A = np.e**theta[0]
    B = -1/theta[1]
    """
    return A, B

# Pomocna funkcia (pre dane, X, A, B vrati predikcie)
# Helper function, for given A, B returns predictions
def predict(X, A, B):
    return A * np.exp(-X/B)

#Avlado, Bvlado = 7, 7
#Avlado, Bvlado = 10, 10
#Avlado, Bvlado = 23, 47
Avlado, Bvlado = 47, 23
#X, Y = load(sys.argv[1])
X, Y = load("test-"+str(Avlado)+"-"+str(Bvlado)+".txt")

A, B = fit(X, Y)

Abest, Bbest = 1,1
err = np.sum((predict(X, 1, 1)-Y)**2)
for Atry in range(1,100):
    for Btry in range(1,100):
        errNew = np.sum((predict(X, Atry, Btry)-Y)**2)
        if errNew < err:
            Abest, Bbest = (Atry, Btry)
            err = errNew
plt.scatter(X,Y)
plt.yscale('log')
plt.plot(np.sort(X),predict(np.sort(X),A,B), label="Regression")
plt.plot(np.sort(X),predict(np.sort(X),Abest,Bbest), label="Grid Search")
plt.plot(np.sort(X),predict(np.sort(X),Avlado, Bvlado), label="Teoretical optimum (from Vlado)")
plt.legend(loc="upper right")

plt.show()

print("A: %.3f, B: %.3f" % (A, B))
print("A: %.3f, B: %.3f" % (Abest, Bbest))

#%%
