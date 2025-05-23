import numpy as np
import sys

def load_data(filename):
  f = open(filename)
  n, m = map(int, next(f).strip().split())

  X = []
  y = []
  for l in f:
    its = l.strip().split()
    k = int(its[0])
    row = {}
    for i in range(k):
      row[int(its[i*2+1])] = float(its[i*2+2])
    y.append(float(its[-1]))
    X.append(row)

  return n, m, X, y

def SGM(X, y, theta0, alfa, thresh):
  i = 0
  err0 = np.sum((np.dot(X, theta0) - y)**2)
  while err0 > thresh:
    #for j in range(y.shape[0]):
    #  theta0 = theta0 - alfa*np.dot(X[j, :], (np.dot(X[j, :], theta0)-y[j]))
    theta1 = theta0 - alfa*np.dot(X.T, (np.dot(X, theta0)-y))
    err1 = np.sum((np.dot(X, theta1) - y)**2)
    if err1 > err0:
      alfa = alfa * 0.9
    else:
      err0 = err1
      theta0 = theta1
    i += 1
    if i%10 == 0:
      print(i, err0)
  return theta0

# Doprogramujte tuto funkciu
# Nech vrati dvojicu (theta, trenovacia chyba)
# X je list dictov, kde kazdy dict ma ako kluc index nenuloveho miesta a hodnotu samotnu hodnotu na
# tom mieste, napr.:
# {0: 4.2, 47: 2.3} znamena, ze na indexe 0 mame hodnotu 4.2 a na indexe 47 mame hodnotu 2.3 a
# ostatne su nulove
#
# Fill out this function
# It should return tuple of (theta, traning error)
# X is list of dictionaries, where each dictionary has key the index of non zero element
# and value the value on that place, for example
# {0: 4.2, 47: 2.3} means, that on index 0 we have value 4.2 and on index 47 we have value 2.3 and
# other elements are zero
def fit(n, m, X, y):
  k = int(np.sqrt(n))
  Xmat = np.zeros((k,m+1))
  yvec = np.zeros(k)
  for i in range(k):
    for j in range(k):
      x = np.zeros(m+1)
      x[list(X[k*i+j].keys())] = list(X[k*i+j].values())
      Xmat[i, :] += x
      yvec[i] += y[k*i+j]
  Xmat[:, m] = np.ones(k)*k
  theta = SGM(Xmat, yvec, np.zeros(m+1), 0.0001, 0.0001)
  print(np.sum((np.dot(Xmat, np.zeros(m+1)) - yvec)**2))

  Xmat2 = np.zeros((n,m+1))
  yvec2 = np.zeros(n)
  for i in range(n):
    x = np.zeros(m+1)
    x[list(X[i].keys())] = list(X[i].values())
    Xmat2[i, :] += x
    yvec2[i] += y[i]
  Xmat2[:, m] = np.ones(n)
  print(np.sum((np.dot(Xmat, theta) - yvec)**2))
  print(np.sum((np.dot(Xmat2, theta) - yvec2)**2))
  return theta, np.sum((np.dot(Xmat, theta) - yvec)**2)


#n, m, X, y = load_data(sys.argv[1])
n, m, X, y = load_data("test_data.txt")
theta, train_error = fit(n, m, X, y)

print("Training error", train_error)



