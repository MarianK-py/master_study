import numpy as np
import sys
from sklearn.model_selection import cross_val_score
from sklearn.linear_model import Lasso, LinearRegression

def load(fn):
  f = open(fn)
  n, m = map(int, next(f).strip().split())
  X = []
  y = []
  for l in f:
    its = l.strip().split()
    X.append([float(x) for x in its[:-1]])
    y.append(float(its[-1]))
  return X, y

# I was not sure if short commentary should be just
# comments in code or separate file, so in the end I decided
# to write it in separate text file
def select_features(X, y):
  X = np.array(X)
  y = np.array(y)
  beta = 0.01
  alfa = 0.2
  bestScore = float("inf")
  bestModel = None
  while abs(beta) > 0.000001:
    clf = Lasso(alpha=alfa+beta)
    clf.fit(X, y)
    #print(np.mean(np.absolute(clf.coef_)*(i/100)))
    #y_pred = clf.predict(X)
    #print(np.mean((np.array(y)-y_pred)**2))
    score = -np.mean(cross_val_score(clf, X, y, cv=5, scoring='neg_mean_squared_error'))
    if score < bestScore:
      alfa += beta
      bestScore = score
      bestModel = clf
    else:
      beta *= -0.9
  sortedCoef = np.sort(np.absolute(bestModel.coef_))[::-1]
  step = 1
  thresh = 12
  bestScore = float("inf")
  counter = 0
  while counter < 2:
    Xx = X[:,np.absolute(bestModel.coef_) > sortedCoef[thresh+step]]
    clf2 = LinearRegression().fit(Xx, y)
    score = -np.mean(cross_val_score(clf2, Xx, y, cv=5, scoring='neg_mean_squared_error'))*(thresh+step)**0.8
    if score < bestScore:
      thresh += step
      bestScore = score
      counter = 0
    else:
      counter += 1
      step *= -1
  return np.array(list(range(len(X[0]))))[np.absolute(bestModel.coef_) > sortedCoef[thresh+step]]

X, y = load(sys.argv[1])
#X, y = load("input.txt")
print("Good features", select_features(X, y) )
