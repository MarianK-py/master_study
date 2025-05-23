import numpy as np
from sklearn.decomposition import PCA
from sklearn.svm import SVC
from sklearn.model_selection import cross_val_score
import sys

def load(fn):
  X = []
  y = []
  f = open(fn)
  for i in range(200):
    its = f.readline().strip().split()
    y.append(int(its[-1]))
    X.append(list(map(float, its[:-1])))
  U = [list(map(float, f.readline().strip().split())) for i in range(10000)]
  T = [list(map(float, f.readline().strip().split())) for i in range(500)]
  return np.array(X), np.array(y), np.array(U), np.array(T)

def magic(X, y, U, T):
  # Firstly we normalize our data
  # because we will use PCA and that algorithm usually
  # require normalized data, I am not sure how it's implemented
  # in sklearn, but even if this step is not essential (because
  # maybe sklearn will do it for me) it will not break anything
  U = (U - np.mean(U, axis=0))/np.std(U, axis=0)
  X = (X - np.mean(X, axis=0))/np.std(X, axis=0)
  T = (T - np.mean(T, axis=0))/np.std(T, axis=0)

  comp = 100
  step = 16
  way = 1
  way_check = 0
  best_pca = None
  best_score = 0
  best_comp = 0
  # Now we number of components that give us best cross-validation score
  # To reduce number of attributes we will use PCA which will transform our data
  # and give us most important components
  while step >= 1:
    pca = PCA(n_components=comp)
    pca.fit(U)
    X_red = pca.transform(X)
    clf = SVC(gamma='auto')
    score = np.mean(cross_val_score(clf, X_red, y, cv=5))
    if score > best_score:
      best_score = score
      best_pca = pca
      best_comp = comp
      comp += step*way
    else:
      if way_check == 1:
        way_check = 0
        step = step//2
        way *= -1
      else:
        way_check = 1
        step = step//2
        way *= -1
      comp += step*way

  # Now we train our model on transformed data
  # I decided to use SVC for this, I am not sure why...
  # well it gave me ~0.9 cross validation score so it looked good
  X_r = best_pca.transform(X)
  clf = SVC(gamma='auto')
  clf.fit(X_r, y)

  # Finally we predict test values
  T_r = best_pca.transform(T)
  y_t = clf.predict(T_r)

  return y_t

X, y, U, T = load(sys.argv[1])
#X, y, U, T = load("data.txt")

output = magic(X, y, U, T)
for o in output:
  print(o)

