import numpy as np
import matplotlib.pyplot as plt
from sklearn import svm
import seaborn as sns

np.random.seed(47)

def gen_data(sample_size, noise, decision_boundary):
  X = []
  y = []
  m = []
  for i in range(sample_size):
    x1 = np.random.uniform(-2, 2)
    x2 = np.random.uniform(-2, 2)
    noise_val = np.random.uniform(0,1)
    X.append([x1, x2])
    if (x1 > decision_boundary[0]) != (x2 > decision_boundary[1]):
      if noise_val > noise:
        y.append(1)
        m.append('+')
      else:
        y.append(0)
        m.append('o')
    else:
      if noise_val > noise:
        y.append(0)
        m.append('o')
      else:
        y.append(1)
        m.append('+')
  return np.array(X), np.array(y), m

# DATASET 1: big (1000 samples), with low noise (1% of samples are wrong)
# and decision boudary in the middle of interval from which data is generated
X1, y1, m1 = gen_data(1000, 0.01, [0,0])
X1_t, y1_t, m1_t = gen_data(100, 0.01, [0,0])

# DATASET 2: small (100 samples), noisy (25% of samples are wrong)
# and decision boudary in the middle of x-axis and 3/4 of y-axis
# of interval from which data is generated
X2, y2, m2 = gen_data(100, 0.25, [0,1])
X2_t, y2_t, m2_t = gen_data(100, 0.25, [0,1])

# DATASET 3: medium (300 saples), medium noise (10% of samples are wrong)
# and decision boudary in the middle of x-axis and 3/4 of y-axis
# of interval from which data is generated
X3, y3, m3 = gen_data(300, 0.1, [1,1])
X3_t, y3_t, m3_t = gen_data(100, 0.1, [1,1])

X, y, m = [X1, X2, X3], [y1, y2, y3], [m1, m2, m3]
X_t, y_t, m_t = [X1_t, X2_t, X3_t], [y1_t, y2_t, y3_t], [m1_t, m2_t, m3_t]

h = .02  # step size in the mesh

# Plot the decision boundary. For that, we will assign a color to each
# point in the mesh [x_min, m_max]x[y_min, y_max].
x_min, x_max = X[0][:, 0].min() - 1, X[0][:, 0].max() + 1
y_min, y_max = X[0][:, 1].min() - 1, X[0][:, 1].max() + 1
xx, yy = np.meshgrid(np.arange(x_min, x_max, h), np.arange(y_min, y_max, h))

for i in range(3):
  print("dataset", i)
  fig, ax = plt.subplots(3,3)
  for a, c in enumerate([0.1, 1, 10]):
    for b, g in enumerate([0.1, 1, 10]):

      # we create an instance of SVM and fit out data.
      clf = svm.SVC(kernel='rbf', gamma=g, C=c)
      clf.fit(X[i], y[i])
      Z = clf.predict(np.c_[xx.ravel(), yy.ravel()])

      # Put the result into a color plot
      Z = Z.reshape(xx.shape)
      ax[a,b].pcolormesh(xx, yy, Z, cmap=plt.cm.Paired)

      # Plot also the training points
      ax[a,b].scatter(X[i][:, 0], X[i][:, 1], c=(3*y[i]+1), s=1)
      ax[a,b].axis('tight')
      if a == 2:
        ax[a,b].set_xlabel("gamma={}".format(g))
      if b == 0:
        ax[a,b].set_ylabel("C={}".format(c))
  fig.savefig("dataset_{}_diff_c_gamma.png".format(i+1))

  plt.show()

  tested_vals = 2.**np.array(list(range(-5, 15)))

  s_train = []
  s_test = []
  for j in tested_vals:
    clf = svm.SVC(kernel='rbf', gamma=j, C=1)
    clf.fit(X[i], y[i])
    #print(clf.score(X[i], y[i]))
    s_train.append(clf.score(X[i], y[i]))
    s_test.append(clf.score(X_t[i], y_t[i]))

  plt.plot(tested_vals, s_train)
  plt.ylim([0,1.1])
  plt.xscale("log")
  plt.xlabel("gamma")
  plt.ylabel("training score")
  plt.savefig("dataset_{}_train_gamma_scores.png".format(i+1))
  plt.show()

  plt.plot(tested_vals, s_test)
  plt.ylim([0,1.1])
  plt.xscale("log")
  plt.xlabel("gamma")
  plt.ylabel("validation score")
  plt.savefig("dataset_{}_validation_gamma_scores.png".format(i+1))
  plt.show()

  s_train = []
  s_test = []
  for j in tested_vals:
    clf = svm.SVC(kernel='rbf', gamma=1, C=j)
    clf.fit(X[i], y[i])
    #print(clf.score(X[i], y[i]))
    s_train.append(clf.score(X[i], y[i]))
    s_test.append(clf.score(X_t[i], y_t[i]))

  plt.plot(tested_vals, s_train)
  plt.ylim([0,1.1])
  plt.xscale("log")
  plt.xlabel("C")
  plt.ylabel("training score")
  plt.savefig("dataset_{}_train_C_scores.png".format(i+1))
  plt.show()

  plt.plot(tested_vals, s_test)
  plt.ylim([0,1.1])
  plt.xscale("log")
  plt.xlabel("C")
  plt.ylabel("validation score")
  plt.savefig("dataset_{}_validation_C_scores.png".format(i+1))
  plt.show()

  scores_train = np.zeros([len(tested_vals),len(tested_vals)])
  scores_test = np.zeros([len(tested_vals),len(tested_vals)])
  for k, u in enumerate(tested_vals):
    for j, v in enumerate(tested_vals):
      clf = svm.SVC(kernel='rbf', gamma=u, C=v)
      clf.fit(X[i], y[i])
      scores_train[k, j] = clf.score(X[i], y[i])
      scores_test[k, j] = clf.score(X_t[i], y_t[i])


  sns.heatmap(scores_train, vmin=0, vmax=1, xticklabels=tested_vals, yticklabels=tested_vals, cmap='autumn')
  plt.xlabel("C")
  plt.ylabel("gamma")
  plt.tight_layout()
  plt.savefig("dataset_{}_train_heatmap_scores.png".format(i+1))
  plt.show()

  sns.heatmap(scores_test, vmin=0, vmax=1, xticklabels=tested_vals, yticklabels=tested_vals, cmap='autumn')
  plt.xlabel("C")
  plt.ylabel("gamma")
  plt.tight_layout()
  plt.savefig("dataset_{}_validation_heatmap_scores.png".format(i+1))
  plt.show()

