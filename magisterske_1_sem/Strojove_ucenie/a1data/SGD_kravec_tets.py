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

# This SGD function is pretty slow because it uses python loops instead
# numpy vectorized calculation, but I was not sure how to write this computation
# in a way I would be able to argue that complexity is O(n*k) which is better than O(n*m)
def SGD(XtrainVals, XtrainKeys, ytrain, XtestVals, XtestKeys, ytest, theta0, alfa, thresh):
    # firstly we prepare our iteration counter, constant contain shape of our data and initial error
    i = 0
    n = ytrain.shape[0]
    m = theta0.shape[0]
    err0 = err(XtestVals, XtestKeys, ytest, theta0)
    VkeyValToVec = np.vectorize(keyValToVec, excluded=['m'])
    Vgrad = np.vectorize(grad, excluded=['theta'], signature='(n),(n),(n),(1)->()')
    # we use threshold in 2 ways
    # first is when error is smaller than threshold we end optimization
    while err0 > thresh:
        # now we compute change in theta in one iteration

        thetachange = np.zeros(m)
        # n-times we compute gradient of each row of X (each train sample)
        """
        for j in range(n):
            # when computing gradient we use only non-zero values of sample
            # thanks to that this computation have complexity only O(k) (where k is average number of non-zero variables)
            # we expect that k << m so complexity O(k) is better than O(m)
            grad = np.dot(XtrainVals[j],
                          (np.dot(XtrainVals[j],
                                  theta0[XtrainKeys[j]])-ytrain[j]))
            thetachange[XtrainKeys[j]] -= alfa*grad
        # so computation of new theta is O(n*k) which should be better than O(n*m)
        """
        gradient = np.sum(VkeyValToVec(Vgrad(XtrainVals, XtrainKeys, ytrain, theta0), XtrainKeys, m), axis=0)
        theta1 = theta0 - alfa*gradient
        print(i)
        #theta1 = theta0+thetachange
        # then we compute new error
        err1 = err(XtestVals, XtestKeys, ytest, theta1)
        # if error increase we still use previous value of theta and make alfa smaller (smaller step)
        # else if difference in errors in smaller than threshold we say that change in error is negligible, and we stop optimalization process
        # else we continue with new theta again
        if err1 > err0:
            alfa = alfa * 0.2
        elif err0-err1 < thresh:
            break
        else:
            err0 = err1
            theta0 = theta1
        i += 1
        # after each 10 iteration we print current state of error
        if i % 10 == 0:
            print("Number of iteration: {:3}  Error: {:.3}".format(i, err0))
    return theta0

def err(XtestVals, XtestKeys, ytest, theta):
    err = 0
    # error is computed in similar fashion as new theta so complexity
    # would be O(n'*k) (in this case n' is different from n because we use different dataset
    # however both n and n' would always smaller than original n because they are positive, and orig. n is their sum)
    for j in range(ytest.shape[0]):
        err += (np.dot(XtestVals[j], theta[XtestKeys[j]])-ytest[j])**2
    return err/ytest.shape[0]


# function to compute vector of length m based on list of keys (positions)
# and list of values on those positions
# it's used to convert computed gradient in which we compute only k positions
# to vector of same length theta
def keyValToVec(vals, keys, m):
    res = np.zeros(m)
    res[keys] = vals
    return(res)

def grad(Xvals, Xkeys, y, theta):
    print(Xkeys)
    print(y)
    return np.dot(Xvals,
                  (np.dot(Xvals,
                          theta[Xkeys])-y))

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
    # Firstly we split input data into training and testing datasets
    # we take 20% samples if we have few samples and 1000 if we have a lot as testing data
    n_test = min(n*0.2, 1_000)
    ytest = np.array(y[0:n_test])
    XtestVals = list()
    XtestKeys = list()
    # X values of training data we split into keys (which parameter they contain)
    # and values (values of those parameters)
    XtrainVals = list()
    XtrainKeys = list()
    ytrain = np.array(y[n_test:n])
    for i in range(n_test):
        # to each testing data point we add intercept
        XtestVals.append(np.array(list(X[i].values())+[1]))
        XtestKeys.append(np.array(list(X[i].keys())+[m]))
    for i in range(n-n_test):
        # to each training data point we add intercept
        XtrainVals.append(np.array(list(X[n_test+i].values())+[1]))
        XtrainKeys.append(np.array(list(X[n_test+i].keys())+[m]))


    print("Starting error: {:.3}".format(err(XtestVals, XtestKeys, ytest, np.zeros(m+1))))

    # we compute theta using SGD method
    theta = SGD(XtrainVals, XtrainKeys, ytrain, XtestVals, XtestKeys, ytest, np.zeros(m+1), 0.01, 0.001)

    print("Final error: {:.3}".format(err(XtestVals, XtestKeys, ytest, theta)))

    return theta, err(XtestVals, XtestKeys, ytest, theta)


#n, m, X, y = load_data(sys.argv[1])
n, m, X, y = load_data("test_data.txt")
theta, train_error = fit(n, m, X, y)

print("Training error", train_error)



