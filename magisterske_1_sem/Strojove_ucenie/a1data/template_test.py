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
def fit(n, m, X_dict, y, max_epochs=1000, tol=0):
    # Initialize theta with zeros and an additional value for the bias
    theta = np.zeros(m + 1)

    # Create a copy to store the previous theta
    theta_prev = np.copy(theta)

    # Initialize learning rate
    alpha = 0.01

    # Initialize previous error to infinity
    prev_error = np.inf

    for epoch in range(max_epochs):
        total_error = 0

        for i in range(n):
            # Calculate prediction using only non-zero values and add bias
            y_hat = sum([X_dict[i][j] * theta[j] for j in X_dict[i].keys()]) + theta[-1]

            # Calculate error
            error = y_hat - y[i]

            # Update the bias term
            theta[-1] -= alpha * error

            # Update only the relevant coefficients in theta
            for j in X_dict[i].keys():
                gradient = error * X_dict[i][j]

                theta[j] -= alpha * gradient

            total_error += error**2

        # Calculate mean squared error
        mse = total_error / n

        # If the error increased, revert to previous theta, and reduce the learning rate
        if mse > prev_error:
            theta = np.copy(theta_prev)
            alpha /= 2
        else:
            # Update previous theta for the next iteration
            theta_prev = np.copy(theta)

        # If the change in error is very small, break
        if abs(mse - prev_error) < tol:
            break

        prev_error = mse

        print(epoch, mse)

    return theta, mse

n, m, X, y = load_data("./test_data.txt")
#n, m, X, y = load_data(sys.argv[1])
theta, train_error = fit(n, m, X, y)
with open('a.txt', 'w') as f:
    for value in theta:
        f.write(str(value) + '\n')

print("Training error", train_error)