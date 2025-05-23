from math import comb
import matplotlib.pyplot as plt

def fuck(n):
    total = 0
    for i in range(1, n):
        c = comb(n, i)
        b = (7*i - 5*(n-i))/n
        temp = c*(0.4**i)*(0.6**(n-i))

        temp2 = i*((b-7)**2)
        temp2 += (n-i)*((b+5)**2)
        total += (temp*temp2)
    return total/n

def fuck2(n):
    total = 0
    for i in range(0, n+1):
        c = comb(n, i)
        b = (7*i - 5*(n-i))/n
        temp = c*(0.4**i)*(0.6**(n-i))

        temp2 = 0.4*((b-7)**2)
        temp2 += 0.6*((b+5)**2)
        total += (temp*temp2)
    return total



n = list(range(1, 100))
err = [fuck2(i) for i in n]

plt.plot(n, err)
plt.xlabel("n")
plt.ylabel("expected testing error")
plt.show()
