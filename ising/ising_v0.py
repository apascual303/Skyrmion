import numpy as np
import matplotlib.pyplot as plt


def Energy(magnetization, external_field, alpha, beta, num):
	return -MU*num*magnetization*external_field*np.cos(alpha)-J*np.power(magnetization,2)*np.cos(beta)*(num/2)


if __name__ == '__main__':
    PI= 3.141592653589793238462643383279
    MU= 4*PI*10**(-7)
    LX = 100; LY = 100
    J = 1.0; m = 1.0; H = -1.0
    N = LX * LY
    T = np.linspace(0.01, 5.0, 10000)
    a = 0.0
    b = np.random.random(10000)
    E = Energy(m,H,a,b,N)

    plt.scatter(T,E)
    plt.show()