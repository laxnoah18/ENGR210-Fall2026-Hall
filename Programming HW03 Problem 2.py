import numpy as np
import matplotlib.pyplot as plt

def der_func(f, x, h):
    return ( ( f(x + h) - f(x - h) )/ (2*h))

def sin_func(x):
    return np.sin(x)

def cos_func(x):
    return np.cos(x)

x_values = np.linspace(0, 2*np.pi, 200)
h_step = float(input("Please input step size: "))
sin_derivative = der_func(sin_func, x_values, h_step)
cos_derivative = der_func(cos_func, x_values, h_step)

if __name__ == "__main__":
    plt.figure(1)
    plt.plot(x_values, sin_derivative)
    plt.xlabel("Derived x values")
    plt.ylabel("Derivatives")
    plt.title("Derived Function of Sine")
    plt.grid(True)
    plt.show()
    
    plt.figure(2)
    plt.plot(x_values, cos_derivative)
    plt.xlabel("Derived x values")
    plt.ylabel("Derivatives")
    plt.title("Derived Function of Cosine")
    plt.grid(True)
    plt.show()
    
    h_arr = np.logspace(0, -6, 7)
    err_arr = []
    for h in h_arr:
        err = der_func(sin_func, 1, h)
        act_der = np.cos(1)
        magn_err = abs(act_der - err)
        err_arr.append(magn_err)
    
    plt.figure(3)
    plt.loglog(h_arr, err_arr)
    plt.show()