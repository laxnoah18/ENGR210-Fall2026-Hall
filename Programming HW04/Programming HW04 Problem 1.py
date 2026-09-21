import numpy as np
import matplotlib.pyplot as plt

def trap_int(f, x_0, x_N, N):
    h = (x_N - x_0)/N
    Sum = 0
    for i in range(1, N):
        x_i = x_0 + i*h
        Sum += f(x_i)
    integral = h * (f(x_0) / 2 + Sum + f(x_N) / 2)
    return integral
def sin_func(x):
    return np.sin(x)
def cos_func(x):
    return np.cos(x)

if __name__ == "__main__":
    x_values = np.linspace(0, 2*np.pi, 100)
    sin_int = []
    cos_int = []
    for x in x_values:
        sin_int.append(trap_int(np.sin, 0, x, 100))
        cos_int.append(trap_int(np.cos, 0, x, 100))
    
    plt.figure(1)
    plt.plot(x_values, sin_int, label = "Numerical", color = "red", ls='--')
    plt.plot(x_values, 1 - np.cos(x_values), label = "Exact", color = "green", ls=':')
    plt.xlabel("Upper Limit x")
    plt.ylabel("Integrals")
    plt.title("Integral of Sine")
    plt.legend()
    plt.grid()
    plt.show()
    
    plt.figure(1)
    plt.plot(x_values, cos_int, label = "Numerical", color = "red", ls='--')
    plt.plot(x_values, np.sin(x_values), label = "Exact", color = "green", ls=':')
    plt.xlabel("Upper Limit x")
    plt.ylabel("Integrals")
    plt.title("Integral of Cosine")
    plt.legend()
    plt.grid()
    plt.show()
    
    sin_zero = trap_int(np.sin, 0, 2*np.pi, 100)
    sin_two = trap_int(np.sin, 0, np.pi, 100)
    cos_one = trap_int(np.cos, 0, np.pi/2, 100)
    cos_zero = trap_int(np.cos, 0, np.pi, 100)
    
    print("Integral of sin from 0 to 2pi equals", sin_zero)
    print("Integral of sin from 0 to pi equals", sin_two)
    print("Integral of cos from 0 to pi/2 equals", cos_one)
    print("Integral of cos from 0 to pi equals", cos_zero)
    
    N_values = [10, 20, 40, 80, 160, 320]
    sin_zero_err = []
    sin_zero_h = []
    sin_two_err = []
    sin_two_h = []
    cos_one_err = []
    cos_one_h = []
    cos_zero_err = []
    cos_zero_h = []
    for N in N_values:
        h_sin_two = (np.pi-0)/N
        numerical_sin_two = trap_int(np.sin, 0, np.pi, N)
        sin_two_exact = 2
        err_sin_two = abs(numerical_sin_two - sin_two_exact)
        sin_two_h.append(h_sin_two)
        sin_two_err.append(err_sin_two)
        
        h_sin_zero = (2*np.pi - 0)/N
        numerical_sin_zero = trap_int(np.sin, 0, 2*np.pi, N)
        sin_zero_exact = 0
        err_sin_zero = abs(numerical_sin_zero - sin_zero_exact)
        sin_zero_h.append(h_sin_zero)
        sin_zero_err.append(err_sin_zero)
        
        h_cos_one = (np.pi/2 - 0)/N
        numerical_cos_one = trap_int(np.cos, 0, np.pi/2, N)
        cos_one_exact = 1
        err_cos_one = abs(numerical_cos_one - cos_one_exact)
        cos_one_h.append(h_cos_one)
        cos_one_err.append(err_cos_one)
        
        h_cos_zero = (np.pi - 0)/N
        numerical_cos_zero = trap_int(np.cos, 0, np.pi, N)
        cos_zero_exact = 0
        err_cos_zero = abs(numerical_cos_zero - cos_zero_exact)
        cos_zero_h.append(h_cos_zero)
        cos_zero_err.append(err_cos_zero)
    
    plt.figure(3)
    plt.loglog(sin_two_h, sin_two_err, color = "blue")
    plt.loglog(sin_zero_h, sin_zero_err, color = "orange")
    plt.loglog(cos_one_h, cos_one_err, color = "red")
    plt.loglog(cos_zero_h, cos_zero_err, color = "green")
    plt.xlabel("Step Size h")
    plt.ylabel("Absolute Error")
    plt.title("Trapezoidal Rule Truncation Error")
    plt.grid()
    plt.show()