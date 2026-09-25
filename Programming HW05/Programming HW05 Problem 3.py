import numpy as np
import matplotlib.pyplot as plt

def f(t,y):
    return -1000*y + 3000 - 2000*np.exp(-t)

def heun_method(f, t_0, y_0, dt, N):
    t_values = np.zeros(N + 1)
    y_values = np.zeros(N + 1)
    t_values[0] = t_0
    y_values[0] = y_0
    t = t_0
    y = y_0
    for i in range(N):
        y_predict = y + f(t, y) * dt
        y = y + (f(t, y) + f(t + dt, y_predict)) / 2 * dt
        t = t + dt
        t_values[i + 1] = t
        y_values[i + 1] = y
    return t_values, y_values

def exact_solution(t):
    return 3 - 0.998*np.exp(-1000*t) - 2.002*np.exp(-t)

if __name__ == "__main__":
    t_0 = 0.0
    y_0 = 0.0
    t_end = 0.1
    dt_values = np.array([1e-3, 5e-4, 2.5e-4, 1.25e-4, 6.25e-5])
    errors = []
    for dt in dt_values:
        N = int(t_end / dt)
        t_values, y_values = heun_method(f, t_0, y_0, dt, N)
        y_exact = exact_solution(t_values)
        error = np.sqrt(dt * np.sum((y_values - y_exact)**2))
        errors.append(error)
    errors = np.array(errors)
    
    plt.loglog(dt_values, errors, "o-")
    plt.xlabel("Time step Δt")
    plt.ylabel("Error")
    plt.title("Heun's Method: Error vs. Time Step")
    plt.grid(True)
    plt.show()

