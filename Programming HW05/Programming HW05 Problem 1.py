import numpy as np
import matplotlib.pyplot as plt

def euler_method(f, t_0, y_0, dt, N):
    t_values = np.zeros(N + 1)
    y_values = np.zeros(N + 1)
    t_values[0] = t_0
    y_values[0] = y_0
    t = t_0
    y = y_0
    for i in range(N):
        y = y + f(t, y) * dt
        t = t + dt
        t_values[i + 1] = t
        y_values[i + 1] = y
    return t_values, y_values

def decay(t, y):
    return -y

def exact_sol(t):
    return np.exp(-t)

def rms_error(y_approx, y_exact, dt):
    error = y_approx[:-1] - y_exact[:-1]
    return np.sqrt(dt * np.sum(error**2))

if __name__ == "__main__":
    t_0 = 0.0
    y_0 = 1.0
    t_end = 5.0
    dt_values = np.array([1e-1, 1e-2, 1e-3, 1e-4, 1e-5])
    rms_errors = []
    for dt in dt_values:
        N = int(t_end / dt)
        t_values, y_values = euler_method(decay, t_0, y_0, dt, N)
        y_exact = exact_sol(t_values)
        rms = rms_error(y_values, y_exact, dt)
        rms_errors.append(rms)
    rms_errors = np.array(rms_errors)
    
    plt.loglog(dt_values, rms_errors, "o-")
    plt.xlabel("Time step Δt")
    plt.ylabel("RMS error")
    plt.title("Euler Method: RMS Error vs. Time Step")
    plt.grid(True)
    plt.show()