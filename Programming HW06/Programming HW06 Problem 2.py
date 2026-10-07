import numpy as np
import matplotlib.pyplot as plt

k_growth = 0.3
k_predation = 0.01
k_death = 0.2
k_hunt = 0.0003
prey_0 = 700
predators_0 = 22
t_0 = 0.0
t_end = 50.0
dt = 1/12     
N = int(t_end / dt)

def f(t, y):
    prey, predators = y
    dprey_dt = k_growth * prey - k_predation * prey * predators
    dpredators_dt = -k_death * predators + k_hunt * prey * predators
    return np.array([dprey_dt, dpredators_dt])

def euler_method(f, t_0, y_0, dt, N):
    t_values = np.zeros(N + 1)
    y_values = np.zeros((N + 1, 2))
    t_values[0] = t_0
    y_values[0] = y_0
    t = t_0
    y = y_0.copy()
    for i in range(N):
        y = y + f(t, y) * dt
        t = t + dt
        t_values[i + 1] = t
        y_values[i + 1] = y
    return t_values, y_values

def heun_method(f, t_0, y_0, dt, N):
    t_values = np.zeros(N + 1)
    y_values = np.zeros((N + 1, 2))
    t_values[0] = t_0
    y_values[0] = y_0
    t = t_0
    y = y_0.copy()
    for i in range(N):
        y_predict = y + f(t, y) * dt
        y = y + (f(t, y) + f(t + dt, y_predict)) / 2 * dt
        t = t + dt
        t_values[i + 1] = t
        y_values[i + 1] = y
    return t_values, y_values

y_0 = np.array([prey_0, predators_0])
t_euler, y_euler = euler_method(f, t_0, y_0, dt, N)
t_heun, y_heun = heun_method(f, t_0, y_0, dt, N)
prey_euler = y_euler[:, 0]
predators_euler = y_euler[:, 1]
prey_heun = y_heun[:, 0]
predators_heun = y_heun[:, 1]

plt.figure(figsize=(8, 5))
plt.plot(t_euler, prey_euler, "-", label="Euler")
plt.plot(t_heun, prey_heun, "-", label="Heun")
plt.xlabel("Time (years)")
plt.ylabel("Prey population")
plt.title("Prey Population vs. Time")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()

plt.figure(figsize=(8, 5))
plt.plot(t_euler, predators_euler, "-", label="Euler")
plt.plot(t_heun, predators_heun, "-", label="Heun")
plt.xlabel("Time (years)")
plt.ylabel("Predator population")
plt.title("Predator Population vs. Time")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()

plt.figure(figsize=(8, 5))
plt.plot(prey_euler, predators_euler, "-", label="Euler")
plt.plot(prey_heun, predators_heun, "-", label="Heun")
plt.xlabel("Prey population")
plt.ylabel("Predator population")
plt.title("Predator-Prey State Space")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()