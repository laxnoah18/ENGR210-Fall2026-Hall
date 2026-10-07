import numpy as np
import matplotlib.pyplot as plt

m = 2.0
k = 20.0
x0 = 0.5
v0 = 0.0
t0 = 0.0
t_end = 200.0
dt = 0.1
N = int(t_end / dt)

def f(t, state):
    x, v = state
    dxdt = v
    dvdt = -(k / m) * x
    return np.array([dxdt, dvdt])

def exact_solution(t):
    omega = np.sqrt(k / m)
    return x0 * np.cos(omega * t)

def euler_method(f, t0, state0, dt, N):
    t_values = np.zeros(N + 1)
    states = np.zeros((N + 1, 2))
    t_values[0] = t0
    states[0] = state0
    t = t0
    state = state0.copy()
    for i in range(N):
        state = state + dt * f(t, state)
        t = t + dt
        t_values[i + 1] = t
        states[i + 1] = state
    return t_values, states

def rk4_step(f, t, state, dt):
    k1 = f(t, state)
    k2 = f(t + dt / 2, state + dt * k1 / 2)
    k3 = f(t + dt / 2, state + dt * k2 / 2)
    k4 = f(t + dt, state + dt * k3)
    return state + dt * (k1 + 2*k2 + 2*k3 + k4) / 6

def rk4_method(f, t0, state0, dt, N):
    t_values = np.zeros(N + 1)
    states = np.zeros((N + 1, 2))
    t_values[0] = t0
    states[0] = state0
    t = t0
    state = state0.copy()
    for i in range(N):
        state = rk4_step(f, t, state, dt)
        t = t + dt
        t_values[i + 1] = t
        states[i + 1] = state
    return t_values, states

def leapfrog_method(t0, x0, v0, dt, N):
    t_values = np.zeros(N + 1)
    x_values = np.zeros(N + 1)
    v_values = np.zeros(N + 1)
    t_values[0] = t0
    x_values[0] = x0
    v_values[0] = v0
    x = x0
    v = v0
    t = t0
    a = -(k / m) * x
    for i in range(N):
        v_half = v + (dt / 2) * a
        x = x + dt * v_half
        a_new = -(k / m) * x
        v = v_half + (dt / 2) * a_new
        t = t + dt
        t_values[i + 1] = t
        x_values[i + 1] = x
        v_values[i + 1] = v
        a = a_new
    return t_values, x_values, v_values

state0 = np.array([x0, v0])

t_euler, state_euler = euler_method(f, t0, state0, dt, N)
x_euler = state_euler[:, 0]

t_rk4, state_rk4 = rk4_method(f, t0, state0, dt, N)
x_rk4 = state_rk4[:, 0]

t_leapfrog, x_leapfrog, v_leapfrog = leapfrog_method(t0, x0, v0, dt, N)

x_exact = exact_solution(t_euler)

error_euler = np.abs(x_euler - x_exact)
error_rk4 = np.abs(x_rk4 - x_exact)
error_leapfrog = np.abs(x_leapfrog - x_exact)

plt.figure(figsize=(10, 6))
plt.plot(t_euler, error_euler, label="Euler")
plt.plot(t_rk4, error_rk4, label="RK4")
plt.plot(t_leapfrog, error_leapfrog, label="Leapfrog")
plt.xlabel("Time (s)")
plt.ylabel("Absolute error (m)")
plt.title("Mass-Spring System: Error vs. Time")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.yscale("log")
plt.show()