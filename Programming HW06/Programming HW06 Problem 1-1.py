import numpy as np
import matplotlib.pyplot as plt

v0 = 22
m = 0.014
C_d = 1.3e-3
g = 9.81
dt = 0.001
angle = np.deg2rad(45)

def rk4_step(f, t, state, dt):
    """Advance one classical RK4 step."""
    k1 = f(t, state)
    k2 = f(t + dt / 2, state + dt * k1 / 2)
    k3 = f(t + dt / 2, state + dt * k2 / 2)
    k4 = f(t + dt, state + dt * k3)
    return state + dt * (k1 + 2 * k2 + 2 * k3 + k4) / 6


def wiffleball_deriv(t,state):
    x, y, vx, vy = state
    speed = np.sqrt(vx**2 + vy**2)
    ax = -(C_d/m)*speed*vx
    ay = -g - (C_d/m)*speed*vy
    return np.array([vx, vy, ax, ay])


state = np.array([0.0, 0.0, v0*np.cos(angle), v0*np.sin(angle)])
t = 0.0
times = [t]
states = [state.copy()]

while True:
    old_state = state.copy()
    old_t = t
    state = rk4_step(wiffleball_deriv, t, state, dt)
    t += dt
    times.append(t)
    states.append(state.copy())
    
    if state[1] < 0:
        fraction = (old_state[1] / (old_state[1] - state[1]))
        landing_state = (old_state + fraction * (state - old_state))
        landing_time = old_t + fraction * dt
        times[-1] = landing_time
        states[-1] = landing_state
        
        break

times = np.array(times)
states = np.array(states)
x = states[:,0]
y = states[:,1]

print(f"Launch angle = {np.rad2deg(angle):.1f} degrees")
print(f"Flight time = {times[-1]:.3f} seconds")
print(f"Range = {x[-1]:.3f} meters")
print(f"Maximum height = {np.max(y):.3f} meters")

plt.figure(figsize=(8, 5))
plt.plot(x, y)
plt.xlabel("Horizontal distance (m)")
plt.ylabel("Height (m)")
plt.title("Wiffle Ball Trajectory")
plt.grid(True)
plt.tight_layout()
plt.show()