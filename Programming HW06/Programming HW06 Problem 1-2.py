import numpy as np
import matplotlib.pyplot as plt

v0 = 22
m = 0.014
C_d = 1.3e-3
g = 9.81
dt = 0.001

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

def simulate(angle):
    theta = np.deg2rad(angle)
    state = np.array([0.0, 0.0, v0*np.cos(theta), v0*np.sin(theta)])
    t = 0.0
    
    while True:
        old_state = state.copy()
        old_t = t
        state = rk4_step(wiffleball_deriv, t, state, dt)
        t += dt
        if state[1] < 0:
            fraction = (old_state[1] /(old_state[1] - state[1]))
            landing_state = (old_state + fraction * (state - old_state))
            return landing_state[0]

angles = np.linspace(0, 90, 91)
ranges = []
for angle in angles:
    range_value = simulate(angle)
    ranges.append(range_value)
ranges = np.array(ranges)

best_index = np.argmax(ranges)
best_angle = angles[best_index]
max_range = ranges[best_index]

print(f"Maximum range = {max_range:.3f} meters")
print(f"Corresponding launch angle = {best_angle:.1f} degrees")

plt.figure(figsize=(8, 5))
plt.plot(angles, ranges, "o-")
plt.xlabel("Launch angle (degrees)")
plt.ylabel("Range (m)")
plt.title("Wiffle Ball Range vs. Launch Angle")
plt.grid(True)
plt.tight_layout()
plt.show()