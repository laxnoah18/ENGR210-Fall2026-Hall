"""Classical fourth-order Runge-Kutta with convergence and solution plots.

Solve scalar y' = f(x, y) initial-value problems on a uniform grid.
Dependencies: NumPy and Matplotlib only.
"""

import numpy as np
import matplotlib.pyplot as plt


def rk4_step(f, x, y, h):
    """Advance one classical RK4 step."""
    k1 = f(x, y)
    k2 = f(x + h / 2, y + h * k1 / 2)
    k3 = f(x + h / 2, y + h * k2 / 2)
    k4 = f(x + h, y + h * k3)
    return y + h * (k1 + 2 * k2 + 2 * k3 + k4) / 6


def integrate(f, y0, x0, x_end, steps):
    """Return NumPy arrays of grid points and solution values."""
    if steps < 1:
        raise ValueError("steps must be a positive integer")
    if x_end == x0:
        raise ValueError("x_end must differ from x0")

    h = (x_end - x0) / steps
    x = np.linspace(x0, x_end, steps + 1)
    y = np.empty(steps + 1, dtype=float)
    y[0] = y0
    for n in range(steps):
        y[n + 1] = rk4_step(f, x[n], y[n], h)
    return x, y


def stability_polynomial(z):
    """RK4 amplification factor R(z) for y' = lambda*y, z=lambda*h."""
    return 1 + z + z**2 / 2 + z**3 / 6 + z**4 / 24


def convergence_demo():
    """Measure endpoint error for y'=y, y(0)=1 on [0,1]."""
    step_counts = np.array([10, 20, 40, 80, 160])
    errors = np.empty(len(step_counts))
    for i, steps in enumerate(step_counts):
        _, y = integrate(lambda x, y: y, 1.0, 0.0, 1.0, int(steps))
        errors[i] = abs(y[-1] - np.e)
    observed_orders = np.full(len(errors), np.nan)
    observed_orders[1:] = np.log(errors[:-1] / errors[1:]) / np.log(2)
    return step_counts, errors, observed_orders


if __name__ == "__main__":
    step_counts, errors, orders = convergence_demo()
    for n, error, order in zip(step_counts, errors, orders):
        print(f"steps={n:3d}  endpoint_error={error:.6e}  observed_order={order:.4f}")

    # Verification plot: endpoint error versus step size.
    h_values = 1.0 / step_counts
    plt.figure(figsize=(7, 5))
    plt.loglog(h_values, errors, "o-", label="RK4 endpoint error")
    plt.loglog(
        h_values,
        errors[0] * (h_values / h_values[0]) ** 4,
        "--",
        label="slope 4 reference",
    )
    plt.xlabel("Step size h")
    plt.ylabel("Absolute endpoint error at x=1")
    plt.title("RK4 fourth-order convergence for y' = y")
    plt.grid(True, which="both", alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.show()

    # Example requested: y'=-y, y(0)=1, whose exact solution is exp(-x).
    x, y = integrate(lambda x, y: -y, 1.0, 0.0, 5.0, 50)
    plt.figure(figsize=(7, 5))
    plt.plot(x, y, "o", markevery=4, label="RK4 numerical solution")
    plt.plot(x, np.exp(-x), "-", label="Exact solution exp(-x)")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.title("RK4 solution of y' = -y, y(0) = 1")
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.show()
