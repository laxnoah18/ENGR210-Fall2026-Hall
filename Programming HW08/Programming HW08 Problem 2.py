import numpy as np
import matplotlib.pyplot as plt

# Coefficients and boundary conditions
hpr = 0.01       # 1/m^2
L = 10           # m
Tleft = 200      # deg(C)
Tambient = 20    # deg(C)

# Grid spacings
dx_values = [1.0, 0.1, 0.01]

# Analytical solution parameter
m = np.sqrt(hpr)

# Analytical solution for an insulated end
def analytical_solution(x):
    return (Tambient
            + (Tleft - Tambient)
            * np.cosh(m * (L - x)) / np.cosh(m * L))

# Store errors
errors = []

# Plot temperature versus distance
plt.figure(figsize=(8, 5))

# Analytical solution
x_exact = np.linspace(0, L, 1001)
T_exact = analytical_solution(x_exact)

plt.plot(x_exact, T_exact,
         label="Analytical solution", linewidth=2)

# Loop over the grid spacings
for dx in dx_values:
    # Number of grid points
    n = int(L / dx) + 1
    # Create A and b
    A = np.zeros((n, n), dtype=float)
    b = np.zeros(n, dtype=float)
    # Fill in A for the interior points
    np.fill_diagonal(A[1:-1, 1:-1], -(2 + hpr * dx**2))
    np.fill_diagonal(A[1:-1, 0:-2], 1)
    np.fill_diagonal(A[1:-1, 2:], 1)
    # Left boundary: fixed temperature
    A[0, 0] = 1
    b[0] = Tleft
    # Right boundary: insulated end
    A[-1, -3] = 1
    A[-1, -2] = -4
    A[-1, -1] = 3
    b[-1] = 0
    # Interior equations
    b[1:-1] = -hpr * dx**2 * Tambient
    # Solve for temperature
    T = np.linalg.solve(A, b)
    # Position array
    x_arr = np.linspace(0, L, n)
    # Analytical temperatures at the same points
    T_analytical = analytical_solution(x_arr)
    # Maximum absolute error
    error = np.max(np.abs(T - T_analytical))
    errors.append(error)

    # Print results
    print(f"dx = {dx} m, n = {n}, "
          f"T at insulated end = {T[-1]:.6f} deg(C), "
          f"maximum error = {error:.6e} deg(C)")

    # Plot numerical solution
    plt.plot(x_arr, T, "o-", markersize=3,
             label=f"Finite difference, dx = {dx} m")


# Finish temperature plot
plt.xlabel("Distance x (m)")
plt.ylabel("Temperature T (deg C)")
plt.title("Temperature Distribution: Insulated Rod End")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()

# Plot error versus grid spacing
plt.figure(figsize=(6, 5))

plt.loglog(dx_values, errors, "o-",
           label="Maximum absolute error")

# Reference line for second-order convergence
reference_errors = (
    errors[0] * (np.array(dx_values) / dx_values[0])**2
)

plt.loglog(dx_values, reference_errors, "--",
           label="Second-order reference")

plt.xlabel("Grid spacing dx (m)")
plt.ylabel("Maximum absolute error (deg C)")
plt.title("Convergence: Insulated Rod End")
plt.grid(True, which="both")
plt.legend()
plt.tight_layout()
plt.show()

# Estimate convergence order
slope = np.polyfit(np.log(dx_values), np.log(errors), 1)[0]

print(f"Estimated convergence order = {slope:.4f}")