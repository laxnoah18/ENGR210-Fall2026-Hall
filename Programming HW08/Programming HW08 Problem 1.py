import numpy as np
import matplotlib.pyplot as plt

# Coefficients and boundary conditions
hpr = 0.01       # 1/m^2
L = 10           # m
Tleft = 200      # deg(C)
Tright = 40      # deg(C)
Tambient = 20    # deg(C)

# Grid spacings to test
dx_values = [1.0, 0.1, 0.01]

# Analytical solution parameter
m = np.sqrt(hpr)

# Analytical solution function
def analytical_solution(x):
    return (Tambient
            + (Tleft - Tambient) * np.sinh(m * (L - x)) / np.sinh(m * L)
            + (Tright - Tambient) * np.sinh(m * x) / np.sinh(m * L))

# Arrays to store results
errors = []

# Plot the numerical and analytical solutions
plt.figure(figsize=(8, 5))

# Smooth curve for the analytical solution
x_exact = np.linspace(0, L, 1001)
T_exact = analytical_solution(x_exact)

plt.plot(x_exact, T_exact, label="Analytical solution", linewidth=2)

# Loop over the three grid spacings
for dx in dx_values:
    # Number of grid points
    n = int(L / dx) + 1
    # Create matrix A and vector b
    A = np.zeros((n, n), dtype=float)
    b = np.zeros(n, dtype=float)
    # Fill in A for the interior points
    np.fill_diagonal(A[1:-1, 1:-1], -(2 + hpr * dx**2))
    np.fill_diagonal(A[1:-1, 0:-2], 1)
    np.fill_diagonal(A[1:-1, 2:], 1)
    # Apply fixed-temperature boundary conditions
    A[0, 0] = 1
    A[-1, -1] = 1
    
    b[1:-1] = -hpr * dx**2 * Tambient
    b[0] = Tleft
    b[-1] = Tright

    # Solve for the temperature distribution
    T = np.linalg.solve(A, b)

    # Create the position array
    x_arr = np.linspace(0, L, n)

    # Calculate the analytical solution at the same points
    T_analytical = analytical_solution(x_arr)

    # Calculate maximum absolute error
    error = np.max(np.abs(T - T_analytical))
    errors.append(error)

    # Print the error for this grid spacing
    print(f"dx = {dx:.2f} m, n = {n}, maximum error = {error:.6e} deg(C)")

    # Plot numerical solution
    plt.plot(x_arr, T, "o-", markersize=3,
             label=f"Finite difference, dx = {dx} m")

# Finish temperature plot
plt.xlabel("Distance x (m)")
plt.ylabel("Temperature T (deg C)")
plt.title("Temperature Distribution Along the Rod")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()

# Plot error versus grid spacing
plt.figure(figsize=(6, 5))
plt.loglog(dx_values, errors, "o-", label="Maximum absolute error")

# Reference line for second-order convergence
reference_errors = errors[0] * (np.array(dx_values) / dx_values[0])**2
plt.loglog(dx_values, reference_errors, "--",
           label="Second-order reference")

plt.xlabel("Grid spacing dx (m)")
plt.ylabel("Maximum absolute error (deg C)")
plt.title("Convergence of the Finite Difference Solution")
plt.grid(True, which="both")
plt.legend()
plt.tight_layout()
plt.show()

# Estimate the convergence order
slope = np.polyfit(np.log(dx_values), np.log(errors), 1)[0]
print(f"Estimated convergence order = {slope:.4f}")