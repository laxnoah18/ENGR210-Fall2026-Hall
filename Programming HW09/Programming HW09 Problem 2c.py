import numpy as np
import matplotlib.pyplot as plt

# Grid dimensions
Nx, Ny = 101, 101

x = np.linspace(0, 1, Nx)
y = np.linspace(0, 1, Ny)

dx = x[1] - x[0]
dy = y[1] - y[0]

# Initialize temperature array
T_arr = np.ones((Nx, Ny), dtype=float) * 50

# Apply boundary conditions
T_arr[1:-1, 0] = 0       # Bottom
T_arr[1:-1, -1] = 100    # Top
T_arr[0, 1:-1] = 75      # Left
T_arr[-1, 1:-1] = 50     # Right

# Corner temperatures
T_arr[0, 0] = 37.5
T_arr[-1, 0] = 25
T_arr[-1, -1] = 75
T_arr[0, -1] = 87.5

# Gauss-Seidel settings
tol = 1e-4
omega = 1.9
max_iterations = 100000

# Gauss-Seidel with over-relaxation
for iters in range(max_iterations):
    max_err = 0
    for i in range(1, Nx-1):
        for j in range(1, Ny-1):
            T_new = (
                T_arr[i-1, j] + T_arr[i+1, j]
                + T_arr[i, j-1] + T_arr[i, j+1]
            ) / 4
            T_new = omega*T_new + (1-omega)*T_arr[i,j]
            max_err = max(max_err, abs(T_new-T_arr[i,j]))
            T_arr[i,j] = T_new
    if max_err < tol:
        print("Converged")
        break

# Analytical solution using Fourier series

X, Y = np.meshgrid(x, y)

def sinh_ratio(a, z):
    return (
        np.exp(a*(z-1))
        * (-np.expm1(-2*a*z))
        / (-np.expm1(-2*a))
    )
T_exact = np.zeros((Ny, Nx))

# Sum odd Fourier terms
for n in range(1, 1000, 2):
    a = n*np.pi
    coefficient = 4/(np.pi*n)
    T_exact += coefficient * (
        100*np.sin(n*np.pi*X)*sinh_ratio(a, Y)
        + 75*np.sin(n*np.pi*Y)*sinh_ratio(a, 1-X)
        + 50*np.sin(n*np.pi*Y)*sinh_ratio(a, X)
    )
# Compare numerical and analytical solutions
T_numerical = T_arr.T

# Compare interior points only
error = np.abs(
    T_numerical[1:-1, 1:-1]
    - T_exact[1:-1, 1:-1]
)
max_error = np.max(error)
rms_error = np.sqrt(np.mean(error**2))

# Print results
print("Iterations:", iters + 1)
print("Maximum numerical error:", max_error, "C")
print("RMS numerical error:", rms_error, "C")

# Numerical solution plot
plt.figure(figsize=(9, 5))
plt.contourf(X, Y, T_numerical, levels=30, cmap="hot")
plt.colorbar(label="Temperature (C)")
plt.xlabel("Horizontal position (m)")
plt.ylabel("Vertical position (m)")
plt.title("Gauss-Seidel Numerical Solution")
plt.tight_layout()
plt.show()
# Analytical solution plot
plt.figure(figsize=(9, 5))
plt.contourf(X, Y, T_exact, levels=30, cmap="hot")
plt.colorbar(label="Temperature (C)")
plt.xlabel("Horizontal position (m)")
plt.ylabel("Vertical position (m)")
plt.title("Analytical Solution")
plt.tight_layout()
plt.show()
# Numerical error plot
plt.figure(figsize=(9, 5))
plt.contourf(X[1:-1, 1:-1],
             Y[1:-1, 1:-1],
             error,
             levels=30,
             cmap="viridis")
plt.colorbar(label="Absolute Error (C)")
plt.xlabel("Horizontal position (m)")
plt.ylabel("Vertical position (m)")
plt.title("Numerical Error")
plt.tight_layout()
plt.show()