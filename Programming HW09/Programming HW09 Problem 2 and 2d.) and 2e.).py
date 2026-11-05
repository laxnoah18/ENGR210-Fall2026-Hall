import numpy as np
import matplotlib.pyplot as plt

# Physical properties
k = 25.0
ho = 1000.0
hi = 200.0

T_hot = 1700.0
T_air = 400.0

# Geometry (meters)
Lx = 0.010
Ly = 0.006

channel_x1 = 0.002
channel_x2 = 0.008

channel_y1 = 0.002
channel_y2 = 0.004

# Grid spacing
dx = 0.00025
dy = 0.00025

nx = round(Lx / dx)
ny = round(Ly / dy)

x_edges = np.linspace(0, Lx, nx + 1)
y_edges = np.linspace(0, Ly, ny + 1)

x = (x_edges[:-1] + x_edges[1:]) / 2
y = (y_edges[:-1] + y_edges[1:]) / 2

# Identify solid cells and channel cells
X, Y = np.meshgrid(x, y)

channel = (
    (X >= channel_x1) & (X < channel_x2) &
    (Y >= channel_y1) & (Y < channel_y2)
)

solid = ~channel

# Initial temperature
T = np.full((ny, nx), 1000.0)

# Gauss-Seidel settings
tolerance = 1e-8
max_iterations = 100000

# Unit length into the page (1 meter)
depth = 1.0

# Gauss-Seidel iteration
for iteration in range(max_iterations):
    max_error = 0.0
    for j in range(ny):
        for i in range(nx):
            if not solid[j, i]:
                continue
            numerator = 0.0
            denominator = 0.0
            # Examine four neighboring faces
            for dj, di in [(0, -1), (0, 1),
                           (-1, 0), (1, 0)]:
                jj = j + dj
                ii = i + di
                # Face area and distance between cell centers
                if di != 0:
                    A = dy * depth
                    distance = dx
                else:
                    A = dx * depth
                    distance = dy
                # Neighbor is inside computational domain
                if 0 <= jj < ny and 0 <= ii < nx:
                    if solid[jj, ii]:
                        # Conduction to another solid cell
                        G = k * A / distance
                        numerator += G * T[jj, ii]
                        denominator += G
                    else:
                        # Convection to cooling air
                        G = A / (distance / (2*k) + 1/hi)
                        numerator += G * T_air
                        denominator += G
                # Outside the computational domain
                else:
                    # Left and right boundaries are symmetry
                    if ii < 0 or ii >= nx:
                        continue
                    # Top and bottom boundaries are hot gas
                    G = A / (distance / (2*k) + 1/ho)
                    numerator += G * T_hot
                    denominator += G
            # Update temperature
            T_old = T[j, i]
            T[j, i] = numerator / denominator
            error = abs(T[j, i] - T_old)
            if error > max_error:
                max_error = error
    # Check convergence
    if max_error < tolerance:
        break

print("Iterations:", iteration + 1)
# Find maximum temperature
T_max = np.max(T[solid])

print("Maximum temperature:", T_max, "K")

# Find location of maximum temperature
j_max, i_max = np.where(solid & (T == T_max))

x_max = x[i_max[0]] * 1000
y_max = y[j_max[0]] * 1000

print("Maximum temperature location:")
print("x =", x_max, "mm")
print("y =", y_max, "mm")

# Plot temperature distribution
T_plot = np.ma.masked_where(~solid, T)

plt.figure(figsize=(10, 5))
plt.pcolormesh(
    x_edges * 1000,
    y_edges * 1000,
    T_plot,
    shading="flat",
    cmap="hot"
)
plt.colorbar(label="Temperature (K)")
plt.xlabel("Horizontal position (mm)")
plt.ylabel("Vertical position (mm)")
plt.title("Turbine Blade Temperature Distribution")
plt.axis("equal")
plt.tight_layout()
plt.show()

# Calculate heat transfer to cooling air
Q_air = 0.0
Q_hot = 0.0

for j in range(ny):
    for i in range(nx):
        if not solid[j, i]:
            continue
        # Four neighboring faces
        for dj, di in [(0, -1), (0, 1),
                       (-1, 0), (1, 0)]:
            jj = j + dj
            ii = i + di
            if di != 0:
                A = dy * depth
                distance = dx
            else:
                A = dx * depth
                distance = dy
            # Internal cooling channel
            if (0 <= jj < ny and 0 <= ii < nx
                    and not solid[jj, ii]):
                G = A / (distance / (2*k) + 1/hi)
                Q_air += G * (T[j, i] - T_air)
            # External hot gas boundaries
            elif (ii >= 0 and ii < nx and
                  (jj < 0 or jj >= ny)):
                G = A / (distance / (2*k) + 1/ho)
                Q_hot += G * (T_hot - T[j, i])

print("Heat transfer to cooling air:", Q_air, "W/m")
print("Heat transfer from hot gases:", Q_hot, "W/m")
print("Energy balance:", Q_hot - Q_air, "W/m")

# Check temperature symmetry
symmetry_x = np.max(np.abs(T - np.fliplr(T)))
symmetry_y = np.max(np.abs(T - np.flipud(T)))

print("Symmetry check:")
print("Horizontal symmetry error:", symmetry_x, "K")
print("Vertical symmetry error:", symmetry_y, "K")

# Check convection boundary conditions
max_convection_error = 0.0
for j in range(ny):
    for i in range(nx):
        if not solid[j, i]:
            continue
        # Examine four neighboring faces
        for dj, di in [(0, -1), (0, 1),
                       (-1, 0), (1, 0)]:
            jj = j + dj
            ii = i + di
            if di != 0:
                A = dy * depth
                distance = dx
            else:
                A = dx * depth
                distance = dy
            # Internal cooling channel
            if (0 <= jj < ny and 0 <= ii < nx
                    and not solid[jj, ii]):
                h = hi
                T_inf = T_air
            # External hot gas surfaces
            elif (0 <= ii < nx and
                  (jj < 0 or jj >= ny)):
                h = ho
                T_inf = T_hot
            else:
                continue
            # Calculate heat transfer
            G = A / (distance / (2*k) + 1/h)
            q = G * abs(T[j, i] - T_inf)
            # Calculate surface temperature
            T_surface = (
                T_inf
                + np.sign(T[j, i] - T_inf) * q / (h*A)
            )
            # Heat transfer by convection
            q_convection = h*A*abs(T_surface - T_inf)
            # Heat transfer by conduction
            q_conduction = (
                k*A*abs(T[j, i] - T_surface)
                / (distance/2)
            )
            error = abs(q_convection - q_conduction)
            if error > max_convection_error:
                max_convection_error = error
print("Maximum convection boundary error:",
      max_convection_error, "W/m")

# Check overall energy balance
energy_error = abs(Q_hot - Q_air)
relative_error = (
    energy_error / Q_hot
) * 100
print("Energy balance check:")
print("Heat entering:", Q_hot, "W/m")
print("Heat leaving:", Q_air, "W/m")
print("Energy balance error:", energy_error, "W/m")
print("Relative energy error:", relative_error, "%")