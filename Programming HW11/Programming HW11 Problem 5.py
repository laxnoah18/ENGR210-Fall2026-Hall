import math
def golden_section_max(f, xL, xU, tol=1e-6):
    # Golden ratio constant
    R = (math.sqrt(5) - 1) / 2

    # Initial interval width
    h = xU - xL

    # Calculate the two interior points
    x1 = xU - R * h
    x2 = xL + R * h

    # Evaluate the function
    fx1 = f(x1)
    fx2 = f(x2)

    # Initialize iteration counter
    iterations = 0

    # Continue until the estimated error is within tolerance
    while h / 2 > tol:
        if fx2 > fx1:
            # Maximum is in the right portion
            xL = x1
            x1 = x2
            fx1 = fx2

            h = xU - xL
            x2 = xL + R * h
            fx2 = f(x2)
        else:
            # Maximum is in the left portion
            xU = x2
            x2 = x1
            fx2 = fx1
            h = xU - xL
            x1 = xU - R * h
            fx1 = f(x1)
        iterations += 1
    # Calculate final optimum
    x_opt = (xL + xU) / 2
    f_opt = f(x_opt)
    error = h / 2
    return x_opt, f_opt, iterations, error

# Given constants
L = 400
E = 52000
I = 32000
w = 4

# Deflection equation
def beam_deflection(x):
    y = (w / (48 * E * I)) * (
        2*x**4 - 3*L*x**3 + L**3*x
    )
    return y

# Golden-section optimization
x_opt, y_opt, iterations, error = golden_section_max(
    beam_deflection, 0, L, tol=1e-6
)

print("Beam Deflection Optimization")
print("Location of maximum deflection:", x_opt, "cm")
print("Maximum deflection:", y_opt, "cm")
print("Iterations:", iterations)
print("Estimated error:", error)