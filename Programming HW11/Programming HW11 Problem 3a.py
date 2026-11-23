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

def parabola(x):
    return -3*x**2 + 9*x + 3
# Search interval
xL = -1
xU = 5
# Golden-section optimization
x_opt, f_opt, iterations, error = golden_section_max(
    parabola, xL, xU
)
# Analytical solution
x_exact = -9 / (2 * -3)

print("Parabola Verification")
print("Calculated optimum:", x_opt)
print("Exact optimum:", x_exact)
print("Maximum function value:", f_opt)
print("Iterations:", iterations)
print("Estimated error:", error)