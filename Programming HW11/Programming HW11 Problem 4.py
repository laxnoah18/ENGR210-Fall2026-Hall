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

K = 2.0
B0 = 100.0

# Calculate equilibrium concentration of C
def equilibrium_C(A0):
    # Physical limits for C
    C_low = 0
    C_high = min(A0 / 2, B0)
    
    # Equilibrium equation
    def equation(C):
        return C - K * (A0 - 2*C)**2 * (B0 - C)

    # Bisection method
    for i in range(200):
        C_mid = (C_low + C_high) / 2
        if equation(C_mid) > 0:
            C_high = C_mid
        else:
            C_low = C_mid
    return (C_low + C_high) / 2

# Calculate profit
def profit(A0):

    C = equilibrium_C(A0)

    return 10*C - A0

# Golden-section optimization
A_opt, P_opt, iterations, error = golden_section_max(
    profit, 0, 500, tol=1e-6
)

# Calculate equilibrium C at optimum
C_opt = equilibrium_C(A_opt)

print("Chemical Reaction Optimization")
print("Optimal initial concentration of A:", A_opt, "M")
print("Equilibrium concentration of C:", C_opt, "M")
print("Maximum profit:", P_opt)
print("Iterations:", iterations)
print("Estimated error:", error)