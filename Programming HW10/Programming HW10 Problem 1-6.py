import math as m
def g(h):
    return m.sqrt(90/(m.pi*(9-h)))

def fixedpoint(g, h_guess, tol):
    err = tol + 1
    while err > tol:
        h_new = g(h_guess)
        err = abs(h_new - h_guess)
        h_guess = h_new
    return h_guess

h_guess = 3
tol = 0.001

h = fixedpoint(g, h_guess, tol)
print(h, 'depth in m')