import math as m
def NewtonRaphson(f, fprime, h0, tol=1e-3, mahiters = 1000):
    hguess =h0
    Converged = False
    for i in range(mahiters):
        hnew = hguess - f(hguess)/fprime(hguess)
        err = abs(hnew - hguess)
        hguess = hnew
        if err < tol:
            Converged = True
            break
    return hguess, Converged

if __name__=='__main__':
    def my_func(h):
        V, R = 30, 3
        return (9*m.pi*h**2) - (m.pi*h**3) - 90
    def my_func_prime(h):
        R = 3
        return 18*m.pi*h - 3*m.pi*h**2
    result = NewtonRaphson(my_func, my_func_prime, 1.875)
    print(result)