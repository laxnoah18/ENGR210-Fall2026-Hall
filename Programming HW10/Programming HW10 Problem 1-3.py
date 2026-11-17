import math as m
def Bisection(f, bracket, abs_tol):
    hL, hU = bracket
    iters = m.ceil(m.log((hU - hL)/abs_tol)/m.log(2) - 1)
    fhL = f(hL)
    fhU = f(hU)
    for i in range(iters):
        hR = (hL + hU)/2
        fhR = f(hR)
        if fhR == 0:
            break
        elif fhR*fhL < 0:
            hU = hR
            fhU = fhR
        else:
            hL = hR
            fhL = fhR
        print(hL, hU)
    return (hL + hU)/2
if __name__=='__main__':
    def my_func(h):
        V, R = 30, 3
        return  (9*m.pi*h**2) - (m.pi*h**3) - 90
    result = Bisection(my_func, (1.5,2.25), 0.001)
    print(result, 'depth in m')