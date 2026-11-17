import math as m
def False_Position(f, bracket, abs_tol):
    hL, hU = bracket
    err = abs_tol + 1
    fhL = f(hL)
    fhU = f(hU)
    while err > abs_tol:
        hR = hU - ((fhU*(hL - hU))/(fhL - fhU))
        fhR = f(hR)
        if fhR == 0:
            break
        elif fhR*fhL < 0:
            hU = hR
            fhU = fhR
        else:
            hL = hR
            fhL = fhR
        err = (abs(hL - hU))/2
        print(hL, hU)
    return hU - ((fhU*(hL - hU))/(fhL - fhU))
if __name__=='__main__':
    def my_func(h):
        V, R = 30, 3
        return  (9*m.pi*h**2) - (m.pi*h**3) - 90
    result = False_Position(my_func, (1.5,2.25), 0.001)
    print(result, 'depth in m')