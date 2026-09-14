import numpy as np
import matplotlib.pyplot as plt

def fwd_deriv(f, x, h):
    '''
    f : callable of a single variable
    x : float
    h: float
    returns the approx derivative of f at x
    '''
    return ( f(x+h) - f(x) )/ h

def PO2_func(depth, total_depth, Thiele_mod, PO2_atm):
    if PO2_atm > 0:
        PO2 = PO2_atm*np.cosh((total_depth - depth)/Thiele_mod)
        PO2 /= np.cosh(total_depth/Thiele_mod)
    else:
        PO2 = np.nan
    return PO2
def my_func(x):
    return PO2_func(x, 30, 10, 0.21)
h_arr = [0.1, 0.01, 0.001, 0.0001, 0.00001]
err_arr = []
for h in h_arr:
    err = fwd_deriv(my_func, 0, h) - (-0.21*0.1*np.tanh(0.1*30))
    err = abs(err)
    err_arr.append(err)
plt.plot(err_arr, h_arr, 'o', color = "magenta")
plt.loglog()
plt.show()