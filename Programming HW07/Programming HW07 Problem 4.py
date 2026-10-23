import numpy as np

A = [[2.0, -2.0, 3.0, -3.0],
    [1.0, 1.0, 1.0, 0.0],
    [1.0, -1.0, 2.0, -1.0],
    [1.0, -1.0, 4.0, 3.0]]

B = [-20.0, -2.0, -8.0, 4.0]

soln = np.linalg.solve(A,B)
print(soln)
print(np.dot(A,soln))