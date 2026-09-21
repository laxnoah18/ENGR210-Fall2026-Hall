import numpy as np
import matplotlib.pyplot as plt

def trap_int_data(strain, strain_i, stress, stress_i):
    integral = ( ((stress_i + stress)/2) * (strain_i - strain) )
    return integral
    
data = np.loadtxt(r'C:/Users/twinl/Downloads/Copy of ENGR210 Fa2026 HW04 Prob2 Data.csv',
                  delimiter = ',',
                  skiprows = 1)

data = np.array(data)
strain = data[:,1]
stress = data[:,0]
print(strain)
print(stress)
n = len(data)

total = 0
for i in range(n - 1):
    integral_sum = ( ((stress[i+1] + stress[i])/2) * (strain[i+1] - strain[i]) )
    total += integral_sum
print(total)