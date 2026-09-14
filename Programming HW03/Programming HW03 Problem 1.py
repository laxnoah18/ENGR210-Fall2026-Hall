import numpy as np
import matplotlib.pyplot as plt

x = np.linspace(0, 48, 100)
y1 = np.sin(np.linspace(0, 48, 100)*(2*np.pi/24))*2
y2 = np.cos(np.linspace(0, 48, 100)*(2*np.pi/24))*2

plt.plot(x, y1, linestyle = '--', color = 'green', label="sin(x)")
plt.plot(x, y2, linestyle = '-.', color = 'red', label="cos(x)")
plt.xlabel("t")
plt.ylabel("f(t)")
plt.legend(loc="upper right")
plt.show()