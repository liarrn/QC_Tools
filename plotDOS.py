import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt('TDOS')
data[:, 0] = data[:, 0] + 5.20172342
x = data[:, 0]
y = data[:, 1]
y = y[np.all(np.array([x>-4, x<4]), axis = 0)]
x = x[np.all(np.array([x>-4, x<4]), axis = 0)]
plt.figure(num=None, figsize=(12, 10), dpi=1200, facecolor='w', edgecolor='k')
plt.xlabel('energy (eV)', fontsize=20)
plt.ylabel('DOS', fontsize=20)
plt.grid()
plt.plot(x, y, color='grey')
plt.fill_between(x, 0, y, where=x<0, color='grey')
plt.show()