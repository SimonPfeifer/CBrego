import numpy as np
import matplotlib.pyplot as plt

density = np.genfromtxt('./output_density.txt', delimiter=',')

fig = plt.figure()
ax = fig.add_subplot(111)
nx = len(density[0])
x = np.arange(nx) / nx + (1. / nx / 2)
ymax = np.nanmax(density)*1.1
for rho in density:
    ax.clear()
    ax.plot(x, rho)
    ax.set_xlim([0,1])
    ax.set_ylim([0,ymax])
    plt.pause(0.01)
plt.show()