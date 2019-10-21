import numpy as np
import matplotlib.pyplot as plt

def func(x):
    return np.exp(-2*abs(x))

r = np.linspace(-5, 5, 1000)
plt.plot(r, func(r))
plt.title("Single Particle Wavefunction", fontsize = 17)
plt.xlabel('x', fontsize = 17)
plt.ylabel('Wavefunction', fontsize = 17)
plt.show()
