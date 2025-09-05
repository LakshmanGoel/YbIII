import numpy as np
import matplotlib.pyplot as plt

# Fundamental constants
planck = 6.626*10**(-34)
boltzmann = 1.38*10**(-23)

# Belt specifications
L = 0.30 # (in m)

# Beam specifications
wavelength = 515*10**(-9) # (in m)
rayleigh = L/2 # (in m)
waist = np.sqrt(wavelength*rayleigh/(np.pi)) # (in m)

# Polarizability
alpha = 280 * (2* 377 * 2.48832 * 10**(-8)) # (in Hz*m^2/W)

# Power calculation from trap depth
depth_K = 100 * 10**(-6) # (in K)
depth_J = depth_K * boltzmann # (in J)
depth_Hz = depth_J / planck # (in Hz)
I = depth_Hz/alpha # (in W m^-2)
P = (np.pi * waist**(2) * I)/8 # (in W)

print(depth_Hz)

def trap_potential(x, w_0):
    """Calculates the trap potential at different positions for a single trap."""
    coefficient = -depth_J
    exponent = -2*(x/w_0)**2
    return coefficient * np.exp(exponent)

# Parameters of the Gaussian distribution
# mean = 0
# std_dev = 1
w_0 = 100*10**(-6)
m = 2.8733965 * 10**(-25)
g = 9.8

print(m*g)

# Generate x values for the plot
# Typically, you'd cover a range like mean - 3*std_dev to mean + 3*std_dev
x = np.linspace(-200*10**(-6), 200*10**(-6), 500)

# Calculate the corresponding y values (PDF values)
y = (trap_potential(x, w_0) + m*g*x)/planck

# Plotting
plt.plot(x, y, label=f'Beam waist={w_0}')
plt.plot(x, trap_potential(x, w_0)/planck, label=f'Beam waist={w_0}')
plt.xlabel('Position')
plt.ylabel('Potential energy')
plt.title('Trap depth + Gravitational potential')
plt.grid(True)
plt.legend()
plt.show()