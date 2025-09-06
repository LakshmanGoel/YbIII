import numpy as np
import matplotlib.pyplot as plt

# Fundamental constants
planck = 6.626*10**(-34) # (in SI units)
boltzmann = 1.38*10**(-23) # (in SI units)

# Polarizability
alpha = 280 * (2* 377 * 2.48832 * 10**(-8)) # (in Hz*m^2/W)

def calculate_waist(L, wavelength):
    # Beam waist calculation from a given conveyor belt length and laser wavelength in m
    # L is conveyor belt length in m
    # wavelength is the laser wavelength in m
    rayleigh = L/2 # (in m)
    waist = np.sqrt(wavelength*rayleigh/(np.pi)) # (in m)
    return waist

def calculate_power(depth_K, L, wavelength):
    # Power calculation from a given trap depth in K and beam waist in m
    # depth_K is trap depth in Kelvin
    # Waist is the beam waist in m
    depth_J = depth_K * boltzmann # (in J)
    depth_Hz = depth_J / planck # (in Hz)
    waist = calculate_waist(L, wavelength)
    I = depth_Hz/alpha # (in W m^-2)
    P = (np.pi * waist**(2) * I)/8 # (in W)
    return P

waist = calculate_waist(0.30, 515*10**(-9))
print(waist)
P = calculate_power(100*10**(-6), 0.30, 515*10**(-9))
print(P)

def trap_potential(x, w_0, depth_K):
    """Calculates the trap potential at different positions for a single trap."""
    coefficient = -depth_K * boltzmann
    exponent = -2*(x/w_0)**2
    return coefficient * np.exp(exponent)

def trap_potential_with_g(x, w_0, depth_K):
    """Calculates the trap potential + gravitational potential at different positions for a single trap."""
    coefficient = -depth_K * boltzmann
    exponent = -2*(x/w_0)**2
    m = 2.8733965 * 10**(-25)
    g = 9.8
    return coefficient * np.exp(exponent) + m*g*x

L = 0.30 # (in m)
wavelength = 515*10**(-9) # in m
w_0 = calculate_waist(L, wavelength)
depth_K_0 = 100*10**(-6) # in uK
depth_K_1 = 50*10**(-6) # in uK
depth_K_2 = 200*10**(-6) # in uK

# Generate x values for the plot
# Typically, you'd cover a range like mean - 3*std_dev to mean + 3*std_dev
x = np.linspace(-300*10**(-6), 300*10**(-6), 500)

# Calculate the corresponding y values (PDF values)
y_0 = trap_potential(x, w_0, depth_K_0)/boltzmann # (in K)
y_1 = trap_potential_with_g(x, w_0, depth_K_0)/boltzmann # (in K)
y_2 = trap_potential_with_g(x, w_0, depth_K_1)/boltzmann # (in K)
y_3 = trap_potential_with_g(x, w_0, depth_K_2)/boltzmann # (in K)

# Plotting
plt.figure(figsize=(10, 6))
plt.plot(x*10**6, y_0*10**6, label=f'No gravity, Power={calculate_power(depth_K_0, L, wavelength):.4g} W')
plt.plot(x*10**6, y_1*10**6, label=f'Power={calculate_power(depth_K_0, L, wavelength):.4g} W')
plt.plot(x*10**6, y_2*10**6, label=f'Power={calculate_power(depth_K_1, L, wavelength):.4g} W')
plt.plot(x*10**6, y_3*10**6, label=f'Power={calculate_power(depth_K_2, L, wavelength):.4g} W')
plt.xlabel('Position (in um)')
plt.ylabel('Potential energy (in uK)')
plt.title('Trap depth + Gravitational potential')
plt.grid(True)
plt.legend(loc='upper right', fontsize=10)
plt.show()
