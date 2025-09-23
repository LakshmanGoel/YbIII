import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks

# Fundamental constants
planck = 6.626*10**(-34) # (in SI units)
boltzmann = 1.38*10**(-23) # (in SI units)

def calculate_waist(L, wavelength):
    # Beam waist calculation from a given conveyor belt length and laser wavelength
    # L is conveyor belt length (in m)
    # wavelength is the laser wavelength (in m)
    rayleigh = L/2 # (in m)
    waist = np.sqrt(wavelength*rayleigh/(np.pi)) # (in m)
    return waist

def calculate_power(depth_K, L, wavelength, a):
    # Power calculation from a given trap depth in K and beam waist
    # depth_K is trap depth in Kelvin
    # Waist is the beam waist in m
    depth_J = depth_K * boltzmann # (in J)
    depth_Hz = depth_J / planck # (in Hz)
    waist = calculate_waist(L, wavelength)
    I = depth_Hz/a # (in W m^-2)
    P = (np.pi * waist**(2) * I)/8 # (in W)
    return P

# Polarizability
alpha = 160 * (2* 377 * 2.48832 * 10**(-8)) # (in Hz*m^2/W) [Depends on laser wavelength and atomic state (In this case, we consider a starting state of 1S0)]
waist = calculate_waist(0.30, 1036*10**(-9))
print(waist)
P = calculate_power(100*10**(-6), 0.30, 1036*10**(-9), alpha)
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
L_1 = 0.40 # (in m)
L_2 = 0.20 # (in m)

wavelength = 515*10**(-9) # in m
wavelength_1 = 532*10**(-9) # in m
wavelength_2 = 1036*10**(-9) # in m

alpha = 280 * (2* 377 * 2.48832 * 10**(-8)) # (in Hz*m^2/W)
alpha_1 = 180 * (2* 377 * 2.48832 * 10**(-8)) # (in Hz*m^2/W)
alpha_2 = 160 * (2* 377 * 2.48832 * 10**(-8)) # (in Hz*m^2/W)

w_0 = calculate_waist(L, wavelength) # in m
w_1 = calculate_waist(L_1, wavelength) # in m
w_2 = calculate_waist(L_2, wavelength) # in m
w_3 = calculate_waist(L, wavelength_1) # in m
w_4 = calculate_waist(L, wavelength_2) # in m

depth_K_0 = 100*10**(-6) # in uK
depth_K_1 = 50*10**(-6) # in uK
depth_K_2 = 200*10**(-6) # in uK

# Generate x values for the plot
# Typically, you'd cover a range like mean - 3*std_dev to mean + 3*std_dev
x = np.linspace(-300*10**(-6), 300*10**(-6), 500)

# Calculate the corresponding y values (PDF values)
y_0 = trap_potential(x, w_0, depth_K_0)/boltzmann # (in K)
y_1 = trap_potential_with_g(x, w_0, depth_K_0)/boltzmann # (in K)
y_2 = trap_potential_with_g(x, w_3, depth_K_0)/boltzmann # (in K)
y_3 = trap_potential_with_g(x, w_4, depth_K_0)/boltzmann # (in K)

# Plotting different trap depths for different input powers
# plt.figure(figsize=(10, 6))
# plt.plot(x*10**6, y_0*10**6, label=f'No gravity, Power={calculate_power(depth_K_0, L, wavelength):.4g} W')
# plt.plot(x*10**6, y_1*10**6, label=f'Power={calculate_power(depth_K_0, L, wavelength):.4g} W')
# plt.plot(x*10**6, y_2*10**6, label=f'Power={calculate_power(depth_K_1, L, wavelength):.4g} W')
# plt.plot(x*10**6, y_3*10**6, label=f'Power={calculate_power(depth_K_2, L, wavelength):.4g} W')
# plt.xlabel('Position (in um)')
# plt.ylabel('Potential energy (in uK)')
# plt.title('Trap depth + Gravitational potential')
# plt.grid(True)
# plt.legend(loc='upper right', fontsize=10)
# plt.show()

# Plotting different trap depths for different conveyor belt lengths
# plt.figure(figsize=(10, 6))
# plt.plot(x*10**6, y_0*10**6, label=f'No gravity, Length={L} m, Power={calculate_power(depth_K_0, L, wavelength, alpha):.4g} W')
# plt.plot(x*10**6, y_1*10**6, label=f'Length={L} m, Power={calculate_power(depth_K_0, L, wavelength, alpha):.4g} W')
# plt.plot(x*10**6, y_2*10**6, label=f'Length={L_1} m, Power={calculate_power(depth_K_0, L_1, wavelength, alpha):.4g} W')
# plt.plot(x*10**6, y_3*10**6, label=f'Length={L_2} m, Power={calculate_power(depth_K_0, L_2, wavelength, alpha):.4g} W')
# # Find local maxima
# peaks1, _ = find_peaks(y_1)
# peaks2, _ = find_peaks(y_2)
# peaks3, _ = find_peaks(y_3)
# minima2, _ = find_peaks(-y_2)
# for i, (x_peaks, y_peaks, color, label) in enumerate([
#     (x[peaks1], y_1[peaks1], "black", None),
#     (x[peaks2], y_2[peaks2], "black", None),
#     (x[peaks3], y_3[peaks3], "black", None),
#     (x[minima2], y_2[minima2], "black", None)
# ]):
#     # Convert to um and uK for plotting
#     x_vals = x_peaks * 1e6
#     y_vals = y_peaks * 1e6
    
#     # Plot the peaks
#     plt.plot(x_vals, y_vals, "o", color=color, label=label)
    
#     # Annotate each peak
#     for x_val, y_val in zip(x_vals, y_vals):
#         plt.text(x_val, y_val + 1, f"{y_val:.2f} µK", fontsize=8, color=color, ha='center')
# plt.xlabel('Position (in um)')
# plt.ylabel('Potential energy (in uK)')
# plt.title('Trap depth + Gravitational potential')
# plt.grid(True)
# plt.legend(loc='upper right', fontsize=10)
# plt.show()

# Plotting different trap depths for different laser wavelengths
plt.figure(figsize=(10, 6))
plt.plot(x*10**6, y_0*10**6, label=f'No gravity, $\lambda$={wavelength*10**9:.4g} nm, Power={calculate_power(depth_K_0, L, wavelength, alpha):.4g} W')
plt.plot(x*10**6, y_1*10**6, label=f'$\lambda$={wavelength*10**9:.4g} nm, Power={calculate_power(depth_K_0, L, wavelength, alpha):.4g} W')
plt.plot(x*10**6, y_2*10**6, label=f'$\lambda$={wavelength_1*10**9:.4g} nm, Power={calculate_power(depth_K_0, L, wavelength_1, alpha_1):.4g} W')
plt.plot(x*10**6, y_3*10**6, label=f'$\lambda$={wavelength_2*10**9:.4g} nm, Power={calculate_power(depth_K_0, L, wavelength_2, alpha_2):.4g} W')
# Find local maxima
peaks1, _ = find_peaks(y_1)
peaks3, _ = find_peaks(y_3)
minima3, _ = find_peaks(-y_3)
for i, (x_peaks, y_peaks, color, label) in enumerate([
    (x[peaks1], y_1[peaks1], "black", None),
    (x[peaks3], y_3[peaks3], "black", None),
    (x[minima3], y_3[minima3], "black", None)
]):
    # Convert to um and uK for plotting
    x_vals = x_peaks * 1e6
    y_vals = y_peaks * 1e6
    
    # Plot the peaks
    plt.plot(x_vals, y_vals, "o", color=color, label=label)
    
    # Annotate each peak
    for x_val, y_val in zip(x_vals, y_vals):
        plt.text(x_val, y_val + 1, f"{y_val:.2f} µK", fontsize=8, color=color, ha='center')
plt.xlabel('Position (in um)')
plt.ylabel('Potential energy (in uK)')
plt.title('Trap depth + Gravitational potential')
plt.grid(True)
plt.legend(loc='upper right', fontsize=10)
plt.show()
