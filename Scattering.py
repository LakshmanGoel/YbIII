import numpy as np
import matplotlib.pyplot as plt

# Using this script, we can calculate the scattering rate for a given transition
# In this case, we want to calculate the scattering rate due to Yb atoms for the 1S0->1P1 and 1S0->3P1 transitions.

# Fundamental constants
planck = 6.626*10**(-34) # (in SI units)
boltzmann = 1.38*10**(-23) # (in SI units)

def scattering_rate(Gamma, depth_K, a, I_sat, wavelength, w_0):
    # This function calculates the scattering rate for a given system
    # Gamma: Natural Linewidth (in Hz)
    # wavelength (in nm)
    # w_0: Resonance frequency for the transition (in Hz)
    # depth_K: Trap depth (in K)
    # a: Polarizability (in amu)
    # I_sat: Saturation intensity for the transition
    w = (2*np.pi*3*10**8)/wavelength
    Delta = w - w_0
    depth_J = depth_K * boltzmann
    depth_Hz = depth_J / planck
    alpha = a * (2* 377 * 2.48832 * 10**(-8)) # in Hz*m^2/W
    I = depth_Hz/alpha # in Wm^-2
    rate = (Gamma/2) * ((I/I_sat)/((I/I_sat) + 4*(Delta/Gamma)**2))
    return rate

# For a 1S0 to 1P1 transition (515 nm):

rate_1 = scattering_rate(2*np.pi*29.13*10**6, 100*10**(-6), 280, 59.97*(10**(-3))/(10**(-4)), 515*10**(-9), 2*np.pi*751.527*10**12)
print(rate_1)

# For a 1S0 to 3P1 transition (515 nm):

rate_2 = scattering_rate(2*np.pi*182.4*10**3, 100*10**(-6), 280, 138.85*(10**(-6))/(10**(-4)), 515*10**(-9), 2*np.pi*539.388*10**12)
print(rate_2)

# Total scattering rate (515 nm):

rate = rate_1 + rate_2
print(rate)