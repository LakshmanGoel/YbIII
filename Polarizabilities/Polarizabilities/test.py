from DynamicPolarizability import *
import matplotlib.pyplot as plt 
import matplotlib.ticker as tick
plt.rcParams["lines.markersize"] = 3


# Initialize atomic states of interest
atom = Yb()
pol_1S0 = DyanmicPolarizability(atom,atom.ground_state,L=0,S=0,J=0)
pol_3P0 = DyanmicPolarizability(atom,atom.ground_state,L=1,S=1,J=0)
pol_3P1 = DyanmicPolarizability(atom,atom.ground_state,L=1,S=1,J=1)
# pol_3P2 = DyanmicPolarizability(atom,atom.ground_state,L=1,S=1,J=2)

# Specify wavelength(s) for polarizability calculations
wavelength = np.linspace(679,680.2,10000)

pol_1S0 = pol_1S0.getPolarizabilities(wavelength, units='a.u.')[0]
pol_3P0 = pol_3P0.getPolarizabilities(wavelength, units='a.u.')[0]

# Calculations for 3P1 states include empirical correction factors
pol_3P1_mJ0 = (pol_3P1.getPolarizabilities(wavelength, units='a.u.')[0]-27) - 2*(pol_3P1.getPolarizabilities(wavelength, units='a.u.')[-1]-61)
pol_3P1_mJ1 = (pol_3P1.getPolarizabilities(wavelength, units='a.u.')[0]-27) + (pol_3P1.getPolarizabilities(wavelength, units='a.u.')[-1]-61)

# Convert from J basis to F basis using CG coefficients
pol_3P1_F_3h_mF_1h = pol_3P1_mJ1/3 + pol_3P1_mJ0*2/3

fig, ax = plt.subplots()

plt.plot(wavelength,pol_1S0, label = '1S0')
# plt.plot(wavelength,pol_3P0, label = '3P0')
# plt.plot(wavelength, pol_3P1_mJ0, label = '3P1, mJ=0')
# plt.plot(wavelength, pol_3P1_F_3h_mF_1h, label = '3P1, |mF|=1/2')
plt.plot(wavelength, pol_3P1_mJ1, label = '3P1, |mF|=3/2')
# plt.legend()
# plt.xlabel("Wavelength (nm)")
# plt.ylabel("Polarizability")
# plt.tight_layout()
# plt.show()


# np.savez('pol_mcr.npz', wavelength=wavelength, pol_1S0=pol_1S0, pol_3P1_mF_1h = pol_3P1_F_3h_mF_1h, pol_3P1_F_3h = pol_3P1_mJ1)
