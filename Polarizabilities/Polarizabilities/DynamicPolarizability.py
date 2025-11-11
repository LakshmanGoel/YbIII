import numpy as np 
import scipy.constants as const 
import lib.wigner_symbols as wig
from lib.Yb import *
import matplotlib.pyplot as plt

class DyanmicPolarizability:

    # important constants
    global c, h, hbar, pi, k, inverse_cm_to_Hz, dip_au, pol_au
    c = const.c 
    h = const.h
    hbar = const.hbar
    pi = const.pi 
    k = const.k
    inverse_cm_to_Hz = c * 100 # convert from cm^-1 to Hz
    # SI to au
    dip_au = const.physical_constants['atomic unit of electric dipole mom.'][0] 
    pol_au = const.physical_constants['atomic unit of electric polarizability'][0] 

    def __init__(self, atom, n, L, S, J, I=0, F=0):
        self.atom = atom
        self.n = n 
        self.L = L 
        self.S = S 
        self.J = J 
        self.I = I
        self.F = F 
        

    
    def alpha_red_nJ(self, K, wavlen): 
        # K is tensor rank
        omega = 2 * pi * c / (wavlen*1e-9)
        target_wavvec = self.atom.getEnergy(self.n, self.L, self.S, self.J)
        target_omega = target_wavvec * 2 * pi * inverse_cm_to_Hz
        if K == 0:
            prefac = 2 / np.sqrt(3 * (2 * self.J + 1)) / hbar

            term = 0

            # look through list of dipole-coupled states
            for key in self.atom.getRDMEKeys(self.n, self.L, self.S, self.J): 
                trans_state = self.atom.decode(key)
                omega_ki = (self.atom.getEnergy(*trans_state) * 2 * pi * inverse_cm_to_Hz) - target_omega
                rdme = self.atom.getRDME(self.n, self.L, self.S, self.J,*trans_state)
                term += omega_ki/ (omega_ki**2 - omega**2) * (rdme * dip_au)**2
            return term * prefac

        elif K == 1:
            prefac = 2 * np.sqrt(3) / hbar
            term = 0

            for key in self.atom.getRDMEKeys(self.n, self.L, self.S, self.J):

                trans_state = self.atom.decode(key)
                trans_state_J = trans_state[-2]
                omega_ki = self.atom.getEnergy(*trans_state) * 2 * pi * inverse_cm_to_Hz - target_omega
                rdme = self.atom.getRDME(self.n, self.L, self.S, self.J,*trans_state)
                term += omega/ (omega_ki**2 - omega**2) * (rdme*dip_au)**2 \
                        * wig.wigner_6j(1, 1, 1, self.J, trans_state_J, self.J) * (-1)**(self.J + trans_state_J)
            return term * prefac


        elif K == 2:
            prefac = -2 * np.sqrt(5) / hbar
            term = 0 

            for key in self.atom.getRDMEKeys(self.n, self.L, self.S, self.J):
                trans_state = self.atom.decode(key)
                trans_state_J = trans_state[-2]
                omega_ki = self.atom.getEnergy(*trans_state) * 2 * pi * inverse_cm_to_Hz - target_omega
                rdme = self.atom.getRDME(self.n, self.L, self.S, self.J,*trans_state)
                term += omega_ki/ (omega_ki**2 - omega**2) * (rdme*dip_au)**2 \
                        * wig.wigner_6j(1, 2, 1, self.J, trans_state_J, self.J) * (-1)**(self.J + trans_state_J)
            return term * prefac



        else:
            print("K must be 0, 1, or 2!")
            return None


    def getPolarizabilities(self, wavlen, units = 'SI'):
        n = self.n 
        J = self.J 
        I = self.I 
        F = self.F 

        if I == 0:
            alpha_s = self.alpha_red_nJ(0, wavlen) / np.sqrt(3 * (2 * J + 1))
            alpha_v = -1 * np.sqrt(2 * J / ( (J+1) * (2 * J + 1) )) * self.alpha_red_nJ(1, wavlen) 
            alpha_t = -1 * np.sqrt( ((2 * J)*(2 * J - 1)) / (3 * (J + 1) * (2 * J + 1) * (2 * J + 3)) ) * self.alpha_red_nJ(2, wavlen) 
        else:
            alpha_s = self.alpha_red_nJ(0, wavlen) / np.sqrt(3 * (2 * J + 1))
            alpha_v = (-1)**(J + I + F) * np.sqrt(2 * F * (2 * F + 1)/ ( (F+1))) * wig.wigner_6j(F, 1, F, J, I, J) * self.alpha_red_nJ(1, wavlen) 
            alpha_t = -(-1)**(J + I + F) * np.sqrt( ((2 * F)*(2 * F - 1)* (2 * F + 1)) / (3 * (F+1)  * (2 * F + 3)) ) * wig.wigner_6j(F, 2, F, J, I, J) * self.alpha_red_nJ(2, wavlen) 

        if units == 'a.u.':
            return np.array([alpha_s, alpha_v, alpha_t])/pol_au
        elif units == 'SI':
            return np.array([alpha_s, alpha_v, alpha_t])



