from Yb import Yb
import scipy.constants as const

global amu
amu = const.physical_constants["atomic mass constant"][0]

class Yb171(Yb):

   
    def __init__(self):
        self.isotope = "171"
        self.I = 1/2
        self.m = Yb.isotope_mass[self.isotope] * amu
        self.abundance = Yb.isotope_abundance[self.isotope]



atom = Yb171()
print(atom.abundance)



    


