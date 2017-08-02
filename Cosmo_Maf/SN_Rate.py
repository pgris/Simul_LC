from astropy import (units as u, constants as const)
from astropy.cosmology import FlatLambdaCDM
import numpy as np

STERADIAN2SQDEG = 180.**2 / np.pi**2
# Mpc^3 -> Mpc^3/sr
norm = 1. / (4. * np.pi) 

class SN_Rate:
    def __init__(self, rate='Ripoche',H0=70,Om0=0.25,survey_area=9.6,effective_duration=0.5):

        self.astropy_cosmo=FlatLambdaCDM(H0=H0, Om0=Om0)
        self.rate=rate
        self.survey_area=survey_area
        self.effective_duration=effective_duration

    def __call__(self,zmin,zmax,dz=0.1):
        
        bins = np.arange(zmin, zmax+dz, dz)
        zz = 0.5 * (bins[1:] + bins[:-1])
        rate = self.sn_rate(zz)
        
        area = self.survey_area / STERADIAN2SQDEG
        dvol = norm * self.astropy_cosmo.comoving_volume(bins).value
        dvol = dvol[1:] - dvol[:-1]
        

        return zz,rate, rate * self.survey_area * dvol * self.effective_duration / (1.+zz)
        


    def ripoche_rate(self, z):
        """The SNLS SNIa rate according to the (unpublished) Ripoche et al study.
        """
        rate = 1.53e-4*0.343
        expn = 2.14
        my_z = np.copy(z)
        my_z[my_z>1.] = 1.
        return rate * np.power((1+my_z)/1.5, expn)

    def perret_rate(self, z):
        """The SNLS SNIa rate according to (Perret et al, 201?)
        """
        rate = 0.17E-4
        expn = 2.11
        my_z = np.copy(z)
        my_z[my_z>1.] = 1.
        return rate * np.power(1+my_z, expn)

    def sn_rate(self,z):
        if self.rate == 'Ripoche':
            return self.ripoche_rate(z)
        if self.rate == 'Perret':
            return self.perret_rate(z)




        
