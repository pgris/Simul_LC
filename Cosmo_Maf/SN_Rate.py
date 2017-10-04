from astropy import (units as u, constants as const)
from astropy.cosmology import FlatLambdaCDM
import numpy as np

STERADIAN2SQDEG = 180.**2 / np.pi**2
# Mpc^3 -> Mpc^3/sr
norm = 1. / (4. * np.pi) 

class SN_Rate:
    def __init__(self, rate='Ripoche',H0=70,Om0=0.25,survey_area=9.6,duration=0.5,min_rf_phase=-15.,max_rf_phase=30.):

        self.astropy_cosmo=FlatLambdaCDM(H0=H0, Om0=Om0)
        self.rate=rate
        self.survey_area=survey_area
        self.duration=duration
        self.min_rf_phase=min_rf_phase
        self.max_rf_phase=max_rf_phase

    def __call__(self,zmin=0.1,zmax=0.2,dz=0.01,bins=None,account_for_edges=False):
        
        if bins is None:
            thebins = np.arange(zmin, zmax+dz, dz)
            zz = 0.5 * (thebins[1:] + thebins[:-1])
            
        else:
            zz=bins
            thebins=bins
            print thebins
            thebins=np.append(thebins,thebins[len(bins)-1]+(bins[1]-bins[0]))

        rate, err_rate = self.sn_rate(zz)
                
        area = self.survey_area / STERADIAN2SQDEG # or area= self.survey_area/41253.
    
        dvol = norm*self.astropy_cosmo.comoving_volume(thebins).value
        #print 'after rate',len(zz),len(rate),len(dvol),area,rate
        dvol = dvol[1:] - dvol[:-1]
        #print 'dvol',dvol

        
        if account_for_edges:
            margin = (1.+zz) * (self.max_rf_phase-self.min_rf_phase) / 365.25
            effective_duration = self.duration - margin
            effective_duration[effective_duration<=0.] = 0.
        else:
            effective_duration = self.duration

        nsn=rate * area * dvol * effective_duration / (1.+zz)
        err_nsn=err_rate* area * dvol * effective_duration / (1.+zz)
        return zz,rate, err_rate,nsn, err_nsn
        

    def ripoche_rate(self, z):
        """The SNLS SNIa rate according to the (unpublished) Ripoche et al study.
        """
        rate = 1.53e-4*0.343
        expn = 2.14
        my_z = np.copy(z)
        my_z[my_z>1.] = 1.
        rate_sn=rate * np.power((1+my_z)/1.5, expn)
        return rate_sn,0.2*rate_sn

    def perrett_rate(self, z):
        """The SNLS SNIa rate according to (Perrett et al, 201?)
        """
        rate = 0.17E-4
        expn = 2.11
        err_rate=0.03E-4
        err_expn=0.28
        my_z = np.copy(z)
        my_z[my_z>1.] = 1.
        rate_sn=rate * np.power(1+my_z, expn)
        
        #err_rate_sn=np.power(err_rate/rate,2.)+np.power(np.log(1+my_z)*err_expn,2.)
        err_rate_sn=np.power(1+my_z, 2.*expn)*np.power(err_rate,2.)+np.power(rate_sn*np.log(1+my_z)*err_expn,2.)
        return rate_sn,np.power(err_rate_sn,0.5)

    def dilday_rate(self,z):

        rate=2.6e-5
        expn=1.5
        err_rate=0.01
        err_expn=0.6
        my_z = np.copy(z)
        my_z[my_z>1.] = 1.
        rate_sn=rate * np.power(1+my_z, expn)
        err_rate_sn=rate_sn*np.log(1+my_z)*err_expn
        return rate_sn,err_rate_sn

    def sn_rate(self,z):
        if self.rate == 'Ripoche':
            return self.ripoche_rate(z)
        if self.rate == 'Perrett':
            return self.perrett_rate(z)
        if self.rate == 'Dilday':
            return self.dilday_rate(z)




        
