import sncosmo
import numpy as np
from astropy import (cosmology, units as u, constants as const)
import matplotlib.pyplot as plt

class Fit_LC:
    def __init__(self,model='salt2-extended',version='1.0',z=0.1,telescope=None,Plot=False,bands='ugrizy'):

        self.Plot=Plot
        transmission=telescope.throughputs
        self.z=z
        #bands=[b[-1:] for b in np.unique(select['band'])]

        #print 'hello',bands,telescope.airmass
        for filtre in bands:
            if telescope.airmass > 0:
                band=sncosmo.Bandpass(transmission.atmosphere[filtre].wavelen,transmission.atmosphere[filtre].sb, name='LSST::'+filtre,wave_unit=u.nm)
            else:
                band=sncosmo.Bandpass(transmission.system[filtre].wavelen,transmission.system[filtre].sb, name='LSST::'+filtre,wave_unit=u.nm) 
            sncosmo.registry.register(band, force=True)

        source=sncosmo.get_source(model,version=version)
        dust = sncosmo.OD94Dust()

          
        self.SN_fit_model=sncosmo.Model(source=source)
        self.SN_fit_model.set(z=self.z)
        #SN_fit_model.set_source_peakabsmag(peakAbsMagBesselB, 'bessellB', 'vega',cosmo=self.astropy_cosmo)
    def __call__(self,meas):
        t=meas
        select=t[np.where(np.logical_and(t['flux']/t['fluxerr']>5.,t['flux']>0.))]
        idx=select['band']!='LSST::u'
        select=select[idx]

        if self.z > 0.35:
           idx=select['band']!='LSST::g'
           select=select[idx] 


        try:
            #print 'trying to fit',len(select)
            res, fitted_model = sncosmo.fit_lc(select, self.SN_fit_model,['t0', 'x0', 'x1', 'c'],bounds={'z':(self.z-0.1, self.z+0.1)})

            #self.sigma_c=res['errors']['c']
            mbfit=fitted_model._source.peakmag('bessellb','vega')
            
            if self.Plot:
                sncosmo.plot_lc(select, model=fitted_model,color='r',pulls=False,errors=res.errors) 
                plt.show()
            return res,fitted_model,mbfit,'ok'
        except (RuntimeError, TypeError, NameError):
                """
                print select
                self.Plot_bands(select)
                plt.show()
                self.sigma_c=0.
                print 'crashed'
                """
                return None,None,-1,'crash'
    @property
    def sigma_c(self):
        return self.sigma_c

    def Plot_bands(self,obj,time_name='time',flux_name='flux',errflux_name='fluxerr',filter_name='band',T0=0):
        
        figa, axa = plt.subplots(ncols=2, nrows=3, figsize=(10,9))
    
        for j,band in enumerate(['u','g','r','i','z','y']):
            if j<2:
                k=0
            if j>= 2 and j < 4:
                k=1
            if j>=4:
                k=2

            selobs=obj[np.where(obj[filter_name]=='LSST::'+band)]

            axa[k][j%2].errorbar(selobs[time_name],selobs[flux_name],yerr=selobs[errflux_name],fmt='.',ecolor='r',color='r')
           
        axa[k][j%2].set_xlim(T0-30,T0+50)
