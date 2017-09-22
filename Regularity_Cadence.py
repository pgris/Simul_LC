import numpy as np
from Generate_Single import Generate_Single_LC
from Telescope import *

class Regularity_Study:
    def __init__(self, cadence=1.,random=True,X1=0.,Color=0.,expTime=dict(zip([b for b in "grizy"],[300, 600., 600., 780.,600.]))):
        
        self.cadence=cadence
        self.X1=X1
        self.Color=Color
        self.DayMax=0.
        self.bands='grizy'
        self.random=random
        self.expTime=expTime

        self.zrange=np.arange(0.,1.,0.05)

        self.m5_lim={'u':23.61,'g':24.83,'r':24.35,'i':23.88,'z':23.30,'y':22.43}
        self.seeing={'u':0.92,'g':0.87,'r':0.83,'i':0.80,'z':0.78,'y':0.76}
        self.sky={'u':22.95,'g':22.24,'r':21.20,'i':20.47,'z':19.60,'y':18.63}

        #this m5_lim values correspond to an exposure time of 30s...needs to be corrected for DDF

        for key, vals in expTime.items():
            self.m5_lim[key]=self.m5_lim[key]+1.25*np.log10(vals/30.)


        restframe_phase_range = (-20., 40.)	
        self.pmin, self.pmax = restframe_phase_range

        self.airmass=1.2
        self.telescope=Telescope(atmos=True,airmass=self.airmass)

        for z in self.zrange:
            self.Simulate(z)

    def Simulate(self,z):

        mjd_min = np.floor(self.pmin * (1.+z) + self.DayMax)
        mjd_max = np.ceil(self.pmax * (1.+z) + self.DayMax)
        mjd_tot = np.arange(mjd_min, mjd_max, self.cadence)
        if self.random:
            mjd=self.Get_mjds(mjd_min,mjd_max,self.cadence)
        else:
            mjd=mjd_tot

        obs=self.Make_Obs(mjd)
        
        lc=Generate_Single_LC(z,self.DayMax,self.X1,self.Color,obs,self.telescope,0,None)



    def Get_mjds(self,mjd_min,mjd_max,cadence):
    
        mjd=[]
        nvals=int((mjd_max-mjd_min)/cadence)
    #print 'hello',mjd_min,mjd_max,cadence,nvals
        for i in range(10*nvals):
            date=int(np.random.uniform(mjd_min,mjd_max))
            if date not in mjd:
                mjd.append(date)
            if len(mjd)==nvals:
                break

        mjd=np.sort(mjd)

        return mjd

    def Make_Obs(self, mjds):
        
        r=[]
        for b in self.bands:
            for mjd in mjds:
                r.append((b,mjd,self.expTime[b],-1,self.seeing[b],0.,self.sky[b],0.,self.airmass,self.m5_lim[b],int(self.expTime[b]/30.),-1,-1))

        return np.rec.fromrecords(r,names=['band','mjd','exptime','rawSeeing ','seeing','moon_frac','sky','kAtm','airmass','m5sigmadepth','Nexp','Ra','Dec'])
            

X1=0.
Color=0.

reg=Regularity_Study(X1=X1,Color=Color)
