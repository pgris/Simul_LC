from Throughputs import Throughputs
import numpy as np
import matplotlib.pyplot as plt
import sncosmo
 
version='1.0'
model ='salt2-extended'
model_min=300.
model_max=180000.
wave_min=3000.
wave_max=11501.

wave= np.arange(wave_min,wave_max,1.)
sn_type='Ia'
source=sncosmo.get_source(model,version=version)
SN=sncosmo.Model(source=source)

z=1.e-8
t0=0.
x1=0.
c=0.

SN.set(z=z)
SN.set(t0=t0)
SN.set(c=c)
SN.set(x1=x1)
#SN.set(x0=X0)

SN.set_source_peakabsmag(-19.6, 'bessellB', 'vega')

fluxes=10.*SN.flux(0.,wave)
        
wavelength=wave/10.

plt.plot(wavelength,fluxes,ls='-',color='b')


throughputs=Throughputs(atmos=False,aerosol=False)

res=[]

lambda_SN=300.

for band in 'ugrizy':
    for thresh in [0.05,0.10,0.15]:
        for z in np.arange(0.,1.5,0.05):
            idx = throughputs.system[band].sb >= thresh
            waves=throughputs.system[band].wavelen[idx]
            lambda_min=np.min(waves)
            #print band,thresh,z,lambda_min/(1.+z)
            res.append((band,thresh,z,lambda_min,lambda_SN*(1.+z)))

rec=np.rec.fromrecords(res,names=['band','thresh','z','lambda_min','lambdaSN_obs'])

print rec
filtercolors = {'u':'b', 'g':'c', 'r':'g', 'i':'y', 'z':'r', 'y':'m'}

figa, axa = plt.subplots(ncols=1, nrows=1, figsize=(10,9))
fontsize=15
myls={0.05 : '-',0.10: '--',0.15 :':'}
for band in 'ugrizy':
    idx = rec['band'] == band
    sela=rec[idx]
    for j,thresh in enumerate([0.05,0.10,0.15]):
        idxb = sela['thresh']==thresh
        selb=sela[idxb]
        axa.plot(selb['z'],selb['lambda_min'],linestyle=myls[thresh],color=filtercolors[band])
        if j == 0:
            axa.text(1.47,np.mean(selb['lambda_min']),band,fontsize=fontsize)

axa.plot(selb['z'],selb['lambdaSN_obs'],linestyle='-',color='k',label='$\lambda_{SN}^{min}$=300 nm')


axa.set_xlabel('z',{'fontsize': fontsize})
axa.set_ylabel('$\lambda^{obs}_{min}$ [nm] ',{'fontsize': fontsize})
axa.legend(loc='best',prop={'size':fontsize})
axa.set_xlim(0.,1.55)
axa.set_xticks([x for x in np.arange(0.,1.5,0.1)])
plt.show()
