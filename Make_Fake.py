import numpy as np
import os

legend=['band','mjd','exptime','rawSeeing','FWHMeff','moon_frac','sky','kAtm','airmass','m5sigmadepth','Nexp','Ra','Dec']

outdir='../Ana_Cadence/OpSimLogs/Mean_Obs_newrefs'
if not os.path.exists(outdir):
    os.makedirs(outdir)

outputfile  = open(outdir+'/Observations_WFD_309.txt','wb')

for leg in legend:
    lego=leg+' : '
    if leg == 'FWHMeff':
        lego='seeing : [was FWHMeff]'
    outputfile.write('# '+lego+'\n')
outputfile.write('# end\n')

#m5_lim={'u':23.61,'g':24.83,'r':24.35,'i':23.88,'z':23.30,'y':22.43}
m5_lim={'u': 23.103921,'g': 24.2292315,'r': 23.8573,'i': 23.4196475,'z': 22.792551,'y': 21.5026765}

expTime=dict(zip([b for b in 'ugrizy'],[30.,30.,30.,30.,30.,30.]))
rawSeeing=-1
seeing={'u':0.92,'g':0.87,'r':0.83,'i':0.80,'z':0.78,'y':0.76}
moonphase=0
sky={'u':22.95,'g':22.24,'r':21.20,'i':20.47,'z':19.60,'y':18.63}
kAtm=-1
airmass=1.2
Nexp=1
Ra=-1
Dec=-1
delay=60/24./3600.
decal={'u':0.,'g':delay,'r':2.*delay,'i':3.*delay,'z':4.*delay,'y':5.*delay}

mjd_min=-20*(1.+0.7)
mjd_max=40*(1.+0.7)
sampling=4

mjds = np.arange(mjd_min, mjd_max, sampling)

bands='ugrizy'

for b in bands:
    for mjd in mjds:
        toprint = 'LSSTPG::'+b+' '
        toprint+=str(format(mjd+decal[b],'3.8f'))+' '
        toprint+=str(expTime[b])+' '
        toprint+=str(format(rawSeeing,'.7f'))+' '
        toprint+=str(format(seeing[b],'.7f'))+' '
        toprint+=str(format(moonphase,'.7f'))+' '
        toprint+=str(format(sky[b],'.7f'))+' '
        toprint+=str(kAtm)+' '
        toprint+=str(format(airmass,'.7f'))+' '
        toprint+=str(format(m5_lim[b],'.7f'))+' '
        toprint+=str(Nexp)+' '
        toprint+=str(format(Ra,'.7f'))+' '
        toprint+=str(format(Dec,'.7f'))
        print toprint
        outputfile.write(toprint+'\n')

outputfile.close()
