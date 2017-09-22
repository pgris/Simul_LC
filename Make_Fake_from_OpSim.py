import numpy as np
from Observations import *

dirmeas='DD'
fieldname='DD'

thedir='../Ana_Cadence/OpSimLogs/'+dirmeas

fieldid=1427

name='Observations_'+fieldname+'_'+str(fieldid)
myobs=Observations(fieldid=fieldid, filename=thedir+'/'+name+'.txt')

bands='ugrizy'
shift_band=dict(zip(bands,[0.,0.,0.05,0.1,0.15,0.20]))


legend=['band','mjd','exptime','rawSeeing','FWHMeff','moon_frac','sky','kAtm','airmass','m5sigmadepth','Nexp','Ra','Dec']

m5_lim={'u': 23.103921,'g': 24.2292315,'r': 23.8573,'i': 23.4196475,'z': 22.792551,'y': 21.5026765}
        
seeing={'u':0.92,'g':0.87,'r':0.83,'i':0.80,'z':0.78,'y':0.76}
sky={'u':22.95,'g':22.24,'r':21.20,'i':20.47,'z':19.60,'y':18.63}


outdir='../Ana_Cadence/OpSimLogs/Mean_Obs_newrefs/'

for season in range(len(myobs.seasons)):
#for season in range(1):
    thename='Observations_'+fieldname
    outputfile  = open(outdir+thename+'_'+str(season)+'_1_'+str(fieldid)+'.txt','wb')
    outputfile_b  = open(outdir+thename+'_'+str(season)+'_2_'+str(fieldid)+'.txt','wb')

    for leg in legend:
        lego=leg+' : '
        if leg == 'FWHMeff':
            lego='seeing : [was FWHMeff]'
        outputfile.write('# '+lego+'\n')
        outputfile_b.write('# '+lego+'\n')
    outputfile.write('# end\n')
    outputfile_b.write('# end\n')

    myseason=myobs.seasons[season]
    #print myseason
    #for b in 'ugrizy':
    for b in bands:
        idx = myseason['band']=='LSSTPG::'+b
        sel=myseason[idx]
        sel.sort(order='mjd')
        print sel
        min_time=np.min(sel['mjd'])
        max_time=np.max(sel['mjd'])
        duration=max_time-min_time

        for vals in sel:
            new_mjd=vals['mjd']-min_time-20.+shift_band[b]
            print >> outputfile,vals['band'],new_mjd,vals['exptime'],vals['rawSeeing'],vals['seeing'],vals['moon_frac'],vals['sky'],vals['kAtm'],vals['airmass'],vals['m5sigmadepth'],vals['Nexp'],vals['Ra'],vals['Dec']
            m5=m5_lim[b]+1.25*np.log10(vals['exptime']/30.)
            print >> outputfile_b,vals['band'],new_mjd,vals['exptime'],vals['rawSeeing'],seeing[b],vals['moon_frac'],sky[b],vals['kAtm'],1.2,m5,vals['Nexp'],vals['Ra'],vals['Dec']

    outputfile.close
    outputfile_b.close
