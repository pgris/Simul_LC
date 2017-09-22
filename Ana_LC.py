import numpy as np
import cPickle as pkl
import matplotlib.pyplot as plt
import glob

dirmeas='Prod_LC/Mean_Obs_newrefs_Obs'
X1=2.
Color=-0.2
zmin=0.95
zmax=0.96

fieldids=[120,121,122,123]

#fieldids=[120]

resu={}

for fieldid in fieldids:
    files = glob.glob(dirmeas+'/'+str(fieldid)+'/Season_0/DD_'+str(fieldid)+'_'+str(zmin)+'_'+str(zmax)+'_X1_'+str(X1)+'_C_'+str(Color)+'*.pkl')

    for fi in files:
        pkl_file = open(fi,'rb')
        print 'loading',fi
    #test=pkl.load(pkl_file)
    #print test,len(test)
        if not fieldid in resu.keys():
            resu[fieldid]=pkl.load(pkl_file)
        else:
            resu[fieldid]=np.concatenate((resu[fieldid],pkl.load(pkl_file)))
        
#band='z'
#print resu[120],resu[120].dtype

print 'blabla',resu.keys()
for fieldid in fieldids:
    
    figa, axa = plt.subplots(ncols=1, nrows=1, figsize=(10,9))
    figa.suptitle('Field '+str(fieldid))
    idx = resu[fieldid]['flux']/resu[fieldid]['fluxerr'] >=5.
    sel=resu[fieldid][idx]
    phase_min=np.min(sel['time'])/(1.+sel['z'][0])
    phase_max=np.max(sel['time'])/(1.+sel['z'][0])
    print sel,phase_min,phase_max
    for band in 'rizy':
        idxb = sel['band']=='LSST::'+band
        selb=sel[idxb]
        if len(selb) > 0:
            idxc = selb['time'] > 0.
            nbef=len(selb[idxc])
            idxc = selb['time'] < 0.
            naft=len(selb[idxc])
            print 'oh yes',fieldid,band,len(selb),nbef,naft
            axa.errorbar(selb['time'],selb['flux'],yerr=selb['fluxerr'],label=band)

        axa.legend(loc='best',prop={'size':12})
    
        
plt.show()
