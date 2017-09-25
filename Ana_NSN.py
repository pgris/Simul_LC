import cPickle as pkl
import numpy as np
import glob
import matplotlib.pyplot as plt

fieldname='DD'
fieldids=[290,744,1427,2412,2786]

tot_nsn=None
for fieldid in fieldids:
    name='N_SN_'+fieldname+'_'+str(fieldid)
    pkl_file = open(name+'.pkl','rb')
    if tot_nsn is None:
        tot_nsn=pkl.load(pkl_file)
    else:
        tot_nsn=np.vstack((tot_nsn,pkl.load(pkl_file)))

print tot_nsn

figa, axa = plt.subplots(ncols=1, nrows=1, figsize=(10,9))
colors = dict(zip(fieldids,['b','k','g','r','m']))

for fieldid in fieldids:
    idx = tot_nsn['fieldid']==fieldid
    sel=tot_nsn[idx]
    nsn=[]
    err=[]
    for i in range(10):
        iid = sel['season']==i
        nsn.append(sel['n_sn_detected'][iid])
        err.append(sel['err_detected'][iid])
    ll='Field '+str(fieldid)
    axa.errorbar([i+1 for i in range(10)],nsn,yerr=err,label=ll,color=colors[fieldid])
minval=220.
interv=15
ntot_nsn=0.
errtot_nsn=0.
for i in range(10):
    idf = tot_nsn['season']==i
    selb= tot_nsn[idf]
    print i, np.sum(selb['n_sn_detected']),np.sqrt(np.sum(np.power(selb['err_detected'],2.)))
    n_sn=int(np.sum(selb['n_sn_detected']))
    err_sn=int(np.sqrt(np.sum(np.power(selb['err_detected'],2.))))
    ntot_nsn+=np.sum(selb['n_sn_detected'])
    errtot_nsn+=np.sum(np.power(selb['err_detected'],2.))
    if i < 5:
        axa.text(1,minval-(i%5)*interv,'Year '+str(i+1)+' : $\mathrm{N_{SN\/ Ia}}$ = '+str(n_sn)+'$\pm$'+str(err_sn))
    else:
       axa.text(5,minval-(i%5)*interv,'Year '+str(i+1)+' : $\mathrm{N_{SN \/Ia}}$ = '+str(n_sn)+'$\pm$'+str(err_sn)) 

fontsize=12
axa.set_xlabel('Year',{'fontsize': fontsize})
axa.set_ylabel('Number of Type Ia Supernovae',{'fontsize': fontsize})
axa.set_xlim([0.5,10.5])
axa.legend(loc='best',prop={'size':fontsize})

figa.suptitle('Perret rate - 10 years - $\mathrm{N_{SN\/Ia}}$ = '+str(int(ntot_nsn))+' $\pm$ '+str(int(np.sqrt(errtot_nsn)))) 
plt.show()
