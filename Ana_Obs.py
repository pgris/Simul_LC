import numpy as np
from Observations import *
import matplotlib.pyplot as plt

bands='ugrizy'

def Ana_Obs(dirmeas,fieldname,fieldid,season):

    dictout={}
    dictout['m5sigmadepth']={}
    dictout['seeing']={}

    filename='../Ana_Cadence/OpSimLogs/'+dirmeas+'/Observations_'+fieldname+'_'+str(fieldid)+'.txt'

    myobs=Observations(fieldid=fieldid, filename=filename)
    obs=myobs.seasons[season]
    for b in bands:
        idx = obs['band'] == 'LSSTPG::'+b
        sel=obs[idx]
        print season+1,b,np.median(sel['m5sigmadepth'])
        dictout['m5sigmadepth'][b]=np.median(sel['m5sigmadepth'])
        dictout['seeing'][b]=np.median(sel['seeing'])
    dictout['airmass']=np.median(obs['airmass'])
    
    return dictout

def All_Seasons(dirmeas,fieldname,fieldid):

    restot=[]
    names=[]
    for season in range(0,10):
        res=[]
        resu=Ana_Obs(dirmeas,fieldname,fieldid,season)
        for b in bands:
            res.append(resu['m5sigmadepth'][b])
            names.append('m5_'+b)
        res.append(resu['airmass'])
        names.append('airmass')
        res.append(season+1)
        names.append('season')
        #print 'hello',res
        restot.append(tuple(res))

#print restot

    return np.rec.fromrecords(restot,names=names[:8])


def Plot(dirmeas,fieldname,fieldid,title):
    
    bcolors = {'u':'b', 'g':'c', 'r':'g', 'i':'y', 'z':'r', 'y':'m'}
    m5_lim={'u':23.61,'g':24.83,'r':24.35,'i':23.88,'z':23.30,'y':22.43}

#fin=All_Seasons('WFD_Rolling_noTwilight','WFD',309)
    fin=All_Seasons(dirmeas,fieldname,fieldid)

    figa, axa = plt.subplots(ncols=1, nrows=2, figsize=(10,9))
    for b in bands:
        ll=b+' band'
        axa[0].plot(fin['season'],fin['m5_'+b],bcolors[b],label=ll)
        refmag=[m5_lim[b]]*10
        axa[0].plot([i for i in range(1,11)],refmag,bcolors[b],ls='--')
        axa[1].plot([i for i in range(1,11)],refmag-fin['m5_'+b],bcolors[b],label=ll)
        print 'hello',b,np.median(fin['m5_'+b]),m5_lim[b]-np.median(fin['m5_'+b])

    axa[0].set_xlim([-1,10])
    axa[0].set_xlabel('Season')
    axa[0].set_ylabel('median m$_{5}$ [mag]')
    axa[0].legend(loc='center left',prop={'size':12})
    
    axa[1].set_xlim([-1,10])
    axa[1].set_xlabel('Season')
    axa[1].set_ylabel('m$_{5}$$^{ref}$ - median m$_{5}$ [mag]')
    axa[1].legend(loc='center left',prop={'size':12})
    
    figa.suptitle(title)

dictplot={}
#dictplot['WFD']=('WFD',309,'WFD - '+str(309))
#dictplot['WFD_Rolling_noTwilight']=('WFD',309,'WFD Rolling - '+str(309))
dictplot['DD']=('DD',290,'DD - '+str(290))

for key,vals in dictplot.items():
    Plot(key,vals[0],vals[1],vals[2])

plt.show()
