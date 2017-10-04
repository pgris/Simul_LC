import numpy as np
import cPickle as pkl
from Observations import *
import pylab as plt
import os
import glob


fieldname='DD'
#fieldids=[116,120,124,128,132]
#fieldids=[120,124,128]
fieldids=[744]
X1=-2.0
Color=0.2

thedir=dict(zip([120,744],['Mean_Obs_newrefs/','DD/']))


seasnum=0
season='Season_'+str(seasnum)

dict_val={}
obs={}

bands='rizy'
zmin=0.6
zmax=0.7
filtercolors = {'u':'b', 'g':'c', 'r':'g', 'i':'y', 'z':'r', 'y':'m'}
res=[]

for fieldid in fieldids:
    dirobs='../Ana_Cadence/OpSimLogs/'+thedir[fieldid]
    dirmeas='Prod_LC/'+thedir[fieldid]


    fichname=fieldname+'_'+str(fieldid)+'_X1_'+str(X1)+'_C_'+str(Color)+'_all.pkl'
    sumfile=dirmeas+str(fieldid)+'/'+season+'/'+fichname
    if os.path.exists(sumfile):
        pkl_file = open(sumfile,'rb') 
        dict_val[fieldid]=pkl.load(pkl_file)
    else:
        files = glob.glob(dirmeas+str(fieldid)+'/'+season+'/'+fieldname+'_'+str(fieldid)+'*_X1_'+str(X1)+'_C_'+str(Color)+'*.pkl')
        for fi in files:
            pkl_file = open(fi,'rb')
            print 'loading',fi
            if not fieldid in dict_val.keys():
                dict_val[fieldid]=pkl.load(pkl_file)
            else:
                dict_val[fieldid]=np.vstack((dict_val[fieldid],pkl.load(pkl_file)))

        pkl_out = open(sumfile,'wb')
                
        pkl.dump(dict_val[fieldid], pkl_out)
        
        pkl_out.close()
       

    name=dirobs+'/Observations_'+fieldname+'_'+str(fieldid)+'.txt'
    
    obs[fieldid]=Observations(fieldid=fieldid, filename=name).seasons[seasnum]

    """
    idx = (dict_val[fieldid]['z']>=zmin)&(dict_val[fieldid]['z']<zmax)
    sel = dict_val[fieldid][idx]

    print sel.dtype
    r=[]
    for b in bands:
        idxb = obs['band']=='LSSTPG::'+b
        selb=obs[idxb]
        cad=np.median([io-jo for jo,io in zip(selb['mjd'][:-1], selb['mjd'][1:])])
        print fieldid,b,np.median(sel['SNR_'+b]),cad
        res.append((b,cad,np.median(sel['SNR_'+b]),sel['date_obs']))
    """
#meas=np.rec.fromrecords(res,names=['band','cadence','SNR'])

figb, axb = plt.subplots(ncols=1, nrows=1, figsize=(10,9))

days=[5.,10.,15.,20.]

days=[1]
lls=dict(zip(fieldids,['-','--',':','-.']))


fontsize=12
tot_label=[]
for iv,daymax in enumerate(days):

    
    for pp,fieldid in enumerate(fieldids):
        meas=dict_val[fieldid]
        meas=meas[np.where(np.logical_and(meas['z']>=zmin,meas['z']<zmax))]
        
    #print meas.dtype
        observ=obs[fieldid]
    #print observ.dtype
        for b in bands:
            idx = (((meas['date_obs']-meas['T0'])>=daymax)&((meas['date_obs']-meas['T0'])<=daymax+0.1))
        
            meas_sel=np.sort(meas,order='date_obs')
        #print 'lll',meas_sel['date_obs']
            if pp==0:
                tot_label.append(axb.errorbar(meas_sel['date_obs'],meas_sel['SNR_'+b],color=filtercolors[b],ls=lls[fieldid],label=b+' band'))
            else:
                axb.errorbar(meas_sel['date_obs'],meas_sel['SNR_'+b],color=filtercolors[b],ls=lls[fieldid],label=b+' band') 
        if iv==0:
            idxb = observ['band']=='LSSTPG::'+b
            obs_sel=observ[idxb]
            axb.plot(obs_sel['mjd'],[8.-2.*pp]*len(obs_sel),'k*')
            #axb.plot(obs_sel['mjd'],[1000.-2.*pp]*len(obs_sel),'k*')


    axb.grid(which='both')
    xlim = axb.get_xlim()
    #axb.set_xlim([xlim[0],180.])
    axb.set_ylabel('SNR',{'fontsize': fontsize})
    axb.set_xlabel('Time [day]',{'fontsize': fontsize})
    labs = [l.get_label() for l in tot_label]
    axb.legend(tot_label, labs, ncol=5,loc='best',prop={'size':12},frameon=False)
    figb.suptitle('faint SN - z = '+str(zmin))

    #plt.plot(meas['date_obs'],meas['T0'],'bo')


plt.show()
