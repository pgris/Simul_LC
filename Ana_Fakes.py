import numpy as np
from Observations import *
import pylab as plt

def Plot_Obs_per_Field(resu):
    fontsize=12.
    
    myls=['-','--']
    colors=dict(zip([i for i in range(10)],['k','k','r','r','b','b','g','g','m','m']))
    
    for fieldid in fieldids:
        figa, axa = plt.subplots(ncols=1, nrows=3, figsize=(10,9))
        idx = resu['fieldid']==fieldid
        figa.suptitle(fieldname+' field - '+str(fieldid))
        tot_label=[]
        sela=resu[idx]
        m5_var={}
        for season in range(10):
            idxb = sela['season']==season
            selb=sela[idxb]
            x=[]
            y=[]
            z=[]
            w=[]
            
            for b in bands:
                if not m5_var.has_key(b):
                    m5_var[b]=[]
                idxc = selb['band']==b
                selc=selb[idxc]
                x.append(selc['ib'][0])
                y.append(selc['airmass'][0])
                z.append(selc['m5'][0])
                w.append(selc['seeing'][0])
                m5_var[b].append(selc['m5'][0])
            print fieldid,season+1,z
            axa[0].plot(x,y,ls=myls[season%2],color=colors[season])
            
            ll='Y'+str(season+1)
            tot_label.append(axa[1].errorbar(x,z,ls=myls[season%2],color=colors[season],label=ll))
                
            axa[2].plot(x,w,ls=myls[season%2],color=colors[season])
        

        axa[0].set_ylabel('Median airmass',{'fontsize': fontsize})
        axa[1].set_ylabel('Median m5 [mag]',{'fontsize': fontsize})
        
        labs = [l.get_label() for l in tot_label]
        axa[1].legend(tot_label, labs, ncol=5,loc='best',prop={'size':12},frameon=False)
        
        axa[2].set_ylabel('Median seeing [\'\']',{'fontsize': fontsize}) 

        print fieldid,'m5 variations'
        for b in bands:
            print b,np.max(m5_var[b])-np.min(m5_var[b])
        
        
        for j in range(3):
            axa[j].set_xlabel('band',{'fontsize': fontsize})
            axa[j].set_xlim([-0.1,4.1])
            axa[j].set_xticks([i for i in range(len(bands))])
            axa[j].set_xticklabels([corresp_inverted[i] for i in range(len(bands))])
        plt.gcf().savefig('Obs_Plots/'+fieldname+'_'+str(fieldid)+'.png')

def Plot_Cadence_per_Field(resu):
    fontsize=12.
    
    myls=['-','--']
    colors=dict(zip([i for i in range(10)],['k','k','r','r','b','b','g','g','m','m']))
    
    for fieldid in fieldids:
        figa, axa = plt.subplots(ncols=1, nrows=3, figsize=(10,9))
        idx = resu['fieldid']==fieldid
        figa.suptitle(fieldname+' field - '+str(fieldid))
        tot_label=[]
        sela=resu[idx]
        for season in range(10):
            idxb = sela['season']==season
            selb=sela[idxb]
            x=[]
            y=[]
            z=[]
            w=[]
            for b in bands:
                idxc = selb['band']==b
                selc=selb[idxc]
                x.append(selc['ib'][0])
                y.append(selc['mean_cadence'][0])
                z.append(selc['rms_cadence'][0])
                w.append(selc['duration'][0])
                print fieldid,season,x,y
            axa[0].plot(x,y,ls=myls[season%2],color=colors[season])
            
            ll='Y'+str(season+1)
            tot_label.append(axa[1].errorbar(x,z,ls=myls[season%2],color=colors[season],label=ll))
                
            axa[2].plot(x,w,ls=myls[season%2],color=colors[season])
        

        axa[0].set_ylabel('Mean cadence [day]',{'fontsize': fontsize})
        axa[1].set_ylabel('RMS cadence [day]',{'fontsize': fontsize})
        
        labs = [l.get_label() for l in tot_label]
        axa[1].legend(tot_label, labs, ncol=5,loc='best',prop={'size':12},frameon=False)
        
        #axa[2].set_ylabel('Duration [days]',{'fontsize': fontsize}) 
        axa[2].set_ylabel('Observation Period [day]',{'fontsize': fontsize}) 
        
        for j in range(3):
            axa[j].set_xlabel('band',{'fontsize': fontsize})
            axa[j].set_xlim([-0.1,4.1])
            axa[j].set_xticks([i for i in range(len(bands))])
            axa[j].set_xticklabels([corresp_inverted[i] for i in range(len(bands))])
        plt.gcf().savefig('Cadence_Plots/'+fieldname+'_'+str(fieldid)+'.png')

def Plot_per_Field(resu,what=('mean_cadence','rms_cadence'),myleg=('Mean cadence [day]','RMS cadence [day]')):
    fontsize=12.
    
    myls=['-','--']
    colors=dict(zip([i for i in range(10)],['k','k','r','r','b','b','g','g','m','m']))
    
    for fieldid in fieldids:
        figa, axa = plt.subplots(ncols=1, nrows=2, figsize=(10,9))
        idx = resu['fieldid']==fieldid
        figa.suptitle(fieldname+' field - '+str(fieldid))
        tot_label=[]
        sela=resu[idx]
        for season in range(10):
            idxb = sela['season']==season
            selb=sela[idxb]
            
            x=[]
            y=[]
            z=[]
            w=[]
            for b in bands:
                idxc = selb['band']==b
                selc=selb[idxc]
                x.append(selc['ib'][0])
                y.append(selc[what[0]][0])
                z.append(selc[what[1]][0])
                
            
            print fieldid,season,what[0],y,what[1],z,np.max(selb[what[0]])-np.min(selb[what[0]])
            axa[0].plot(x,y,ls=myls[season%2],color=colors[season])
            
            ll='Y'+str(season+1)
            tot_label.append(axa[1].errorbar(x,z,ls=myls[season%2],color=colors[season],label=ll))
                
            
        idxc = sela['band']=='a'
        print fieldid,season,sela['duration'][idxc]
        
        ax2 = axa[0].twiny()
        ax2.plot(sela['season'][idxc]+1,sela['duration'][idxc],ls='-',color='k',marker='s',label='grizy')
        
        axa[0].set_ylabel(myleg[0],{'fontsize': fontsize})
        axa[1].set_ylabel(myleg[1],{'fontsize': fontsize})
        ax2.set_xlabel('Year',{'fontsize': fontsize})

        labs = [l.get_label() for l in tot_label]
        axa[1].legend(tot_label, labs, ncol=5,loc='best',prop={'size':12},frameon=False)
        ax2.legend(loc='best',prop={'size':12})
        for j in range(2):
            axa[j].set_xlabel('band',{'fontsize': fontsize})
            axa[j].set_xlim([-0.1,4.1])
            axa[j].set_xticks([i for i in range(len(bands))])
            axa[j].set_xticklabels([corresp_inverted[i] for i in range(len(bands))])

        
        #plt.gcf().savefig('Cadence_Plots/'+fieldname+'_'+str(fieldid)+'.png')

def Plot_Cadence_per_Year(resu):

    #myls=['-','--']
    colors=dict(zip(fieldids,['k','r','b','g','m']))
    fontsize=12

    for season in range(10):
        idx = resu['season']==season
        sela=resu[idx]
        figa, axa = plt.subplots(ncols=1, nrows=3, figsize=(10,9))
        figa.suptitle(fieldname+' - Year '+str(season+1))
        tot_label=[]
        for fieldid in fieldids:
            idxb = sela['fieldid'] == fieldid
            selb=sela[idxb]
            x=[]
            y=[]
            z=[]
            w=[]
            for b in bands:
                idxc = selb['band']==b
                selc=selb[idxc]
                x.append(selc['ib'][0])
                y.append(selc['mean_cadence'][0])
                z.append(selc['rms_cadence'][0])
                w.append(selc['duration'][0])
                #print fieldid,season,x,y
            axa[0].plot(x,y,ls='-',color=colors[fieldid])
            axa[0].set_xlabel('band',{'fontsize': fontsize})
            axa[0].set_ylabel('Mean cadence [day]',{'fontsize': fontsize})
            axa[0].set_xlim([-0.1,4.1])
            axa[0].set_xticks([i for i in range(len(bands))])
            axa[0].set_xticklabels([corresp_inverted[i] for i in range(len(bands))])
            ll=str(fieldid)
            tot_label.append(axa[1].errorbar(x,z,ls='-',color=colors[fieldid],label=ll))
                
            axa[2].plot(x,w,ls='-',color=colors[fieldid])


        axa[0].set_ylabel('Mean cadence [day]',{'fontsize': fontsize})
        axa[1].set_ylabel('RMS cadence [day]',{'fontsize': fontsize})
        
        labs = [l.get_label() for l in tot_label]
        axa[1].legend(tot_label, labs, ncol=5,loc='best',prop={'size':12},frameon=False)
        
        #axa[2].set_ylabel('Duration [days]',{'fontsize': fontsize}) 
        axa[2].set_ylabel('Observation period [day]',{'fontsize': fontsize})
        
        for j in range(3):
            axa[j].set_xlabel('band',{'fontsize': fontsize})
            axa[j].set_xlim([-0.1,4.1])
            axa[j].set_xticks([i for i in range(len(bands))])
            axa[j].set_xticklabels([corresp_inverted[i] for i in range(len(bands))])
        
        plt.gcf().savefig('Cadence_Plots/'+fieldname+'_Y'+str(season+1)+'.png')


def Plot_Nobs(Nobs):

    for season in range(1):
        Plot_Nobs_Indiv(Nobs[Nobs['season']==season],season)

def Plot_Nobs_Indiv(Nobs,season):
 
   for band in ['g','r','i','z','y','all']:    
        figb, axb = plt.subplots(ncols=1, nrows=2, figsize=(10,9))
            #figb.suptitle(band+' band')
        figb.suptitle(fieldname+' - '+str(fieldid)+' - Year '+str(season+1)+' '+band+' band') 
        idx = Nobs['band']==band
            #sel=Nobs[np.where(np.logical_and(Nobs['T0']>=62.,Nobs['T0']<63.))]
        sel=Nobs[idx]
        ll=''
        if band == 'all':
            ll='grizy'
            
        axb[0].plot(sel['T0'],sel['Nbef'],label=ll)
        axb[1].plot(sel['T0'],sel['Naft'],label=ll)

        if band == 'all':
            idx = Nobs['band']==band+'_no_g'
            #sel=Nobs[np.where(np.logical_and(Nobs['T0']>=62.,Nobs['T0']<63.))]
            sel=Nobs[idx]
            ll='rizy'
            axb[0].plot(sel['T0'],sel['Nbef'],label=ll,color='r')
            axb[1].plot(sel['T0'],sel['Naft'],label=ll,color='r')
            axb[0].plot(sel['T0'],[4.]*len(sel['T0']),color='k')
            axb[1].plot(sel['T0'],[10.]*len(sel['T0']),color='k')

        axb[0].set_xlabel('T0 [day]',{'fontsize': fontsize})
        axb[0].set_ylabel('Nobs in [T0-20, T0]',{'fontsize': fontsize})
        axb[1].set_xlabel('T0 [day]',{'fontsize': fontsize})
        axb[1].set_ylabel('Nobs in [T0, T0+40]',{'fontsize': fontsize})

        if band =='all':
            axb[0].legend(loc='upper left',prop={'size':fontsize})
            axb[1].legend(loc='upper left',prop={'size':fontsize})

def Plot_Diffs(diffs,fieldname,fieldids,colors):
    
        fontsize=12
        r=[]
        for fieldid in fieldids:
            for key,diff in diffs[fieldid].items():
            #print diff.keys()
                figa, axa = plt.subplots(ncols=1, nrows=1, figsize=(10,9))
                figa.suptitle(fieldname+' - '+str(fieldid)+' - Year '+str(key+1))
                for band in 'grizy':
                    ll=band+' band'
                    axa.errorbar([val for val in diff[band][0]],[val for val in diff[band][1]],color=filtercolors[band],marker='o',label=ll)
                    axa.errorbar([val for val in diff[band][0]],[np.median([val for val in diff[band][1]])]*len(diff[band][0]),color=filtercolors[band],ls='--')
                    r.append((fieldid,key+1,band,np.median([val for val in diff[band][1]])))
                axa.legend(loc='best',prop={'size':fontsize})
                axa.set_ylabel('$\Delta$T = T$_{obs}$-T$_{obs-1}$ [day]',{'fontsize': fontsize})
                axa.set_xlabel('MJD [day]',{'fontsize': fontsize})
                plt.gcf().savefig('Obs_Plots/DeltaT_'+fieldname+'_'+str(fieldid)+'_Y'+str(key+1)+'.png')
                plt.close(figa)
        med_diffs=np.rec.fromrecords(r,names=['fieldid','season','band','median_diff'])

        myband='z'
        figb, axb = plt.subplots(ncols=1, nrows=1, figsize=(10,9))
        figb.suptitle(myband+' band')
        for fieldid in fieldids:
            idx=(med_diffs['fieldid']==fieldid)&(med_diffs['band']==myband)
            sela=med_diffs[idx]
            ll='Field '+str(fieldid)
            axb.plot(sela['season'],sela['median_diff'],color=colors[fieldid],label=ll)
        axb.set_xlabel('Year',{'fontsize': fontsize})
        axb.set_ylabel('Median $\Delta$T [day]',{'fontsize': fontsize})
        axb.legend(loc='best',prop={'size':fontsize})
        ylim = axb.get_ylim()
        axb.set_ylim([ylim[0]-0.1,ylim[1]])
        xlim = axb.get_xlim()
        axb.set_xlim([xlim[0]-0.1,xlim[1]+0.1])
        major_ticks = np.arange(ylim[0],ylim[1], 1)
        axb.set_yticks(major_ticks)
        major_ticks = np.arange(xlim[0],xlim[1]+1,1)
        axb.set_xticks(major_ticks)
        axb.grid(which='both')

def Plot_median(res,fieldids,fieldcolors):

    figa, axa = plt.subplots(ncols=1, nrows=2, figsize=(10,9))
    fontsize=12

    for fieldid in fieldids:
        idx = res['fieldid']==fieldid
        sela=res[idx]
        ras=[]
        for band in 'grizy':
            idxb = sela['band']==band
            selb=sela[idxb]
            print 'hello',fieldid,selb
            ras.append((np.mean(selb['ib']),np.median(selb['duration']),np.median(selb['obstime'])))
            
        ll='Fieldid '+str(fieldid)
        axa[0].plot([vv[0] for vv in ras],[vv[1] for vv in ras],color=fieldcolors[fieldid],label=ll)
        axa[1].plot([vv[0] for vv in ras],[vv[2] for vv in ras],color=fieldcolors[fieldid],label=ll)


    #axa[0].set_ylabel('Median duration [day] ',{'fontsize': fontsize})
    axa[0].set_ylabel('Median observation period [day] ',{'fontsize': fontsize})
    axa[1].set_ylabel('Median expTime [s]',{'fontsize': fontsize})
    axa[0].legend(loc='best',prop={'size':fontsize},frameon=False)
    axa[1].legend(loc='best',prop={'size':fontsize},frameon=False)
    for j in range(2):
        axa[j].set_xlabel('band',{'fontsize': fontsize})
        axa[j].set_xlim([-0.1,4.1])
        axa[j].set_xticks([i for i in range(len(bands))])
        axa[j].set_xticklabels([corresp_inverted[i] for i in range(len(bands))])

dirmeas='Mean_Obs_newrefs'
dirmeas='DD'
fieldname='DD'

thedir='../Ana_Cadence/OpSimLogs/'+dirmeas

myobs={}

#fieldids=[120,121,122,123]
fieldids=[124,125,126,127]
fieldids=[128,129,130,131]
fieldids=[field+4 for field in fieldids]
print 'alors',fieldids
#fieldids=[123]
#fieldids=[123]
fieldids=[290,744,1427,2412,2786]
#fieldids=[744,1427]
#fieldids=[120]
for fieldid in fieldids:
    name='Observations_'+fieldname+'_'+str(fieldid)+'.txt'
    myobs[fieldid]=Observations(fieldid=fieldid, filename=thedir+'/'+name)

r=[]
bands='grizy'
corresp=dict(zip(bands+'a',[i for i in range(len(bands)+1)]))
corresp_inverted=dict(zip([i for i in range(len(bands)+1)],bands+'a'))

filtercolors = {'u':'c', 'g':'b', 'r':'g', 'i':'y', 'z':'r', 'y':'m'}
fieldcolors=dict(zip([290,744,1427,2412,2786],'bgyrm'))
ra=[]
all_diff={}
for key, vals in myobs.items():
    print key,len(vals.seasons)
    iseason=-1
    all_diff[key]={}
    for season in range(len(vals.seasons)):
        
        myseason=vals.seasons[season]
        all_diff[key][season]={}
        full_season=myseason.copy()
        full_season.sort(order='mjd')
        idx = full_season['band'] != 'LSSTPG::u'
        full_season=full_season[idx]
        print full_season[full_season['band'] == 'LSSTPG::i']
        min_season=np.min(full_season['mjd'])
        max_season=np.max(full_season['mjd'])
        
        for val in np.arange(min_season,max_season,1.):
        #for val in np.arange(62.,63.,1.):
            idxa = np.logical_and(full_season['mjd']> val -20. ,full_season['mjd'] < val)
            idxb = np.logical_and(full_season['mjd']> val,full_season['mjd'] < val+40.)
            sela=full_season[idxa]
            selb=full_season[idxb]
            ra.append((iseason+1,'all',val,len(sela),len(selb),len(sela)+len(selb)))
            ppa=sela[sela['band'] != 'LSSTPG::g']
            ppb=selb[selb['band'] != 'LSSTPG::g']
            ra.append((iseason+1,'all_no_g',val,len(ppa),len(ppb),len(ppa)+len(ppb)))
            print val,len(sela),len(selb)
            for b in 'grizy':
                selbf=sela[np.where(sela['band']=='LSSTPG::'+b)]
                selaft=selb[np.where(selb['band']=='LSSTPG::'+b)]
                print 'hello',b,selbf
                ra.append((iseason+1,b,val,len(selbf),len(selaft),len(selbf)+len(selaft)))
       
        for b in 'grizy':
            idx = myseason['band']=='LSSTPG::'+b
            sel = myseason[idx]
            
            sel.sort(order='mjd')
            """
            if len(sel) >=2:
                if b == 'g':
                    iseason+=1
            """ 
            if b == 'g':
                iseason+=1
            diff=[io-jo for jo,io in zip(sel['mjd'][:-1], sel['mjd'][1:])]
            all_diff[key][season][b]=(sel['mjd'][1:],diff)
            themin=np.min(diff)
            themax=np.max(diff)
            """
            #plt.hist(diff,range=[int(themin),int(themax)],bins=int(themax)-int(themin))
            plt.plot(sel['mjd'][1:],diff,'bo')
            plt.show()
            """
                #print key, np.mean(diff),np.std(diff),np.max(sel['mjd'])-np.min(sel['mjd'])
            r.append((key,iseason,b,np.mean(diff),np.std(diff),np.max(sel['mjd'])-np.min(sel['mjd']),corresp[b],np.median(sel['airmass']),np.median(sel['m5sigmadepth']),np.median(sel['seeing']),np.sum(sel['exptime'])))

        idxc = myseason['band']!='LSSTPG::u'
        selcb=myseason[idxc]
            
        selcb.sort(order='mjd')
        diff=[io-jo for jo,io in zip(selcb['mjd'][:-1], selcb['mjd'][1:])]
        r.append((key,iseason,'a',np.mean(diff),np.std(diff),np.max(selcb['mjd'])-np.min(selcb['mjd']),corresp['a'],np.median(selcb['airmass']),np.median(selcb['m5sigmadepth']),np.median(selcb['seeing']),np.sum(selcb['exptime'])))

        
        """
        plt.hist(sel['mjd'],bins=int(np.max(sel['mjd']))-int(np.min(sel['mjd'])))

        plt.show()
        """

resu=np.rec.fromrecords(r,names=['fieldid','season','band','mean_cadence','rms_cadence','duration','ib','airmass','m5','seeing','obstime'])
Nobs=np.rec.fromrecords(ra,names=['season','band','T0','Nbef','Naft','Nmeas'])

print resu
iseason=0
idx= resu['season'] == iseason
sela=resu[idx]
for fieldid in fieldids:
    idxb=sela['fieldid']==fieldid
    selb=sela[idxb]
    for band in 'grizy':
        idxc = selb['band']==band
        selc = selb[idxc]
        print fieldid, band, selc['mean_cadence'],selc['rms_cadence'],selc['duration']

#Plot_per_Field(resu)
#Plot_per_Field(resu,what=('duration','obstime'),myleg=('Duration [day]','Observing Time [s]'))
Plot_Diffs(all_diff,fieldname,fieldids,fieldcolors)
#Plot_Obs_per_Field(resu)
#Plot_Cadence_per_Year(resu)
#Plot_Nobs(Nobs)
#Plot_median(resu,fieldids,fieldcolors)
plt.show()
