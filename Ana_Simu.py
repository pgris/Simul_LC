import cPickle as pkl
import numpy as np
import glob
import matplotlib.pyplot as plt
from SN_Rate import *
from scipy import interpolate
from scipy.interpolate import UnivariateSpline
from scipy.integrate import quad
import os
import collections
from optparse import OptionParser

class Id:
    def __init__(self,thedir,fieldname,fieldid,X1,Color,season,colorfig):
        self.thedir=thedir
        self.fieldname=fieldname
        self.fieldid=fieldid
        self.X1=X1
        self.Color=Color
        self.season=season
        self.colorfig=colorfig

    @property
    def thedir(self):
        return self.thedir
    @property
    def fieldname(self):
        return self.fieldname
    @property
    def fieldid(self):
        return self.fieldid
    @property
    def X1(self):
        return self.X1
    @property
    def Color(self):
        return self.Color
    @property
    def season(self):
        return self.season
    @property
    def colorfig(self):
        return self.colorfig

class Ana_Simu:
    def __init__(self, dict_ana,zmin=0.,zmax=1.2,bin_z=0.01):

        thedir='Prod_LC'

        tot_resu={}

        self.zmin=zmin
        self.zmax=zmax
        self.bin_z=bin_z

        for key,val in dict_ana.items():
            #key=vals[0]
            #val=vals[1]
            print 'hee',val.thedir
            dirmeas=thedir+'/'+val.thedir+'/'+str(val.fieldid)+'/Season_'+str(val.season)
            sum_file=dirmeas+'/'+val.fieldname+'_'+str(val.fieldid)+'_X1_'+str(val.X1)+'_C_'+str(val.Color)+'_all.pkl'

            if os.path.exists(sum_file):
                pkl_file = open(sum_file,'rb')
                tot_resu[key]=pkl.load(pkl_file)
            else:

                files = glob.glob(dirmeas+'/'+val.fieldname+'_'+str(val.fieldid)+'*_X1_'+str(val.X1)+'_C_'+str(val.Color)+'*.pkl')
                print 'hello',dirmeas
                for fi in files:
                    pkl_file = open(fi,'rb')
                    print 'loading',fi
                    if not key in tot_resu.keys():
                        tot_resu[key]=pkl.load(pkl_file)
                    else:
                        tot_resu[key]=np.vstack((tot_resu[key],pkl.load(pkl_file)))

                pkl_out = open(sum_file,'wb')
                
                pkl.dump(tot_resu[key], pkl_out)
                
                pkl_out.close()
       

            print 'there',key,len(tot_resu[key]),tot_resu[key].dtype
         
        # this is for simulated parameters
    
        figb, axb = plt.subplots(ncols=2, nrows=2, figsize=(10,9))
        for key in dict_ana.keys():
            self.Plot_Sims(axb,tot_resu[key])
        

        self.nsn_tot=0.
        self.err_tot=0.
        self.ms=['o','o','s','s','.','.','^','^','<','<']
        self.color=['b','r','b','r','b','r','b','r','b','r']
        figa, axa = plt.subplots(ncols=1, nrows=1, figsize=(10,9))
        axca = axa.twinx()
        figc, axc = plt.subplots(ncols=1, nrows=1, figsize=(10,9))
        title=val.fieldname+' - '+str(val.fieldid)
        

        #figbb, axbb = plt.subplots(ncols=2, nrows=2, figsize=(10,9))
        col=['k','r','b']
        idraw=-1
        tot_resu_o=collections.OrderedDict(sorted(tot_resu.items()))
        
        for key, val in tot_resu_o.items():
            #self.Plot_Sims(axbb,val)
            """
            for band in 'y':
                fig, ax = plt.subplots(ncols=1, nrows=1, figsize=(10,9))
                self.Plot(ax,val,'N_bef_'+band,'z',col[0])
            """
            
            #idx=val['status']=='go_fit'
            #print np.unique(val[idx]['fit_status'])

            idraw+=1
            sel=val.copy()
            selb=val.copy()
            selc=val.copy()
            sel=sel[np.where(np.logical_and(sel['N_bef']>=4,sel['N_aft']>=10))]
            #sel=val[np.where(val['N_bef']>=4)]
            #sel=val[np.where(val['N_aft']>=10)]
            #sel=sel[np.where(sel['status']=='go_fit')]
            sel=sel[np.where(np.logical_and(sel['status']=='go_fit',sel['fit_status']=='fit_ok'))]
            sel=sel[np.where(np.logical_and(sel['phase_first']<=-5,sel['phase_last']>=20))]
            
            selb=selb[np.where(selb['status']=='no_obs')]
            #selc=selc[np.where(selc['fit_status']=='crashd')]
            #selc=selc[np.where(selc['status']=='unknown')]
            #selc=selc[np.where(np.logical_or(selc['phase_first']>=-5,selc['phase_last']<=20))]
            #selc=selc[np.where(np.logical_and(selc['phase_first']>=-500.,selc['phase_last']>=-500.))]
            selc=selc[np.where(np.logical_and(selc['status']=='go_fit',selc['fit_status']=='fit_ok'))]
            selc=selc[np.where(np.logical_or(selc['phase_first']>-5,selc['phase_last']<20))]
            
            """
        #print sel['phase_first'],sel['phase_last']
        for i,val in enumerate([(0.,0.),(2.0,-0.2),(-2.0,0.2)]):
        #for i,val in enumerate([(-2.0,0.2)]):    
            sela=tot_resu[np.where(np.logical_and(tot_resu['X1']==val[0],tot_resu['Color']==val[1]))]
            selc=sel[np.where(np.logical_and(sel['X1']==val[0],sel['Color']==val[1]))]
            """
            #Efficiency per season vs z plot
            self.Plot_Eff(axa,val,sel,'z',dict_ana[key].colorfig,dict_ana[key].season,key,0)
            self.Plot_Eff(axca,val,selb,'z',dict_ana[key].colorfig,dict_ana[key].season,key,1)
            self.Plot_Eff(axca,val,selc,'z',dict_ana[key].colorfig,dict_ana[key].season,key,2)
            # number of supernovae per season
            self.Plot_N_SN(axc,val,sel,'z',dict_ana[key].colorfig,dict_ana[key].season,key,idraw)

            #print 'boouh',key,dict_ana[key].colorfig,dict_ana[key].season

            #self.Plot(axc,sel,'salt2.CovColorColor','z',dict_ana[key].colorfig)
            #print 'Number of Events',val[0],val[1],len(val)
            for vval in np.arange(self.zmin,self.zmax,0.1):
                ssel=val[np.where(np.logical_and(val['z']>=vval,val['z']<vval+0.1))]
                #print vval,vval+0.1,len(ssel)

        #print val['T0']
        print 'in total :',self.nsn_tot,'+-',np.power(self.err_tot,0.5)
        axa.set_title(title)
        #axa.set_xlim(self.zmin,self.zmax+0.01)
        title+=' - N$_{SN Ia}$ ='+str(int(self.nsn_tot))+'$\pm$'+str(int(np.power(self.err_tot,0.5)))
        axc.set_title(title)
        #axc.set_xlim(self.zmin,self.zmax+0.01)

        self.Plot_m5(tot_resu)
        #self.Plot_Cadences(tot_resu)
        #self.Plot_Nobs(tot_resu)
        #self.Plot_Phases(tot_resu)
        #self.Plot_Color(tot_resu)
        self.Plot_Effi_vs_Cadence(tot_resu)
        plt.show()

    def Plot_Sims(self,axbb,tab_resu):
        
        
        fontsize=10

        for (j,vals) in enumerate(['T0','z','X1','Color']):
            if j==0 or j==2:
                k=0
            else:
                k=1

            axbb[j/2][k].hist(tab_resu[vals],bins=10,histtype='step')
            
        
            axbb[j/2][k].set_xlabel(r''+vals+'_sim',{'fontsize': fontsize})
            axbb[j/2][k].set_ylabel(r'Number of Entries',{'fontsize': fontsize})
            print vals,np.mean(tab_resu[vals]),np.std(tab_resu[vals])

    def Histo_ratio(self,sela,selb,varname):

        range=[self.zmin,self.zmax]
        bin_width=self.bin_z
        num_bins=int((range[1]-range[0])/bin_width)
        #range=[0.0,1.2]
        
        hista, bin_edgesa = np.histogram(sela[varname],bins=num_bins,range=range)
        histb, bin_edgesb = np.histogram(selb[varname],bins=num_bins,range=range)
        bin_center = (bin_edgesa[:-1] + bin_edgesa[1:]) / 2

        ratio=[]
        ratio_err=[]
        norm=[]
        norm_err=[]
        
        for a,b in zip(histb,hista):
            if b==0:
                ratio.append(a)
                ratio_err.append(0)
            else:
                effi=float(a)/float(b)
                ratio.append(effi)
                ratio_err.append(np.sqrt(1.-effi)*effi/np.sqrt(float(b)))
            eff=float(a)/float(np.sum(hista))
            norm.append(eff)
            norm_err.append(np.sqrt(1.-eff)*eff/np.sqrt(float(np.sum(hista))))
    
        return bin_center, ratio, ratio_err,norm,norm_err

    def Plot_Eff(self,axc,sela,selb,varname,color,season,ll,idraw):
        tot_label=[]
        #figc, axc = plt.subplots(ncols=1, nrows=1, figsize=(10,9))
        self.Plot_Eff_Indiv(axc,sela,selb,varname,tot_label,ll,color,season,idraw)
        

    def Plot_Eff_Indiv(self,axc,sela,selb,varname,tot_label,ll,color,season,idraw):

        bin_center, ratio, ratio_err,norm,norm_err= self.Histo_ratio(sela,selb,varname)
        #tot_label.append(axc.errorbar(bin_center,ratio, yerr=ratio_err,marker=marker, mfc=colors[key], mec=colors[key], ms=8, linestyle=myfmt[i],color='k',label=ll))
        ll='Y'+str(season+1)
        
       
        axc.set_xlabel('z')
        if idraw == 0:
            axc.errorbar(bin_center,ratio, yerr=ratio_err,marker=self.ms[season], mfc=self.color[season], mec=self.color[season], ms=8, linestyle='-',color=color,label=ll)
            axc.set_ylabel('Efficiency')
            axc.legend(loc='best',prop={'size':12})
            axc.set_xlim(self.zmin,np.max(bin_center)+0.01)
        else:
            if idraw == 1:
                myls='--'
            else:
                myls = ':'
            axc.errorbar(bin_center,ratio, yerr=ratio_err,marker=self.ms[season], mfc='k', mec='k', ms=8, linestyle=myls,color='k',label=ll)
            axc.set_ylabel('Fraction of Events') 


    def Plot_N_SN(self,axb,sela,selb,varname,color,season,ll,idraw):
        tot_label=[]
        #figc, axc = plt.subplots(ncols=1, nrows=1, figsize=(10,9))
        self.Plot_Eff_Indiv_N_SN(axb,sela,selb,varname,tot_label,ll,color,season,idraw)

    def Plot_Eff_Indiv_N_SN(self,axb,sela,selb,varname,tot_label,ll,color,season,idraw):

        bin_center, ratio, ratio_err,norm,norm_err= self.Histo_ratio(sela,selb,varname)

        effi = interpolate.interp1d(bin_center,ratio)

        rate_name='Perret'
        sn_rate=SN_Rate(rate=rate_name)

        #zz,rate,err_rate,nsn,err_nsn=sn_rate(self.zmin,self.zmax,self.bin_z)
        zz,rate,err_rate,nsn,err_nsn=sn_rate(bins=bin_center)
        
        print 'Nsn',np.sum(nsn),zz,bin_center,nsn,err_nsn,np.power(np.sum(np.power(err_nsn,2.)),0.5)
        
        nsn_season = interpolate.interp1d(zz,nsn)
        err_nsn_season = interpolate.interp1d(zz,err_nsn)

        combi=nsn_season(zz)*effi(zz)
        #yerr_combi=np.power(np.power(ratio_err*nsn_season(zz),2.)+np.power(0.1*nsn_season(zz)*effi(zz),2.),0.5)
        yerr_combi=np.power(np.power(ratio_err*nsn_season(zz),2.)+np.power(err_nsn_season(zz)*effi(zz),2.),0.5)

        print 'Number of SN Ia',season,np.sum(combi),np.power(np.sum(np.power(yerr_combi,2.)),0.5)
        N_sn=np.sum(combi)
        err_N_sn=np.power(np.sum(np.power(yerr_combi,2.)),0.5)
        self.nsn_tot+=np.sum(combi)
        self.err_tot+=np.power(err_N_sn,2.)
        #axc[1].plot(zz,nsn,'b+')
        #axc[1].plot(zz,nsn_season(zz),marker='.', mfc='black', mec='black', ms=8, linestyle='-',color='black')
        if idraw == 0:
            axb.errorbar(zz,nsn_season(zz),yerr=err_nsn_season(zz),marker='.', mfc='black', mec='black', ms=8, linestyle='-',color='black',label='N$_{SN Ia}$ expected ('+rate_name+' rate) : '+str(int(np.sum(nsn_season(zz))))+'$\pm$'+str(int(np.power(np.sum(np.power(err_nsn_season(zz),2.)),0.5))))
        #axc[1].plot(zz,nsn_season(zz)*effi(zz),marker='.', mfc='red', mec='red', ms=8, linestyle='-',color=color)
        if season < 9:
            ll='Y'+str(season+1)+'   - N$_{SN Ia}$ = '+str(int(N_sn))+' $\pm$ '+str(int(err_N_sn))
        else:
            ll='Y'+str(season+1)+' - N$_{SN Ia}$ = '+str(int(N_sn))+' $\pm$ '+str(int(err_N_sn))
        axb.errorbar(zz,combi,yerr=yerr_combi,marker=self.ms[season], mfc=self.color[season], mec=self.color[season], ms=8, linestyle='-',color=color,label=ll)
        #axc[1].set_yscale('log')
        axb.set_xlabel('z')
        axb.set_ylabel('Number of SN Ia')
        axb.legend(loc='best',prop={'size':12})
        axb.set_xlim(self.zmin,np.max(bin_center)+0.01)
        

    def Plot(self, axc,sel,vary,varx,color):
        print 'there',sel[varx],sel[vary]
        axc.plot(sel[varx],sel[vary],color+'.')

    def Plot_m5(self,tab_resu):
        for band in 'grizy':
            self.Plot_m5_Indiv(tab_resu,band)

    def Plot_m5_Indiv(self,tab_resu,band):
        figa, axa = plt.subplots(ncols=1, nrows=1, figsize=(10,9))
        what='m5sigma'
        res=[]
        for key, val in tab_resu.items():
            print key,val.dtype
            season=int(key.split('_')[2])-1
            spl=key.split('_')
            fieldname=spl[0]
            fieldid=spl[1]
            sel=val[np.where(np.logical_and(val['status']=='go_fit',val['fit_status']=='fit_ok'))]
            sel = sel[~np.isnan(sel['m5sigma_'+band])]
            res.append((season+1,np.median(sel['m5sigma_'+band]),np.std(sel['m5sigma_'+band])))

        tot=np.rec.fromrecords(res,names=['season','m5_mean','m5_rms'])
        print band,tot['m5_mean']
        axa.errorbar(tot['season'],tot['m5_mean'], yerr=tot['m5_rms'],marker='.', mfc='k', mec='k',ms=8,color='k')
        axa.set_xlabel('year')
        axa.set_ylabel('<m5> [mag]')
        title=fieldname+' - '+fieldid+' - band '+band
        axa.set_title(title)
        axa.set_xlim(0.8,10.2)

    def Plot_Cadences(self,tab_resu):
        
        for band in 'grizy':
            self.Plot_Cadences_Indiv(tab_resu,band)
        self.Plot_Cadences_Indiv(tab_resu,'') 

    def Plot_Cadences_Indiv(self,tab_resu,band):

        figa, axa = plt.subplots(ncols=1, nrows=1, figsize=(10,9))
        what='cad'
        if band != '':
            what+='_'+band
        
        zstep=0.05
        for key, val in tab_resu.items():
            res=[]
            print key,val.dtype
            for z in np.arange(0.,1.,zstep):
                sel=val[np.where(np.logical_and(val['status']=='go_fit',val['fit_status']=='fit_ok'))]
                #sel=sel[np.where(np.logical_and(sel['phase_first']<=-5,sel['phase_last']>=20))]
                sel=sel[np.where(np.logical_or(sel['phase_first']>-5,sel['phase_last']<20))]
                sel=sel[np.where(np.logical_and(sel['z']>=z,sel['z']<z+zstep))]
                #print z,np.median(sel[bef]),np.median(sel[aft]),np.median(sel[bef]+sel[aft])
                idx = sel[what] > 0.
                sel=sel[idx]
                res.append((z+zstep/2.,np.median(sel[what]),np.std(sel[what])))

            season=int(key.split('_')[2])-1
            spl=key.split('_')
            fieldname=spl[0]
            fieldid=spl[1]
            tot=np.rec.fromrecords(res,names=['z','cad_mean','cad_rms'])
            ll='Y'+str(season+1)
            axa.errorbar(tot['z'],tot['cad_mean'], yerr=tot['cad_rms'],marker=self.ms[season], mfc=self.color[season], mec=self.color[season],ms=8, linestyle='-',color='k',label=ll)          
           
 
        axa.set_xlabel('z')
        axa.set_ylabel('Mean cadence (observer frame) [days$^{-1}$]')
        
        title=fieldname+' - '+fieldid
        if band != '':
            title +=' - '+band+' band'
        axa.set_title(title)
        axa.legend(loc='best',prop={'size':12})
        #axa.plot(tot['z'],tot['Nbef'],'bo')
    
    def Plot_Nobs(self,tab_resu):
        
        for band in 'grizy':
            self.Plot_Nobs_Indiv(tab_resu,band)
        self.Plot_Nobs_Indiv(tab_resu,'') 

    def Plot_Nobs_Indiv(self,tab_resu,band):

        figa, axa = plt.subplots(ncols=1, nrows=3, figsize=(10,9))
        bef='N_bef'
        aft='N_aft'
        
        if band != '':
            bef+='_'+band
            aft+='_'+band
        
        zstep=0.05
        zrange=np.arange(0.,1.,zstep)
        for key, val in tab_resu.items():
            res=[]
            #print key,val.dtype
            for z in zrange:
                sel=val[np.where(np.logical_and(val['status']=='go_fit',val['fit_status']=='fit_ok'))]
                sel=sel[np.where(np.logical_and(sel['z']>=z,sel['z']<z+zstep))]
                #print z,np.median(sel[bef]),np.median(sel[aft]),np.median(sel[bef]+sel[aft])
                res.append((z+zstep/2.,np.median(sel[bef]),np.std(sel[bef]),np.median(sel[aft]),np.std(sel[aft]),np.median(sel[bef]+sel[aft]),np.std(sel[bef]+sel[aft])))

            season=int(key.split('_')[2])-1
            spl=key.split('_')
            fieldname=spl[0]
            fieldid=spl[1]
            
            tot=np.rec.fromrecords(res,names=['z','Nbef_mean','Nbef_rms','Naft_mean','Naft_rms','Ntot_mean','Ntot_rms'])
            ll='Y'+str(season+1)
            axa[0].errorbar(tot['z'],tot['Nbef_mean'], yerr=tot['Nbef_rms'],marker=self.ms[season], mfc=self.color[season], mec=self.color[season],ms=8, linestyle='-',color='k',label=ll)          
            axa[1].errorbar(tot['z'],tot['Naft_mean'], yerr=tot['Naft_rms'],marker=self.ms[season], mfc=self.color[season], mec=self.color[season],ms=8, linestyle='-',color='k',label=ll)
            axa[2].errorbar(tot['z'],tot['Ntot_mean'], yerr=tot['Ntot_rms'],marker=self.ms[season], mfc=self.color[season], mec=self.color[season],ms=8, linestyle='-',color='k',label=ll)
 
        axa[2].legend(loc='upper right',prop={'size':12})
        axa[0].set_xlabel('z')
        axa[0].set_ylabel('N$_{obs}$ before T0')
        axa[1].set_xlabel('z')
        axa[1].set_ylabel('N$_{obs}$ after T0')
        axa[2].set_xlabel('z')
        axa[2].set_ylabel('N$_{obs}$ in [T0-20,T0+40]')
        
        title=fieldname+' - '+fieldid
        if band != '':
            title +=' - '+band+' band'
        else:
            axa[0].plot(zrange,[4]*len(zrange),ls='-',color='k')
            axa[1].plot(zrange,[10]*len(zrange),ls='-',color='k') 

        axa[0].set_title(title)
        #axa.plot(tot['z'],tot['Nbef'],'bo')
    
    def Plot_Phases(self,tab_resu):

        figa, axa = plt.subplots(ncols=1, nrows=2, figsize=(10,9))

        phasef='phase_first'
        phasel='phase_last'
        zstep=0.05
        zrange=np.arange(0.,1.,zstep)
        for key, val in tab_resu.items():
            res=[]
            sel=val[np.where(np.logical_and(val['status']=='go_fit',val['fit_status']=='fit_ok'))]
            #print sel['z']
            #sel=sel[np.where(np.logical_or(sel['phase_first']>-5,sel['phase_last']<20))]
            for z in zrange:
                sol=sel[np.where(np.logical_and(sel['z']>=z,sel['z']<z+zstep))]
                #print z,np.median(sel[phasef]),np.median(sel[phasel])
                res.append((z+zstep/2.,np.median(sol[phasef]),np.std(sol[phasef]),np.median(sol[phasel]),np.std(sol[phasel])))
            
            tot=np.rec.fromrecords(res,names=['z','phase_first','phase_first_rms','phase_last','phase_last_rms']) 
            season=int(key.split('_')[2])-1
            spl=key.split('_')
            fieldname=spl[0]
            fieldid=spl[1]
            ll='Y'+str(season+1)
            #print 'legend',key,season,self.color[season]
            axa[0].errorbar(tot['z'],tot['phase_first'], yerr=tot['phase_first_rms'],marker=self.ms[season], mfc=self.color[season], mec=self.color[season],ms=8, linestyle='-',color='k',label=ll)
            axa[1].errorbar(tot['z'],tot['phase_last'], yerr=tot['phase_last_rms'],marker=self.ms[season], mfc=self.color[season], mec=self.color[season],ms=8, linestyle='-',color='k',label=ll)

        axa[1].legend(loc='upper right',prop={'size':12})
        axa[0].set_xlabel('z')
        axa[0].set_ylabel('<phase first point> (restframe) [days]')
        axa[1].set_xlabel('z')
        axa[1].set_ylabel('<phase last point> (restframe) [days]')

        axa[0].plot(zrange,[-5.]*len(zrange),ls='-',color='k')
        axa[1].plot(zrange,[20.]*len(zrange),ls='-',color='k')
        title=fieldname+' - '+fieldid
        axa[0].set_title(title)

    def Plot_Color(self,tab_resu):

        figa, axa = plt.subplots(ncols=1, nrows=1, figsize=(10,9))
        for key, val in tab_resu.items():
            sel=val[np.where(np.logical_and(val['status']=='go_fit',val['fit_status']=='fit_ok'))]
            sel=sel[np.where(np.logical_and(sel['phase_first']<=-5,sel['phase_last']>=20))]
            axa.plot(sel['z'],np.sqrt(sel['salt2.CovColorColor']),'bo')
            """
            axa[1].plot(sel['phase_first'],np.sqrt(sel['salt2.CovColorColor']),'bo')
            axa[2].plot(sel['phase_last'],np.sqrt(sel['salt2.CovColorColor']),'bo')
            """


    def Plot_Effi_vs_Cadence(self,tab_resu):
        for band in 'grizy':
            self.Plot_Effi_vs_Cadence_Indiv(tab_resu,band)
        self.Plot_Effi_vs_Cadence_Indiv(tab_resu,'')

    def Plot_Effi_vs_Cadence_Indiv(self,tab_resu,band):
        
        figa, axa = plt.subplots(ncols=1, nrows=1, figsize=(10,9))
        #band='g'
        what='cad'
        if band != '':
            what+='_'+band
        zstep=0.1
        zrange=np.arange(0.,0.3,zstep)
        for key, val in tab_resu.items():
            res=[]
            sel=val.copy()
            #sel=sel[np.where(np.logical_and(sel['N_bef']>=4,sel['N_aft']>=10))]
            sel=sel[np.where(sel['N_bef']>=4)]
            sel=sel[np.where(np.logical_and(sel['status']=='go_fit',sel['fit_status']=='fit_ok'))]
            #sel=sel[np.where(np.logical_and(sel['phase_first']<=-5,sel['phase_last']>=20))]
            for z in zrange:
                sol=sel[np.where(np.logical_and(sel['z']>=z,sel['z']<z+zstep))]
                ref=val[np.where(np.logical_and(val['z']>=z,val['z']<z+zstep))]
                res.append((z+zstep/2.,float(len(sol))/float(len(ref)),np.median(sol[what]),np.std(sol[what])))
            season=int(key.split('_')[2])-1
            spl=key.split('_')
            fieldname=spl[0]
            fieldid=spl[1]
            
           
            tot=np.rec.fromrecords(res,names=['z','effi','cad_mean','cad_rms'])
            print tot
            ll='Y'+str(season+1)
            #axa.errorbar(tot['effi'],tot['cad_mean'], yerr=tot['cad_rms'],marker=self.ms[season], mfc=self.color[season], mec=self.color[season],ms=8, linestyle='-',color='k',label=ll) 
            axa.plot(tot['cad_mean'],tot['effi'],marker=self.ms[season], mfc=self.color[season], mec=self.color[season],ms=8, linestyle='-',color='k',label=ll) 
            
        axa.set_xlabel('Median cadence [days$^{-1}$]')
        axa.set_xlabel('Detection efficiency')
        title=fieldname+' - '+fieldid+' - band '+band
        axa.set_title(title)

parser = OptionParser()
parser.add_option("-f", "--fieldname", type="string", default='DD', help="filter [%default]")
parser.add_option("-i", "--fieldid", type="int", default=290, help="filter [%default]")
parser.add_option("-x", "--stretch", type="float", default=-999., help="filter [%default]")
parser.add_option("-c", "--color", type="float", default=-999., help="filter [%default]")
parser.add_option("-d", "--dirmeas", type="string", default="DD", help="filter [%default]")

opts, args = parser.parse_args()

dict_ana={}
list_ana=[]
#dict_ana['Mean_Obs_Faint_T0_0']=Id(thedir='Mean_Obs_T0_0',fieldname='WFD',fieldid=309,X1=-2.0,Color=0.2,season=0,colorfig='k')
#dict_ana['Mean_Obs']=Id(thedir='Mean_Obs',fieldname='WFD',fieldid=309,X1=-999.,Color=-999.,season=0,colorfig='r')

fieldid=opts.fieldid
fieldname=opts.fieldname
thedir=opts.dirmeas

#X1=0.0
#Color=0.0

X1=opts.stretch
Color=opts.color
#for seas in [1]:
for seas in [i for i in range(10)]:
    #dict_ana['Rolling_Faint_Seas'+str(seas+1)]=Id(thedir='WFD_Rolling_noTwilight',fieldname='WFD',fieldid=309,X1=-2.0,Color=0.2,season=seas,colorfig='k')
    dict_ana['DD_'+str(fieldid)+'_'+str(seas+1)]=Id(thedir=thedir,fieldname=fieldname,fieldid=fieldid,X1=X1,Color=Color,season=seas,colorfig='k')
    #list_ana.append(('DD_290_'+str(seas+1),Id(thedir='DD',fieldname='DD',fieldid=290,X1=-1.,Color=-1.,season=seas,colorfig='k')))
    
#dict_ana['Mean_Obs_Faint_newrefs']=Id(thedir='Mean_Obs_newrefs',fieldname='WFD',fieldid=309,X1=-2.0,Color=0.2,season=0,colorfig='k')

#dict_ana_ordered=collections.OrderedDict(sorted(dict_ana.items()))
Ana_Simu(dict_ana,zmin=0.,zmax=1.0,bin_z=0.1)
