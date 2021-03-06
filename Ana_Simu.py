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
from astropy.table import vstack,Table
from Observations import *
from ID import *

class Ana_Simu:
    def __init__(self, dict_ana,zmin=0.,zmax=1.2,bin_z=0.01):

        thedir='Prod_LC'

        tot_resu={}

        self.zmin=zmin
        self.zmax=zmax
        self.bin_z=bin_z
        dir_observations_OpSim='../Ana_Cadence/OpSimLogs/'
        dict_obs={}
        myobs={}

        for key,val in dict_ana.items():
            #key=vals[0]
            #val=vals[1]
            print 'hee',val.thedir
            dirmeas=thedir+'/'+val.thedir+'/'+str(val.fieldid)+'/Season_'+str(val.season)
            sum_file=dirmeas+'/'+val.fieldname+'_'+str(val.fieldid)+'_X1_'+str(val.X1)+'_C_'+str(val.Color)+'_all.pkl'
            if not myobs.has_key(val.fieldname):
                myobs[val.fieldname]={}
                if not myobs[val.fieldname].has_key(val.fieldid):
                   name=val.thedir+'/Observations_'+val.fieldname+'_'+str(val.fieldid)+'.txt'
                   myobs[val.fieldname][val.fieldid]=Observations(fieldid=val.fieldid, filename=dir_observations_OpSim+'/'+name) 
                   print 'loading observations',name,key,val.season

            if not dict_obs.has_key(key):
                dict_obs[key]={}
           
            dict_obs[key][val.season]=myobs[val.fieldname][val.fieldid].seasons[val.season]
            

            if os.path.exists(sum_file):
                pkl_file = open(sum_file,'rb')
                loaded_file=pkl.load(pkl_file)
                #print 'done',key,tot_resu[key],tot_resu[key].dtype,sum_file
                #idx = loaded_file['z']>=self.zmin and loaded_file['z']<self.zmax
                tot_resu[key]=loaded_file
                #print 'done',key,sum_file,tot_resu[key].dtype
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
       
            """
            print 'there',key,tot_resu[key],tot_resu[key].dtype
            for band in 'grizy':
                print band,tot_resu[key]['N_bef_'+band],tot_resu[key]['N_aft_'+band]
            print 'all bands',tot_resu[key]['N_bef'],tot_resu[key]['N_aft']
            """

        # this is for simulated parameters

        
        figb, axb = plt.subplots(ncols=2, nrows=2, figsize=(10,9))
        for key in dict_ana.keys():
            self.Plot_Sims(axb,tot_resu[key])
        

        self.nsn_tot=0.
        self.err_tot=0.
        self.nsn_theo=0
        self.err_theo=0
        self.res_nsn=[]
        self.ms=['o','o','s','s','.','.','^','^','<','<']
        self.color=['b','r','b','r','b','r','b','r','b','r']
        figa, axa = plt.subplots(ncols=1, nrows=1, figsize=(10,9))
        axca = axa.twinx()

        figaa, axaa = plt.subplots(ncols=1, nrows=1, figsize=(10,9))
        axcaa = axaa.twinx()

        figc, axc = plt.subplots(ncols=1, nrows=1, figsize=(10,9))
        title=val.fieldname+' - '+str(val.fieldid)
        

        #figbb, axbb = plt.subplots(ncols=2, nrows=2, figsize=(10,9))
        col=['k','r','b']
        idraw=dict(zip([0,1,2,3,4],[0,0,0,0,0]))
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
            for keyb in range(4):
                idraw[keyb]+=1
            print 'Nentries',len(val),np.unique(val['z'])

            sel=val.copy()
            selb=val.copy()
            selc=val.copy()
            seld=val.copy()
            sel=sel[np.where(np.logical_and(sel['N_bef']>=4,sel['N_aft']>=10))]
            logical_test={}
            
            for band in 'grizy':
                logical_test[band]=np.logical_and(sel['N_bef_'+band]>=1,sel['N_aft_'+band]>=1)

            logical_and_g_r=np.logical_and(logical_test['g'],logical_test['r'])
            logical_and_r_i=np.logical_and(logical_test['r'],logical_test['i'])
            logical_and_i_z=np.logical_and(logical_test['i'],logical_test['z'])
            logical_and_z_y=np.logical_and(logical_test['z'],logical_test['y'])
                    
            #sel=sel[np.where(np.logical_or(np.logical_or(np.logical_or(logical_and_g_r,logical_and_r_i),logical_and_i_z),logical_and_z_y))]

            #sel=val[np.where(val['N_bef']>=4)]
            #sel=val[np.where(val['N_aft']>=10)]
            #sel=sel[np.where(sel['status']=='go_fit')]
            sel=sel[np.where(np.logical_and(sel['status']=='go_fit',sel['fit_status']=='fit_ok'))]
            sel=sel[np.where(np.logical_and(sel['phase_first']<=-5,sel['phase_last']>=20))]
            #sel=sel[np.where(sel['phase_first']<=-5.)]
            #sel=sel[np.where(np.sqrt(sel['salt2.CovColorColor'])<0.04)]

            selb=selb[np.where(selb['status']=='no_obs')]
            #selc=selc[np.where(selc['fit_status']=='crashd')]
            #selc=selc[np.where(selc['status']=='unknown')]
            #selc=selc[np.where(np.logical_or(selc['phase_first']>=-5,selc['phase_last']<=20))]
            #selc=selc[np.where(np.logical_and(selc['phase_first']>=-500.,selc['phase_last']>=-500.))]
            selc=selc[np.where(np.logical_and(selc['status']=='go_fit',selc['fit_status']=='fit_ok'))]
            #selc=selc[np.where(np.logical_or(selc['phase_first']>-5,selc['phase_last']<20))]
            selc=selc[np.where(~np.logical_and(selc['phase_first']<=-5,selc['phase_last']>=20))]
            #seld=seld[np.where(np.logical_and(seld['status']=='go_fit',seld['fit_status']=='fit_ok'))]
            
            #seld=seld[np.where(np.logical_and(seld['status']=='go_fit',seld['fit_status']=='crashd'))]
            #seld=seld[np.where(seld['status']=='no_pha')]
            seld=seld[np.where(~np.logical_and(seld['N_bef']>=4,seld['N_aft']>=10))]
            """
        #print sel['phase_first'],sel['phase_last']
        for i,val in enumerate([(0.,0.),(2.0,-0.2),(-2.0,0.2)]):
        #for i,val in enumerate([(-2.0,0.2)]):    
            sela=tot_resu[np.where(np.logical_and(tot_resu['X1']==val[0],tot_resu['Color']==val[1]))]
            selc=sel[np.where(np.logical_and(sel['X1']==val[0],sel['Color']==val[1]))]
            """
            #Efficiency per season vs z plot
            
            self.Plot_Eff(axa,val,sel,'z',dict_ana[key].colorfig,dict_ana[key].season,key,0,idraw[0],self.zmin,self.zmax,self.bin_z)
            #self.Plot_Eff(axca,val,selb,'z',dict_ana[key].colorfig,dict_ana[key].season,key,1,idraw[1],self.zmin,self.zmax,self.bin_z)
            #self.Plot_Eff(axca,val,selc,'z',dict_ana[key].colorfig,dict_ana[key].season,key,2,idraw[2],self.zmin,self.zmax,self.bin_z)
            #self.Plot_Eff(axca,val,seld,'z',dict_ana[key].colorfig,dict_ana[key].season,key,3,idraw[3],self.zmin,self.zmax,self.bin_z)
            
            #Efficiency per season vs T0 plot
            
            idraw[4]+=1
            min_T0=np.min(val['T0'])
            max_T0=np.max(val['T0'])
            bin_T0=2.
            min_z=0.3
            max_z=0.4
            sel_val=self.Select('z',val,min_z,max_z)
            sel_sel=self.Select('z',sel,min_z,max_z)
            sel_selb=self.Select('z',selb,min_z,max_z)
            sel_selc=self.Select('z',selc,min_z,max_z)
            sel_seld=self.Select('z',seld,min_z,max_z)
            """
            mmin_T0=61640
            mmax_T0=61650
            
            sel_val=self.Select('T0',sel_val,mmin_T0,mmax_T0)
            sel_sel=self.Select('T0',sel_sel,mmin_T0,mmax_T0)
            sel_selb=self.Select('T0',sel_selb,mmin_T0,mmax_T0)
            sel_selc=self.Select('T0',sel_selc,mmin_T0,mmax_T0)
            sel_seld=self.Select('T0',sel_seld,mmin_T0,mmax_T0)
            
            test=self.Select('T0',sel_val,mmin_T0,mmax_T0)
            testb=self.Select('phase_first',test,-5.,1000.)
            #testb=self.Select('phase_last',test,-999,20.)
            testb=testb[np.where(testb['status']=='go_fit')]
            print len(test),len(testb),testb['phase_first'],testb['T0'],testb['X1'],testb['Color'],testb['z'],testb['phase_last'],test['phase_last']
            """
            
            print 'before T0',len(sel_val),len(sel_sel)
            self.Plot_Eff(axaa,sel_val,sel_sel,'T0',dict_ana[key].colorfig,dict_ana[key].season,key,4,idraw[4],min_T0,max_T0,bin_T0,dict_obs[key][dict_ana[key].season])
            self.Plot_Eff(axcaa,sel_val,sel_sel,'T0',dict_ana[key].colorfig,dict_ana[key].season,key,5,idraw[1],min_T0,max_T0,bin_T0)
                #self.Plot_Eff(axcaa,sel_val,sel_selb,'T0',dict_ana[key].colorfig,dict_ana[key].season,key,1,idraw[1],min_T0,max_T0,bin_T0)
                #self.Plot_Eff(axcaa,sel_val,sel_selc,'T0',dict_ana[key].colorfig,dict_ana[key].season,key,2,idraw[2],min_T0,max_T0,bin_T0)
                #self.Plot_Eff(axcaa,sel_val,sel_seld,'T0',dict_ana[key].colorfig,dict_ana[key].season,key,3,idraw[3],min_T0,max_T0,bin_T0)


            # number of supernovae per season
            self.Plot_N_SN(axc,val,sel,'z',dict_ana[key].colorfig,dict_ana[key].season,key,idraw,dict_ana[key].fieldname,dict_ana[key].fieldid,self.zmin,self.zmax,self.bin_z,cumul=False)

            #print 'boouh',key,dict_ana[key].colorfig,dict_ana[key].season

            #self.Plot(axc,sel,'salt2.CovColorColor','z',dict_ana[key].colorfig)
            #print 'Number of Events',val[0],val[1],len(val)
            for vval in np.arange(self.zmin,self.zmax,0.1):
                ssel=val[np.where(np.logical_and(val['z']>=vval,val['z']<vval+0.1))]
                #print vval,vval+0.1,len(ssel)

        #print val['T0']
        print 'in total :',self.nsn_tot,'+-',np.power(self.err_tot,0.5)
        
        axa.set_title(title)
        axa.set_xlabel('z')
        axa.set_ylim([0.,1.05])
        axca.set_ylim([0.,1.05])
        
        axaa.set_title(title)
        axaa.set_xlabel('DayMax [day]')
        axaa.set_ylim([0.,1.05])
        #axcaa.set_ylim([0.,10000.])

        #axa.set_xlim(self.zmin,self.zmax+0.01)
        nsntot_str=str(int(self.nsn_tot))+'$\pm$'+str(int(np.power(self.err_tot,0.5)))
        nsntheo_str=str(int(self.nsn_theo))+'$\pm$'+str(int(np.power(self.err_theo,0.5)))
        title+=' - N$_{SN Ia}$ = '+nsntot_str +' / '+nsntheo_str
        axc.set_title(title)

        if len(self.res_nsn) > 0:
            rate_nsn=np.rec.fromrecords(self.res_nsn,names=['fieldname','fieldid','season','n_sn_detected','err_detected','n_sn_expected','err_expected'])
            fieldname=rate_nsn['fieldname'][0]
            fieldid=rate_nsn['fieldid'][0]
            figc.savefig('Plots_NSN/NSN_'+fieldname+'_'+str(fieldid)+'.png')
            figa.savefig('Plots_NSN/Eff_'+fieldname+'_'+str(fieldid)+'.png')
            """
            pkl_out=open('N_SN_'+fieldname+'_'+str(fieldid)+'_colorcut.pkl','wb')
            pkl.dump(rate_nsn,pkl_out)
            pkl_out.close()
            """
        #axc.set_xlim(self.zmin,self.zmax+0.01)

        #self.Plot_m5(tot_resu)
        #self.Plot_Cadences(tot_resu)
        #self.Plot_Nobs(tot_resu)
        #self.Plot_Phases(tot_resu)
        #self.Plot_Color(tot_resu)
        #self.Plot_Effi_vs_Cadence(tot_resu)
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

    def Histo_ratio(self,sela,selb,varname,zmin,zmax,bin_z):

        range=[zmin,zmax]
        bin_width=bin_z
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

    def Plot_Eff(self,axc,sela,selb,varname,color,season,ll,idraw,icount,min_bin,max_bin,delta_bin,obs=None):
        tot_label=[]
        #figc, axc = plt.subplots(ncols=1, nrows=1, figsize=(10,9))
        self.Plot_Eff_Indiv(axc,sela,selb,varname,tot_label,ll,color,season,idraw,icount,min_bin,max_bin,delta_bin,obs)
        

    def Plot_Eff_Indiv(self,axc,sela,selb,varname,tot_label,ll,color,season,idraw,icount,min_bin,max_bin,delta_bin,obs=None):

        filtercolors = {'u':'b', 'g':'c', 'r':'g', 'i':'y', 'z':'r', 'y':'m'}

        fontsize=12
        bin_center, ratio, ratio_err,norm,norm_err= self.Histo_ratio(sela,selb,varname,min_bin,max_bin,delta_bin)
        #tot_label.append(axc.errorbar(bin_center,ratio, yerr=ratio_err,marker=marker, mfc=colors[key], mec=colors[key], ms=8, linestyle=myfmt[i],color='k',label=ll))
        ll='Y'+str(season+1)
        
        
        if idraw == 0 or idraw==4:
            if icount <=10:
                axc.errorbar(bin_center,ratio, yerr=ratio_err,marker=self.ms[season], mfc=self.color[season], mec=self.color[season], ms=8, linestyle='-',color=color,label=ll)
            else:
                axc.errorbar(bin_center,ratio, yerr=ratio_err,marker=self.ms[season], mfc=self.color[season], mec=self.color[season], ms=8, linestyle='-',color=color) 
            #axc.errorbar(bin_center,ratio, yerr=ratio_err,marker=self.ms[season], mfc=color, mec=color, ms=8, linestyle='-',color=color,label=ll)
            if idraw == 0:
                axc.legend(loc='center right',prop={'size':fontsize})
                #axc.legend(loc='center left',prop={'size':fontsize})
            if icount == 1:
                axc.set_ylabel('Efficiency')
            if obs is not None:
                   #print 'yyyyyyyyy',np.min(obs['mjd']),np.max(obs['mjd'])
                   idx=obs['band']!= 'LSSTPG::u'
                   print 'season obs ',season,len(sela),len(selb)
                   axc.plot(obs['mjd'][idx],[0.6]*len(obs['mjd'][idx]),'k*')
                   for icol,f in enumerate('grizy'):
                       idxb=obs['band']== 'LSSTPG::'+f
                       axc.plot(obs['mjd'][idxb],[0.5-0.1*icol]*len(obs['mjd'][idxb]),filtercolors[f]+'*')
            if idraw == 0:
                axc.set_xlim([min_bin,max_bin])
        else:
            if idraw == 1:
                myls='--'
                ll='No obs - SNR >5.'
                if icount == 1:
                    axc.plot([min_bin],[0.],linestyle=myls,color='k',label=ll)
            if idraw == 2:
                myls = ':'
                ll='Phase cut'
                if icount == 1:
                    axc.plot([min_bin],[0.],linestyle=myls,color='k',label=ll)
            if idraw == 3:
                myls = '-.'
                ll='No Obs in phase [-30,50]' 
                if icount == 1:
                    axc.plot([min_bin],[0.],linestyle=myls,color='k',label=ll)

            axc.legend(loc='best',prop={'size':fontsize},frameon=False)
            #axc.legend(bbox_to_anchor=(0.5,0.42), loc=2, borderaxespad=0.,prop={'size':fontsize},frameon=False)
            if idraw <=3:
                axc.errorbar(bin_center,ratio, yerr=ratio_err,marker=self.ms[season],  mfc=self.color[season], mec=self.color[season], ms=8, linestyle=myls,color='k')
            else:
                #print 'aloo',sela.dtype
                print selb.dtype,len(selb),len(sela)
                selb.sort(order='T0')
                #curve = interpolate.interp1d(selb['T0'],selb['SNR_tot'])
                #axc.errorbar(bin_center,curve(bin_center),color='k')
                axc.errorbar(selb['T0'],selb['SNR_tot'],color='k')
                print selb.dtype
                for f in 'grizy':
                    #curve = interpolate.interp1d(sela['T0'],sela['SNR_'+f]/(sela['N_aft_'+f]+sela['N_bef_'+f]))
                    #curve = interpolate.interp1d(selb['T0'],selb['SNR_'+f])
                    #print curve(bin_center)
                    #axc.errorbar(bin_center,curve(bin_center),color=filtercolors[f])
                    axc.errorbar(selb['T0'],selb['SNR_'+f],color=filtercolors[f])
            #axc.errorbar(bin_center,ratio, yerr=ratio_err,marker=self.ms[season],  mfc=color, mec=color, ms=8, linestyle=myls,color='k',label=ll)
            axc.set_ylabel('Fraction of Events') 
            

    def Plot_N_SN(self,axb,sela,selb,varname,color,season,ll,idraw,fieldname,fieldid,min_bin,max_bin,delta_bin,cumul):
        tot_label=[]
        #figc, axc = plt.subplots(ncols=1, nrows=1, figsize=(10,9))
        self.Plot_Eff_Indiv_N_SN(axb,sela,selb,varname,tot_label,ll,color,season,idraw,fieldname,fieldid,min_bin,max_bin,delta_bin,cumul)

    def Plot_Eff_Indiv_N_SN(self,axb,sela,selb,varname,tot_label,ll,color,season,idraw,fieldname,fieldid,min_bin,max_bin,delta_bin,cumul=False):

        bin_center, ratio, ratio_err,norm,norm_err= self.Histo_ratio(sela,selb,varname,min_bin,max_bin,delta_bin)

        effi = interpolate.interp1d(bin_center,ratio)

        rate_name='Perrett'
        sn_rate=SN_Rate(rate=rate_name,duration=(selb['duration'][0]+0.)/365.25)

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
        self.nsn_theo+=np.sum(nsn_season(zz))
        self.err_theo+=np.sum(np.power(err_nsn_season(zz),2.))
        self.res_nsn.append((fieldname,fieldid,season,np.sum(combi),err_N_sn,np.sum(nsn_season(zz)),np.sum(np.power(err_nsn_season(zz),2.))))
        #axc[1].plot(zz,nsn,'b+')
        #axc[1].plot(zz,nsn_season(zz),marker='.', mfc='black', mec='black', ms=8, linestyle='-',color='black')
        if idraw >= 0:
            if not cumul:
                axb.errorbar(zz,nsn_season(zz),yerr=err_nsn_season(zz),marker=self.ms[season], mfc=self.color[season], mec=self.color[season], ms=8, linestyle='--',color='black')
            else:
                axb.errorbar(zz,np.cumsum(nsn_season(zz)),marker=self.ms[season], mfc=self.color[season], mec=self.color[season], ms=8, linestyle='--',color='black')
            
        #axc[1].plot(zz,nsn_season(zz)*effi(zz),marker='.', mfc='red', mec='red', ms=8, linestyle='-',color=color)
        nsn_str= str(int(np.sum(nsn_season(zz))))+'$\pm$'+str(int(np.power(np.sum(np.power(err_nsn_season(zz),2.)),0.5)))
        if season < 9:
            ll='Y'+str(season+1)+'   - N$_{SN Ia}$ = '+str(int(N_sn))+' $\pm$ '+str(int(err_N_sn))+' / '+nsn_str
        else:
            ll='Y'+str(season+1)+' - N$_{SN Ia}$ = '+str(int(N_sn))+' $\pm$ '+str(int(err_N_sn))+ ' / '+nsn_str

        if not cumul:
            axb.errorbar(zz,combi,yerr=yerr_combi,marker=self.ms[season], mfc=self.color[season], mec=self.color[season], ms=8, linestyle='-',color=color,label=ll)
        else:
            axb.errorbar(zz,np.cumsum(combi),marker=self.ms[season], mfc=self.color[season], mec=self.color[season], ms=8, linestyle='-',color=color,label=ll)
        #axc[1].set_yscale('log')
        axb.set_xlabel('z')
        if cumul:
            axb.set_ylabel('Number of SN Ia < z')
        else:
           axb.set_ylabel('Number of SN Ia per z bin') 
        axb.legend(loc='best',prop={'size':12})
        axb.set_xlim(self.zmin,np.max(bin_center)+0.01)
        

    def Plot(self, axc,sel,vary,varx,color):
        print 'there',sel[varx],sel[vary]
        axc.plot(sel[varx],sel[vary],color+'.')

    def Plot_m5(self,tab_resu):
        
        figa, axa = plt.subplots(ncols=1, nrows=1, figsize=(10,9))
        what='m5sigma'
        res=[]
        filtercolors = {'u':'b', 'g':'c', 'r':'g', 'i':'y', 'z':'r', 'y':'m'}
        for key, val in tab_resu.items():
            print key,val.dtype
            season=int(key.split('_')[2])-1
            spl=key.split('_')
            fieldname=spl[0]
            fieldid=spl[1]
            sel=val[np.where(np.logical_and(val['status']=='go_fit',val['fit_status']=='fit_ok'))]
            
            for band in 'grizy':
                sel = sel[~np.isnan(sel['m5sigma_'+band])]
                res.append((season+1,np.median(sel['m5sigma_'+band]),np.std(sel['m5sigma_'+band]),band))

        tot=np.rec.fromrecords(res,names=['season','m5_mean','m5_rms','band'])
        tot.sort(order='season')
        print band,tot['m5_mean']
        for band in 'grizy':
            idx = tot['band'] == band
            sel=tot[idx]
            axa.errorbar(sel['season'],sel['m5_mean'], marker='.', mfc='k', mec='k',ms=8,color=filtercolors[band],label=band+' band')
            axa.set_xlabel('year')
            axa.set_ylabel('Median m5 [mag]')
            title=fieldname+' - '+fieldid
        axa.set_title(title)
        axa.legend(loc='best',prop={'size':12})
        #axa.set_xlim(0.8,10.2)

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

        figa, axa = plt.subplots(ncols=1, nrows=2, figsize=(10,9))
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
            axa[0].errorbar(tot['z'],tot['cad_mean'],marker=self.ms[season], mfc=self.color[season], mec=self.color[season],ms=8, linestyle='-',color='k',label=ll) 
            axa[1].errorbar(tot['z'],tot['cad_rms'],marker=self.ms[season], mfc=self.color[season], mec=self.color[season],ms=8, linestyle='-',color='k',label=ll)
            
 
        axa[0].set_xlabel('z')
        axa[0].set_ylabel('Mean cadence (observer frame) [days$^{-1}$]')
        
        axa[1].set_xlabel('z')
        axa[1].set_ylabel('RMS cadence (observer frame) [days$^{-1}$]')


        title=fieldname+' - '+fieldid
        if band != '':
            title +=' - '+band+' band'
        axa[0].set_title(title)
        axa[0].legend(loc='best',prop={'size':12})
        #axa.plot(tot['z'],tot['Nbef'],'bo')
    
    def Plot_Nobs(self,tab_resu):
        
        for band in 'grizy':
            self.Plot_Nobs_Indiv_all(tab_resu,band)
        self.Plot_Nobs_Indiv_all(tab_resu,'') 

    def Plot_Nobs_Indiv_all(self,tab_resu,band):

        figa, axa = plt.subplots(ncols=1, nrows=3, figsize=(10,9))
        figa.suptitle(band+' band')
        bef='N_bef'
        aft='N_aft'
        
        if band != '':
            bef+='_'+band
            aft+='_'+band
        
        for key, val in tab_resu.items():
            sel=val[np.where(np.logical_and(val['status']=='go_fit',val['fit_status']=='fit_ok'))]
            sel=sel[np.where(np.logical_and(sel['T0']>=60.,sel['T0']<61.))]
            sel=sel[np.where(sel['z']<0.001)]
            axa[0].plot(sel['T0'],sel[bef],'bo')
            axa[1].plot(sel['T0'],sel[aft],'bo')
            axa[2].plot(sel['T0'],sel[bef]+sel[aft],'bo')

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
        zrange=np.arange(0.,1.+zstep,zstep)
        for key, val in tab_resu.items():
            res=[]
            sel=val[np.where(np.logical_and(val['status']=='go_fit',val['fit_status']=='fit_ok'))]
            sel=sel[np.where(np.logical_and(sel['N_bef']>=4,sel['N_aft']>=10))]
            #print sel['z']
            #sel=sel[np.where(np.logical_or(sel['phase_first']>-5,sel['phase_last']<20))]
            for z in zrange:
                sol=sel[np.where(np.logical_and(sel['z']>=z,sel['z']<z+zstep))]
                #print z,np.median(sel[phasef]),np.median(sel[phasel])
                res.append((z+zstep/2.,np.mean(sol[phasef]),np.std(sol[phasef]),np.mean(sol[phasel]),np.std(sol[phasel])))
            
            tot=np.rec.fromrecords(res,names=['z','phase_first','phase_first_rms','phase_last','phase_last_rms']) 
            season=int(key.split('_')[2])-1
            spl=key.split('_')
            fieldname=spl[0]
            fieldid=spl[1]
            season=int(fieldid)-120
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
        axa[0].set_xlim([np.min(zrange),np.max(zrange)])
        axa[1].set_xlim([np.min(zrange),np.max(zrange)])

    def Plot_Color(self,tab_resu):

        figa, axa = plt.subplots(ncols=1, nrows=1, figsize=(10,9))
        bcolors='bcgyrm'
        icol=-1
        fontsize=12
        for key, val in tab_resu.items():
            icol+=1
            sel=val[np.where(np.logical_and(val['status']=='go_fit',val['fit_status']=='fit_ok'))]
            sel=sel[np.where(np.logical_and(sel['N_bef']>=4,sel['N_aft']>=10))]
            selc=sel[np.where(~np.logical_and(sel['phase_first']<=-5,sel['phase_last']>=20))]
            sel=sel[np.where(np.logical_and(sel['phase_first']<=-5,sel['phase_last']>=20))]
            logical_or_r_i=np.logical_or(sel['N_bef_r']>=2,sel['N_bef_i']>=2)
            logical_or_z_y=np.logical_or(sel['N_bef_z']>=2,sel['N_bef_y']>=2)
            
                    
            sel=sel[np.where(np.logical_or(logical_or_r_i,logical_or_z_y))]


            #axa.plot(sel['z'],np.sqrt(sel['salt2.CovColorColor']),'.',label=key)
            axa.plot(sel['T0'],np.sqrt(sel['salt2.CovColorColor']),'.',label=key)
            #axa.plot(selc['z'],np.sqrt(selc['salt2.CovColorColor']),'k*',label=key)
            print 'ahahah',key
            """
            axa[1].plot(sel['phase_first'],np.sqrt(sel['salt2.CovColorColor']),'bo')
            axa[2].plot(sel['phase_last'],np.sqrt(sel['salt2.CovColorColor']),'bo')
            """
        axa.legend(loc='best',prop={'size':fontsize})
        axa.set_xlabel('z',{'fontsize': fontsize})
        axa.set_ylabel('$\sigma_{C}$',{'fontsize': fontsize})

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

    def Select(self,varname,val,valmin,valmax):
        idx=(val[varname]>=valmin)&(val[varname]<valmax)
        return val[idx]

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

Corresp=dict(zip([120,124,128,132,136,140,144,148],[3,4,5,6,7,8,9,10]))

filtercolors = dict(zip([120,124,128,132,136,140,144,148],['b','c','g','y','r','m','k','#eeefff']))
Corresp=dict(zip([123,127,131,135,139,143,147,151],[3,4,5,6,7,8,9,10]))
filtercolors = dict(zip([123,127,131,135,139,143,147,151],['b','c','g','y','r','m','k','#eeefff']))
filtercolors = dict(zip(range(10),['b','c','g','y','r','m','k','#ffa400','#fe871c','#1e4033']))
#for seas in [1]:
#for seas in [i for i in range(10)]:
for seas in [0]:
    #dict_ana['DD_'+str(fieldid)+'_'+str(seas+1)]=Id(thedir=thedir+'_last',fieldname=fieldname,fieldid=fieldid,X1=X1,Color=Color,season=seas,colorfig='k')
    dict_ana['DD_'+str(fieldid)+'_'+str(seas+1)]=Id(thedir=thedir,fieldname=fieldname,fieldid=fieldid,X1=X1,Color=Color,season=seas,colorfig='k')
    
    #dict_ana['DD1_'+str(fieldid)+'_'+str(seas)]=Id(thedir=thedir,fieldname=fieldname,fieldid=fieldid,X1=0.0,Color=0.0,season=seas,colorfig='k')
    #dict_ana['DD2_'+str(fieldid)+'_'+str(seas)]=Id(thedir=thedir,fieldname=fieldname,fieldid=fieldid,X1=2.0,Color=-0.2,season=seas,colorfig='k')
    #dict_ana['DD3_'+str(fieldid)+'_'+str(seas)]=Id(thedir=thedir,fieldname=fieldname,fieldid=fieldid,X1=-2.0,Color=0.2,season=seas,colorfig='k')
    
    #dict_ana['Rolling_Faint_Seas'+str(seas+1)]=Id(thedir='WFD_Rolling_noTwilight',fieldname='WFD',fieldid=309,X1=-2.0,Color=0.2,season=seas,colorfig='k')
    """
    for pseason in range(10):
        dict_ana['DD_'+str(pseason)+'_1_'+str(fieldid)+'_'+str(seas+1)]=Id(thedir=thedir,fieldname=fieldname+'_'+str(pseason)+'_1',fieldid=fieldid,X1=X1,Color=Color,season=seas,colorfig=filtercolors[pseason])
    """
    

    #dict_ana['DD_1_2_'+str(fieldid)+'_'+str(seas+1)]=Id(thedir=thedir,fieldname=fieldname+'_1_2',fieldid=fieldid,X1=X1,Color=Color,season=seas,colorfig='k')
    #dict_ana['DD_'+str(fieldid+1)+'_'+str(seas+1)]=Id(thedir=thedir,fieldname=fieldname,fieldid=fieldid+1,X1=X1,Color=Color,season=seas,colorfig='r')
    #dict_ana['DD_'+str(fieldid+2)+'_'+str(seas+1)]=Id(thedir=thedir,fieldname=fieldname,fieldid=fieldid+2,X1=X1,Color=Color,season=seas,colorfig='g')
    #list_ana.append(('DD_290_'+str(seas+1),Id(thedir='DD',fieldname='DD',fieldid=290,X1=-1.,Color=-1.,season=seas,colorfig='k')))
    
    """
    for X1,Color in zip([0.,2.,-2.],[0.,-0.2,0.2]):
    #for X1,Color in zip([2.],[-0.2]):
        print 'hello',X1,Color
        for i in range(1):
            val=i*4
            dict_ana['DD_'+str(Corresp[fieldid+val])+'_'+str(seas+1)+'_'+str(X1)+'_'+str(Color)]=Id(thedir=thedir,fieldname=fieldname,fieldid=fieldid+val,X1=X1,Color=Color,season=seas,colorfig=filtercolors[fieldid+val])
     """
       
#dict_ana['Mean_Obs_Faint_newrefs']=Id(thedir='Mean_Obs_newrefs',fieldname='WFD',fieldid=309,X1=-2.0,Color=0.2,season=0,colorfig='k')

#dict_ana_ordered=collections.OrderedDict(sorted(dict_ana.items()))
Ana_Simu(dict_ana,zmin=0.,zmax=1.,bin_z=0.05)
