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
            dirmeas=thedir+'/'+val.thedir+'/Season_'+str(val.season)
            sum_file=dirmeas+'/'+val.fieldname+'_'+str(val.fieldid)+'_X1_'+str(val.X1)+'_C_'+str(val.Color)+'_all.pkl'

            if os.path.exists(sum_file):
                pkl_file = open(sum_file,'rb')
                tot_resu[key]=pkl.load(pkl_file)
            else:

                files = glob.glob(dirmeas+'/'+val.fieldname+'_'+str(val.fieldid)+'*_X1_'+str(val.X1)+'_C_'+str(val.Color)+'*.pkl')
        
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
       

            print 'there',len(tot_resu[key]),tot_resu[key].dtype
        

        
        #self.Plot_Sims(tot_resu)

        
        figa, axa = plt.subplots(ncols=1, nrows=2, figsize=(10,9))
        #figc, axc = plt.subplots(ncols=1, nrows=1, figsize=(10,9))
        #figbb, axbb = plt.subplots(ncols=2, nrows=2, figsize=(10,9))
        col=['k','r','b']
        idraw=-1
        tot_resu_o=collections.OrderedDict(sorted(tot_resu.items()))
        for key, val in tot_resu_o.items():
            #self.Plot_Sims(axbb,val)
            idraw+=1
            sel=val.copy()
            sel=val[np.where(np.logical_and(val['N_bef']>=4,val['N_aft']>=10))]
            sel=sel[np.where(np.logical_and(sel['status']=='go_fit',sel['fit_status']=='fit_ok'))]
            sel=sel[np.where(np.logical_and(sel['phase_first']<=-5,sel['phase_last']>=20))]
            
            """
        #print sel['phase_first'],sel['phase_last']
        for i,val in enumerate([(0.,0.),(2.0,-0.2),(-2.0,0.2)]):
        #for i,val in enumerate([(-2.0,0.2)]):    
            sela=tot_resu[np.where(np.logical_and(tot_resu['X1']==val[0],tot_resu['Color']==val[1]))]
            selc=sel[np.where(np.logical_and(sel['X1']==val[0],sel['Color']==val[1]))]
            """
            self.Plot_Eff(axa,val,sel,'z',dict_ana[key].colorfig,dict_ana[key].season,key,idraw)
            #self.Plot(axc,sel,'salt2.CovColorColor','z',dict_ana[key].colorfig)
            #print 'Number of Events',val[0],val[1],len(val)
            for vval in np.arange(self.zmin,self.zmax,0.1):
                ssel=val[np.where(np.logical_and(val['z']>=vval,val['z']<vval+0.1))]
                print vval,vval+0.1,len(ssel)

        #print val['T0']
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
        axc[0].errorbar(bin_center,ratio, yerr=ratio_err,marker='.', mfc='red', mec='red', ms=8, linestyle='-',color=color,label=ll)
        axc[0].set_xlabel('z')
        axc[0].set_ylabel('Efficiency')
        axc[0].legend(loc='best',prop={'size':12})
        effi = interpolate.interp1d(bin_center,ratio)

        rate_name='Perret'
        sn_rate=SN_Rate(rate=rate_name)

        #zz,rate,err_rate,nsn,err_nsn=sn_rate(self.zmin,self.zmax,self.bin_z)
        zz,rate,err_rate,nsn,err_nsn=sn_rate(bins=bin_center)
        
        print 'Nsn',np.sum(nsn),zz,bin_center

        nsn_season = interpolate.interp1d(zz,nsn)
        err_nsn_season = interpolate.interp1d(zz,err_nsn)

        combi=nsn_season(zz)*effi(zz)
        yerr_combi=np.power(np.power(ratio_err*nsn_season(zz),2.)+np.power(0.1*nsn_season(zz)*effi(zz),2.),0.5)

        print 'Number of SN Ia',np.sum(combi),np.power(np.sum(np.power(yerr_combi,2.)),0.5)
        N_sn=np.sum(combi)
        err_N_sn=np.power(np.sum(np.power(yerr_combi,2.)),0.5)
        #axc[1].plot(zz,nsn,'b+')
        #axc[1].plot(zz,nsn_season(zz),marker='.', mfc='black', mec='black', ms=8, linestyle='-',color='black')
        if idraw == 0:
            axc[1].errorbar(zz,nsn_season(zz),yerr=err_nsn_season(zz),marker='.', mfc='black', mec='black', ms=8, linestyle='-',color='black',label='N$_{SN Ia}$ expected ('+rate_name+' rate)')
        #axc[1].plot(zz,nsn_season(zz)*effi(zz),marker='.', mfc='red', mec='red', ms=8, linestyle='-',color=color)
        
        ll='Y'+str(season+1)+' - N$_{SN Ia}$ = '+str(int(N_sn))+'$\pm$'+str(int(err_N_sn))
        axc[1].errorbar(zz,combi,yerr=yerr_combi,marker='.', mfc='red', mec='red', ms=8, linestyle='-',color=color,label=ll)
        #axc[1].set_yscale('log')
        axc[1].set_xlabel('z')
        axc[1].set_ylabel('Number of SN Ia')
        axc[1].legend(loc='best',prop={'size':8})

    def Plot(self, axc,sel,vary,varx,color):
        print 'there',sel[varx],sel[vary]
        axc.plot(sel[varx],np.sqrt(sel[vary]),color+'.')

dict_ana={}
list_ana=[]
#dict_ana['Mean_Obs_Faint_T0_0']=Id(thedir='Mean_Obs_T0_0',fieldname='WFD',fieldid=309,X1=-2.0,Color=0.2,season=0,colorfig='k')
#dict_ana['Mean_Obs']=Id(thedir='Mean_Obs',fieldname='WFD',fieldid=309,X1=-999.,Color=-999.,season=0,colorfig='r')

#for seas in [0,1,4,7]:
for seas in [i for i in range(10)]:
    #dict_ana['Rolling_Faint_Seas'+str(seas+1)]=Id(thedir='WFD_Rolling_noTwilight',fieldname='WFD',fieldid=309,X1=-2.0,Color=0.2,season=seas,colorfig='k')
    dict_ana['DD_290_'+str(seas+1)]=Id(thedir='DD',fieldname='DD',fieldid=290,X1=-1.,Color=-1.,season=seas,colorfig='k')
    #list_ana.append(('DD_290_'+str(seas+1),Id(thedir='DD',fieldname='DD',fieldid=290,X1=-1.,Color=-1.,season=seas,colorfig='k')))
    
#dict_ana['Mean_Obs_Faint_newrefs']=Id(thedir='Mean_Obs_newrefs',fieldname='WFD',fieldid=309,X1=-2.0,Color=0.2,season=0,colorfig='k')

#dict_ana_ordered=collections.OrderedDict(sorted(dict_ana.items()))
Ana_Simu(dict_ana,zmin=0.,zmax=1.4,bin_z=0.1)
