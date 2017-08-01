import cPickle as pkl
import numpy as np
import glob
import matplotlib.pyplot as plt

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
    def __init__(self, dict_ana):

        thedir='Prod_LC'

        tot_resu={}

        for key, val in dict_ana.items():
            print 'hee',val.thedir
            dirmeas=thedir+'/'+val.thedir+'/Season_'+str(val.season)
            files = glob.glob(dirmeas+'/'+val.fieldname+'_'+str(val.fieldid)+'*_X1_'+str(val.X1)+'_C_'+str(val.Color)+'*.pkl')
        
            for fi in files:
                pkl_file = open(fi,'rb')
                print 'loading',fi
                if not key in tot_resu.keys():
                    tot_resu[key]=pkl.load(pkl_file)
                else:
                    tot_resu[key]=np.vstack((tot_resu[key],pkl.load(pkl_file)))

            print 'there',len(tot_resu[key]),tot_resu[key].dtype
        

        
        #self.Plot_Sims(tot_resu)

        
        figa, axa = plt.subplots(ncols=1, nrows=1, figsize=(10,9))
        figc, axc = plt.subplots(ncols=1, nrows=1, figsize=(10,9))
        #figbb, axbb = plt.subplots(ncols=2, nrows=2, figsize=(10,9))
        col=['k','r','b']
        for key, val in tot_resu.items():
            #self.Plot_Sims(axbb,val)
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
            self.Plot_Eff(axa,val,sel,'z',dict_ana[key].colorfig)
            self.Plot(axc,sel,'salt2.CovColorColor','z',dict_ana[key].colorfig)
            print 'Number of Events',val[0],val[1],len(val)
            for vval in np.arange(0.,0.7,0.1):
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

        num_bins=36
        range=[0.0,1.2]
        
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

    def Plot_Eff(self,axc,sela,selb,varname,color):
        tot_label=[]
        #figc, axc = plt.subplots(ncols=1, nrows=1, figsize=(10,9))
        self.Plot_Eff_Indiv(axc,sela,selb,varname,tot_label,'test',color)
        


    def Plot_Eff_Indiv(self,axc,sela,selb,varname,tot_label,ll,color):

        bin_center, ratio, ratio_err,norm,norm_err= self.Histo_ratio(sela,selb,varname)
        #tot_label.append(axc.errorbar(bin_center,ratio, yerr=ratio_err,marker=marker, mfc=colors[key], mec=colors[key], ms=8, linestyle=myfmt[i],color='k',label=ll))
        tot_label.append(axc.errorbar(bin_center,ratio, yerr=ratio_err,marker='.', mfc='red', mec='red', ms=8, linestyle='-',color=color,label=ll))

    def Plot(self, axc,sel,vary,varx,color):
        print 'there',sel[varx],sel[vary]
        axc.plot(sel[varx],np.sqrt(sel[vary]),color+'.')

dict_ana={}

dict_ana['Mean_Obs_Faint_T0_0']=Id(thedir='Mean_Obs_T0_0',fieldname='WFD',fieldid=309,X1=-2.0,Color=0.2,season=0,colorfig='k')
dict_ana['Mean_Obs']=Id(thedir='Mean_Obs',fieldname='WFD',fieldid=309,X1=-999.,Color=-999.,season=0,colorfig='r')

dict_ana['Rolling_Faint']=Id(thedir='WFD_Rolling_noTwilight',fieldname='WFD',fieldid=309,X1=-2.0,Color=0.2,season=1,colorfig='b')
dict_ana['Mean_Obs_Faint_newrefs']=Id(thedir='Mean_Obs_newrefs',fieldname='WFD',fieldid=309,X1=-2.0,Color=0.2,season=0,colorfig='k')

Ana_Simu(dict_ana)
