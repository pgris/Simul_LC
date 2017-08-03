from Generate_LC import *
import multiprocessing
from astropy.table import vstack,Table
from Fit_LC import *
import cPickle as pkl

class Generate_Single_LC:
    def __init__(self,z,T0,X1,Color,obs,telescope,inum,output_q):
        
        params={}
        params['z']=z
        params['DayMax']=T0
        params['X1']=X1
        params['Color']=Color

        
        self.telescope=telescope
        self.z=z
        self.T0=T0
        self.X1=X1
        self.Color=Color
        self.outdict={}
        self.outdict=params
        self.outdict['status']='unkown'
        self.outdict['fit']=None
        self.outdict['mbsim']=-999.
        self.outdict['observations']=None

        self.mysn=Generate_LC(params,telescope=self.telescope)
       
        self.obs=obs
        self.bands=[b[-1:] for b in np.unique(self.obs['band'])]

        """
        print 'hello',self.bands,self.z
        for b in self.bands:
            idx = obs['band'] == 'LSSTPG::'+b
            sel=obs[idx]
            print 'ttttt',b,np.median(sel['m5sigmadepth'])
        """

        self.bands_rest = 'grizy'
        self.Gen_LC()
        self.dict_quality={}
        self.Get_Quality_LC()
        Nmeas=self.dict_quality['all'][0]+self.dict_quality['all'][1]
        if Nmeas >= 5:
            self.outdict['status']='go_fit'
            self.outdict['fit_status']='unknow'
            self.Fit_LC()
        else:
           self.outdict['status']='no_obs'
           self.outdict['fit_status']='unknow' 
        #out_q.put({irun :(params,self.tot_obs)})
        

        """
        pkl_file = open('Prod_LC/'+name_for_output+'.pkl','wb')
        
        pkl.dump(self.outdict, pkl_file)
            
        pkl_file.close()
        """
        output_q.put({inum : self.Summary()})
        
    def Gen_LC(self):

        result_queue = multiprocessing.Queue()
        process=[]

        
        for b in self.bands_rest:
            idx= self.obs['band']=='LSSTPG::'+b
            sel=self.obs[idx]
            p=multiprocessing.Process(name='Subprocess-'+b,target=self.mysn,args=(sel['mjd'],sel['airmass'],sel['m5sigmadepth'],b,sel['exptime'],result_queue))
            process.append(p)
            p.start()

        resultdict = {}
        for b in self.bands_rest:
            resultdict.update(result_queue.get())

        for p in multiprocessing.active_children():
            p.join()
    
        self.tot_obs=None
        for b in self.bands_rest:
            if resultdict[b][1] is not None:
                if self.tot_obs is None:
                    self.tot_obs=resultdict[b][1]
                else:
                    self.tot_obs=vstack([self.tot_obs,resultdict[b][1]])

        self.outdict['observations']=self.tot_obs

    def Fit_LC(self):
        
        
        myfit=Fit_LC(z=self.z,telescope=self.telescope,Plot=False)
        res,fitted_model,mbfit,fit_status=myfit(self.tot_obs)

        #print 'hello',res,fitted_model,mbfit,fit_status
        if fit_status == 'ok':
            #print 'hello',res,fitted_model
            self.outdict['sncosmo_res']=res
        
            self.outdict['sncosmo_fitted']={}
            for i,par in enumerate(fitted_model.param_names):
                self.outdict['sncosmo_fitted'][par]=fitted_model.parameters[i]
                #self.outdict['fitted_model']=fitted_model
        
            self.outdict['mbfit']=mbfit
            self.outdict['fit_status']='fit_ok'
        
        #print 'fit',myfit.sigma_c

        #now get number of observations before and after peak
        """
        for b in self.bands:
            idx= self.obs['band']=='LSSTPG::'+b
            sel=self.obs[idx]
            idb=np.logical_and(sel['mjd'] <= self.T0,sel['mjd'] > self.T0-20)
            ida=np.logical_and(sel['mjd']> self.T0,sel['mjd'] <= self.T0+40)
            print b,len(sel[idb]),len(sel[ida]),len(sel[idb])+len(sel[ida])

        """
    def Get_Quality_LC(self):

        #estimate the number of LC points (5 sigma) before and after T0 - observer frame
        for band in self.bands_rest:
            self.dict_quality[band]=(0,0)
        self.dict_quality['all']=(0,0)
        self.dict_quality['phase']=(0.0,0.0)
        
        if self.tot_obs is not None:
            obs_sel=self.tot_obs[np.where(np.logical_and(self.tot_obs['flux']/self.tot_obs['fluxerr']>5.,self.tot_obs['flux']>0.))]
        
        #print self.Get_nums(obs_sel)
            if len(obs_sel) > 0:
                n_bef_tot=0
                n_aft_tot=0
                for band in self.bands_rest:
                    idx=obs_sel['band']=='LSST::'+band
                    n_bef, n_aft=self.Get_nums(obs_sel[idx])
            #print 'eheh',band,n_bef,n_aft
                    n_bef_tot+=n_bef
                    n_aft_tot+=n_aft
                    self.dict_quality[band]=(n_bef,n_aft)
                self.dict_quality['all']=(n_bef_tot,n_aft_tot) 

                obs_sel.sort('time')
                phase_first=(obs_sel[0]['time']-self.T0)/(1.+self.z)
                phase_last=(obs_sel[len(obs_sel)-1]['time']-self.T0)/(1.+self.z)
                self.dict_quality['phase']=(phase_first,phase_last)
           

        #print phase_first,phase_last

    def Get_nums(self, sel):
        
        idxa=np.logical_and(sel['time'] <= self.T0,sel['time'] > self.T0-20)
        idxb=np.logical_and(sel['time']> self.T0,sel['time'] <= self.T0+40)

        return len(sel[idxa]),len(sel[idxb])


    def Summary(self):
        
        resu={}
 
        resu['mbsim']=self.mysn.mbsim
        resu['sn_type']=self.mysn.sn_type
        resu['X0']=self.mysn.X0
        resu['z']=self.z
        resu['T0']=self.T0
        resu['X1']=self.X1
        resu['Color']=self.Color
        for band in self.bands_rest:
            resu['N_bef_'+band]=self.dict_quality[band][0]
            resu['n_aft_'+band]=self.dict_quality[band][1]
        resu['phase_first']=self.dict_quality['phase'][0]
        resu['phase_last']=self.dict_quality['phase'][1]
        resu['N_bef']=self.dict_quality['all'][0]
        resu['N_aft']=self.dict_quality['all'][1]
        resu['status']=self.outdict['status']
        resu['fit_status']=self.outdict['fit_status']

        resu['salt2.X0']=-999.
        resu['salt2.X1']=-999.
        resu['salt2.Color']=-999.
        resu['salt2.CovX1X1']=-999.
        resu['salt2.CovColorColor']=-999.
        resu['salt2.CovX0X1']=-999.
        resu['salt2.CovColorX0']=-999.
        resu['salt2.CovColorX1']=-999.
        resu['mbfit']=-999.

        if self.outdict['status']=='go_fit':
            #print 'rg',self.outdict['sncosmo_fitted']
            if self.outdict['fit_status']=='fit_ok':
                corr={}
                for i,pal in enumerate(self.outdict['sncosmo_res']['vparam_names']):
                    corr[pal]=i
                #print 'hhe',i,pal

                resu['salt2.X0']=self.outdict['sncosmo_fitted']['x0']
                resu['salt2.X1']=self.outdict['sncosmo_fitted']['x1']
                resu['salt2.Color']=self.outdict['sncosmo_fitted']['c']
                resu['salt2.CovX1X1']=self.outdict['sncosmo_res']['covariance'][corr['x1']][corr['x1']]
                resu['salt2.CovColorColor']=self.outdict['sncosmo_res']['covariance'][corr['c']][corr['c']]
                resu['salt2.CovX0X1']=self.outdict['sncosmo_res']['covariance'][corr['x0']][corr['x1']]
                resu['salt2.CovColorX0']=self.outdict['sncosmo_res']['covariance'][corr['c']][corr['x0']]
                resu['salt2.CovColorX1']=self.outdict['sncosmo_res']['covariance'][corr['c']][corr['x1']]
                resu['mbfit']=self.outdict['mbfit']

        
        #print 'hello',resu.values(),resu.keys()
        return np.rec.fromarrays(tuple([res for res in resu.values()]),names=[key for key in resu.keys()])
