import cPickle as pkl
import numpy as np
import glob
import matplotlib.pyplot as plt
from optparse import OptionParser
from ID import *
import os
from astropy.table import vstack,Table


def Select(varname,val,valmin,valmax):
    idx=(val[varname]>=valmin)&(val[varname]<valmax)
    return val[idx]

class Read_Observations: 
    def __init__(self, dirobs, fieldname, fieldid, season,X1, Color,zmin,zmax,fichnum):
        
        thedir='Prod_LC'
        self.observations=None
        self.zmin=zmin
        self.zmax=zmax
        dirmeas=thedir+'/'+dirobs+'_Obs/'+str(fieldid)+'/Season_'+str(season)
        dir_SN=thedir+'/'+dirobs+'/'+str(fieldid)+'/Season_'+str(season)
        filename=dirmeas+'/'+fieldname+'_'+str(fieldid)+'_'+str(self.zmin)+'_'+str(self.zmax)+'_X1_'+str(X1)+'_C_'+str(Color)+'_'+str(fichnum)+'.pkl'
        filesn=dir_SN+'/'+fieldname+'_'+str(fieldid)+'_'+str(self.zmin)+'_'+str(self.zmax)+'_X1_'+str(X1)+'_C_'+str(Color)+'_'+str(fichnum)+'.pkl'
        pkl_file = open(filename,'rb')
        self.observations=pkl.load(pkl_file)
        pkl_file = open(filesn,'rb')
        self.SN=pkl.load(pkl_file)


    """
    def __init__(self, dirobs, fieldname, fieldid, season,X1, Color,zmin=0.,zmax=0.05):
        
        thedir='Prod_LC'
        self.observations=None
        self.zmin=zmin
        self.zmax=zmax

        
        print 'hee',dirobs
        dirmeas=thedir+'/'+dirobs+'/'+str(fieldid)+'/Season_'+str(season)
        sum_file=dirmeas+'/'+fieldname+'_'+str(fieldid)+'_'+str(self.zmin)+'_'+str(self.zmax)+'_X1_'+str(X1)+'_C_'+str(Color)+'_all.pkl'

        if os.path.exists(sum_file):
            pkl_file = open(sum_file,'rb')
            loaded_file=pkl.load(pkl_file)
                #print 'done',key,tot_resu[key],tot_resu[key].dtype,sum_file
                #idx = loaded_file['z']>=self.zmin and loaded_file['z']<self.zmax
            self.observations=loaded_file
            
        else:
            
            files = glob.glob(dirmeas+'/'+fieldname+'_'+str(fieldid)+'_'+str(self.zmin)+'_'+str(self.zmax)+'_X1_'+str(X1)+'_C_'+str(Color)+'*.pkl')
            print 'hello',dirmeas
            for fi in files:
                pkl_file = open(fi,'rb')
                print 'loading',fi
                theload=pkl.load(pkl_file)
                if self.observations is None:
                    self.observations=theload
                else:
                        #print tot_resu[key].dtype,
                    self.observations=vstack([self.observations,theload])

                pkl_out = open(sum_file,'wb')
                
                pkl.dump(self.observations, pkl_out)
                
                pkl_out.close()


       
        #print self.observations
       """



parser = OptionParser()
parser.add_option("-f", "--fieldname", type="string", default='DD', help="filter [%default]")
parser.add_option("-i", "--fieldid", type="int", default=290, help="filter [%default]")
parser.add_option("-x", "--stretch", type="float", default=-999., help="filter [%default]")
parser.add_option("-c", "--color", type="float", default=-999., help="filter [%default]")
parser.add_option("-d", "--dirmeas", type="string", default="DD_Obs", help="filter [%default]")

opts, args = parser.parse_args()

fieldid=opts.fieldid
fieldname=opts.fieldname
dirmeas=opts.dirmeas

#X1=0.0
#Color=0.0

X1=opts.stretch
Color=opts.color
season=4
num=10
print 'alors',opts.dirmeas
T0_ref=51074.5070051
for num in range(200):
    read=Read_Observations(dirmeas,fieldname, fieldid, season,X1, Color,0.,0.05,num)

    obs=read.observations
    sn=read.SN

    print obs.dtype
    """
    idx = (obs['T0']>=51077.)&(obs['T0']<=51078.)
    idx = obs['x1'] == -3.7864694
    """
    
    print len(obs),np.min(obs['T0']),np.max(obs['T0']),set(obs['T0'])
    print len(sn),np.min(sn['T0']),np.max(sn['T0']),sn['T0']
    idx=((sn['T0']/T0_ref)>0.999999)&((sn['T0']/T0_ref)<1.000001)
    if len(sn[idx])>0:
        print 'yeyeyeyeyeyey',sn['T0'][idx],sn['phase_first'][idx],sn['phase_last'][idx],sn['X1'][idx],sn['Color'][idx]
        break
