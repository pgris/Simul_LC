from Generate_LC import *
from Observations import *
import numpy as np
from Telescope import *
#from Fit_LC import *
from Generate_Single import Generate_Single_LC
import time
import multiprocessing
import cPickle as pkl
from optparse import OptionParser
import os

parser = OptionParser()

parser.add_option("-N", "--nevts", type="int", default=100, help="filter [%default]")
parser.add_option("-f", "--fieldname", type="string", default="WFD", help="filter [%default]")
parser.add_option("-n", "--fieldid", type="int", default=309, help="filter [%default]")
parser.add_option("-s", "--season", type="int", default=1, help="filter [%default]")
parser.add_option("-z", "--zmin", type="float", default=0.0, help="filter [%default]")
parser.add_option("-Z", "--zmax", type="float", default=0.1, help="filter [%default]")
parser.add_option("-x", "--stretch", type="float", default=2.0, help="filter [%default]")
parser.add_option("-c", "--color", type="float", default=-0.2, help="filter [%default]")
parser.add_option("-t", "--sntype", type="string", default="Ia", help="filter [%default]")
parser.add_option("-d", "--dirmeas", type="string", default="None", help="filter [%default]")
parser.add_option("-r", "--T0random", type="string", default="yes", help="filter [%default]")

#parser.add_option("-r", "--rolling", type="int", default=0, help="filter [%default]")
#parser.add_option("--nrolling", type="int", default=0, help="filter [%default]")
#parser.add_option("--merge_factor", type="int", default=0, help="filter [%default]")

opts, args = parser.parse_args()

telescope=Telescope(atmos=True,airmass=1.2)

fieldid=opts.fieldid
num_season=opts.season
N_sn=opts.nevts
zmin=opts.zmin
zmax=opts.zmax
T0random=opts.T0random
X1=opts.stretch
Color=opts.color

filename='../Ana_Cadence/OpSimLogs/'+opts.dirmeas+'/Observations_'+opts.fieldname+'_'+str(fieldid)+'.txt'

#filename='OpSimLogs/WFD/Observations_WFD_'+str(fieldid)+'.txt'

myobs=Observations(fieldid=fieldid, filename=filename)


outdir='Prod_LC/'+opts.dirmeas+'/Season_'+str(num_season)
if not os.path.exists(outdir):
    os.makedirs(outdir)

print(len(myobs.seasons))

if X1==-999. or Color==-999.:
    if zmax < 0.1:
        X1_Color_npzfile = np.load('Dist_X1_Color_low_z.npz','r')
    else:
        X1_Color_npzfile = np.load('Dist_X1_Color_high_z.npz','r')
        


#num_season=1

myseason=myobs.seasons[num_season]

print myseason.dtype,np.min(myseason['mjd']),np.max(myseason['mjd'])

iddx=myseason['band']=='LSSTPG::g'
mysel=myseason[iddx]

min_season=np.min(mysel['mjd'])
max_season=np.max(mysel['mjd'])

diff=max_season-min_season



n_multi=5
n_batch=N_sn/n_multi

time_begin=time.time()

for i in range(0,n_batch):
    result_queue = multiprocessing.Queue()
#process=[]
    #print 'processing main',i
    name_for_output=opts.fieldname+'_'+str(fieldid)+'_'+str(zmin)+'_'+str(zmax)+'_X1_'+str(X1)+'_C_'+str(Color)+'_'+str(i)
    for j in range(0,n_multi):
        
        if T0random == 'yes':
            T0 = np.random.uniform(min_season,max_season)
        else:
            T0=0
            
        z=np.random.uniform(zmin,zmax)
        X1_val=X1
        Color_val=Color
        if X1 == -999.:
            X1_val=np.random.choice(X1_Color_npzfile['x1_vals'],1,p=X1_Color_npzfile['x1_weights'])[0]
        if Color==-999.:
            Color_val=np.random.choice(X1_Color_npzfile['c_vals'],1,p=X1_Color_npzfile['c_weights'])[0]

        #print 'Processing',j,T0,z
        p=multiprocessing.Process(name='Subprocess-'+str(i),target=Generate_Single_LC,args=(z,T0,X1_val,Color_val,myseason,telescope,j,result_queue))
    #process.append(p)
        p.start()
    
    resultdict = {}
    for j in range(0,n_multi):
        resultdict.update(result_queue.get())

    for p in multiprocessing.active_children():
        p.join()

    tot_obs=None
    for j in range(0,n_multi):
        #print 'hello there',resultdict[j].dtype
        if tot_obs is None:
            tot_obs=resultdict[j]
            #print 'there',tot_obs
        else:
            #print 'careful',tot_obs.dtype
            #print 'carefulb',resultdict[j].dtype
            tot_obs=np.vstack((tot_obs,resultdict[j]))

    pkl_file = open(outdir+'/'+name_for_output+'.pkl','wb')
    
    pkl.dump(tot_obs, pkl_file)
    
    pkl_file.close()
print 'total elapse time',time.time()-time_begin
