import numpy as np
import os
from Observations import *
import pylab as plt
from scipy.stats import truncnorm
from optparse import OptionParser

class Make_Obs:
    def __init__(self, cadence=1.,random=False,expTime=dict(zip([b for b in 'grizy'],[300.,600.,600.,600.,780.,600.])),Ra=-1.,Dec=-1.,fieldtype='DD',fieldid=290):

        self.cadence=cadence
        self.random=random
        self.bands=expTime.keys()

        outdir='../Ana_Cadence/OpSimLogs/Mean_Obs_newrefs'
        if not os.path.exists(outdir):
            os.makedirs(outdir)

        self.fieldid=fieldid
        self.name_output=outdir+'/Observations_'+fieldtype+'_'+str(fieldid)+'.txt'

#m5_lim={'u':23.61,'g':24.83,'r':24.35,'i':23.88,'z':23.30,'y':22.43}
        self.m5_lim={'u': 23.103921,'g': 24.2292315,'r': 23.8573,'i': 23.4196475,'z': 22.792551,'y': 21.5026765}
        
        self.expTime=expTime
        self.rawSeeing=-1
        
        
        self.seeing={'u':0.92,'g':0.87,'r':0.83,'i':0.80,'z':0.78,'y':0.76}
        self.moonphase=0
        self.sky={'u':22.95,'g':22.24,'r':21.20,'i':20.47,'z':19.60,'y':18.63}
        self.kAtm=-1
        self.airmass=1.2
        self.Nexp={}
        
        for key, vals in self.expTime.items():
            self.m5_lim[key]=self.m5_lim[key]+1.25*np.log10(vals/30.)
            self.Nexp[key]=int(vals/30.)
            
        self.Ra=Ra
        self.Dec=Dec
        delay=60/24./3600.
        self.decal={'u':0.,'g':delay,'r':2.*delay,'i':3.*delay,'z':4.*delay,'y':5.*delay}

        self.outputfile  = open(self.name_output,'wb')
        legend=['band','mjd','exptime','rawSeeing','FWHMeff','moon_frac','sky','kAtm','airmass','m5sigmadepth','Nexp','Ra','Dec']
        for leg in legend:
            lego=leg+' : '
            if leg == 'FWHMeff':
                lego='seeing : [was FWHMeff]'
            self.outputfile.write('# '+lego+'\n')
        self.outputfile.write('# end\n')

        

        """
        scale = 2.
        range = 4.5
        size = 1000

        X = truncnorm(a=-range/scale, b=+range/scale, scale=scale).rvs(size=size)
        X = X.round().astype(int)
        
        #bins = 2 * range + 1
        self.nvals=(int(self.mjd_max)-int(self.mjd_max))/self.cadence
        self.deltaT =[x + self.cadence for x in X]
        self.Numobs=[self.nvals + self.cadence for x in X]
        bins=np.max(X)-np.min(X)
        plt.hist(X, bins,range=[np.min(X),np.max(X)])
        print X

        self.Numobs, self.Binsobs = np.histogram(X,bins=bins,range=[np.min(X),np.max(X)])
        #self.Numobs=hista
        #print hista,bin_edgesa
        plt.show()
        """
        
    @property
    def cadences(self):
        return self.cadences

    def __call__(self,shift=0.):
            
        
        #sampling=self.cadence
        z=1.
        self.mjd_min=-20*(1.+z)+shift
        self.mjd_max=60*(1.+z)+shift


        mjd_tot = np.arange(self.mjd_min, self.mjd_max, self.cadence)

        if self.random:
            mjds=self.Get_mjds(self.mjd_min,self.mjd_max,self.cadence)
            #print 'yes pal',mjds,len(mjds)
        else:
            mjds=mjd_tot

        for b in self.bands:
            for mjd in mjds:
                toprint = 'LSSTPG::'+b+' '
                toprint+=str(format(mjd+self.decal[b],'3.8f'))+' '
                toprint+=str(self.expTime[b])+' '
                toprint+=str(format(self.rawSeeing,'.7f'))+' '
                toprint+=str(format(self.seeing[b],'.7f'))+' '
                toprint+=str(format(self.moonphase,'.7f'))+' '
                toprint+=str(format(self.sky[b],'.7f'))+' '
                toprint+=str(self.kAtm)+' '
                toprint+=str(format(self.airmass,'.7f'))+' '
                toprint+=str(format(self.m5_lim[b],'.7f'))+' '
                toprint+=str(self.Nexp[b])+' '
                toprint+=str(format(self.Ra,'.7f'))+' '
                toprint+=str(format(self.Dec,'.7f'))
                #print toprint
                self.outputfile.write(toprint+'\n')

        #self.Ana_Sim(self.fieldid,self.name_output)

    def End(self):
        self.outputfile.close()

    def Get_mjds(self,mjd_min,mjd_max,cadence_ref):
    
        mjd=[]
        cadence=np.random.randint(cadence_ref,cadence_ref+2)
        cadence=cadence_ref
        ndays=int(mjd_max)-int(mjd_min)
        nbins=ndays/cadence
        #print 'hello pal',mjd_max-mjd_min,cadence
    #print 'hello',mjd_min,mjd_max,cadence,nvals
        #for i in range(10*nvals):
        #mjd.append(mjd_min)
        
        
        """
        r = truncnorm.rvs(3., 2, size=1000,scale=3.)

        plt.hist(r, normed=True, histtype='stepfilled', alpha=0.2)
        plt.show()
        """

        mjd.append(mjd_min)
        #mjd.append(int(mjd_max))
        nobs=(int(mjd_max)-int(mjd_min))/cadence_ref
        """
        while len(mjd) <= nbins:
            #date=np.random.uniform(mjd_min,mjd_max)
            date=np.random.randint(int(mjd_min),int(mjd_max))
            
            #if date not in mjd and self.Diff(date,mjd) > cadence-jitter and self.Diff(date,mjd) <= cadence+jitter+1:
            #if date not in mjd and self.Diff(date,mjd) >= cadence-jitter and self.Diff(date,mjd) <= cadence+jitter+1:
                #print 'hello',self.Diff(date,mjd),cadence,date,mjd
            #if date not in mjd and self.Diff(date,mjd) == cadence:
            if date not in mjd and np.mod(float(date-mjd_min),float(cadence)) <=0.:
                #print 'add',date
                mjd.append(date)
            #if len(mjd)==nvals:
                #break
        """
        #nobs=np.random.choice(self.Binsobs[:-1],1,p=self.Numobs)[0]
        #while len(mjd) <= nobs-1:
        #120->123 : cadence=3
        #myweights=dict(zip([i for i in range(1,int(cadence)+5)],[0.,0.,1.,0.0,0.0,0.0,0.0])) #120
        #myweights=dict(zip([i for i in range(1,int(cadence)+5)],[0.,0.1,0.8,0.1,0.0,0.0,0.0])) #121
        #myweights=dict(zip([i for i in range(1,int(cadence)+5)],[0.,0.1,0.7,0.05,0.05,0.05,0.05]))#122
        #myweights=dict(zip([i for i in range(1,int(cadence)+5)],[0.,0.,0.7,0.05,0.05,0.1,0.1])) #123
        #124->127 : cadence=4
        if self.fieldid in [120,124,128,132,136,140,144,148]:
            myweights=dict(zip([i for i in range(int(cadence)-2,int(cadence)+5)],[0.,0.,1.,0.0,0.0,0.0,0.0])) #124
        if self.fieldid in [121,125,129,133,137,141,145,149]:
            myweights=dict(zip([i for i in range(int(cadence)-2,int(cadence)+5)],[0.,0.1,0.8,0.1,0.0,0.0,0.0])) #121
        if self.fieldid in [122,126,130,134,138,142,146,150]:   
            myweights=dict(zip([i for i in range(int(cadence)-2,int(cadence)+5)],[0.,0.1,0.7,0.05,0.05,0.05,0.05]))#122
        if self.fieldid in [123,127,131,135,139,143,147,151]:
            myweights=dict(zip([i for i in range(int(cadence)-2,int(cadence)+5)],[0.,0.,0.7,0.05,0.05,0.1,0.1])) #123
        
        

        #print 'sum',np.sum([vals for key,vals in myweights.items()])
        """
        for val in range(1,cadence+2):
            xvals.append(val)
            weight.append(myweights[val])
        """
        while mjd_max-mjd[-1]>= -1:
            #delta_time=np.random.choice(np.range(0,2),1,p=self.Delta_T)[0]
            #delta_time = np.random.normal(self.cadence,0.001, 10)[5]
            delta_time=np.random.choice([key for key,vals in myweights.items()],1,p=[vals for key,vals in myweights.items()])[0]
            diff=mjd_max-mjd[-1]-delta_time
            #print 'hhh',delta_time,mjd[-1],diff
            #if delta_time>=1. and self.Check_Outliers(cadence,delta_time,mjd):
            #print 'hhh', mjd_max-mjd[-1]-delta_time,delta_time
            if mjd_max-mjd[-1]-delta_time >= 0.:
                if delta_time>=1. :
                    mjd.append(mjd[-1]+delta_time)
            else:
                break

        mjd.append(mjd_max)

        mjd=np.sort(mjd)
        #mjd=mjd[:-1]
        
        """
        print mjd,nbins,int(mjd_min),int(mjd_max)
        plt.hist(mjd,bins=ndays,range=[np.min(mjd),np.max(mjd)])
        plt.show()
        
        print 'pool',np.mean([io-jo for jo,io in zip(mjd[:-1], mjd[1:])]),nbins,cadence,len(mjd)
        
        jitter=1.
        mjd=[float(val)+np.random.randint(-jitter,jitter) for val in mjd]
        #mjd=[val for val in mjd]
        diffs=[io-jo for jo,io in zip(mjd[:-1], mjd[1:])]
        if [x <= 0 for x in diffs].count(True):
            print 'Big Problem'
        """
        diffs=[io-jo for jo,io in zip(mjd[:-1], mjd[1:])]
        mean_cadence=np.mean(diffs)
        #print 'hello',mjd,len(mjd)
        
        #if mean_cadence < cadence_ref-0.1 or mean_cadence> cadence_ref+0.1:
        #if not self.Check_Outliers(cadence,-1,mjd):
            #return self.Get_mjds(mjd_min,mjd_max,cadence_ref)
        
        #else:
            #print 'go',mjd
        
        return mjd

    def Check_Outliers(self,cadence, deltaT, mjd):

        diffs=[io-jo for jo,io in zip(mjd[:-1], mjd[1:])]
        if deltaT > 0.:
            diffs.append(deltaT)

        sela=[dd for dd in diffs if dd >= cadence+1. and dd < cadence+2. ]
        #print 'hello',sela,len(sela)
        selb=[dd for dd in diffs if dd >= cadence+2. ]
        #print 'hello',selb,len(selb)
        rata=float(len(sela))/float(len(diffs))
        ratb=float(len(selb))/float(len(diffs))
                                    

        return rata<0.05 and ratb<0.20

    def Diff(self,val, listv):
        ls=list(listv)
        ls.append(val)
        ls=np.sort(ls)
        diff=[io-jo for jo,io in zip(ls[:-1], ls[1:])]
        #print 'hello',np.min(diff)
        return np.min(diff)

    def Ana_Sim(self, fieldid,filename):

        myobs=Observations(fieldid=fieldid, filename=filename)
        season=0

        r=[]
        for b in self.bands:
            idx = myobs.seasons[season]['band']=='LSSTPG::'+b
            sel=myobs.seasons[season][idx]
            diff=[io-jo for jo,io in zip(sel['mjd'][:-1], sel['mjd'][1:])]
            print b,np.mean(diff),np.std(diff)
            r.append((b,np.mean(diff),np.std(diff)))

        self.cadences=np.rec.fromrecords(r,names=['band','mean_cadence','rms_cadence'])


parser = OptionParser()

parser.add_option("-n", "--fieldid", type="int", default=309, help="filter [%default]")
parser.add_option("-c", "--cadence", type="float", default=1, help="filter [%default]")
opts, args = parser.parse_args()

fieldid=opts.fieldid
cadence=opts.cadence

for field in range(fieldid,fieldid+4):
    my_run=Make_Obs(cadence=cadence,random=True,fieldid=field)

    tot_cad=None
    for i in range(1):
        my_run(float(i)*250)
    #print 'eeee',i
    #my_run(0.)
    """
    print my_run.cadences
    if tot_cad is None:
        tot_cad=my_run.cadences
    else:
        tot_cad=np.vstack([tot_cad,my_run.cadences])
    """
    my_run.End()



    myobs=Observations(fieldid=fieldid, filename=my_run.name_output)

    print len(myobs.seasons)

    res=[]
    for i in range(len(myobs.seasons)):
        season=myobs.seasons[i]
        for b in 'grizy':
            idx = season['band'] == 'LSSTPG::'+b
            sel= season[idx]
            diff=[io-jo for jo,io in zip(sel['mjd'][:-1], sel['mjd'][1:])]
            print i,b,np.mean(diff),np.std(diff)
            res.append((b,np.mean(diff),np.std(diff)))

        tot_cad=np.rec.fromrecords(res,names=['band','mean_cadence','rms_cadence'])

        print 'allo',len(tot_cad)
        for b in 'g':
            idx = tot_cad['band']==b
            print b,len(tot_cad[idx])
            sel=tot_cad[idx]
            figa, axa = plt.subplots(ncols=1, nrows=2, figsize=(10,9))

            axa[0].hist(sel['mean_cadence'])
            axa[1].hist(sel['rms_cadence'])

    #print sel['mean_cadence']
    plt.show()

