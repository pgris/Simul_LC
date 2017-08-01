import numpy as np

class Observations:
    def __init__(self,fieldid,ra=-1,dec=-1,filename='',nseasons=10):

        self.fieldid=fieldid

        if filename != '':
            data=self.Load(filename)
            data.sort(order='mjd')
            self.seasons=self.Get_Seasons(data)
            self.Ra=np.unique(data['Ra'])[0]
            self.Dec=np.unique(data['Dec'])[0]
        else:
            self.Ra=ra
            self.Dec=dec
            self.seasons={}
            for i in range(nseasons):
                self.seasons[i]=None

        #print 'Nseasons',len(self.seasons)
        #print data

    def Load(self,filename):

        sfile=open(filename,'r')
        varname=[]
        r=[]
        for line in sfile.readlines():
            if line[0] == '#':
                varname.append(line.split(' ')[1])
            else:
                tofill=[]
                thesplit=line.strip().split(' ')
                tofill.append(thesplit[0])
                for i in range(1,len(thesplit)):
                    tofill.append(float(thesplit[i]))
                r.append(tuple(tofill))
        return np.rec.fromrecords(r,names=varname)
            
    def Get_Seasons(self,data,season_length=110):
        
        thediff=data['mjd'][1:]-data['mjd'][:-1]
        idx,=np.where(thediff > season_length)
        lidx=[val+1 for val in list(idx)]
        
        lidx.insert(0,0)
        lidx.append(len(data['mjd']))
        #print 'ay',lidx
        
        seasons={}
        
        for i in range(len(lidx)-1):
            #print lidx[i],lidx[i+1]
            seasons[i]=data[lidx[i]:lidx[i+1]]

        return seasons
            



    
