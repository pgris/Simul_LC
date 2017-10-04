from optparse import OptionParser
import numpy as np
import os

parser = OptionParser()
parser.add_option("-z", "--zmin", type="float", default=0.1, help="filter [%default]")
parser.add_option("-Z", "--zmax", type="float", default=0.2, help="filter [%default]")
parser.add_option("--zstep", type="float", default=0.1, help="filter [%default]")
parser.add_option("-N", "--nevts", type="int", default=10, help="filter [%default]")
parser.add_option("-m", "--model", type="string", default='salt2-extended', help="filter [%default]")
parser.add_option("-v", "--version", type="string", default='', help="filter [%default]")
parser.add_option("-f", "--fieldname", type="string", default='WFD', help="filter [%default]")
parser.add_option("-i", "--fieldid", type="int", default=309, help="filter [%default]")
parser.add_option("-s", "--season", type="int", default=-1, help="filter [%default]")
parser.add_option("-t", "--sntype", type="string", default='Ia', help="filter [%default]")
#parser.add_option("-d", "--dbFile", type="string",default='None', help="dbFile to process [%default]")
parser.add_option("-x", "--stretch", type="float", default=2.0, help="filter [%default]")
parser.add_option("-c", "--color", type="float", default=-0.2, help="filter [%default]")
parser.add_option("-d", "--dirmeas", type="string", default="None", help="filter [%default]")
parser.add_option("-r", "--T0random", type="string", default="yes", help="filter [%default]")
parser.add_option("--zrandom", type="string", default="yes", help="filter [%default]")

opts, args = parser.parse_args()

nevts=opts.nevts
fieldname=opts.fieldname
fieldid=opts.fieldid
season=opts.season
sntype=opts.sntype
stretch=opts.stretch
color=opts.color
dirmeas=opts.dirmeas
T0random=opts.T0random
zstep=opts.zstep
zrandom=opts.zrandom

range_loop=np.arange(opts.zmin,opts.zmax,zstep)
#print range_loop
nvals=int(round((opts.zmax-opts.zmin)/zstep))

if len(range_loop) > nvals:
    range_loop=range_loop[:-1]

for z in range_loop:
    #print 'here',z
    cmd='python batch.py --zmin '+str(z)+' --zmax '+str(z+zstep)+' --nevts '+str(nevts)+' --fieldname '+fieldname+' --fieldid '+str(fieldid)+' --season '+str(season)+' --sntype '+sntype+' --stretch '+str(stretch)+' --color '+str(color)+' --dirmeas '+dirmeas+' --T0random '+T0random+' --zrandom '+zrandom
    #print cmd
    os.system(cmd)
