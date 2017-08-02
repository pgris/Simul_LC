from SN_Rate import *
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter

rate=SN_Rate()

zz,rate_one,nsn=rate(0.001,1.2,0.01)

rateb=SN_Rate(rate='Perret')
zzb,rate_two,nsnb=rateb(0.001,1.2,0.01)

print zz, nsn

figbb, axbb = plt.subplots(ncols=1, nrows=1, figsize=(10,9))
axbb.plot(zz,1.e4*rate_one,label='Ripoche rate')
axbb.plot(zzb,1.e4*rate_two,label='Perret rate')
axbb.set_yscale('log')
axbb.set_xlabel('Redshift')
axbb.set_ylabel('10$^{-4}$ SNe Ia yr$^{-1}$ Mpc$^{-3}$ h$^{3}$$_{70}$')
axbb.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
plt.legend(loc='best')

figc, axc = plt.subplots(ncols=1, nrows=1, figsize=(10,9))
axc.plot(zz,nsn,label='Ripoche rate')
axc.plot(zzb,nsnb,label='Perret rate')
axc.set_yscale('log')
axc.set_xlabel('Redshift')
axc.set_ylabel('N(SNe Ia) - one year (effective search duration=0.5)')
#axc.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
plt.legend(loc='best')

plt.show()
