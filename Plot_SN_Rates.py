from SN_Rate import *
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import numpy as np

survey_area= 9.6*2293
duration = 0.5
zmin=0.0
zmax=0.4
zstep=0.05

rate=SN_Rate(rate='Dilday',survey_area=survey_area,duration=duration)

zz,rate_one,err_rate,nsn,err_nsn=rate(zmin,zmax,zstep)

rateb=SN_Rate(rate='Perrett',survey_area=survey_area,duration=duration)
zzb,rate_two,err_rateb,nsnb,err_nsnb=rateb(zmin,zmax,zstep)

print zz, nsn,nsnb

figbb, axbb = plt.subplots(ncols=1, nrows=1, figsize=(10,9))
axbb.errorbar(zz,1.e4*rate_one,yerr=1.e4*err_rate,marker='.', mfc='black', mec='black', ms=8, linestyle='-',color='blue',label='Dilday rate')
axbb.errorbar(zzb,1.e4*rate_two,yerr=1.e4*err_rateb,marker='.', mfc='black', mec='black', ms=8, linestyle='-',color='red',label='Perrett rate')

axbb.set_yscale('log')
axbb.set_xlabel('z')
axbb.set_ylabel('10$^{-4}$ SNe Ia yr$^{-1}$ Mpc$^{-3}$ h$^{3}$$_{70}$')
axbb.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
axbb.set_ylim(top=1.1)
plt.legend(loc='best')

print 'Tot nsn',np.sum(nsn),np.sum(nsnb)
tot_sn=np.sum(nsn)
err_tot_sn=np.power(np.sum(np.power(err_nsn,2.)),0.5)
tot_snb=np.sum(nsnb)
err_tot_snb=np.power(np.sum(np.power(err_nsnb,2.)),0.5)

figc, axc = plt.subplots(ncols=1, nrows=1, figsize=(10,9))
ll='N$_{SN Ia}$ = '+str(int(tot_sn))+'$\pm$'+str(int(err_tot_sn))
figc.suptitle('Survey area: '+str(survey_area)+' deg$^{2}$ - Duration: '+str(duration)+' year')
axc.errorbar(zz,nsn,yerr=err_nsn,marker='.', mfc='black', mec='black', ms=8, linestyle='-',color='blue',label='Dilday rate - '+ll)
ll='N$_{SN Ia}$ = '+str(int(tot_snb))+'$\pm$'+str(int(err_tot_snb))
axc.errorbar(zzb,nsnb,yerr=err_nsnb,marker='.', mfc='black', mec='black', ms=8, linestyle='-',color='red',label='Perrett rate - '+ll)
#axc.set_yscale('log')
axc.set_xlabel('z')
axc.set_ylabel('N(SNe Ia) / z bin')
#axc.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
axc.set_xlim([0.,1.])
axc.legend(loc='best',prop={'size':12})

figd, axd = plt.subplots(ncols=1, nrows=1, figsize=(10,9))
figd.suptitle('Survey area: '+str(survey_area)+' deg$^{2}$ - Duration: '+str(duration)+' year')
ll='N$_{SN Ia}$ = '+str(int(tot_sn))+'$\pm$'+str(int(err_tot_sn))
axd.errorbar(zz,np.cumsum(nsn),yerr=err_nsn,marker='.', mfc='black', mec='black', ms=8, linestyle='-',color='blue',label='Dilday rate - '+ll)
ll='N$_{SN Ia}$ = '+str(int(tot_snb))+'$\pm$'+str(int(err_tot_snb))
axd.errorbar(zzb,np.cumsum(nsnb),yerr=err_nsnb,marker='.', mfc='black', mec='black', ms=8, linestyle='-',color='red',label='Perrett rate - '+ll)
#axd.set_yscale('log')
axd.set_xlabel('z')
axd.set_ylabel('N(SNe Ia) < z')
axd.set_xlim([zmin,zmax])
axd.legend(loc='best',prop={'size':12})


plt.show()
