from SN_Rate import *
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter

rate=SN_Rate(survey_area=10.)

zz,rate_one,err_rate,nsn,err_nsn=rate(0.1,1.,0.01)

rateb=SN_Rate(rate='Perret',survey_area=10.)
zzb,rate_two,err_rateb,nsnb,err_nsnb=rateb(0.1,1.,0.01,account_for_edges=True)

print zz, nsn,nsnb

figbb, axbb = plt.subplots(ncols=1, nrows=1, figsize=(10,9))
axbb.errorbar(zz,1.e4*rate_one,yerr=1.e4*err_rate,marker='.', mfc='black', mec='black', ms=8, linestyle='-',color='blue',label='Ripoche rate (20% uncert.)')
axbb.errorbar(zzb,1.e4*rate_two,yerr=1.e4*err_rateb,marker='.', mfc='black', mec='black', ms=8, linestyle='-',color='red',label='Perret rate')
axbb.set_yscale('log')
axbb.set_xlabel('Redshift')
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
axc.errorbar(zz,nsn,yerr=err_nsn,marker='.', mfc='black', mec='black', ms=8, linestyle='-',color='blue',label='Ripoche rate - '+ll)
ll='N$_{SN Ia}$ = '+str(int(tot_snb))+'$\pm$'+str(int(err_tot_snb))
axc.errorbar(zzb,nsnb,yerr=err_nsnb,marker='.', mfc='black', mec='black', ms=8, linestyle='-',color='red',label='Perret rate - '+ll)
#axc.set_yscale('log')
axc.set_xlabel('Redshift')
axc.set_ylabel('N(SNe Ia) - one year (survey aera=9.6 deg$^2$, effective search duration=0.5)')
#axc.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
plt.legend(loc='best')



plt.show()
