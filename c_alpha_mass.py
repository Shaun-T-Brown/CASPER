from casper import casper
import numpy as np 
import matplotlib.pyplot as plt 

#define cosmologival parameters
omega_m=0.2793
h=0.7

#load the power spectra
pk_0=np.loadtxt('./Pk_z_0.dat')
pk_1=np.loadtxt('./Pk_z_1.dat')
pk_2=np.loadtxt('./Pk_z_2.dat')

#define mass range (in M_200c) to calucalte the density paramters
M=np.logspace(12,14,100)


#use the fucntion casper to predict the concentration and shape parameter for the various redshifts
c0,alpha0,nu_c0,nu_alpha0=casper(M,pk_0,omega_m,h,return_peak_height=True)
c1,alpha1,nu_c1,nu_alpha1=casper(M,pk_1,omega_m,h,return_peak_height=True)
c2,alpha2,nu_c2,nu_alpha2=casper(M,pk_2,omega_m,h,return_peak_height=True)

#plot the data

plt.figure()
plt.plot(M,c0,label='z=0')
plt.plot(M,c1,label='z=1')
plt.plot(M,c2,label='z=2')

plt.xscale('log')
plt.yscale('log')
plt.xlabel('$M_{200c} [h^{-1} M_{\\odot}]$')
plt.ylabel('$c$')
plt.legend()


plt.figure()
plt.plot(M,alpha0,label='z=0')
plt.plot(M,alpha1,label='z=1')
plt.plot(M,alpha2,label='z=2')

plt.xscale('log')
plt.yscale('log')
plt.xlabel('$M_{200c} [h^{-1} M_{\\odot}]$')
plt.ylabel('$\\alpha$')
plt.legend()




plt.show()
