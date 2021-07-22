import numpy as np 
import matplotlib.pyplot as plt 
from casper import casper
from casper import Einasto_density
from casper import Einasto_mass
from casper import Critical_density
#define cosmologival parameters
omega_m=0.2793
h=0.7
rho_crit=Critical_density()

#load the power spectra
pk_0=np.loadtxt('./Power_spectra/WMAP9_z_0.dat')
pk_1=np.loadtxt('./Power_spectra/WMAP9_z_1.dat')
pk_2=np.loadtxt('./Power_spectra/WMAP9_z_2.dat')

#define mass range (in M_200c [h^-1 M_sun]) to calculate the density paramters
M=np.logspace(12,14,100)


#use the fucntion casper to predict the concentration and shape parameter for the various redshifts
c0,alpha0=casper(M,pk_0,omega_m,return_peak_height=False)
c1,alpha1=casper(M,pk_1,omega_m,return_peak_height=False)
c2,alpha2=casper(M,pk_2,omega_m,return_peak_height=False)


#plot the c and alpha mass relations
plt.figure()
plt.plot(M,c0,label='z=0')
plt.plot(M,c1,label='z=1')
plt.plot(M,c2,label='z=2')

plt.xscale('log')
plt.yscale('log')
plt.xlabel('$M_{200c} [h^{-1} M_{\\odot}]$')
plt.ylabel('$c$')
plt.legend(frameon=False)


plt.figure()
plt.plot(M,alpha0,label='z=0')
plt.plot(M,alpha1,label='z=1')
plt.plot(M,alpha2,label='z=2')

plt.xscale('log')

plt.xlabel('$M_{200c} [h^{-1} M_{\\odot}]$')
plt.ylabel('$\\alpha$')
plt.legend(frameon=False)


#plot the resulting density and circular velocity profiles and  for a mass of 10^12 and 10^15 at z=0
M200=np.array([10**12,10**15])
R200=(M200/(200*rho_crit*4/3*np.pi))**(1/3)
c,alpha=casper(M200,pk_0,omega_m,return_peak_height=False)

r=np.logspace(-1.5,0,100) #radii to plot the density profiles

plt.figure()
ax1=plt.subplot()
plt.figure()
ax2=plt.subplot()

for i in range(len(M200)):
    rho_2=M200[i]/Einasto_mass(1,1,c[i],alpha[i],R200[i]) #normalise the denisty profiles
    
    rho=Einasto_density(r,rho_2,c[i],alpha0[i])
    v_circ=np.sqrt(Einasto_mass(r,rho_2,c[i],alpha0[i],R200[i])/r)

    ax1.plot(r,rho,label='$M_{200c}=%.0e$, $z=0$'%M200[i])
    ax2.plot(r,v_circ/np.max(v_circ),label='$M_{200c}=%.0e$, $z=0$'%M200[i])

ax1.set_xlabel('$r/R_{200c}$')
ax1.set_ylabel('$\\rho$ $[h^2 M_{\\odot} Mpc^{-3}]$')
ax1.legend(frameon=False)
ax1.set_xscale('log')
ax1.set_yscale('log')

ax2.set_xlabel('$r/R_{200c}$')
ax2.set_ylabel('$v_{circ}/v_{circ,max}$')
ax2.legend(frameon=False)


plt.show()
