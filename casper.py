import numpy as np 
from scipy import integrate 


def Critical_density(scale_factor,omega_matt,hubble_const):
    #'hubble constant is little h

    G=4.301*10**(-9) #big g in Mpc M_sun (km s^-1)^2
    
    rho_crit=3/(8*np.pi*G)*(100)**2*omega_matt
    
    return(rho_crit)

def smooth_kspace(k,R,mu,beta):
    y=1/(1+(10**mu*k*R/2.50)**(beta*3.12))
    return(y)

def FT_spherical_top_hat(k,R):
    y=3*(np.sin(k*R)-k*R*np.cos(k*R))/(k*R)**3
    return(y)

def density_rms(R,pk, window_function,*filter_args):

    rms=np.empty(len(R))
    for i in range(len(R)):
        rms[i]=integrate.simps(pk[:,0]**2*pk[:,1]*window_function(pk[:,0],R[i],*filter_args)**2/(2*np.pi**2) , pk[:,0])

    return(rms)


def peak_height(pk,M,h,omega_m,window_function,*filter_args):

    rho_crit=Critical_density(1.0,omega_m,h)
    rho_mean=rho_crit*omega_m
    R=(M/(4/3*np.pi*rho_mean))**(1/3)

    rms=density_rms(R,pk,window_function,*filter_args)

    peak_height=1.68/rms**0.5

    return(peak_height)

def casper(M,pk,omega_m,h,return_peak_height=False):

    #best fit window function for concentration
    beta_c=2; mu_c=-0.55

    #best fit window function for the shape parameter
    beta_alpha=2; mu_alpha=0.0

    #calcate peak heights for the two window functions 
    nu_c=peak_height(pk,M,h,omega_m,smooth_kspace,mu_c,beta_c)
    nu_alpha=peak_height(pk,M,h,omega_m,smooth_kspace,mu_alpha,beta_alpha)

    c=4.39*nu_c**(-0.96)
    alpha=0.0019*nu_alpha**4+0.167

    if return_peak_height==False:
        return(c,alpha)

    elif return_peak_height==True:
        return(c,alpha,nu_c,nu_alpha)

