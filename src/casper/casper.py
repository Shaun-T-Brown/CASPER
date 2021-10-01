import numpy as np 
from scipy import integrate 
from scipy import special
import warnings as wn

def Einasto_density(r,rho_s,c,alpha):
    """

    Function to generate an Einasto density profile.

    Parameters:

    r (array_like): The radial range over which to calculate the density, in units of R200c (or similar).

    rho_s (float): Density normalisation 

    c (float): Concentration

    alpha (float): Shape parameter

    Returns:

    rho (array_like): The density at the specified radii for the given parameters. 

    """

    rho=rho_s*np.exp(-2/alpha*((c*r)**alpha-1))

    return(rho)

def Einasto_mass(r,rho_s,c,alpha,R200):
    """

    Function to generate the integrated mass (i.e. the mass within a given radius) for an Einasto profile.

    Parameters:

    r (array_like): The radial range over which to calculate the density, in units of R200c (or similar).

    rho_s (float): Density normalisation

    c (float): Concentration

    alpha (float): Shape parameter

    Returns:

    rho (array_like): The total mass below the specified radii for the given parameters.


    """

    n=1/alpha
    d_n=2*n

    M=4*np.pi*rho_s*np.exp(2*n)*n/((2*n)**(3*n))*(R200/c)**3*special.gamma(3*n)

    s=d_n**n*r*c
    M_less=M*(1-special.gammaincc(3*n,s**(1/n)))

    return(M_less)

def Critical_density():
    """Function to calculate the critical density, in units of h^2 M_sun Mpc^-3."""
    G=4.301*10**(-9) #big g in Mpc M_sun (km s^-1)^2
    rho_crit=3/(8*np.pi*G)*(100)**2
    
    return(rho_crit)



def smooth_kspace(k,R,mu,beta):
    """Smooth k-space window function."""

    y=1/(1+(mu*k*R/2.50)**(beta*3.12))
    return(y)

def spherical_top_hat(k,R):
    """Spherical top hat window function."""

    y=3*(np.sin(k*R)-k*R*np.cos(k*R))/(k*R)**3
    return(y)

def spherical_top_hat_generalised(k,R,mu_g):
    """Generalised spherical top hat window function."""

    y=3*(np.sin(mu_g*k*R)-mu_g*k*R*np.cos(mu_g*k*R))/(mu_g*k*R)**3
    return(y)

def density_rms(R,pk, window_function,*filter_args):

    """

    Function to calculate the density rms, for a specified window function.

    Parameters:
    
    M (array_like): The radius over which to calculate the density rms, this
    will mot often be the lagrangian radius.

    pk (2d array): The linear power spectrum for the given cosmology. Assumes
    pk[:,0] is k while pk[:,1] is P(k). Units are assumed to be h Mpc^-1 for
    k while h^3Mpc^-3 for P(k)

    omega_m (float): The mass density at z=0.
    
    window_function (function): function for the window function, i.e. W(kR).

    *filter_args: any free parameters associated with the window function
    

    Returns:

    rms (array_like): Peak height values as an array the same length as M.

    """

    rms=np.empty(len(R))
    for i in range(len(R)):
            rms[i]=integrate.simps(pk[:,0]**2*pk[:,1]*window_function(pk[:,0],R[i],*filter_args)**2/(2*np.pi**2) , pk[:,0])

    return(rms)


def peak_height(pk,M,omega_m,window_function,*filter_args):
    """

    Function to calculate the peak height.

    Parameters:

    pk (2d array): The linear power spectrum for the given cosmology. Assumes
    pk[:,0] is k while pk[:,1]  is P(k). Units are assumed to be h Mpc^-1 for
    k while h^3Mpc^-3 for P(k)

    M (array_like): The masses at which to calculate the shape parameters. Any
    mass definition can be used, however our model assumes M200c.

    omega_m (float): The mass density at z=0.
    
    window_function (function): function for the window function, i.e. W(kR)

    *filter_args: any free parameters associated with the window function
    

    Returns:

    peak_height (array_like): Peak height values as an array the same length as M.

    """

    rho_crit=Critical_density()
    rho_mean=rho_crit*omega_m

    R=(M/(4/3*np.pi*rho_mean))**(1/3)

    rms=density_rms(R,pk,window_function,*filter_args)

    peak_height=1.68/rms**0.5

    return(peak_height)

def casper(M,pk,omega_m,return_peak_height=False):
    """

    Function to calculate the concentration and shape parameter for a given
    mass, redshift and cosmology.

    Parameters:


    M (array_like): The masses at which to calculate the shape parameters.
    Assumes the mass is M_200c in  units of h^-1 M_sun.

    pk (2d array): The linear power spectrum for the given cosmology. Assumes
    pk[:,0] is k while pk[:,1]  is P(k). Units are assumed to be h Mpc^-1 for
    k while h^3Mpc^-3 for P(k)

    omega_m (float): The mass density at z=0.
    

    Keyword arguments:

    return_peak_height (boolean): Option to return the peak height values,
    using the two different window function definitions, for the given masses
    and cosmology.
    

    Returns:

    c (array_like): Halo concentration as an array the same length as M.

    alpha (array_like): Halo shape parameter as an array the same length as M.
    

    if return_peak_height==True also returns

    nu_c (array_like): Peak height using window function for concentration as
    an array the same length as M.
    
    nu_alpha (array_like): Peak height using window function for shape parameter as
    an array the same length as M.


    """


    #best fit window function for concentration
    mu_c=10**(-0.67)

    #best fit window function for the shape parameter
    mu_alpha=10**(-0.01)

    nu_c=peak_height(pk, M, omega_m, spherical_top_hat_generalised, mu_c)
    nu_alpha=peak_height(pk, M, omega_m, spherical_top_hat_generalised, mu_alpha)

    #maximum and minumum values of nu_c and nu_alpha that were used for the
    #empirical fits. Outside these ranges is an extroplation and may result in
    #unreliable/inaccurate predictions.
    
    nu_c_min=0.299
    nu_c_max=1.891
    nu_alpha_min=1.029
    nu_alpha_max=3.445

    if np.min(nu_c)<nu_c_min:
        wn.warn("The minimum value of nu_c used for calibration was %.3f. You have requested a minimum value of %.3f, the reliablity/accuracy of the model may be reduced in this regime."%(nu_c_min,np.min(nu_c)), stacklevel = 1)
        
    if np.max(nu_c)>nu_c_max:
        wn.warn("The maximum value of nu_c used for calibration was %.3f. You have requested a maximum value of %.3f, the reliablity/accuracy of the model may be reduced in this regime."%(nu_c_max,np.max(nu_c)), stacklevel = 1)
        
    if np.min(nu_alpha)<nu_alpha_min:
        wn.warn("The minimum value of nu_alpha used for calibration was %.3f. You have requested a minimum value of %.3f, the reliablity/accuracy of the model may be reduced in this regime."%(nu_alpha_min,np.min(nu_alpha)), stacklevel = 1)
        
    if np.max(nu_alpha)>nu_alpha_max:
        wn.warn("The maximum value of nu_alpha used for calibration was %.3f. You have requested a maximum value of %.3f, the reliablity/accuracy of the model may be reduced in this regime."%(nu_alpha_max,np.max(nu_alpha)), stacklevel = 1)
        
    #empirical fits
    c=4.39*nu_c**(-0.87)
    alpha=8.52*10**(-4)*nu_alpha**4+0.166

    if return_peak_height==False:
        return(c,alpha)

    elif return_peak_height==True:
        return(c,alpha,nu_c,nu_alpha)




