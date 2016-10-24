#Module containing frequently used physics constants and functions

import math
import cmath
import numpy as np

#%%
#==============================================================================
# Fundamental constants
#==============================================================================
me = 9.10938*10**-31 #Electron mass in kg
mp = 1.67262158*10**-27 #Proton mass in kg
q = 1.602177*10**-19 #Elementary charge in C
G = 6.674*10**-11 #Gravitational constant in N*m^2/kg^2
c = 2.998*10**8 #Speed of light in m/s
epsilon0 = 8.854*10**-12 #Electric constant in F/m
mu0 = 4*math.pi*10**-7 #Magnetic constant in T*m/A
hbar = 1.054572*10**-34 #Reduced Planck constant in J*s
kb = 1.38065*10**-23 #Boltzmann constant in J/K
NA = 6.022*10**23 #Avagadro's number


#%%
#==============================================================================
# Functions
#==============================================================================

def savitzky_golay(y, window_size, order, deriv=0, rate=1):
    """Smooth (and optionally differentiate) data with a Savitzky-Golay filter.
    
    y : array_like, shape (N,)
        the values of the time history of the signal.
    window_size : int
        the length of the window. Must be an odd integer number.
    order : int
        the order of the polynomial used in the filtering.
        Must be less then `window_size` - 1.
    deriv: int
        the order of the derivative to compute (default = 0 means only smoothing)
    Returns
    -------
    ys : ndarray, shape (N)
        the smoothed signal (or it's n-th derivative)."""
    
    try:
        window_size = np.abs(np.int(window_size))
        order = np.abs(np.int(order))
    except ValueError:#, msg:
        raise ValueError("window_size and order have to be of type int")
    if window_size % 2 != 1 or window_size < 1:
        raise TypeError("window_size size must be a positive odd number")
    if window_size < order + 2:
        raise TypeError("window_size is too small for the polynomials order")
    order_range = range(order+1)
    half_window = (window_size -1) // 2
    # precompute coefficients
    b = np.mat([[k**i for i in order_range] for k in range(-half_window, half_window+1)])
    m = np.linalg.pinv(b).A[deriv] * rate**deriv * math.factorial(deriv)
    # pad the signal at the extremes with
    # values taken from the signal itself
    firstvals = y[0] - np.abs( y[1:half_window+1][::-1] - y[0] )
    lastvals = y[-1] + np.abs(y[-half_window-1:-1][::-1] - y[-1])
    y = np.concatenate((firstvals, y, lastvals))
    return np.convolve( m[::-1], y, mode='valid')


def Sellmeier(wavelength,B1,C1,B2,C2,B3,C3):
    """Calculates the refractive index of a material given the wavelength in
    microns and the Sellmeier coefficients."""
    nsquared = 1 + B1*wavelength**2/(wavelength**2-C1) + \
    B2*wavelength**2/(wavelength**2-C2) + \
    B3*wavelength**2/(wavelength**2-C3)
    
    return math.sqrt(nsquared)
    
class Drude:
    
    def omegap(ne,meff):
        """Returns the plasma frequency in Hz.
        ne: electron density (m^-3)
        meff: effective mass (kg)"""

        return np.sqrt(ne*q**2/(meff*epsilon0))
        
    def epsilon(epsilonc,ne,meff,tau,wavelength):
        """Returns the laser-excited permittivity.
        epsilonc: non-excited permittivity
        ne: electron density (m^-3)
        meff: effective mass (kg)
        tau: average electron-phonon collision time (s)
        wavelength: laser wavelength (m)."""
        
        omega = 2*math.pi*c/wavelength
        return epsilonc - Drude.omegap(ne,meff)**2/(omega*(omega + 1/tau*1.0j))
        
    def RefractiveIndex(epsilonc,ne,meff,tau,wavelength):
        """Returns the laser-excited refractive index.
        epsilonc: non-excited permittivity
        ne: electron density (m^-3)
        meff: effective mass (kg)
        tau: average electron-phonon collision time (s)
        wavelength: laser wavelength (m)."""
        
        return np.sqrt(Drude.epsilon(epsilonc,ne,meff,tau,wavelength))
        

class Optics:
    
    def RefractiveIndex(epsilon):
        """Returns the refractive index given the complex permittivity"""
        n = cmath.sqrt(epsilon)
        return n
        
    def Reflectivity(n,theta,polarization):
        """Returns the reflectivity given the complex refractive index n, the 
        angle of incidence theta, and the polarization"""
        
        if polarization=='p':
            return abs((cmath.sqrt(1-(1/n*math.sin(theta))**2)-n*math.cos(theta))/\
            (cmath.sqrt(1-(1/n*math.sin(theta))**2)+n*math.cos(theta)))**2
            
        elif polarization=='s':
            return abs((math.cos(theta)-n*cmath.sqrt(1-(1/n*math.sin(theta))**2))/\
            (math.cos(theta)+n*cmath.sqrt(1-(1/n*math.sin(theta))**2)))**2
            
        else: #unpolarized
            Rs = abs((math.cos(theta)-n*cmath.sqrt(1-(1/n*math.sin(theta))**2))/\
            (math.cos(theta)+n*cmath.sqrt(1-(1/n*math.sin(theta))**2)))**2
        
            Rp = abs((cmath.sqrt(1-(1/n*math.sin(theta))**2)-n*math.cos(theta))/\
            (cmath.sqrt(1-(1/n*math.sin(theta))**2)+n*math.cos(theta)))**2
            
            return (Rs+Rp)/2
            
    def BrewsterAngle(n):
        """Returns the Brewster angle given the refractive index n"""
        return math.tan(n.real)*360/(2*math.pi)
        

class Greek:
    
    greek_letters=[chr(i) for i in range(945,970)]
    
    textalpha = greek_letters[0]
    textbeta = greek_letters[1]
    textgamma = greek_letters[2]
    textdelta = greek_letters[3]
    textepsilon = greek_letters[4]
    textzeta = greek_letters[5]
    texteta = greek_letters[6]
    texttheta = greek_letters[7]
    textiota = greek_letters[8]
    textkappa = greek_letters[9]
    textlambda = greek_letters[10]
    textmu = greek_letters[11]
    textnu = greek_letters[12]
    textxi = greek_letters[13]
    textomicron = greek_letters[14]
    textpi = greek_letters[15]
    textrho = greek_letters[16]
    textvarsigma = greek_letters[17]
    textsigma = greek_letters[18]
    texttau = greek_letters[19]
    textupsilon = greek_letters[20]
    textphi = greek_letters[21]
    textchi = greek_letters[22]
    textpsi = greek_letters[23]
    textomega = greek_letters[24]

#%%
#==============================================================================
# Colormaps
#==============================================================================
