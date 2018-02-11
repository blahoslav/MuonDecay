# -*- coding: utf-8 -*-
"""
Created on Tue Dec 13 09:30:30 2016

@author: br616
"""
import scipy as sp
import pylab as pl
import scipy.optimize as spo
# Constants given
r_s = 0.14
M_z = 91.0
j_a = 0.81
Gamma_z = -2.5
j_s = -0.033
r_a = 0.0027

# Task 2 - distribution of particles - with constant s
# Choose value of s
def calculate_muons(s):
    sigma_s, sigma_o = sigma_equations(s)
    # Integrate over variable angle (limits (-0.95, 0.95)) with constant s
    number_of_muons = []
    uperr_limit = -0.90
    lower_limit = -0.95
    while (uperr_limit<0.95): 
        number_of_muons_error = sp.integrate.quad(number_of_muons_gradient_cos, lower_limit, uperr_limit,args=(sigma_s, sigma_o,1000000.0))
        number_of_muons.append( number_of_muons_error[0])
        uperr_limit += 0.05
        lower_limit += 0.05
        print lower_limit
    # Print number of muons per 24h at set s and different angle
    print number_of_muons
    # Plot number of muons per 24h at set s and different angle - dodelat spravne uhly -  cosθ bins
    pl.title('n/theta')
    angles = sp.arange(-0.95, 0.90, 0.05)
    pl.subplot(4,2,5)
    pl.plot(angles, number_of_muons, 'r.')
    #poiss = sp.random.poisson(lam=sp.mean(number_of_muons), size=(3., 4.))
    pl.subplot(4,2,6)
    poiss = []
    for y in range(0, len(angles)):
        poiss.append(sp.random.poisson(number_of_muons[y]))
    pl.plot(angles, poiss, 'b+')
    
    #Error of muons at every bin
    muons_error = sp.sqrt(poiss)
    print len(muons_error)
    
    # Fitting, reference manual
    def function_for_fit(x, a, ratio):
        return a*(1+x**2)+a*(ratio*x)
    
    pl.subplot(4,2,7)
    # Random initial guess of constants
    initial_guess=[0.5, 1]
    po,po_cov = spo.curve_fit(function_for_fit, angles, poiss, initial_guess, muons_error)
    pl.errorbar(angles,poiss,yerr=muons_error, fmt='ro',label="Measured data")
    xx=sp.arange(-0.95,0.95,0.05)
    pl.plot(xx,function_for_fit(xx,po[0], po[1]),'b-',label='Fit results')
    pl.xlabel("x")
    pl.ylabel("y")
    pl.grid(1)
    pl.show()
    # Print ratio determined by fitting, for single s
    return po[1] 
    
def sigma_equations(s_square):
    # Sigma_s calculation
    a = (s_square**2)*r_s+(((s_square**2)-(M_z)**2)*j_s)
    b = ((s_square**2)-(M_z**2))**2+((M_z**2)*(Gamma_z**2))
    sigma_s = ((4.0/3.0)*sp.pi)*((1/(s_square**2))+(a/b))
    

    
    # Sigma_a calculation    
    sigma_o = sp.pi*((((s_square**2)*r_a)+((s_square**2)-M_z**2)*j_a)/(((s_square**2)-(M_z**2))**2+(M_z**2)*(Gamma_z**2)))
    return(sigma_s, sigma_o) #return 2 separate rrays as tuple

def number_of_muons_gradient_cos(fita, sigma_s, sigma_o, k):
    # Function to integrate- return number of muons as a function of s
     return(k*(sigma_s*(1+fita*fita)+(sigma_o*fita)))