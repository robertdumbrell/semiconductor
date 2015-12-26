
import numpy as np
import matplotlib.pylab as plt
import sys
import os
import ConfigParser


def add_mobilities(self, mobility_list):
    imobility = 0

    for i in mobility_list:
        imobility += 1. / i

    return 1. / imobility



def CaugheyThomas( vals,  impurity, min_car_den, **kwargs):
    '''
    emperical form for one temperature taken from:
    D. M. Caughey and R. E. Thomas, Proc. U.E.E., pp. 2192,
    2193 (Dec. 1977).

    inputs:
        impurty: the number of impurities (cm^-3)
        min_carr_den: the number of minoirty carrier densities (cm^-3)
        maj_car:  the majority carrier type
        temp: temperature (K)
    output:
        mobility (cm^2 V^-1 s^-1)

    '''

    mu = vals['mu_min'] + (vals['mu_max'] - vals['mu_min']
                           ) / (1. + (impurity / vals['nr'])**vals['alpha'])
    return mu 

def dorkel(vals, impurity, min_car_den, maj_car_den, temp, carrier):
    '''
    not consistent with PVlihthouse at high injection

    inputs:
        impurty: the number of impurities (cm^-3)
        min_carr_den: the number of minoirty carrier densities (cm^-3)
        maj_car_den: the number of majority carrier densities (cm^-3)
        temp: temperature (K)
    output:
         electron mobility (cm^2 V^-1 s^-1)
         hole mobility (cm^2 V^-1 s^-1)
    '''

    # this relatves the carrier to the extension in the variable name
    if carrier == 'electron':
        carrier = 'e'
    elif carrier == 'hole':
        carrier = 'h'

    # determine hole dependent carrier partial mobilities
    mu_L = lattice_mobility(vals, temp, carrier)
    mu_i = impurity_mobility(vals, impurity, temp, carrier)

    # determine both carrier scattering mobilities
    mu_css = carrier_scattering_mobility(
        vals, min_car_den, maj_car_den, temp)

    # determine sudo function
    X = np.sqrt(6. * mu_L * (mu_i + mu_css) / (mu_i * mu_css))

    # combine partial moblities into total
    mu = mu_L * (1.025 / (1. + (X / 1.68)**(1.43)) - 0.025)

    return mu

def lattice_mobility(vals, temp, carrier):
    ''' 
    due to scattering of acoustic phonons
    '''
    mu_L = vals['mul0' + carrier] * \
        (temp / vals['temp0'])**(-vals['alpha' + carrier])
    return mu_L

def impurity_mobility( vals, impurity, temp, carrier):
    '''
    interactions between the carriers and the ionized impurities.
    This partial mobility increases as the temperature
    increases or the doping concentration
    decreases. The relationship which we use in the calculatipn
    of the pr component is that of Brooks and
    Herring
    '''
    A = np.log(1. + vals['b' + carrier] * temp**2 / impurity)
    B = (vals['b' + carrier] * temp ** 2.) / \
        (impurity + vals['b' + carrier] * temp**2)
    mu_i = vals['a' + carrier] * temp**(3. / 2) / impurity / (A - B)
    return mu_i

def impurity_neutral():
    pass

def carrier_scattering_mobility( vals, min_car_den, maj_car_den, temp):
    '''
    The coefficient in B is 2e17 (equation 3) and not 2e7 (Equation 7) as presented in the paper

    '''

    A = np.log(1. + 8.28e8 * temp**2 / (min_car_den * maj_car_den)**(1. / 3))
    B = 2e17 * temp**(3. / 2) / np.sqrt(min_car_den * maj_car_den)
    mu_css = B / A
    return mu_css


