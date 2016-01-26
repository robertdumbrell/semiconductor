#!/usr/local/bin/python
# UTF-8

import numpy as np
import scipy.constants as const


def Couderc(vals, temp, **kwargs):
    '''
    Parameterisation of fit based on change in band gap with temperature
    of ml (temp in dependent) and mt (band gap dependent). This provides
    how mdc and mdv change with temp, and have been fitted with a 3rd order
    polynomial. 
    '''

    dc_value = vals['ac'] * temp**3 + vals['bc'] * \
        temp**2 + vals['cc'] * temp + vals['dc']
    dv_value = vals['av'] * temp**3 + vals['bv'] * \
        temp**2 + vals['cv'] * temp + vals['dv']

    # The constant is
    coef = (2. * const.pi / const.h / const.h * const.k * const.m_e)
    # adjust Plancs constant to be in cm
    coef /= 10000
    coef = 2 * np.power(coef, 3. / 2)

    Nc = coef * dc_value * np.power(temp, (3. / 2.))
    Nv = coef * dv_value * np.power(temp, (3. / 2.))
    return Nc, Nv


def Green_full(vals, temp, Egratio):

    coef = (2. * const.pi / const.h / const.h * const.k * const.m_e * 300)
    # adjust Plancs constant to be in cm
    coef /= 10000
    coef = 2 * np.power(coef, 3. / 2)

    # print Egratio

    mdcdme = 6.**(2. / 3.) * (
        (vals['mt']*Egratio)**2. * (vals['ml'])) ** (1. / 3.)

    mdvdme = (vals['a'] +
              vals['b'] * temp +
              vals['c'] * temp**2. +
              vals['d'] * temp**3. +
              vals['e'] * temp**4.
              ) / (
        1. +
        vals['f'] * temp +
        vals['g'] * temp**2. +
        vals['h'] * temp**3. +
        vals['i'] * temp**4.)

    # print vals

    mdvdme = mdvdme**(2. / 3.)

    Nc = coef * mdcdme**(3. / 2) * (temp / 300.) ** (3. / 2.)
    Nv = coef * mdvdme**(3. / 2) * (temp / 300.) ** (3. / 2.)

    return Nc, Nv


def Green_param(vals, temp, **kwargs):

    Nc = 2.86e19 * np.power(temp / 300., (1.58))
    Nv = 3.10e19 * np.power(temp / 300., (1.85))

    return Nc, Nv

# print const.physical_constants['Planck constant']
# Green(1, 300)
