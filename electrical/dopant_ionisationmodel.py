#!/usr/local/bin/python
# UTF-8

import numpy as np
import scipy.constants as C
import matplotlib.pylab as plt
# from semiconductor.helper.helper import HelperFunctions
# import dopant_ionisation_models
# import semiconductor.matterial.bandgap_narrowing_models as Bgn
# import semiconductor.general_functions.carrierfunctions as GF
from glob import glob

SiP_table3 = {
    'dopant': 'p',
    'E_dop0': 45.5e-3,
    'N_ref': 2.2e18,
    'c': 2,
    'N_b': 6e18,
    'd': 2.3,
    'g': 0.5,
    'tpe': 'donor'}

SiP = {
    'dopant': 'p',
    'E_dop0': 45.5e-3,
    'N_ref': 3e18,
    'c': 2,
    'N_b': 6e18,
    'd': 2.3,
    'g': 0.5,
    'tpe': 'donor'}

SiB = {
    'dopant': 'b',
    'E_dop0': 44.39e-3,
    'N_ref': 1.7e18,
    'c': 1.4,
    'N_b': 6e18,
    'd': 2.4,
    'g': 0.25,
    'tpe': 'acceptor'}


SiB_table3 = {
    'dopant': 'b',
    'E_dop0': 44.39e-3,
    'N_ref': 1.3e18,
    'c': 1.4,
    'N_b': 4.5e18,
    'd': 2.4,
    'g': 1./4,
    'tpe': 'acceptor'}

SiAs_table3 = {
    'dopant': 'as',
    'E_dop0': 53.7e-3,
    'N_ref': 3e18,
    'c': 1.5,
    'N_b': 9e18,
    'd': 1.8,
    'g': 0.5,
    'tpe': 'donor'}

SiAs = {
    'dopant': 'as',
    'E_dop0': 53.7e-3,
    'N_ref': 4e18,
    'c': 1.5,
    'N_b': 1.4e19,
    'd': 3,
    'g': 0.5,
    'tpe': 'donor'}
# r = 4.2e-12
# s = 1e19
# T = np.linspace(30,300)
# E = C.k * T / C.e


def altermatt(dopant):
    pass


# def DOS_gaussian(E):
#     ''' density of states assuming a gaussian function'''
#     print b(), E_dop(), delta()
#     return (N_dop * b()) / np.sqrt(s * C.pi) / delta() *\
#         np.exp(-(E - E_dop())**2 / (2. * delta()**2))


def E_dop(values, Ni):
    '''retuns the Dopant energy level in eV'''
    return values['E_dop0'] / (1. + (Ni / values['N_ref'])**values['c'])


# def delta():
#     '''half width of dopants'''
#     return r * np.sqrt(Ni) * (1. - np.exp(-s / Ni))


def b(values, Ni):
    '''fration of carriers in localised states'''
    return 1. / (1. + (Ni / values['N_b'])**values['d'])


def altermatt2006(values, N_impurity, ne, nh, T):
    '''
    Dopant ionisation of single doped material
    One exists for co doped material, its just not
    put together
    '''

    # something still not correct

    # TO DO: get real values of Nc and Nv,
    # these are from Hu_ch01v4.fm Page 21 Thursday, February 12, 2009 10:14 AM
    Nc = 2.8e19
    Nv = 1.04e19

    Av = 2.525e-9
    Bv = -4.689e-6
    Cv = 3.376e-3
    Dv = 3.426e-1

    Ac = -4.609e-10
    Bc = 6.7523e-7
    Cc = -1.312e-5
    Dc = 1.094 

    dc_value = Ac*T**3 + Bc*T**2 + Cc*T + Dc
    dv_value = Av*T**3 + Bv*T**2 + Cv*T + Dv

    #couders values
    Nc = 4.83e15 * dc_value * T**(3./2.)
    Nv = 4.83e15 * dv_value * T**(3./2.)

    # sentarous
    Nc = 2.89e19
    Nv = 3.14e19
    print Nc, Nv

    # green 1990 values
    # Nc = 2.86e19
    # Nv = 3.1e19
    vt = C.k * T / C.e
    if values['tpe'] == 'donor':
        ne1 = Nc * np.exp(-E_dop(values, N_impurity) / vt)
        ratio = 1. - b(values, N_impurity) * ne / (ne + values['g'] * ne1)

    elif values['tpe'] == 'acceptor':

        nh1 = Nv * np.exp(-E_dop(values, N_impurity) / vt)
        ratio = 1. - b(values, N_impurity) * nh / (nh + values['g'] * nh1)

    return ratio


# plot dos
T = 300.
ne = nh = 1e10
N = np.logspace(15, 20)

for i in [SiP, SiB, SiAs]:
    nh = ne = N
    for j in range(5):
        ne = altermatt2006(i, N, ne, ne, T)
        ne *= N
    ne /= N
    plt.plot(N, ne, '.', label=i['dopant'])

plt.semilogx()
plt.legend(loc=0)
plt.ylim(.65, 1.05)



for fname in glob('.\Si\check data\*.csv'):
    data = np.genfromtxt(fname, delimiter=',', names = ['dopants', 'ii'])
    plt.plot(data['dopants'], data['ii'], label=fname)




plt.show()
