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


# SiP = {
#     'dopant': 'p',
#     'E_dop0': 45.5e-3,
#     'N_ref': 3e18,
#     'c': 2,
#     'N_b': 6e18,
#     'd': 2.3,
#     'g': 0.5,
#     'tpe': 'donor'}

# SiB = {
#     'dopant': 'b',
#     'E_dop0': 44.39e-3,
#     'N_ref': 1.7e18,
#     'c': 1.4,
#     'N_b': 6e18,
#     'd': 2.4,
#     'g': 0.25,
#     'tpe': 'acceptor'}

# SiAs = {
#     'dopant': 'as',
#     'E_dop0': 53.7e-3,
#     'N_ref': 4e18,
#     'c': 1.5,
#     'N_b': 1.4e19,
#     'd': 3,
#     'g': 0.5,
#     'tpe': 'donor'}

# SiP_table3 = {
#     'dopant': 'p',
#     'E_dop0': 45.5e-3,
#     'N_ref': 2.2e18,
#     'c': 2,
#     'N_b': 6e18,
#     'd': 2.3,
#     'g': 0.5,
#     'tpe': 'donor'}


# SiB_table3 = {
#     'dopant': 'b',
#     'E_dop0': 44.39e-3,
#     'N_ref': 1.3e18,
#     'c': 1.4,
#     'N_b': 4.5e18,
#     'd': 2.4,
#     'g': 1./4,
#     'tpe': 'acceptor'}

# SiAs_table3 = {
#     'dopant': 'as',
#     'E_dop0': 53.7e-3,
#     'N_ref': 3e18,
#     'c': 1.5,
#     'N_b': 9e18,
#     'd': 1.8,
#     'g': 0.5,
#     'tpe': 'donor'}


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


def E_dop(values, Ni, dopant):
    '''retuns the Dopant energy level in eV'''
    return values['e_dop0_' + dopant] / (
        1. + (Ni / values['n_ref_' + dopant])**values['c_' + dopant])


# def delta():
#     '''half width of dopants'''
#     return r * np.sqrt(Ni) * (1. - np.exp(-s / Ni))


def b(values, Ni, dopant):
    '''fration of carriers in localised states'''
    return 1. / (1. + (Ni / values['n_b_' + dopant])**values['d_' + dopant])


def altermatt2006(values, N_impurity, ne, nh, T, Nc, Nv, dopant):
    '''
    This function returns the fraction of ionisated dopants.
    Dopant ionisation of single doped material
    One exists for co doped material, its just not
    put together
    '''

    vt = C.k * T / C.e
    if values['tpe_' + dopant] == 'donor':
        ne1 = Nc * np.exp(-E_dop(values, N_impurity, dopant) / vt)
        ratio = 1. - b(values, N_impurity, dopant) * ne / \
            (ne + values['g_' + dopant] * ne1)

    elif values['tpe_' + dopant] == 'acceptor':

        nh1 = Nv * np.exp(-E_dop(values, N_impurity, dopant) / vt)
        ratio = 1. - b(values, N_impurity, dopant) * nh / \
            (nh + values['g_' + dopant] * nh1)

    return ratio


# plot dos
# T = 300.
# ne = nh = 1e10
# N = np.logspace(15, 20)

# for i in [SiP, SiB, SiAs]:
#     nh = ne = N
#     for j in range(5):
#         ne = altermatt2006(i, N, ne, ne, T)
#         ne *= N
#     ne /= N
#     plt.plot(N, ne, '.', label=i['dopant'])

# plt.semilogx()
# plt.legend(loc=0)
# plt.ylim(.65, 1.05)


# for fname in glob('.\Si\check data\*.csv'):
#     data = np.genfromtxt(fname, delimiter=',', names=['dopants', 'ii'])
#     plt.plot(data['dopants'], data['ii'], label=fname)


# plt.show()
