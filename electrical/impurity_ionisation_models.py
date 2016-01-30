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

