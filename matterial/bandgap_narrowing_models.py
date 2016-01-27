#!/usr/local/bin/python
# UTF-8

import numpy as np
import scipy.constants as Const


def apparent_BGN(vals, doping, **kargs):
    '''
    It returns the 'apparent BGN'. This estimates the real band gap
    narrowing, but uses the accurate boltzman statastics. N is the
    net dopant concentration.
    This model incorporates degeneracy as well as band gap narrowing,
    which is why it determines an 'apparent BGN' abd not an actual BGN.
    The apparent occurs as Boltzmann statistics are used rather than
    Fermi-Dirac satastic.
    '''

    if type(doping).__module__ != np.__name__:
        doping = np.array(doping)

    BGN = np.zeros(np.array(doping).shape)

    index = doping > vals['n_onset']
    BGN[index] = (
        vals['de_slope'] * np.log(doping[index] / vals['n_onset']))

    return BGN


def not_implimented(vals, doping, **kargs):
    '''
    model not implimented, returning 0 values
    '''
    return np.zeros(doping.shape)


def BGN(vals, doping, **kargs):
    '''
    Returns the BGN when applied for carriers with fermi distribution.
    This estimates the real band gap narrowing,
    but uses boltzman stats
    where N is the net dopant concentration.
    '''
    # BGN = np.zeros(doping.shape)

    bgn = vals['de_slope']\
        * np.power(np.log(doping / vals['n_onset']), vals['b'])\
        + vals['de_offset']

    # ensures no negitive values
    if bgn.size > 1:
        bgn[bgn < 0] = 0

    return bgn


def Schenk(vals, Nd, Na, ne, nh, temp, **args):
    '''
    Based on the two principles:
    1. The rigid quasi-particle shifts of the conduction and
    valence band edges, which depend from the free-carrier concentrations
    2. the ionic quasi-particle shifts of the conduction and
    valence band edges resulting from ionised dopants concentrations
    '''

    # makes the values unitless
    ne *= vals['aex']**3.
    nh *= vals['aex']**3.
    Na *= vals['aex']**3.
    Nd *= vals['aex']**3.

    n_sum = ne + nh
    n_ionic = Na + Nd

    n_p = vals['alphae'] * ne + vals['alphah'] * nh

    # cacualtes curly T, to give linear dependence with temp
    vals['t'] = Const.k * temp / vals['ryex'] / Const.e

    delta_Ec = ridged_shift(vals, n_sum, n_p, ne, 'e')\
        + ionic_shift(vals, n_sum, n_p, n_ionic, 'e')
    delta_EV = ridged_shift(vals, n_sum, n_p, nh, 'h')\
        + ionic_shift(vals, n_sum, n_p, n_ionic, 'h')

    return delta_Ec + delta_EV


def ridged_shift(vals, n_sum, n_p, num_carrier, carrier):
    '''
    The rigid quasi-particle shift for a band
    the subscript a represents a carrier value (electron or hole)
    '''

    delta = -(
        (4. * Const.pi)**3. * n_sum**2. *
        (
            (48. * num_carrier / Const.pi / vals['g' + carrier])**(1. / 3.)
            + vals['c' + carrier] * np.log(1. + vals['d' + carrier] * n_p**vals['p' + carrier]))
        + (8. * Const.pi * vals['alpha' + carrier] / vals['g' + carrier])
        * num_carrier * vals['t']**2.
        + np.sqrt(8. * Const.pi * n_sum) * vals['t']**(5. / 2.)
    ) / (
        (4. * Const.pi)**3. * n_sum**2. + vals['t']**3. + vals['b' + carrier] *
        np.sqrt(n_sum) * vals['t']**2. + 40. * n_sum**1.5 * vals['t'])
    return -vals['ryex'] * delta


def ionic_shift(vals, n_sum, n_p, n_ionic, carrier):
    '''
    The ionic quasi-particle shift for a band
    '''

    U = n_sum**2. / vals['t']**3.

    delta = -n_ionic * (1. + U) / (
        np.sqrt(vals['t'] * n_sum / 2. / Const.pi) *
        (
            1. + vals['h' + carrier] * np.log(1. + np.sqrt(n_sum) / vals['t'])
        )
        + vals['j' + carrier] * U * n_p**.75 *
        (1. + vals['k' + carrier] * n_p**vals['q' + carrier])
    )

    return -vals['ryex'] * delta
