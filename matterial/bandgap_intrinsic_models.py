
import numpy as np


def Passler(vals, temp):
    """
    taken from the Couderc2014 paper
    depent on temperature

    This model different from proceeding ones by( taken
    from 10.1063/1.4867776):

    1. The model applies to all physically relevant
    regimes of dispersion comprising the whole interval
    from the Bose-Einstein limit Delta -> 0 up to sufficiently
    large magnitudes of the dispersion coefficient i.e., the
    experimentally relevant range of 0 < Delta < 3/4, at least!.

    2. The theoretical E(T) function should tend in the cryogenic
    region to a quadratic asymptote, i.e., E(0)
    [E(0) - E(T)] prop T^2 for T < Theta , in accordance with basic theoretical
    expectations as well as experimental observations.

    3. The new model should be suited for the construction of a
    practicable analytical expression that is capable of providing
    accurate self-consistent values for the
    dispersion-related parameters Theta and Delta.

    4. An important additional aim of the model is to be able to
    represent the final E(T) expression explicitly in terms of
    the model-independent parameters Theta and Delta, instead of
    requiring preliminary determinations of various model specific
    auxiliary quantities.

    returns Eg in eV
    """

    gamma = (1. - 3. * vals['delta']**2) / \
        (np.exp(vals['theta'] / temp) - 1)
    xi = 2. * temp / vals['theta']

    # Values for each sum component
    No2 = np.pi**2. * xi**2. / (3. * (1 + vals['delta']**2))
    No3 = (3. * vals['delta']**2 - 1) / 4. * xi**3
    No4 = 8. / 3. * xi**4.
    No5 = xi**6.

    E = vals['e0'] - vals['alpha'] * vals['theta'] * \
        (gamma + 3. * vals['delta']**2 / 2 *
         ((1. + No2 + No3 + No4 + No5)**(1. / 6.) - 1))
    return E


def Varshni(vals, temp):
    '''
    Passler's paper suggests that this model is for very
    high dispersion relations  Delta  = 5/4
    '''
    return vals['e0'] - vals['alpha'] * temp**2 / (temp + vals['beta'])


def Cubic_partial(vals, temp):
    '''
    Is a cublic paramterisation for several given temp range, spliced together.
    The first paper where this is seen for silicon is believed to be Bludau1974a.

    inputs:
        vals a dictionary containing the coefs for a fit in the from
            Eg = \sum_{i=0}^3 ai + bi \times temp + ci \times temp^2
        and the temp range for each coeffieinct given by "ti". It is assumed that the ith
        values apply up to this temperature value.

    output:
        returns the band gap in eV
    '''

    # this line is to make sure the number is a numpy 1D array
    # really this should be somewhere else, as a general function and not here
    temp = np.asarray([temp * 1.]).flatten()

    Eg = np.copy(temp)

    for i in [2, 1, 0]:
        index = temp < vals['t' + str(i)]

        Eg[index] = vals['a' + str(i)] + \
            vals['b' + str(i)] * temp[index] + \
            vals['c' + str(i)] * temp[index]**2.
    if np.any(temp > vals['t2']):
        print 'Intrinsic bandgap does not cover this temperature range'
        index = temp > vals['t2']
        Eg[index] = vals['a2'] + \
            vals['b2'] * temp[index] + \
            vals['c2'] * temp[index]**2.
    return Eg
