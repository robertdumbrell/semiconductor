
import numpy as np


def Passler(vals, temp):
    """
    taken from the Couderc2014 paper
    depent on temperature

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
