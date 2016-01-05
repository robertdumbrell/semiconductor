
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