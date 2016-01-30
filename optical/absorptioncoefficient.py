import numpy as np
import matplotlib.pylab as plt
import scipy.constants as const
import ConfigParser
import os
import sys

import Si.opticalproperties as opticalproperties

#constaints the stuff for tabulated data
OP = opticalproperties.OpticalProperties()

sys.path.append(
    os.path.abspath(os.path.join(os.getcwd(), os.pardir, os.pardir)))
from semiconductor.helper.helper import HelperFunctions


class absorptioncoefficient(HelperFunctions):

    ''' This purpose of this it to provide access if the absorption
        coefficient have a model'''
    beta = 0
    gamma = 0

    def __init__(self, matterial='Si', author=None):
        self.Models = ConfigParser.ConfigParser()
        self.matterial = matterial

        constants_file = os.path.join(
            os.path.dirname(os.path.realpath(__file__)),
            matterial,
            r'absorptioncoefficient.const')

        self.Models.read(constants_file)

        self.vals, self.model = self.change_model(author)


    def update_absorptioncoefficients(self, f=None, Input=None):
        '''
        updates the absorption coefficients
        inputs:
            f: can input a frequency range that you want alpha returned for
                if not supplied will use self.f
            Input: (optional) set to wavelength to input wavelength in nm

        returns:
            absorption coefficients in cm^-1
        '''
        if f is None:
            f = self.f
        else:
            if Input == 'wavelength':
                f = self._wavelength2frequency(wavelength)

        self.alpha = getattr(self, self.model)(self.vals, f)

        return self.alpha

    def _EgwithT(self, Eg, T, gamma=None, beta=None):
        if beta is None:
            beta = self.beta

        if gamma is None:
            gamma = self.gamma
        return Eg - beta * T**2 / (T + gamma)

    def _checkf(self, f):
        if f is None:
            f = self.f
        return f

    def _wavelength2frequency(self, Lambda):
        '''
        Takes lambda in nm
        provides f in Hz
        '''
        self.f = const.c / Lambda / 1e-9
        return self.f

    def _alpha_function(self, f, Eth, A, power, T=300.):
        '''
        Generic function to determine alpha
        inputs:
            f is the photon energy in hz
            Eth is the threshold energy in eV
            A is the probability coefficient 
            power is the power its applied to

        '''

        # Use provided f or self.f
        f = self._checkf(f)

        # change E with T, assumes that the
        # band narrowing is not a function of Eg
        Eth = self._EgwithT(Eth, T)

        # calculate alpha
        alpha = A * (const.h * f / const.e - Eth)**power

        # make sure we don't have any impossible values
        alpha[const.h * f / const.e - Eth < 0] = 0

        return alpha

    def alpha_p_absorption(self, Eg, Ep, A, T, f=None):
        '''Phonon absorption, with E as power 2'''
        A /= (np.exp(Ep / const.k / T * const.e) - 1)
        return self._alpha_function(f, Eg - Ep, A, 2, T)

    def alpha_p_emission(self, Eg, Ep, A, T, f=None):
        '''Phonon emission, with E as power 2'''
        A /= (1 - np.exp(-Ep / const.k / T * const.e))
        return self._alpha_function(f, Eg + Ep, A, 2, T)

    def alpha_exciton_emission(self, Eg, Ep, Ee, A, T, f=None):
        '''Phonon and exciton emission, E power of 0.5'''
        A /= (1 - np.exp(-Ep / const.k / T * const.e))
        return self._alpha_function(f, Eg + Ep + Ee, A, 0.5, T)

    def alpha_exciton_absorption(self, Eg, Ep, Ee, A, T, f=None):
        '''Phonon and absorption emission, E power of 0.5
        This process may not be possible is high quality semiconductors'''
        A /= (1 - np.exp(-Ep / const.k / T * const.e))
        return self._alpha_function(f, Eg + Ep - Ee, A, 0.5, T)

    def alpha_direct(self, Eg, A, T=300, f=None, power=0.5):
        '''Direct band gap'''
        return self._alpha_function(f, Eg, A, power, T)

    def alpha_indirect(self, Eg, Ephonon, A, T):
        """
        The change in phonon absorption probability between emission 
        and absorption is just e^{E_p / kt}, taken from MacFarlane
        """

        alpha_in = np.zeros(self.f.shape)
        for Ep, Aa in zip(Ephonon, A):
            Ae = Aa * np.exp(Ep * const.e / const.k / T)
            alpha_in += self.alpha_p_absorption(Eg, Ep, Aa, T)
            alpha_in += self.alpha_p_emission(Eg, Ep, Ae, T)
        return alpha_in

    def MacFarlane(self, vals, f):

        # Eg = vals['egi']

        Egi_list = [vals[s] for s in vals if 'egi' in s]
        Ep_list = [vals[s] / 1000 for s in vals if 'ep' in s]
        Ap_list = [vals[s.replace('e', 'a')] for s in vals if 'ep' in s]

        Egd_list = [vals[s] for s in vals if 'egd' in s]
        Ad_list = [vals[s.replace('e', 'a')] for s in vals if 'ad' in s]

        # print Ep_list, Ap_list
        alpha = np.zeros(self.f.shape)

        for Eg in Egi_list:
            # print Eg
            alpha += self.alpha_indirect(Eg, Ep_list, Ap_list, 300.)

        for Eg in Egd_list:
            alpha += self.alpha_direct(Eg, Ad_list)

        return alpha

    def Rajkanan(self, vals, f):
        # Based on indirect theory from Elliot

        # print vals
        self.gamma = vals['gamma']
        self.beta = vals['beta']

        Egi_list = [vals[s] for s in vals if 'egi' in s]
        Ep_list = [vals[s] / 1000 for s in vals if 'ep' in s]
        Ap_list = [vals[s.replace('e', 'a')] for s in vals if 'ap' in s]
        C_list = [vals[s.replace('e', 'a')] for s in vals if 'c' in s]

        Egd_list = [vals[s] for s in vals if 'egd' in s]
        Ad_list = [vals[s.replace('e', 'a')] for s in vals if 'ad' in s]

        alpha = np.zeros(self.f.shape)
        T = 300.
        for Eg, Ap in zip(Egi_list, Ap_list):
            for Ep, C in zip(Ep_list,  C_list):

                alpha += self.alpha_p_absorption(Eg, Ep, Ap * C, T)
                alpha += self.alpha_p_emission(Eg, Ep, Ap * C, T)

        for Eg in Egd_list:
            alpha += self.alpha_direct(Eg, Ad_list)

        return alpha

    def Bucher(self, vals, f):
        # Based on indirect theory from Elliot

        # print vals
        self.gamma = vals['gamma']
        self.beta = vals['beta']

        Egi_list = [vals[s] for s in vals if 'egi' in s]
        Ep_list = [vals[s] / 1000 for s in vals if 'ep' in s]
        Ap_list = [vals[s.replace('e', 'a')] for s in vals if 'ap' in s]
        C_list = [vals[s.replace('e', 'a')] for s in vals if 'c' in s]

        Egd_list = [vals[s] for s in vals if 'egd' in s]
        Ad_list = [vals[s.replace('e', 'a')] for s in vals if 'ad' in s]

        alpha = np.zeros(self.f.shape)
        T = 300.
        for Eg, Ap in zip(Egi_list, Ap_list):
            for Ep, C in zip(Ep_list,  C_list):
                alpha += self.alpha_p_absorption(Eg, Ep, Ap * C, T)
                alpha += self.alpha_p_emission(Eg, Ep, Ap * C, T)

        power = vals['power']

        for Eg in Egd_list:
            alpha += self.alpha_direct(Eg,
                                       Ad_list[0] / const.h / f * const.e, 1.5)

        return alpha



