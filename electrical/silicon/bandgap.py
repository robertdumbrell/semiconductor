import numpy
from pylab import *
import BlackBody
from glob import glob
import sys
import os
import scipy.constants as Const


class BandGap_Si():

    """
    Contains analtical models for silicons band gap.
    A comparison is made in Couderc2014
    """

    E0 = 'Passler'

    def E0(self, temp=False, E0=False):
        """
        just a function that provide E0 based on the E0 provided
        If no E0 is provided it defults to the self.E0 value
        """

        if not E0:
            E0 = self.E0
        else:
            self.E0 = E0

        self.E = getattr(self, 'E0_' + E0)()
        return self.E

    def E0_Thurmond(self, temp=False):
        """
        Doesn't work, gives too high numbers
        Taken from Couderc2014
        Think its a mistake in the thta paper
        """

        if not numpy.all(temp):
            temp = self.Temp

        E0 = 1.17

        alpha = 4.73e-4
        beta = 636.

        E = E0 - alpha * temp**2 / (temp + beta)
        # print self.E0_Passler(temp)/1.602e-19,E
        return E * Const.e

    def E0_Passler(self, temp=False):
        """
        taken from Couderc2014
        """

        if not numpy.all(temp):
            temp = self.Temp

        E0 = 1.17
        alpha = 3.23e-4
        theta = 446.
        delta = 0.51
        gamma = (1. - 3. * delta**2) / (numpy.exp(theta / temp) - 1)
        xi = 2. * temp / theta

        # Values for each sum component
        No2 = numpy.pi**2. * xi**2. / (3. * (1 + delta**2))
        No3 = (3. * delta**2 - 1) / 4. * xi**3
        No4 = 8. / 3. * xi**4.
        No5 = xi**6.

        E = E0 - alpha * theta * \
            (gamma + 3. * delta**2 / 2 *
             ((1. + No2 + No3 + No4 + No5)**(1. / 6.) - 1))
        return E * Const.e