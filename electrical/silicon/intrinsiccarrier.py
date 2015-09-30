import numpy
from pylab import *
import BlackBody
from glob import glob
import sys
import os
import scipy.constants as Const



class IntrinsicCarrierDensity_Si(BandGap_Si):

    """
    Place that holds all the intensic carrier densties formulatsion for silicon
    """
    Temp = 300
    ni_version = 'Sproul_Couderc2014'

    def __init__(self):

        self.load_ni()

    def load_ni(self, ni_version=False, temp=False):
        """
        Obtains ni at a given termpature. Different ni_version can be provided
        If no ni_version is provided it defults to the self.ni_version value
        """

        if not ni_version:
            ni_version = self.ni_version

        if not numpy.all(temp):
            temp = self.Temp
        data = getattr(self, 'ni_' + ni_version)(temp)

    def ni_Misiakos1993(self, temp):
        """
        Capactiance measuremnts of ni. Not impact by bandgap narrowing
        Comment on in Couderc2014
        Couderc2014 paper doesn't normalise to 300
        """

        self.ni = 5.29e19 * (temp / 300)**2.54 * numpy.exp(-6726. / temp)
        # self.ni = 50.27e14*(temp/300)**2.54*numpy.exp(-6726/temp)
        return self.ni

    def ni_Sproul_1991(self, temp):
        """
        ni from spectral response. from 275 K to 375 K
        Did not account for bandgap narrowing
        """
        # taken from Couderc2014 paper

        if not numpy.all(temp):
            temp = self.Temp

        # print self.E0_Passler(temp)/1.6e-19
        # print self.E0_Thurmond(temp)/1.6e-19
        self.ni = 10.2e14 * (temp)**(2) * numpy.exp(-6880. / Const.k / temp)
        return self.ni

    def ni_Sproul_1993(self, temp):
        """
        ni from spectral response from 77 -300
        Did not account for bandgap narrowing
        """
        # taken from Couderc2014 paper

        self.ni = 16.4e14 * (temp)**(1.706)
        self.ni *= numpy.exp(-self.E0_Passler(temp) / 2. / Const.k / temp)
        return self.ni

    def ni_Sproul_Couderc2014(self, temp):
        """
        Revaluation of sproul account for ni accounting for band bap narrowing
        Band gap model was Schenk1998
        Same work as  Altermatt2003a but over wide temp range
        claims validity for 50 K to 400 K
        """
        # print self.E0_Passler(temp)/1.6e-19
        # print self.E0_Thurmond(temp)/1.6e-19
        self.ni = 1.541e15 * (temp)**(1.712)
        self.ni *= numpy.exp(-self.E0_Passler(temp) / 2. / Const.k / temp)
        return self.ni
