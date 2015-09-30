import numpy as np
from pylab import *
import sys
import os
import scipy.constants as Const
import ConfigParser

class BandGap():
    Eg = 1.12
    temp = 300

    def __init__(self, matterial='Si', model_author=None):
        self.Models = ConfigParser.ConfigParser()
        self.matterial = matterial

        constants_file = os.path.join(
            os.path.dirname(os.path.realpath(__file__)),
            matterial,
            'bandgap.const')

        self.Models.read(constants_file)

        self.change_model(model_author)

    def change_model(self, model_author=None):

        if model_author is None:
            self.model_author = self.Models.get('default', 'model')
        else:
            # Need a check to make sure craNhcan't be passed
            self.model_author = model_author
        # print self.model_author
        self.model = self.Models.get(self.model_author, 'model')

        self.vals = dict(self.Models.items(self.model_author))
        # List.remove('model')
        del self.vals['model']
        self.vals = {k: float(v) for k, v in self.vals.iteritems()}

    # def E0(self, temp=False, E0=False):
    #     """
    #     just a function that provide E0 based on the E0 provided
    #     If no E0 is provided it defults to the self.E0 value
    #     """

    #     if not E0:
    #         E0 = self.E0
    #     else:
    #         self.E0 = E0

    #     self.E = getattr(self, 'E0_' + E0)()
    #     return self.E

    # def E0_Thurmond(self, temp=False):
    #     """
    #     Doesn't work, gives too high numbers
    #     Taken from Couderc2014
    #     Think its a mistake in the thta paper
    #     """

    #     if not np.all(temp):
    #         temp = self.Temp

    #     E0 = 1.17

    #     alpha = 4.73e-4
    #     beta = 636.

    #     E = E0 - alpha * temp**2 / (temp + beta)
    # print self.E0_Passler(temp)/1.602e-19,E
    #     return E * Const.e

    def update_Eg(self, temp=None):
        if temp is None:
            temp = self.temp
        self.Eg = getattr(self, self.model)(self.vals, temp)
        
        return self.Eg

    def Eg_Passler(self, vals, temp):
        """
        taken from Couderc2014
        """

        gamma = (1. - 3. * vals['delta']**2) / (np.exp(vals['theta'] / temp) - 1)
        xi = 2. * temp / vals['theta']

        # Values for each sum component
        No2 = np.pi**2. * xi**2. / (3. * (1 + vals['delta']**2))
        No3 = (3. * vals['delta']**2 - 1) / 4. * xi**3
        No4 = 8. / 3. * xi**4.
        No5 = xi**6.

        E = vals['e0'] - vals['alpha'] * vals['theta'] * \
            (gamma + 3. * vals['delta']**2 / 2 *
             ((1. + No2 + No3 + No4 + No5)**(1. / 6.) - 1))
        return E * Const.e


