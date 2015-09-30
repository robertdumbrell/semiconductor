import numpy as np
import matplotlib.pylab as plt
import sys
import os
import ConfigParser
import scipy.constants as Const


class absorptioncoefficient():

    ''' This purpose of this it to provide acess if the absorption
        coefficient have a model'''

        def __init__(self, matterial='Si', model_author=None):
        self.Models = ConfigParser.ConfigParser()
        self.matterial = matterial

        constants_file = os.path.join(
            os.path.dirname(os.path.realpath(__file__)),
            matterial,
            'absorptioncoefficient.const')

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

    def update_absorptioncoefficients(self, wavelength=None):
        if wavelength is None:
            wavelength = self.wavelength
        self.Eg = getattr(self, self.model)(self.vals, wavelength)

        return self.Eg

    def alpha_Macfarlane(self, vals, wavelength):
        '''This has not been added properly, 
        but will be taken from Macfarlane1958, as has been used by
        mcdonald at eu pvsec 2015'''

        pass
