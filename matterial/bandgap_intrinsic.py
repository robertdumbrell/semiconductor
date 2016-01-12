
import matplotlib.pylab as plt
import os
import ConfigParser
import numpy as np
import semiconductor.matterial.intrinsic_bandgap_models as iBg
from semiconductor.helper.helper import HelperFunctions


class IntrinsicBandGap(HelperFunctions):

    '''
    The intrinsic band-gap as a function of temperature
        it changes as a result of:
             different effective carrier mass (band strucutre)
    '''

    model_file = 'bandgap.models'

    def __init__(self, matterial='Si', model=None):
        self.Models = ConfigParser.ConfigParser()
        self.matterial = matterial

        constants_file = os.path.join(
            os.path.dirname(os.path.realpath(__file__)),
            matterial,
            self.model_file)

        self.Models.read(constants_file)

        self.change_model(model)

    def update_iEg(self, temp=None, model=None, multiplier=1.01):
        '''
        a function to update the intrinsic BandGap

        inputs:
            temperature in kelvin
            model: (optional)
                  the model used.
                  If not provided the last provided model is used
                  If no model has been provided Passler model is used
            multiplier: A band gap multipler. 1.01 is suggested.

        output:
            the intrinsic bandgap in eV
        '''

        if temp is None:
            temp = self.temp
        if model is not None:
            self.change_model(model)

        Eg = getattr(iBg, self.model)(self.vals, temp)

        return Eg * multiplier

    def check_models(self):
        '''
        Displays a plot of the models against that taken from a
        respected website (https://www.pvlighthouse.com.au/)
        '''
        plt.figure('Intrinsic bandgap')
        t = np.linspace(1, 500)

        for model in self.available_models():

            Eg = self.update_iEg(t, model=model, multiplier=1.0)
            plt.plot(t, Eg, label=model)

        test_file = os.path.join(
            os.path.dirname(os.path.realpath(__file__)),
            'Si', 'check data', 'iBg.csv')

        data = np.genfromtxt(test_file, delimiter=',', names=True)

        for temp, name in zip(data.dtype.names[0::2], data.dtype.names[1::2]):
            plt.plot(
                data[temp], data[name], '--', label=name)

        plt.xlabel('Temperature (K)')
        plt.ylabel('Intrinsic Bandgap (eV)')

        plt.legend(loc=0)
        self.update_iEg(0, model=model, multiplier=1.01)
