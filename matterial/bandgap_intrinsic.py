
import matplotlib.pylab as plt
import os
import ConfigParser
import numpy as np
import semiconductor.matterial.bandgap_intrinsic_models as iBg
from semiconductor.helper.helper import HelperFunctions


class IntrinsicBandGap(HelperFunctions):

    '''
    The intrinsic band-gap as a function of temperature
        it changes as a result of:
             different effective carrier mass (band strucutre)
    '''
    cal_dts = {
        'matterial': 'Si',
        'temp': 300.,
        'author': None,
        'multiplier' : 1.,
    }
    author_list = 'bandgap.models'

    def __init__(self, **kwargs):

        # update any values in cal_dts
        # that are passed
        temp = locals().copy()
        del temp['self']
        self._update_dts(**temp)

        # get the address of the authors list
        author_file = os.path.join(
            os.path.dirname(os.path.realpath(__file__)),
            self.cal_dts['matterial'],
            self.author_list)

        # get the models ready
        self._int_model(author_file)

        # initiate the first model
        self.change_model(self.cal_dts['author'])

    def update(self, **kwargs):
        '''
        a function to update the intrinsic BandGap

        inputs:
            temperature in kelvin
            author: (optional)
                  the author used.
                  If not provided the last provided author is used
                  If no author has then the author Passler's is used
            multiplier: A band gap multipler. 1.01 is suggested.

        output:
            the intrinsic bandgap in eV
        '''

        self._update_dts(**kwargs)

        if 'author' in kwargs.keys():
            self.change_model(self.cal_dts['author'])
        Eg = getattr(iBg, self.model)(self.vals, temp=self.cal_dts['temp'])

        return Eg * self.cal_dts['multiplier']

    def check_models(self):
        '''
        Displays a plot of the models against that taken from a
        respected website (https://www.pvlighthouse.com.au/)
        '''
        plt.figure('Intrinsic bandgap')
        t = np.linspace(1, 500)

        for author in self.available_models():

            Eg = self.update(temp=t, author=author, multiplier=1.0)
            plt.plot(t, Eg, label=author)

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
        self.update(temp=0, author=author, multiplier=1.01)
