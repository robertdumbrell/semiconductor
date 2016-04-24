#!/usr/local/bin/python
# UTF-8

import numpy as np
from pylab import *
import os
import scipy.constants as Const
from bandgap_intrinsic import IntrinsicBandGap
from semiconductor.helper.helper import HelperFunctions
from semiconductor.matterial import ni_models


class IntrinsicCarrierDensity(HelperFunctions):

    '''
    The intrinisc carrier density is the number of carriers
    that exist the a matterial at thermal equlibrium.
    It is impacted by the band gap

    The effective intrinsic carrier density refers to a modification
    of the intrinsic to account for band gap narrowing.
    '''

    cal_dts = {
        'matterial': 'Si',
        'temp': 300.,
        'author': None,
    }

    author_list = 'ni.models'

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
            temperature: (optional)
                         in kelvin
            author:  (optional)
                    the author used.
                    If not provided the last provided author is used
                    If no author has been provided,  Couderc's model is used
        output:
            The intrinsic carrier concentration in cm^-3.
            Note this is not the effectice intrinsic carrier concentration
        '''
        self._update_dts(**kwargs)

        # a check to make sure the model hasn't changed
        if 'author' in kwargs.keys():
            self.change_model(self.cal_dts['author'])

        # if the model required the energy gap, caculate it
        if self.model == 'ni_temp_eg':
            Eg = IntrinsicBandGap(
                matterial=self.cal_dts['matterial'],
                author=self.vals['eg_model']
            ).update(
                temp=self.cal_dts['temp'],
                multiplier=1
            )
        else:
            Eg = 0

        self.ni = getattr(ni_models, self.model)(
            self.vals, temp=self.cal_dts['temp'], Eg=Eg)

        return self.ni

    def check_models(self):
        '''
        Displays a plot of all the models against experimental data
        '''
        # fig = plt.figure('Intrinsic carriers')
        fig, ax = plt.subplots(1)
        fig.suptitle('Intrinsic carrier concentration')
        # fig.set_title('Intrinsic carriers')
        temp = np.linspace(100, 500)

        Eg = IntrinsicBandGap(matterial='Si', author='Passler2002'
                              ).update(temp=1, multiplier=1)
        Eg = 1.17

        for author in self.available_models():
            ni = self.update(temp=temp, author=author)
            ax.plot(np.log(temp),
                    np.log(ni * np.exp(Eg / 2. * Const.e / Const.k / temp)),
                    label=author)
            print author, '\t {0:.2e}'.format(self.update(temp=300, author=author))

        test_file = os.path.join(
            os.path.dirname(os.path.realpath(__file__)),
            'Si', 'check data', 'ni.csv')

        data = np.genfromtxt(
            test_file,
            delimiter=',', names=True, skip_header=1, filling_values=np.inf)
        for name in data.dtype.names[1:]:
            ax.plot(np.log(data['Temp']),
                    np.log(data[name] *
                           np.exp(Eg / 2. * Const.e / Const.k / data['Temp'])),
                    'o', label='experimental values\'s: ' + name)

        ax.set_xlabel('log(Temperature (K))')
        ax.set_ylabel(r'$log(n_i \times e^{Eg_0(0)/kT}  )$')

        ax.legend(loc=0)
        ax.set_xlim(4, 6)
        ax.set_ylim(42, 47)
        plt.show()
