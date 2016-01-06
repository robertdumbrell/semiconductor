
import numpy as np
from pylab import *
import sys
import os
import scipy.constants as Const
import ConfigParser
from bandgap import IntrinsicBandGap

from semiconductor.helper.helper import HelperFunctions
from semiconductor.matterial import ni_models


class IntrinsicCarrierDensity(HelperFunctions):

    '''
    The intrinisc carrier density is the number of carriers
    that exist the a matterial at thermal equlibrium.
    It is impacted by the band gap (and bandgap narrowing)
    '''

    # ni = 1e10
    temp = 300.
    model_file = 'ni.models'

    def __init__(self, matterial='Si', model_author=None, temp=300.):
        self.Models = ConfigParser.ConfigParser()
        self.matterial = matterial

        constants_file = os.path.join(
            os.path.dirname(os.path.realpath(__file__)),
            matterial,
            self.model_file)

        self.Models.read(constants_file)

        self.change_model(model_author)
        self.temp = temp

    def update_ni(self, temp=None, model=None):
        '''
        a function to update the intrinsic BandGap

        inputs:
            temperature: (optional)
                         in kelvin
            model:  (optional)
                    the model used.
                    If not provided the last provided model is used
                    If no model has been provided Passler model is used
        output:
            the intrinsic carrier concentration in cm^-3. This is not the effectice
            carrier concentration
        '''
        # able to input a temperature to change
        if temp is None:
            temp = self.temp

        # a check to make sure the model hasn't changed
        if model is not None:
            self.change_model(model)

        # if the model required the energy gap, caculate it
        if self.model == 'ni_temp_eg':
            Eg = IntrinsicBandGap(self.matterial, self.vals['eg_model']
                                  ).update_iEg(temp=temp, multiplier=1)
        else:
            Eg = 0

        self.ni = getattr(ni_models, self.model)(
            self.vals, temp=temp, Eg=Eg)

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

        Eg = IntrinsicBandGap('Si', 'Passler2002'
                              ).update_iEg(temp=1, multiplier=1)
        Eg = 1.17

        for model in self.available_models():
            ni = self.update_ni(temp, model=model)
            ax.plot(np.log(temp),
                    np.log(ni * np.exp(Eg / 2. * Const.e / Const.k / temp)),
                    label=model)
            print model, '\t {0:.2e}'.format(self.update_ni(300, model=model))

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

