import numpy as np
import matplotlib.pylab as plt
import os
import ConfigParser

from semiconductor.helper.helper import HelperFunctions
import semiconductor.matterial.bandgap_narrowing_models as Bgn
from semiconductor.matterial.ni import IntrinsicCarrierDensity as ni


class BandGapNarrowing(HelperFunctions):

    '''
    Bang gap narrowing accounts for a reduction in bandgap that
    occurs as a result from no thermal effects. These include:
        doping 
        excess carrier density (non thermal distribution)

    As it depends ont eh excess carriers, it also depends on the 
    intrinsic carrier density. 
    Note: I currently believed that the impact of dopants
        is much larger than the impact of the carrier distribution
    '''

    model_file = 'bandgap_narrowing.models'

    def __init__(self, matterial='Si', model=None, ni_model=None):
        self.Models = ConfigParser.ConfigParser()
        self.matterial = matterial

        constants_file = os.path.join(
            os.path.dirname(os.path.realpath(__file__)),
            matterial,
            self.model_file)

        self.Models.read(constants_file)
        self.change_model(model)

        self.ni = ni(matterial, model=ni_model, temp=300)

    def update_BGN(self, Na, Nd, min_car_den=None, model=None, temp=300, ni_model=None):
        '''
        Calculates the band gap narrowing

        Inputs:
        Na, Nd, delta n, temp, ni

        output:
            band gap narrowing in eV
        '''

        # this should be change an outside function alter

        ne, nh = self.get_carriers(Na,
                                   Nd,
                                   min_car_den,
                                   self.ni.update_ni(temp=temp,
                                                     model=ni_model))
        print self.ni.update_ni(temp=temp, model=ni_model)
        doping = np.abs(Na - Nd)

        if model is not None:
            self.change_model(model)

        return getattr(Bgn, self.model)(
            self.vals,
            Na=np.copy(Na),
            Nd=np.copy(Nd),
            ne=ne,
            nh=nh,
            temp=temp,
            doping=doping)

    def check_models(self):
        plt.figure('Bandgap narrowing')
        Na = np.logspace(12, 20)
        Nd = 0
        dn = 1e14
        temp = 550

        for model in self.available_models():
            BGN = self.update_BGN(Na=Na, Nd=Nd, min_car_den=dn,
                                  model=model,
                                  temp=temp)

            if not np.all(BGN == 0):
                plt.plot(Na, BGN, label=model)

        test_file = os.path.join(
            os.path.dirname(os.path.realpath(__file__)),
            'Si', 'check data', 'Bgn.csv')

        data = np.genfromtxt(test_file, delimiter=',', names=True)

        for name in data.dtype.names[1:]:
            plt.plot(
                data['N'], data[name], 'r--',
                label='PV-lighthouse\'s: ' + name)

        plt.semilogx()
        plt.xlabel('Doping (cm$^{-3}$)')
        plt.ylabel('Bandgap narrowing (K)')

        plt.legend(loc=0)
