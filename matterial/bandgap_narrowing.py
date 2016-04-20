#!/usr/local/bin/python
# UTF-8

import numpy as np
import matplotlib.pylab as plt
import os
import ConfigParser

from semiconductor.helper.helper import HelperFunctions
import semiconductor.matterial.bandgap_narrowing_models as Bgn
import semiconductor.general_functions.carrierfunctions as GF


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

    author_list = 'bandgap_narrowing.models'

    def __init__(self, matterial='Si', author=None, ni_author=None):
        self.Models = ConfigParser.ConfigParser()
        self.matterial = matterial

        constants_file = os.path.join(
            os.path.dirname(os.path.realpath(__file__)),
            matterial,
            self.author_list)

        self.Models.read(constants_file)
        self.vals, self.model = self.change_model(author)

    def update_BGN(self, Na, Nd, min_car_den=None, author=None, temp=300, ni_author=None):
        '''
        Calculates the band gap narrowing

        Inputs:
        Na, Nd, delta n, temp, ni

        output:
            band gap narrowing in eV
        '''

        # this should be change an outside function alter
        ne, nh = GF.get_carriers(Na,
                                 Nd,
                                 min_car_den,
                                 temp=temp)

        doping = np.abs(Na - Nd)

        if author is not None:
            self.vals, self.model = self.change_model(author)

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
        temp = 300

        for author in self.available_models():
            BGN = self.update_BGN(Na=Na, Nd=Nd, min_car_den=dn,
                                  author=author,
                                  temp=temp)

            if not np.all(BGN == 0):
                plt.plot(Na, BGN, label=author)

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


def check_Schenk(fig, ax):
    '''compared to values taken from XXX'''
    BGN = BandGapNarrowing(matterial='Si', author='Schenk1988fer')

    folder = os.path.join(
        os.path.dirname(__file__), 'Si', r'check data')

    fnames = ['BGN_Schenk_asN-dn-1e14.csv']
    min_car_den = 1e14

    ax.set_color_cycle(['c', 'c', 'm', 'm', 'b', 'b', 'r', 'r', 'g', 'g'])

    for f_name in fnames:
        data = np.genfromtxt(os.path.join(folder, f_name),
                             names=True,
                             delimiter=',',
                             skip_header=1)
        ND = np.zeros(data.shape)
        for temp in data.dtype.names[1::2]:
            bgn = BGN.update_BGN(data['N'], ND, min_car_den, temp=float(temp))
            ax.plot(data['N'], bgn,
                    '.')
            ax.plot(data['N'], data[temp],
                    '--',
                    label=temp)

        ax.legend(loc=0, title='Temperature (K)')

    ax.set_title('BGN comparison to PV-lighthouse: $\Delta$n=1e14:')
    ax.set_ylabel('Bang gap narrowing (eV)')
    ax.set_xlabel('Ionised Doping (cm$^{-3}$)')
    ax.semilogx()
