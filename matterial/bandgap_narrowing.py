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
    cal_dts = {
        'matterial': 'Si',
        'temp': 300.,
        'author': None,
    }

    author_list = 'bandgap_narrowing.models'

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

    def update(self, Na, Nd, nxc=None, **kwargs):
        '''
        Calculates the band gap narrowing

        Inputs:
        Na, Nd, delta n, temp, ni

        output:
            band gap narrowing in eV
        '''

        self._update_dts(**kwargs)

        # a check to make sure the model hasn't changed
        if 'author' in kwargs.keys():
            self.change_model(self.cal_dts['author'])

        # this should be change an outside function alter
        ne, nh = GF.get_carriers(Na,
                                 Nd,
                                 nxc,
                                 temp=self.cal_dts['temp'])

        doping = np.abs(Na - Nd)

        return getattr(Bgn, self.model)(
            self.vals,
            Na=np.copy(Na),
            Nd=np.copy(Nd),
            ne=ne,
            nh=nh,
            temp=self.cal_dts['temp'],
            doping=doping)

    def check_models(self):
        plt.figure('Bandgap narrowing')
        Na = np.logspace(12, 20)
        Nd = 0
        dn = 1e14
        temp = 300.

        for author in self.available_models():
            BGN = self.update(Na=Na, Nd=Nd, nxc=dn,
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
    nxc = 1e14

    ax.set_color_cycle(['c', 'c', 'm', 'm', 'b', 'b', 'r', 'r', 'g', 'g'])

    for f_name in fnames:
        data = np.genfromtxt(os.path.join(folder, f_name),
                             names=True,
                             delimiter=',',
                             skip_header=1)
        ND = np.zeros(data.shape)
        for temp in data.dtype.names[1::2]:
            bgn = BGN.update(data['N'], ND, nxc, temp=float(temp))
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
