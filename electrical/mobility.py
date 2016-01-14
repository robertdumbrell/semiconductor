#!/usr/local/bin/python
# UTF-8

import numpy as np
import matplotlib.pylab as plt
import os
import ConfigParser
import mobilitymodels as model

from semiconductor.helper.helper import HelperFunctions


class Mobility(HelperFunctions):
    model_file = 'mobility.models'
    ni = 1e10
    temp = 300

    def __init__(self, matterial='Si', model_author=None, temp=300.):
        self.Models = ConfigParser.ConfigParser()
        self.matterial = matterial

        constants_file = os.path.join(
            os.path.dirname(os.path.realpath(__file__)),
            matterial,
            self.model_file)

        self.Models.read(constants_file)

        self.change_model(model_author)

    def electron_mobility(self, min_car_den, Na, Nd, **kwargs):

        self.Nh_0, self.Ne_0 = self.check_doping(Na, Nd)

        return getattr(model, self.model)(
            self.vals, Na, Nd, min_car_den, carrier='electron', **kwargs)

    def hole_mobility(self, min_car_den, Na, Nd, **kwargs):

        self.Nh_0, self.Ne_0 = self.check_doping(Na, Nd)

        return getattr(model, self.model)(
            self.vals, Na, Nd, min_car_den, carrier='hole', **kwargs)

    def mobility_sum(self, min_car_den, Na, Nd, **kwargs):

        return self.hole_mobility(min_car_den, Na, Nd, **kwargs) +\
            self.electron_mobility(min_car_den, Na, Nd, **kwargs)

    def check_models(self):
        check_klaassen()
        check_dorkel()


def check_klaassen():
    '''compares to values taken from www.PVlighthouse.com.au'''

    a = Mobility('Si')
    a.change_model('klaassen1992')

    dn = np.logspace(10, 20)
    # dn = np.array([1e14])
    Nd = 1e14
    Na = 0

    folder = os.path.join(
        os.path.dirname(__file__), 'Si', 'test_mobility_files')
    fnames = ['Klassen_1e14_dopants.dat',
              'Klassen_1e14_temp-450.dat']

    for temp, f_name in zip([300, 450], fnames):

        plt.figure('Mobility - Klaassen: Deltan at ' + str(temp))

        plt.plot(dn, a.hole_mobility(dn, Na, Nd, temp=temp),
                 'r-',
                 label='hole-here')
        plt.plot(dn, a.electron_mobility(dn, Na, Nd, temp=temp),
                 'b-',
                 label='electron-here')

        data = np.genfromtxt(os.path.join(folder, f_name), names=True)

        plt.plot(data['deltan'], data['uh'], 'b--',
                 label='hole - PV-lighthouse')
        plt.plot(data['deltan'], data['ue'], 'r--',
                 label='electron - PV-lighthouse')
        plt.legend(loc=0, title='Mobility from')

        plt.semilogx()
        plt.xlabel(r'$\Delta$n (cm$^{-3}$)')
        plt.xlabel(r'Moblity  (cm$^2$V$^{-1}$s$^{-1}$)')


def check_dorkel():
    '''compares to values taken from www.PVlighthouse.com.au'''

    a = Mobility('Si')
    a.change_model('dorkel1981')

    dn = np.logspace(10, 20)
    # dn = np.array([1e14])
    Nd = 1e14
    Na = 0

    folder = os.path.join(
        os.path.dirname(__file__), 'Si', 'test_mobility_files')
    fnames = ['dorkel_1e14_carriers.dat',
              'dorkel_1e14_temp-450.dat']

    for temp, f_name in zip([300, 450], fnames):

        plt.figure('Mobility - Dorkel: Deltan at ' + str(temp))

        plt.plot(dn, a.hole_mobility(dn, Na, Nd, temp=temp),
                 'r-',
                 label='hole-here')
        plt.plot(dn, a.electron_mobility(dn, Na, Nd, temp=temp),
                 'b-',
                 label='electron-here')

        data = np.genfromtxt(os.path.join(folder, f_name), names=True)

        plt.plot(data['deltan'], data['uh'], 'b--',
                 label='hole - PV-lighthouse')
        plt.plot(data['deltan'], data['ue'], 'r--',
                 label='electron - PV-lighthouse')
        plt.legend(loc=0, title='Mobility from')

        plt.semilogx()
        plt.xlabel(r'$\Delta$n (cm$^{-3}$)')
        plt.xlabel(r'Moblity  (cm$^2$V$^{-1}$s$^{-1}$)')

if __name__ == "__main__":
    check_klaassen()
    check_dorkel()
    plt.show()
