#!/usr/local/bin/python
# UTF-8
import numpy as np
import matplotlib.pylab as plt
import sys
import os
import ConfigParser
import mobilitymodels as model

sys.path.append(
    os.path.abspath(os.path.join(os.getcwd(), os.pardir, os.pardir)))

from semiconductor.helper.helper import HelperFunctions


class Mobility(HelperFunctions):
    model_file = 'mobility.models'
    ni = 1e10
    temp = 300

    def __init__(self, matterial='Si', model_author=None, temp=300., ni=9.65e9):
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
        impurity = Na + Nd

        return getattr(model, self.model)(
            self.vals, impurity, min_car_den,  ni=self.ni, carrier='electron', **kwargs)

    def hole_mobility(self, min_car_den, Na, Nd, **kwargs):

        self.Nh_0, self.Ne_0 = self.check_doping(Na, Nd)
        impurity = Na + Nd

        return getattr(model, self.model)(
            self.vals, impurity, min_car_den, ni=self.ni, carrier='hole', **kwargs)

    def mobility_sum(self, min_car_den, Na, Nd, **kwargs):

        return self.hole_mobility(min_car_den, Na, Nd, **kwargs) +\
            self.electron_mobility(min_car_den, Na, Nd, **kwargs)


if __name__ == "__main__":

    a = Mobility('Si')

    # a.change_model('dorkel1981')

    dn = np.logspace(10, 20)
    # dn = np.array([1e14])
    Nd = 1e14
    Na = 0
    # mob_e = model.dorkel(a.vals, Nd, dn, dn + Nd, 300., carrier='electron')

    # plt.plot(dn, mob_e)
    # # plt.plot(dn, mob_h)

    # a.change_model('dorkel1981')

    # print dn.shape, a.hole_mobility(dn, Nd, 0, maj_car_den=1e16, temp=300).shape
    # plt.plot(dn, a.hole_mobility(
    #     dn, Nd, 0, maj_car_den=1e16 + dn, temp=300), '--')
    # plt.plot(dn, a.electron_mobility(
    #     dn, Nd, 0, maj_car_den=1e16 + dn, temp=300), '--')
    # plt.ylim(bottom=0)

    a.change_model('klaassen1992')
    # (self, vals, p0, n0, dn, temp, carrier)
    plt.plot(dn, a.hole_mobility(dn, Nd, Na, N_d=Nd, N_a=Na, temp=300))
    plt.plot(dn, a.electron_mobility(dn, Nd, Na, N_d = Nd, N_a = Na, temp=300))

    folder = os.path.join(os.getcwd(), 'Si')
    fname = 'Klassen_1e14_dopants.dat'
    data = np.genfromtxt(os.path.join(folder, fname), names=True)
    # print data.dtype.names
    plt.plot(data['deltan'], data['ue'], '--',
             label='electron - PV-lighthouse')
    plt.plot(data['deltan'], data['uh'], '--', label='hole - PV-lighthouse')
    plt.legend(loc=0)


    print model.uLS('e', a.vals, 300), model.uLS('h', a.vals, 300)
    print model.un('e', a.vals, 300), model.un('h', a.vals, 300)
    print model.uc('e', a.vals), model.uc('h', a.vals)


    # plt.figure('components')
    # plt.plot(dn, model.Nsc('e', a.vals,  Nd+dn, dn, 0, Nd))
    # plt.plot(dn, model.Nsc('h', a.vals,  Nd+dn, dn, 0, Nd))

    # plt.plot(dn, model.Nsceff('e', a.vals,Nd+dn , dn, Nd, 0, 300))
    # plt.plot(dn, model.Nsceff('h', a.vals,Nd+dn , dn, Nd, 0, 300))

    # plt.loglog()

    plt.semilogx()
    plt.show()
