
import numpy as np
import matplotlib.pylab as plt
import os
import ConfigParser

from semiconductor.helper.helper import HelperFunctions
from semiconductor.general_functions.carrierfunctions import get_carriers
import radiative_models as radmdls
import auger_models as augmdls


class Intrinsic():

    def __init__(self, matterial='Si', rad_author=None, aug_author=None, **kwargs):

        self.Radiative = Radiative(matterial, rad_author, **kwargs)
        self.Auger = Auger(matterial, aug_author, **kwargs)

    def intrisic_carrier_lifetime(self, nxc, Na, Nd, inverse=False):
        itau = 1. / \
            self.Radiative.tau(nxc, Na, Nd) + 1. / \
            self.Auger.tau(nxc, Na, Nd)
        if inverse is False:
            itau = 1. / itau
        return itau


class Radiative(HelperFunctions):
    model_file = 'radiative.model'

    def __init__(self, matterial='Si', author=None, temp=300., ni=9.65e9):
        self.Models = ConfigParser.ConfigParser()
        self.matterial = matterial

        constants_file = os.path.join(
            os.path.dirname(os.path.realpath(__file__)),
            matterial,
            self.model_file)

        self.Models.read(constants_file)

        self.vals, self.model = self.change_model(author)

        'sets the temp for the thing'
        self.temp = temp

        'This is ni not the effective ni'
        self.ni = ni

    def tau(self, nxc, Na, Nd, temp=None):
        self.Nh_0, self.Ne_0 = self.check_doping(Na, Nd)

        if 'blow_model' in self.vals.keys():
            B = getattr(radmdls, self.blow_model)(self.vals, temp)
        else:
            B = self.vals['b']

        doping = np.amax([Na, Nd])
        return getattr(radmdls, self.model)(self.vals, nxc, self.Nh_0, self.Ne_0, temp, B=B)
#

    def itau(self, nxc, Na, Nd, temp=None):
        return 1. / self.tau(nxc, Na, Nd, temp)

    def B(self, nxc, doping, temp):
        if 'blow_model' in self.vals.keys():
            vals, model = self.change_model(self.vals['blow_model'])
            B = getattr(radmdls, model)(vals, temp)
            temp_tt = np.array([77, 90, 112, 170, 195, 249, 300])
            plt.figure()
            B_tt = np.array(
                [8.01e-14, 4.57e-14, 2.14e-14, 8.84e-15, 7.35e-15, 5.48e-15, 4.73e-15])
            plt.semilogy()
            plt.plot(temp_tt, B_tt, 'o')
            plt.plot(temp, B)
            plt.plot()
            # print np.vstack((B, temp)).T
            B = getattr(
                radmdls, self.model + '_B')(self.vals, nxc, doping, temp, B)
            # print B.shape
        else:
            B = self.vals['B']

        return B


class Auger(HelperFunctions):
    model_file = 'auger.model'

    def __init__(self, matterial='Si', author=None, temp=300, ni=9.65e9):

        self.Models = ConfigParser.ConfigParser()
        self.matterial = matterial

        constants_file = os.path.join(
            os.path.dirname(os.path.realpath(__file__)),
            matterial,
            self.model_file)

        self.Models.read(constants_file)

        self.vals, self.model = self.change_model(author)
        self.temp = temp
        self.ni = ni

    def tau(self, nxc, Na, Nd, temp=300):

        ne0, nh0 = get_carriers(
            Na, Nd, nxc=0, ni_author=None, temp=temp)
        print self.vals
        return getattr(augmdls, self.model)(self.vals, nxc, ne0, nh0)

    def itau_aug(self, nxc, Na, Nd):
        return 1 / self.tau(nxc, Na, Nd)

    def check(self, author):
        fig, ax = plt.subplots(1)
        vals, model = self.change_model(author)
        func = getattr(augmdls, self.model)

        getattr(augmdls, author + '_check')(vals, func, fig, ax)


