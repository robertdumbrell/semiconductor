
import numpy as np
import matplotlib.pylab as plt
import sys
import os
import ConfigParser


sys.path.append(
    os.path.abspath(os.path.join(os.getcwd(), os.pardir, os.pardir)))
import semiconductor.defults
import semiconductor.matterial.ni as niclass
from semiconductor.helper.helper import HelperFunctions

class Helper():

    def tau():
        pass

    def _PlotAll(self):
        fig, ax = plt.subplots(1)
        # ax = plt.add_subplot(111)
        for model in self.AvailableModels():
        # ax.plot(np.inf,np.inf,'k-',label = 'Auger')
            self.change_model(model)
            tau = self.tau(self.min_car_den, 1e16, 0)
            if tau is not None:
                ax.plot(self.min_car_den, tau * 1e6, label=model)

        ax.legend(loc=0)
        ax.loglog()
        ax.set_xlabel('$\Delta$ n (cm$^{-3}$)')
        ax.set_ylabel('Lifetime (us)')

        # Helper routiens



class Intrinsic():

    def __init__(self, matterial='Si', rad_model_author=None, aug_model_author=None, **kwargs):

        self.Radiative = Radiative(matterial, rad_model_author, **kwargs)

        self.Auger = Auger(matterial, aug_model_author, **kwargs)

    def intrisic_carrier_lifetime(self, min_car_den, Na, Nd):

        return 1. / (1. / self.Radiative.tau(min_car_den, Na, Nd) + 1. / self.Auger.tau(min_car_den, Na, Nd))


class Radiative(HelperFunctions):
    model_file = 'radiative.model'

    def __init__(self, matterial='Si', model_author=None, temp=300., ni=9.65e9):
        self.Models = ConfigParser.ConfigParser()
        self.matterial = matterial

        constants_file = os.path.join(
            os.path.dirname(os.path.realpath(__file__)),
            matterial,
            self.model_file)

        self.Models.read(constants_file)

        self.change_model(model_author)
        
        'sets the temp for the thing'
        self.temp = temp

        'This is ni not the effective ni'
        self.ni = ni


    def tau(self, min_car_den, Na, Nd):
        self.Nh_0, self.Ne_0 = self.check_doping(Na, Nd)
        doping = np.max(Na, Nd)
        return getattr(self, self.model)(min_car_den, doping)

    def itau(self, min_car_den, Na, Nd):
        return 1. / self.tau(min_car_den, Na, Nd)

    def Roosbroeck(self,  min_car_den, doping, B=None):
        if B is None:
            B = self.Models.getfloat(self.model_author, 'B')

        Nh = self.Nh_0 + min_car_den
        Ne = self.Ne_0 + min_car_den

        R = B * (Ne * Nh - self.Ne_0 * self.Nh_0)
        return min_car_den / R

    def Roosbroeck_with_screening(self, min_car_den, doping):
        """ 
        This is the roosbroeck model that accounts for many things
        It needs temperature, min_car_den, doping and blow to be defined
        """

        bmin = self.vals['rmax'] + (self.vals['rmin'] - self.vals['rmax']) / (
            1. + (self.temp / self.vals['r1'])**self.vals['r2'])
        b1 = (self.vals['smax'] + (self.vals['smin'] - self.vals['smax']) / (
            1. + (self.temp / self.vals['s1'])**self.vals['s2'])) * 2
        b3 = (self.vals['wmax'] + (self.vals['wmin'] - self.vals['wmax']) / (1. + (self.temp / self.vals['w1'])**self.vals['w2'])) * 2

        # print bmin
        B = self.vals['blow'] * (bmin + (self.vals['bmax'] - bmin) / (
            1. + ((2. * min_car_den + doping) / b1
                 )**self.vals['b2'] + ((2. * min_car_den + doping) / b3)**self.vals['b4']))

        return self.Roosbroeck(min_car_den, doping, B)


class Auger(HelperFunctions):
    model_file = 'auger.model'

    def __init__(self, matterial, model_author=None, temp=300, ni=9.65e9):

        self.Models = ConfigParser.ConfigParser()
        self.matterial = matterial

        constants_file = os.path.join(
            os.path.dirname(os.path.realpath(__file__)),
            matterial,
            self.model_file)

        self.Models.read(constants_file)

        self.change_model(model_author)
        self.temp = temp
        self.ni = ni

    # def change_model(self, model_author=None):

    #     if model_author is None:
    #         self.model_author = self.Models.get('default', 'constants')
    #     else:
    #         # Need a check to make sure craNhcan't be passed
    #         self.model_author = model_author
    #     # print self.constants
    #     self.model = self.Models.get(self.model_author, 'model')
    #     self.model_values = self.Models._sections[self.model_author]

    def tau(self, min_car_den, Na, Nd):

        self.Nh_0, self.Ne_0 = self.check_doping(Na, Nd)
        doping = np.max(Na, Nd)
        return getattr(self, self.model)(min_car_den, doping)

    def itau_aug(self, min_car_den, Na, Nd):
        return 1 / self.tau(min_car_den, Na, Nd)

    def auger_dopants(self, min_car_den, doping=None):
        '''This is the classic auger model that only depends on doping

        The model requires the inputs
        C1, C2, and C3.

        It also requires the dark carrier concentrations
        I'm not sure where C1 is being used... need to check this
        '''
        Ce = self.Models.getfloat(self.model_author, 'cn')
        Ch = self.Models.getfloat(self.model_author, 'cp')

        Nh = self.Nh_0 + min_car_den
        Ne = self.Ne_0 + min_car_den

        R = Ce * Ne**2 * Nh + Ch * Ne * Nh**2

        return min_car_den / R

    def auger(self, min_car_den, doping=None):
        '''This is the classic auger model that has the impact of many carriers

        The model requires the inputs
        C1, C2, and C3.

        It also requires the  number of dark carrier concentrations
        I'm not sure where C1 is being used... need to check this
        '''

        Nh = self.Nh_0 + min_car_den
        Ne = self.Ne_0 + min_car_den

        Ce = self.vals['ced'] * self.Ne_0 / \
            (self.Ne_0 + Nh) + self.vals['ccc'] / 2 * Nh / (Nh + self.Ne_0)
        Ch = self.vals['chd'] * self.Nh_0 / \
            (self.Nh_0 + Ne) + self.vals['ccc'] / 2 * Ne / (Ne + self.Nh_0)

        R = (Ce * Ne + Ch * Nh) * (Ne * Nh - self.ni**2)

        return min_car_den / R

    def auger_part(self, *args):
        """ A dummy class for models with incomplete parameters"""
        pass

    def coulomb_enhanced_auger(self, min_car_den, doping):
        '''The coulomb enhanced auger model, as proposed by Kerr2002

        The model requires the inputs
        K_n, K_p, K_Delta, delta, L_eeh, N_eeh, K_eeh, P_ehh, L_ehh, K_ehh

        It also requires the doping and the dark carrier concentrations

        There may be an older form from schmit. 
        '''
        # then need to get doping and delta n
        # p, n, Nh_0, Ne_0 = self.n_and_p(self.min_car_den)

        Nh_0, Ne_0 = doping, self.ni**2 / doping
        Nh = Nh_0 + min_car_den
        Ne = Ne_0 + min_car_den
        # ni = 1e10

        g_eeh = 1 + self.vals['l_eeh'] * (1 - np.tanh(
            (Ne_0 / self.vals['k_eeh'])**(self.vals['n_eeh'])))
        g_ehh = 1 + self.vals['l_ehh'] * (1 - np.tanh(
            (Nh_0 / self.vals['k_ehh'])**(self.vals['p_ehh'])))

        # Then the lifetime is provided by this
        return min_car_den / (Ne * Nh - self.ni**2) /\
            (self.vals['k_n'] * g_eeh * Ne_0 +
             self.vals['k_p'] * g_ehh *
             Nh_0 +
             self.vals['k_delta'] *
             min_car_den**self.vals['delta']
             )

if __name__ == "__main__":
    # a = IntrinsicCarrierDensity()
    # a._PlotAll()
    # plt.show()
    a = Intrinsic('Si')
    deltan = np.logspace(12, 17)
    # print a.Radiative.ni, a.Auger.ni

    print plt.plot(deltan, a.intrisic_carrier_lifetime(deltan, 1e16, 0) * 1e6)
    plt.loglog()
    plt.show()
