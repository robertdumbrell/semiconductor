
import numpy as np
import matplotlib.pylab as plt
import sys
import os
import ConfigParser


sys.path.append(
    os.path.abspath(os.path.join(os.getcwd(), os.pardir, os.pardir)))
import semiconductor.defults


class Helper():

    def tau():
        pass

    def _PlotAll(self):
        fig, ax = plt.subplots(1)
        # ax = plt.add_subplot(111)
        for model in self.AvailableModels():
        # ax.plot(np.inf,np.inf,'k-',label = 'Auger')
            self.change_model(model)
            tau = self.tau()
            if tau is not None:
                ax.plot(self.Deltan, tau * 1e6, label=model)

        ax.legend(loc=0)
        ax.loglog()
        ax.set_xlabel('$\Delta$ n (cm$^{-3}$)')
        ax.set_ylabel('Lifetime (us)')

        # Helper routiens

    def AvailableModels(self):
        a = self.Models.sections()
        a.remove('default')
        return a


class Radiative(Helper):

    Deltan = np.logspace(12, 19)
    Doping = 1e16
    Ne_0 = 1e3
    Nh_0 = 1e16
    ni = 9.65e9
    T = 300

    def __init__(self, matterial='Si', model_author=None):
        self.Models = ConfigParser.ConfigParser()
        self.matterial = matterial

        constants_file = os.path.join(
            os.path.dirname(os.path.realpath(__file__)),
            matterial,
            'radiative.model')

        self.Models.read(constants_file)

        self.change_model(model_author)

    def change_model(self, model_author=None):

        if model_author is None:
            self.model_author = self.Models.get('default', 'constants')
        else:
            # Need a check to make sure craNhcan't be passed
            self.model_author = model_author
        # print self.constants
        self.model = self.Models.get(self.model_author, 'model')
        self.model_values = self.Models._sections[self.model_author]

    def tau(self):
        return getattr(self, self.model)()

    def itau(self):
        return getattr(self, self.model)()

    def Roosbroeck(self, B=None):
        if B is None:
            B = self.Models.getfloat(self.model_author, 'B')

        Nh = self.Nh_0 + self.Deltan
        Ne = self.Ne_0 + self.Deltan

        R = B * (Ne * Nh - self.Ne_0 * self.Nh_0)
        return self.Deltan / R

    def Roosbroeck_with_screening(self):
        """ 
        This is the roosbroeck model that accounts for many things
        It needs temperature, deltan, doping and blow to be defined

        """

        # self.Models.items(self.model_author)


        Blow = self.Models.getfloat(self.model_author, 'blow')
        bmax = self.Models.getfloat(self.model_author, 'bmax')
        rmax = self.Models.getfloat(self.model_author, 'rmax')
        smax = self.Models.getfloat(self.model_author, 'smax')
        wmax = self.Models.getfloat(self.model_author, 'wmax')
        rmin = self.Models.getfloat(self.model_author, 'rmin')
        smin = self.Models.getfloat(self.model_author, 'smin')
        wmin = self.Models.getfloat(self.model_author, 'wmin')
        b2 = self.Models.getfloat(self.model_author, 'b2')
        r1 = self.Models.getfloat(self.model_author, 'r1')
        s1 = self.Models.getfloat(self.model_author, 's1')
        w1 = self.Models.getfloat(self.model_author, 'w1')
        b4 = self.Models.getfloat(self.model_author, 'b4')
        r2 = self.Models.getfloat(self.model_author, 'r2')
        s2 = self.Models.getfloat(self.model_author, 's2')
        w2 = self.Models.getfloat(self.model_author, 'w2')

        bmin = rmax + (rmin - rmax) / (1. + (self.T / r1)**r2)
        b1 = (smax + (smin - smax) / (1. + (self.T / s1)**s2)) * 2
        b3 = (wmax + (wmin - wmax) / (1. + (self.T / w1)**w2)) * 2

        # print bmin
        B = Blow * (bmin + (bmax - bmin) / (
            1 + ((2 * self.Deltan + self.Doping) / b1
                 )**b2 + ((2 * self.Deltan + self.Doping) / b3)**b4))

        return self.Roosbroeck(B)


class Auger(Helper):
    # needs to be dealt with
    Deltan = np.logspace(12, 19)
    Doping = 1e16
    Ne_0 = 1e3
    Nh_0 = 1e16
    ni = 9.65e9

    def __init__(self, matterial, model_author=None):
        self.Models = ConfigParser.ConfigParser()
        self.matterial = matterial

        constants_file = os.path.join(
            os.path.dirname(os.path.realpath(__file__)),
            matterial,
            'auger.model')

        self.Models.read(constants_file)

        self.change_model(model_author)

    def change_model(self, model_author=None):

        if model_author is None:
            self.model_author = self.Models.get('default', 'constants')
        else:
            # Need a check to make sure craNhcan't be passed
            self.model_author = model_author
        # print self.constants
        self.model = self.Models.get(self.model_author, 'model')
        self.model_values = self.Models._sections[self.model_author]

    def tau(self):
        return getattr(self, self.model)()

    def itau_aug(self):
        return 1 / self.tau()

    def auger_dopants(self):
        '''This is the classic auger model that only depends on doping

        The model requires the inputs
        C1, C2, and C3.

        It also requires the doping and the dark carrier concentrations
        I'm not sure where C1 is being used... need to check this
        '''
        Ce = self.Models.getfloat(self.model_author, 'cn')
        Ch = self.Models.getfloat(self.model_author, 'cp')

        Nh = self.Nh_0 + self.Deltan
        Ne = self.Ne_0 + self.Deltan

        R = Ce * Ne**2 * Nh + Ch * Ne * Nh**2

        return self.Deltan / R

    def auger(self):
        '''This is the classic auger model that has the impact of many carriers

        The model requires the inputs
        C1, C2, and C3.

        It also requires the doping and the dark carrier concentrations
        I'm not sure where C1 is being used... need to check this
        '''
        Ccc = self.Models.getfloat(self.model_author, 'ccc')
        Chd = self.Models.getfloat(self.model_author, 'chd')
        Ced = self.Models.getfloat(self.model_author, 'ced')

        Nh = self.Nh_0 + self.Deltan
        Ne = self.Ne_0 + self.Deltan

        Ce = Ced * self.Ne_0 / \
            (self.Ne_0 + Nh) + Ccc / 2 * Nh / (Nh + self.Ne_0)
        Ch = Chd * self.Nh_0 / \
            (self.Nh_0 + Ne) + Ccc / 2 * Ne / (Ne + self.Nh_0)

        R = (Ce * Ne + Ch * Nh) * (Ne * Nh - self.ni**2)

        return self.Deltan / R

    def auger_part(self):
        """ A dummy class for models with incomplete parameters"""
        pass

    def coulomb_enhanced_auger(self):
        '''The coulomb enhanced auger model, as proposed by Kerr2002

        The model requires the inputs
        K_n, K_p, K_Delta, delta, L_eeh, N_eeh, K_eeh, P_ehh, L_ehh, K_ehh

        It also requires the doping and the dark carrier concentrations

        There may be an older form from schmit. 
        '''
        K_n = self.Models.getfloat(self.model_author, 'k_n')
        K_Nh = self.Models.getfloat(self.model_author, 'k_p')
        K_Delta = self.Models.getfloat(self.model_author, 'k_delta')
        Delta = self.Models.getfloat(self.model_author, 'delta')
        L_eeh = self.Models.getfloat(self.model_author, 'l_eeh')
        N_eeh = self.Models.getfloat(self.model_author, 'n_eeh')
        K_eeh = self.Models.getfloat(self.model_author, 'k_eeh')
        P_ehh = self.Models.getfloat(self.model_author, 'p_ehh')
        L_ehh = self.Models.getfloat(self.model_author, 'l_ehh')
        K_ehh = self.Models.getfloat(self.model_author, 'k_ehh')
        # then need to get doping and delta n
        # p, n, Nh_0, Ne_0 = self.n_and_p(self.Deltan)
        print 'Have not impimented doping and deltan yet or ni'
        # Deltan = 1e14
        Nh_0, Ne_0 = self.Doping, self.ni**2 / self.Doping
        Nh = Nh_0 + self.Deltan
        Ne = Ne_0 + self.Deltan
        # ni = 1e10

        g_eeh = 1 + L_eeh * (1 - np.tanh(
            (Ne_0 / K_eeh)**(N_eeh)))
        g_ehh = 1 + L_ehh * (1 - np.tanh(
            (Nh_0 / K_ehh)**(P_ehh)))

        # Then the lifetime is provided by this
        return self.Deltan / (Ne * Nh - self.ni**2) /\
            (K_n * g_eeh * Ne_0 +
             K_Nh * g_ehh *
             Nh_0 +
             K_Delta *
             self.Deltan**Delta
             )


a = Radiative('Si')
print  dict(a.Models.items(a.model_author))['model']
# print a._PlotAll()
# print a.tau()
plt.show()
