
import numpy as np
import sys
import os
import matplotlib.pylab as plt
import ConfigParser

sys.path.append(
    os.path.abspath(os.path.join(os.getcwd(), os.pardir, os.pardir)))
import semiconductor.defults


class SRH(semiconductor.defults.Constants):
    vth = 2.05E+7  # taken from PV light house

    def __init__(self, matterial, defect=None, N=None):

        self.Defects = ConfigParser.ConfigParser()
        self.matterial = matterial

        constants_file = os.path.join(
            os.path.dirname(os.path.realpath(__file__)),
            matterial,
            'SRH.defects')

        self.Defects.read(constants_file)

        self.select_defect(defect, N)

    def select_defect(self, defect=None, N=None):

        if defect is None:
            self.defect = 'none'
        else:
            # Need a check to make sure craNhcan't be passed
            self.defect = defect

        self.Cal_tau(N)

    def Cal_tau(self, N=None):
        if N is None:
            N = 1e10

        # print N,self.sigma_n,self.vth, self.Et
        # print self.defect
        sigma_n = self.Defects.getfloat(self.defect, 'sigma_n')
        sigma_p = self.Defects.getfloat(self.defect, 'sigma_p')
        self.Et = self.Defects.getfloat(self.defect, 'et')

        self.tau_n = 1. / N / sigma_n / self.vth  # s
        self.tau_p = 1. / N / sigma_p / self.vth  # s

    def tau(self, Deltan=[0, 0], tau_n=0, tau_p=0, Et=None):
        if all(Deltan) != 0:
            self.Deltan = Deltan
            # print 'inputed Deltan values used'
        if(tau_n != 0):
            self.tau_n = tau_n
            # print 'inputed tau_n values used',
        if(tau_p != 0):
            self.tau_p = tau_p

        if Et is None:
            Et = self.Et
            # print 'inputed tau_p values used'
        # print self.tau_n,self.tau_p,self.Et

        return self.SchroderTextbook_SRH(self.tau_n, self.tau_p, Et)

    def SchroderTextbook_SRH(self, tau_n, tau_p, Et=0):
        """
        this is been confirmed to agree with SZE as well as PVlighthouse
        It is assumed that Nt << p, n, and that the semiconductor is not
        degenerate.
        """

        p1 = self.ni * np.exp(-(Et) / self.Vth())
        n1 = self.ni * np.exp((Et) / self.Vth())

        p, n, p_0, n_0 = self.n_and_p(self.Deltan)

        return (tau_p * (n + n1) + tau_n * (p + p1)) / (p)

    def itau(self, Deltan=0, tau_n=0, tau_p=0, Et=0):
        # print tau_n,tau_p,Et
        return 1. / self.tau(Deltan, tau_n, tau_p, Et)

    def AvailableModels(self):
        a = self.Defects.sections()
        print self.Defects.defaults()
        a.remove('none')
        return a

    def _PlotAll(self):
        fig, ax = plt.subplots(1, 2, figsize=(16, 6))
        # ax = plt.add_subplot(111)
        counter = 0
        defects = self.AvailableModels()
        deltan = np.logspace(12, 17)
        for defect in defects:
        # ax.plot(np.inf,np.inf,'k-',label = 'Auger')
            self.select_defect(defect)
            self.Cal_tau(1e10)
            ax[0].plot(counter, self.Et, 'o', label=defect)
            ax[1].plot(deltan, self.tau(deltan) * 1e6, label=defect)
            counter += 1
        # ax.legend(loc=0)
        # ax.loglog()
        x = np.arange(counter)
        print defects
        ax[0].set_xticks(x)
        ax[0].set_xticklabels(defects, rotation=90)
        ax[0].set_xlabel('Defect')
        ax[0].set_ylabel('E$_t$ from Ei (eV)')
        ax[1].loglog()
        ax[1].set_xlabel('$\Delta$ n (cm$^{-3}$)')
        ax[1].set_ylabel('Lifetime (us)')


# class SHR_defects(semiconductor.defults.Constants):
#     Ev = 0

#     def __init__(self, SRH_Class):
#         """used to pass the SRH class"""
#         self.SHR = SRH_Class
#         self.Ec = self.Eg

#     def Fe(self, N):
# Taken from Iron detection in crystalline silicon by carrier lifetime
# measurements for arbitrary injection and doping 2004
#         Et = self.Ev + .38
#         self.SHR.Et = Et - self.Eg / 2

# self.SHR.sigma_n = 5e-14  # s
# self.SHR.sigma_p = 7e-17  # s

#     def FeB(self, N):
# Taken from Iron detection in crystalline silicon by carrier lifetime
# measurements for arbitrary injection and doping 2004

#         Et = self.Ec - 0.23
# self.SHR.Et = Et - self.Eg / 2  # eV  defined as Et = E_trap - E_i

# self.SHR.sigma_n = 3e-14  # s
# self.SHR.sigma_p = 2e-15  # s

# deltan = np.logspace(12, 17, 100)
# a = SRH('Si', 'Fei_d')
# a._PlotAll()
# tau = a.tau(deltan)
# plt.plot(deltan, tau)
# print tau
# a = SRH('Si', 'FeB')
# tau = a.tau(deltan)
# plt.plot(deltan, tau)
# plt.loglog()
# plt.sh/ow()
