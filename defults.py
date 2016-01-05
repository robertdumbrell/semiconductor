

import scipy.constants as const
import numpy as np
import semiconductor

class Constants(object):

    # def __init__(self):
    # constantly constant contants, lol
    kb = 1.3806488e-23  # J/K
    q = 1.6e-19  # C
    # not so constant constants
    T = 300  # K
    ni = 9.65e9  # /cm3 from memory need reference, make a function of T?
    Doping = 5e15  # /cm3
    Deltan = np.logspace(11, 19, 1000)
    MajorityCarrier = 'p'
    Width = 0.018
    Eg = 1.12  # eV
    # cm/s Thermal velocity from W. M. Bullis and H. R. Huff, J. Electrochem.
    # Soc 143 1399 1996
    # vth = 1.1e7

    def Vth(self):
        return const.k * self.T / const.e

    def n_and_p(self, Deltan):
        if (self.MajorityCarrier == 'n'):
            n = self.Doping + Deltan
            p = Deltan
            n_0 = self.Doping
            p_0 = self.ni**2 / self.Doping
        elif(self.MajorityCarrier == 'p'):
            p = self.Doping + Deltan
            n = Deltan
            p_0 = self.Doping
            n_0 = self.ni**2 / self.Doping
        return p, n, p_0, n_0


if __name__ == "__main__":
    a = semiconductor.matterial.ni.IntrinsicCarrierDensity()
    a.check_models()
    # temp = np.linspace(0, 600)
    # a.plot_all_models('update_ni', temp=temp)

    plt.show()