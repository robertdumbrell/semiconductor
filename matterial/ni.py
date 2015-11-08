
import numpy as np
from pylab import *
import sys
import os
import scipy.constants as Const
import ConfigParser
from bandgap import IntrinsicBandGap as BandGap

# TODO:
# Need to make bandgap class, that includes that impact of BNG

sys.path.append(
    os.path.abspath(os.path.join(os.getcwd(), os.pardir, os.pardir)))
from semiconductor.helper.helper import HelperFunctions


class IntrinsicCarrierDensity(HelperFunctions):
    '''
    The intrinisc carrier density is the number of carriers 
    that exist the a matterial at thermal equlibrium.
    It is impacted by the band gap (and bandgap narrowing)
    '''

    # ni = 1e10
    temp = 300.
    model_file = 'ni.models'
    def __init__(self, matterial='Si', model_author=None, temp = 300.):
        self.Models = ConfigParser.ConfigParser()
        self.matterial = matterial

        constants_file = os.path.join(
            os.path.dirname(os.path.realpath(__file__)),
            matterial,
            self.model_file)

        self.Models.read(constants_file)

        self.change_model(model_author)
        self.temp = temp


    def update_ni(self, temp=None):
        if temp is None:
            temp = self.temp

        self.ni = getattr(self, self.model)(self.vals, temp)
        # print self.ni
        return self.ni

    def ni_temp(self, vals, temp):
        """
         This form comes from Bludau, Onton, and
         Heinke3 and Macfarlane et a1.31 as cited by Green,3 is
         given by
        """
        ni = vals['a'] * (temp)**vals['power'] * \
            np.exp(- vals['eg'] / temp)
        # print self.ni, - vals['eg'], 2., Const.k, self.temp

        return ni

    def ni_temp_eg(self, vals, temp):
        """
         This form comes from Bludau, Onton, and
         Heinke3 and Macfarlane et a1.31 as cited by Green,3 
        """
        Eg = BandGap(self.matterial, vals['eg_model']).update_Eg(temp=temp)

        # print vals, Eg
        ni = vals['a'] * temp**vals['power'] * \
            np.exp(- Eg / 2. / Const.k / temp)

        return ni

if __name__ == "__main__":
    a = IntrinsicCarrierDensity()
    temp = np.linspace(0,600)
    a.plot_all_models('update_ni', temp=temp)
    plt.semilogy()
    plt.show()
