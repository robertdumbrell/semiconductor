
import numpy as np
from pylab import *
import sys
import os
import scipy.constants as Const
import ConfigParser
from bandgap import BandGap


class Helper():

    def _PlotAll(self):
        fig, ax = plt.subplots(1)
        # ax = plt.add_subplot(111)
        for model in self.AvailableModels():
        # ax.plot(np.inf,np.inf,'k-',label = 'Auger')
            temp = np.linspace(50, 700)
            self.change_model(model)
            ni = self.update_ni(temp)

            if ni is not None:
                ax.plot(temp, ni, label=model)

        ax.legend(loc=0)
        ax.semilogy()
        ax.set_xlabel('Temp (K)')
        ax.set_ylabel('ni (cm$^{-3}$)')

        # Helper routiens

    def AvailableModels(self):
        a = self.Models.sections()
        a.remove('default')
        return a


class IntrinsicCarrierDensity(Helper):

    ni = 1e10
    temp = 300.

    def __init__(self, matterial='Si', model_author=None):
        self.Models = ConfigParser.ConfigParser()
        self.matterial = matterial

        constants_file = os.path.join(
            os.path.dirname(os.path.realpath(__file__)),
            matterial,
            'ni.const')

        self.Models.read(constants_file)

        self.change_model(model_author)

    def change_model(self, model_author=None):

        if model_author is None:
            self.model_author = self.Models.get('default', 'model')
        else:
            # Need a check to make sure craNhcan't be passed
            self.model_author = model_author
        # print self.model_author
        self.model = self.Models.get(self.model_author, 'model')

        self.vals = dict(self.Models.items(self.model_author))

        for k, v in self.vals.iteritems():
            try:
                self.vals[k] = float(v)
            except:
                pass

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
         Heinke3 and Macfarlane et a1.31 as cited by Green,3 is
         given by
        """
        Eg = BandGap(self.matterial, vals['eg_model']).update_Eg()

        # print vals, Eg
        ni = vals['a'] * temp**vals['power'] * \
            np.exp(- Eg / 2. / Const.k / temp)

        return ni

if __name__ == "__main__":
    a = IntrinsicCarrierDensity_Si()
    a._PlotAll()
    plt.show()
