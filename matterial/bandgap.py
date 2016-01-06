import numpy as np
import matplotlib.pylab as plt
import os
import ConfigParser

from semiconductor.helper.helper import HelperFunctions
import semiconductor.matterial.intrinsicbandgapmodels as iBg
import semiconductor.matterial.BandgapNarrowingModels as Bgn


class BandGap():

    '''
    A simple class to combine the intrinsic band gap and
    band gap narrowing classes for easy access
    '''

    def __init__(self,
                 matterial='Si', Egi_model=None, BNG_model=None, dopant=None):

        self.Egi = IntrinsicBandGap(matterial, model=Egi_model)
        self.BGN = BandGapNarrowing(matterial, model=BNG_model)
        self.dopant = dopant

    def plot_all_models(self):
        print 'See IntrinsicBandGap class and BGN class for models'

    def caculate_Eg(self, temp, doping, min_car_den=None, dopant=None):
        '''
        Calculates the band gap
        '''

        if dopant is None:
            dopant = self.dopant

        # just prints a warning if the model is for the incorrect dopants
        dopant_model_list = self.BGN.available_models('dopant', dopant)
        if self.BGN.model not in dopant_model_list:
            print 'You have the incorrect model for your dopant'

        # print 'The band gaps are:', self.Egi.update_Eg(temp),
        # self.BGN.update_BGN(doping, min_car_den)
        Eg = self.Egi.update_iEg(
            temp) - self.BGN.update_BGN(doping, min_car_den)
        return Eg


class IntrinsicBandGap(HelperFunctions):

    '''
    The intrinsic band-gap as a function of temperature
        it changes as a result of:
             different effective carrier mass (band strucutre)
    '''

    model_file = 'bandgap.models'

    def __init__(self, matterial='Si', model=None):
        self.Models = ConfigParser.ConfigParser()
        self.matterial = matterial

        constants_file = os.path.join(
            os.path.dirname(os.path.realpath(__file__)),
            matterial,
            self.model_file)

        self.Models.read(constants_file)

        self.change_model(model)

    def update_iEg(self, temp=None, model=None, multiplier=1.01):
        '''
        a function to update the intrinsic BandGap

        inputs:
            temperature in kelvin
            model: (optional)
                  the model used.
                  If not provided the last provided model is used
                  If no model has been provided Passler model is used
            multiplier: A band gap multipler. 1.01 is suggested.

        output:
            the intrinsic bandgap in eV
        '''

        if temp is None:
            temp = self.temp
        if model is not None:
            self.change_model(model)

        Eg = getattr(iBg, self.model)(self.vals, temp)

        return Eg * multiplier

    def check_models(self):
        '''
        Displays a plot of the models against that taken from a
        respected website (https://www.pvlighthouse.com.au/)
        '''
        plt.figure('Intrinsic bandgap')
        t = np.linspace(1, 500)

        for model in self.available_models():

            Eg = self.update_iEg(t, model=model, multiplier=1.0)
            plt.plot(t, Eg, label=model)

        test_file = os.path.join(
            os.path.dirname(os.path.realpath(__file__)),
            'Si', 'check data', 'iBg.csv')

        data = np.genfromtxt(test_file, delimiter=',', names=True)

        for temp,name in zip(data.dtype.names[0::2],data.dtype.names[1::2]):
            plt.plot(
                data[temp], data[name], '--', label=name)

        plt.xlabel('Temperature (K)')
        plt.ylabel('Intrinsic Bandgap (eV)')

        plt.legend(loc=0)
        self.update_iEg(0, model=model, multiplier=1.01)


class BandGapNarrowing(HelperFunctions):

    '''
    Bang gap narrowing accounts for a reduction in bandgap that 
    occurs as a result from no thermal effects. These include:
        doping. It is dopant dependent
        excess carrier density (non thermal distribution)
    Note: I currently believed that the impact of dopants 
        is much larger than the impact of the carrier distribution
    '''

    model_file = 'bandgap_narrowing.models'

    def __init__(self, matterial='Si', model=None):
        self.Models = ConfigParser.ConfigParser()
        self.matterial = matterial

        constants_file = os.path.join(
            os.path.dirname(os.path.realpath(__file__)),
            matterial,
            self.model_file)

        self.Models.read(constants_file)
        self.change_model(model)

    def update_BGN(self, doping, min_car_den=None, model=None):

        if model is not None:
            self.change_model(model)

        return getattr(Bgn, self.model)(self.vals, doping, min_car_den)

    def check_models(self):
        plt.figure('Bandgap narrowing')
        doping = np.logspace(12, 20)

        for model in self.available_models():
            BGN = self.update_BGN(doping, model=model)
            plt.plot(doping, BGN, label=model)

        test_file = os.path.join(
            os.path.dirname(os.path.realpath(__file__)),
            'Si', 'check data', 'Bgn.csv')

        data = np.genfromtxt(test_file, delimiter=',', names=True)
        print data.dtype.names
        for name in data.dtype.names[1:]:
            plt.plot(
                data['N'], data[name], 'r--', label='PV-lighthouse\'s: ' + name)

        plt.semilogx()
        plt.xlabel('Doping (cm$^{-3}$)')
        plt.ylabel('Bandgap narrowing (K)')

        plt.legend(loc=0)


if __name__ == "__main__":
    # a = IntrinsicCarrierDensity()
    # a._PlotAll()
    # plt.show()

    a = BandGap('Si')

    IntrinsicBandGap().check_models()
    BandGapNarrowing().check_models()
    plt.show()
    # a.print_model_notes()
    # deltan = np.logspace(12, 17)
    # print a.Radiative.ni, a.Auger.ni
    # doping = np.logspace(12, 20, 100)
    # a.plot_all_models('update_BGN', xvalues=doping, doping=doping)

    # a.caculate_Eg(300., 1e16, 1e16, 'boron')
    # print plt.plot(doping, a.update_BGN(doping))
    # plt.semilogx()
    # plt.show()
