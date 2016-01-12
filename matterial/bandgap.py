

from bandgap_intrinsic import IntrinsicBandGap
from bandgap_narrowing import BandGapNarrowing


class BandGap():

    '''
    A simple class to combine the intrinsic band gap and
    band gap narrowing classes for easy access
    '''

    def __init__(self,
                 matterial='Si', iEg_model=None, BNG_model=None, dopant=None):

        self.iEg = IntrinsicBandGap(matterial, model=iEg_model)
        self.BGN = BandGapNarrowing(matterial, model=BNG_model)
        self.dopant = dopant

    def plot_all_models(self):
        self.iEg.plot_all_models()
        self.iEg.plot_all_models()

    def caculate_Eg(self, temp, doping, min_car_den=None, dopant=None):
        '''
        Calculates the band gap
        '''

        if dopant is None:
            dopant = self.dopant

        # just prints a warning if the model is for the incorrect dopants
        dopant_model_list = self.Egi.available_models('dopant', dopant)
        if self.BGN.model not in dopant_model_list:
            print 'You have the incorrect model for your dopant'

        # print 'The band gaps are:', self.iEg.update_Eg(temp),
        # self.BGN.update_BGN(doping, min_car_den)
        Eg = self.iEg.update_iEg(
            temp) - self.BGN.update_BGN(doping, min_car_den)
        return Eg

    def check_models(self):
        self.iEg.check_models()
        self.BGN.check_models()
