#!/usr/local/bin/python
# UTF-8

import numpy as np
import matplotlib.pylab as plt
import os
import ConfigParser

from semiconductor.helper.helper import HelperFunctions
import dopant_ionisation_models 
# import semiconductor.matterial.bandgap_narrowing_models as Bgn
# import semiconductor.general_functions.carrierfunctions as GF


class DopantIonisation(HelperFunctions):

    '''
    Depending on a dopant level from a band, and the thermal 
    energy available, a dopant is electrical active (donating or
    accepting an electron to the band)  or inactive. 
    '''

    author_list = 'ionisation.models'

    def __init__(self, matterial='Si', author=None):
        self.Models = ConfigParser.ConfigParser()
        self.matterial = matterial

        constants_file = os.path.join(
            os.path.dirname(os.path.realpath(__file__)),
            matterial,
            self.author_list)

        self.Models.read(constants_file)
        self.change_model(author)

    def update_DI(self, Na, Nd, temp, author=None,  ni_author=None):
        '''
        Calculates the dopant ionisation fraction

        Inputs:
            ???

        output:
            ???
        '''

        # this should be change an outside function alter
        pass 
        
    def check_models(self):
        plt.figure('Bandgap narrowing')
        Na = np.logspace(12, 20)
        Nd = 0
        dn = 1e14
        temp = 300

        for author in self.available_models():
            BGN = self.update_BGN(Na=Na, Nd=Nd, min_car_den=dn,
                                  author=author,
                                  temp=temp)

            if not np.all(BGN == 0):
                plt.plot(Na, BGN, label=author)

        test_file = os.path.join(
            os.path.dirname(os.path.realpath(__file__)),
            'Si', 'check data', 'Bgn.csv')

        data = np.genfromtxt(test_file, delimiter=',', names=True)

        for name in data.dtype.names[1:]:
            plt.plot(
                data['N'], data[name], 'r--',
                label='PV-lighthouse\'s: ' + name)

        plt.semilogx()
        plt.xlabel('Doping (cm$^{-3}$)')
        plt.ylabel('Bandgap narrowing (K)')

        plt.legend(loc=0)
