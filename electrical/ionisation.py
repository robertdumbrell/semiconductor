#!/usr/local/bin/python
# UTF-8

import numpy as np
import matplotlib.pylab as plt
import os
import ConfigParser

from semiconductor.helper.helper import HelperFunctions
import impurity_ionisation_models as IIm
import semiconductor.matterial.densityofstates as DOS
import semiconductor.general_functions.carrierfunctions as CF


class Ionisation(HelperFunctions):

    '''
    Depending on a dopant level from a band, and the thermal 
    energy available, a dopant is electrical active (donating or
    accepting an electron to the band)  or inactive. 
    '''

    author_list = 'ionisation.models'

    def __init__(self, matterial='Si', author=None, temp=300):
        self.Models = ConfigParser.ConfigParser()
        self.matterial = matterial

        constants_file = os.path.join(
            os.path.dirname(os.path.realpath(__file__)),
            matterial,
            self.author_list)

        self.Models.read(constants_file)
        self.vals, self.model = self.change_model(author)
        self.temp = temp

    def update(self, N_imp, ne, nh, impurity, temp=None, author=None):
        '''
        Calculates the number of ionisied impurities

        Inputs:
            N_imp: (float or numpy array)
                number of the impurity in the material
            ne: (float or numpy array)
                number of electrons
            nh: (float or numpy array)
               number of holes
            impurity: (str)
                the element name of the impurity
            temp: (optional)
                the  temperature in Kelvin to be evaluated
            author: (optional str)
                the author of the impurity model to use

        output:
            the number of ionised impurities
        '''
        if temp is None:
            temp = self.temp

        # a check to make sure the model hasn't changed
        if author is not None:
            self.vals, self.model = self.change_model(author)
        # this should be change an outside function alter

        # checks if and get the required density of states model
        if 'dos_author' in self.vals.keys():
            Nc, Nv = DOS.DOS('Si').update(
                temp=temp, author=self.vals['dos_author'])
        else:
            Nc, Nv = 0, 0

        # 
        if impurity in self.vals.keys():
            iN_imp = getattr(IIm, self.model)(
                self.vals, N_imp, ne, nh, temp, Nc, Nv, self.vals[impurity])
            iN_imp *= N_imp
        else:
            print 'No such impurity, please check your model and spelling'
            print 'Returning zero array'
            iN_imp = np.zeros(np.asarray(N_imp).flatten().shape[0])

        return iN_imp

    def update_dopant_ionisation(self, N_dop, dn, impurity,
                                 temp=None, author=None):
        '''
        This is a special function used to determine the number of ionised dopants
        given a number of excess carriers, and a single dopant type.    
        '''

        if temp is None:
            temp = self.temp

        if author is not None:
            self.vals, self.model = self.change_model(author)

        iN_dop = N_dop

        if impurity in self.vals.keys():
            # TO DO, change this from just running 10 times to a proper check
            for i in range(10):
                if self.vals['tpe_' + self.vals[impurity]] == 'donor':
                    Nd = iN_dop
                    Na = 0
                elif self.vals['tpe_' + self.vals[impurity]] == 'acceptor':
                    Na = iN_dop
                    Nd = 0

                ne, nh = CF.get_carriers(Na, Nd, dn, temp=temp, matterial='Si')

                iN_dop = self.update(
                    N_dop, ne, nh, impurity, temp=temp, author=None)
        else:
            print 'Not a valid impurity'
        return iN_dop

    def check_models(self):
        '''
        Plots a check of the modeled data against Digitised data from either
        papers or from other implementations of the model.
        '''
        plt.figure('Ionised impurities')

        iN_imp = N_imp = np.logspace(15, 20)

        dn = 1e10
        temp = 300

        for impurity in ['phosphorous', 'arsenic']:

            iN_imp = self.update_dopant_ionisation(N_imp,
                                                   dn,
                                                   impurity,
                                                   temp, author=None)

            if not np.all(iN_imp == 0):
                plt.plot(
                    N_imp, iN_imp / N_imp * 100, label='Altermatt: ' + impurity)

        test_file = os.path.join(
            os.path.dirname(os.path.realpath(__file__)),
            'Si', 'check data', 'donors.csv')

        data = np.genfromtxt(test_file, delimiter=',', skip_header=1)

        for i in range(0, (data.shape[1] + 2) / 2, 2):
            print i
            plt.plot(
                data[:, i], data[:, i + 1] * 100, 'r.',
                label='Digitised data')

        plt.semilogx()
        plt.xlabel('Impurity (cm$^{-3}$)')
        plt.ylabel('Fraction of Ionised impurities (%)')

        plt.legend(loc=0)
