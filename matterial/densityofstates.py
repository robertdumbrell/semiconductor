#!/usr/local/bin/python
# UTF-8

import numpy as np
import matplotlib.pylab as plt
import sys
import os
import ConfigParser

from semiconductor.helper.helper import HelperFunctions
import densityofstates_models as dos_models
from bandgap_intrinsic import IntrinsicBandGap as Egi


class DOS(HelperFunctions):

    '''
    The density of states is a value that determines the 
    number of free states for electrons and holes in the conduction
    and valance band
    '''

    author_list = 'DOS.models'

    def __init__(self, matterial='Si', author=None, temp=300.):
        self.Models = ConfigParser.ConfigParser()
        self.matterial = matterial

        constants_file = os.path.join(
            os.path.dirname(os.path.realpath(__file__)),
            matterial,
            self.author_list)

        self.Models.read(constants_file)

        self.change_model(author)
        self.temp = temp

    def update(self, temp=None, author=None):
        '''
        a function to update the density of states

        inputs:
            temperature: (optional)
                         in kelvin
            author:  (optional)
                    the author used.
                    If not provided the last provided author is used
                    If no author has been provided,  Couderc's model is used
        output:
            the density of states of the valance band
            the density of states of the conduction band
        '''
        # able to input a temperature to change
        if temp is None:
            temp = self.temp

        # a check to make sure the model hasn't changed
        if author is not None:
            self.change_model(author)

        if 'egi_author' in self.vals.keys():

            Eg0 = Egi(matterial=self.matterial).update_iEg(
                temp=0, author=self.vals['egi_author'])
            Egratio = Eg0 / Egi(matterial=self.matterial).update_iEg(
                temp=temp, author=self.vals['egi_author'])
        else:
            Egratio = None

        self.Nc, self.Nv = getattr(dos_models, self.model)(
            self.vals, temp=temp, Egratio=Egratio)

        return self.Nc, self.Nv

    def check_models(self):
        temp = np.logspace(0, np.log10(600))
        num = len(self.available_models())

        fig, ax = plt.subplots(1)
        self.plotting_colours(num, fig, ax, repeats=2)

        for author in self.available_models():
            Nc, Nv = self.update(temp, author)
            # print Nc.shape, Nv.shape, temp.shape
            ax.plot(temp, Nc, '--')
            ax.plot(temp, Nv, '.', label=author)

        ax.loglog()
        leg1 = ax.legend(loc=0, title='colour legend')

        Nc, = ax.plot(np.inf, np.inf, 'k--', label='Nc')
        Nv, = ax.plot(np.inf, np.inf, 'k.', label='Nv')

        plt.legend([Nc, Nv], ['Nc', 'Nv'], loc=4, title='Line legend')
        plt.gca().add_artist(leg1)

        ax.set_xlabel('Temperature (K)')
        ax.set_ylabel('Density of states (cm$^{-3}$)')
        plt.show()
