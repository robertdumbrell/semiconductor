
import numpy as np
import matplotlib.pylab as plt
import sys
import os
import ConfigParser


sys.path.append(
    os.path.abspath(os.path.join(os.getcwd(), os.pardir, os.pardir)))

from semiconductor.helper.helper import HelperFunctions


class Mobility(HelperFunctions):
    model_file = 'mobility.models'

    def __init__(self, matterial='Si', model_author=None, temp=300., ni=9.65e9):
        self.Models = ConfigParser.ConfigParser()
        self.matterial = matterial

        constants_file = os.path.join(
            os.path.dirname(os.path.realpath(__file__)),
            matterial,
            self.model_file)

        self.Models.read(constants_file)

        self.change_model(model_author)

    def mobility(self, min_car_den, Na, Nd, carrier_type):
        self.Nh_0, self.Ne_0 = self.check_doping(Na, Nd)
        impurity = Na + Nd
        return getattr(self, self.model)(min_car_den, impurity, carrier_type)

    def add_mobilities(self, mobility_list):
        imobility = 0

        for i in mobility_list:
            imobility += 1. / i

        return 1. / imobility

    def CaugheyThomas(self, vals,  min_car_den, impurity, carrier_type):
        '''
        emperical form
        D. M. Caughey and R. E. Thomas, Proc. U.E.E., pp. 2192,
        2193 (Dec. 1977).
        '''
        mu = vals['mu_min'] + (vals['mu_max'] - vals['mu_min']
                               ) / (1. + (impurity / vals['Nr'])**vals['alpha'])
        return mu

    def lattice_mobility(self, vals, temp, carrier):
        ''' 
        due to scattering of acoustic phonons
        '''
        mu_L = vals['mul0' + carrier] * \
            (temp / vals['temp0'])**(-vals['alpha' + carrier])
        return mu_L

    def impurity_mobility(self, vals, impurity, temp, carrier):
        '''
        interactions between the carriers and the ionized impurities.
        This partial mobility increases as the temperature
        increases or the doping concentration
        decreases. The relationship which we use in the calculatipn
        of the pr component is that of Brooks and
        Herring
        '''
        A = np.log(1. + vals['b' + carrier] * temp**2 / impurity)
        B = (vals['b' + carrier] * temp ** 2.) / \
            (impurity + vals['b' + carrier] * temp**2)
        mu_i = vals['a' + carrier] * temp**(3. / 2) / impurity / (A - B)
        return mu_i

    def impurity_neutral(self):
        pass

    def carrier_scattering_mobility(self, vals, min_car_den, maj_car_den, temp):
        '''
        The coefficient in B is 2e17 (equation 3) and not 2e7 (Equation 7) as presented in the paper

        '''

        A = np.log(1. + 8.28e8 * temp**2 / (min_car_den * maj_car_den)**(1. / 3))
        B = 2e17 * temp**(3. / 2) / np.sqrt(min_car_den * maj_car_den)
        mu_css = B / A
        return mu_css

    def dorkel(self, vals, impurity, min_car_den, maj_car_den, temp):
        '''
        not consistent with PVlihthouse at high injection
        '''

        # determine hole dependent carrier partial mobilities
        mu_Lh = self.lattice_mobility(vals, temp, 'h')
        mu_ih = self.impurity_mobility(vals, impurity, temp, 'h')

        # determine electron dependent carrier partial mobilities
        mu_Le = self.lattice_mobility(vals, temp, 'e')
        mu_ie = self.impurity_mobility(vals, impurity, temp, 'e')

        # determine both carrier scattering mobilities
        mu_css = self.carrier_scattering_mobility(
            vals, min_car_den, maj_car_den, temp)

        # determine sudo function
        Xh = np.sqrt(6. * mu_Lh * (mu_ih + mu_css) / (mu_ih * mu_css))
        Xe = np.sqrt(6. * mu_Le * (mu_ie + mu_css) / (mu_ie * mu_css))

        # combine partial moblities into total
        mu_h = mu_Lh * (1.025 / (1. + (Xh / 1.68)**(1.43)) - 0.025)
        mu_e = mu_Le * (1.025 / (1. + (Xe / 1.68)**(1.43)) - 0.025)

        print mu_ie, mu_Le, mu_css, mu_e
        return mu_e, mu_h

if __name__ == "__main__":

    a = Mobility('Si')

    dn = np.logspace(10, 20)
    Nd = 1e16
    mob_e, mob_h = a.dorkel(a.vals, Nd, dn, dn + Nd, 300.)
    # print a.Radiative.ni, a.Auger.ni

    plt.plot(dn, mob_e)
    plt.plot(dn, mob_h)
    plt.ylim(bottom = 0)

    plt.semilogx()
    plt.show()
