
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

        A = np.log(
            1. + 8.28e8 * temp**2 / (min_car_den * maj_car_den)**(1. / 3))
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

        # print mu_ie, mu_Le, mu_css, mu_e
        return mu_e, mu_h

    def klassen(self, vals, impurity, min_car_den, maj_car_den, temp):
        
        return



class Mobility_Klassen():

    """
    Thaken from: 

    [1] D. B. M. Klaassen, 
    "A unified mobility model for device simulation-I. Model equations and concentration dependence"
     Solid. State. Electron., vol. 35, no. 7, pp. 953-959, Jul. 1992. 

    [2] D. B. M. Klaassen,
    "A unified mobility model for device simulation-II. Temperature dependence of carrier mobility and lifetime,"
    Solid. State. Electron., vol. 35, no. 7, pp. 961-967, Jul. 1992.

    additional comments taken from https://www.pvlighthouse.com.au/calculators/Mobility%20calculator/Mobility%20calculator.aspx

    This is the Klaassen's mobility model, for which the calculations  with two exceptions: 
        (i) r5 is set to -0.8552 rather than -0.01552 (see Table 2 of [1]), 
        (ii) Eq. A3 of [1] is adjusted such that PCWe is determined with Ne,sc rather than (Z^3 Ni) 
         and PCWh is determined with Nh,sc rather than (Z^3 Ni);

    these changes give a better fit to the solid calculated lines in Figures 6 and 7 of [1], which better fits the experimental data. 
    These modifications are also contained in Sentaurus's version of Klaassen's model.
    Klaassen's mobility model fits reasonably with experimental data over an estimated temperature range of 100 - 450 K.
    Its accuracy is greatest at 300 K (see [1,2]).
    """

    # these are the values for phosphorous and boron respectively.
    umax = np.array([1414, 470.5])
    umin = np.array([68.5, 44.9])
    theta = np.array([2.285, 2.247])

    ni = 9.66e9

    # Nref        = array([9.68e16, 2.23e17]) #sentarous - arsnic
    Nref = np.array([9.2e16, 2.23e17])
    alpha = np.array([.711, .719])

    c = np.array([0.21, 0.5])
    Nref2 = np.array([4e20, 7.2e20])

    fCW, fBH = 2.459, 3.828

    T = 300
    mr = np.array([1., 1.258])
    # other values of mr?
    # mr = [1./1.258,1.258]  #value taken from
    # users.cecs.anu.edu.au/~u5096045/QSSModel52.xls is m1/m2

    s1, s2, s3, = .89233, .41372, .19778
    s4, s5, s6, s7 = .28227, .005978, 1.80618, 0.72169

    r1, r2, r3, r4, r5, r6 = .7643, 2.2999, 6.5502, 2.367, -0.8552, .6478
    # Original value
    # r5 = -0.01552, changing this means changing 2 equations as well

    # a switch used for different types
    # change to hle and electron for clarity
    type_dic = {'hole': 1, 'electron': 0}

    def update_carriers(self, deltan):

        # finding the majority carriers
        self.p0 = self.return_dopant('hole') - self.return_dopant('electron')
        if self.p0.all() > 0:
            self.n0 = self.ni**2 / self.p0
            # print 'p-type'
        else:
            # print 'n-type'
            self.n0 = -self.p0
            self.p0 = self.ni**2 / self.p0

        self.p = deltan + self.p0
        self.n = deltan + self.n0

    def mobility_sum(self, deltan):

        return self.mobility_hole(deltan) + self.mobility_electron(deltan)

    def mobility_hole(self, deltan):

        self.update_carriers(deltan)

        return 1. / (1. / self.uDCS('hole') + 1. / self.uLS('hole'))

    def mobility_electron(self, deltan):

        self.update_carriers(deltan)

        return 1 / (1 / self.uDCS('electron') + 1 / self.uLS('electron'))

    def uLS(self, Type):
        i = self.type_dic[Type]
        return self.umax[i] * (300. / self.T)**self.theta[i]

    def uDCS(self, Type):

        i = self.type_dic[Type]
        print 'here'
        # print  self.Nsc(Type) ,'\n', self.Nsceff(Type) ,'\n', (
        #     self.Nref[i] / self.Nsc(Type))**(self.alpha[i]), (
        #     self.uc(Type) * self.carrier_sum() / self.Nsceff(Type))

        return self.un(Type) * self.Nsc(Type) / self.Nsceff(Type) * (
            self.Nref[i] / self.Nsc(Type))**(self.alpha[i]) + (
            self.uc(Type) * self.carrier_sum() / self.Nsceff(Type))

    def un(self, Type):
        """
        majority dopant scattering (with screening)
        """
        # Done
        i = self.type_dic[Type]

        return self.umax[i] * self.umax[i] / (self.umax[i] - self.umin[i])

    def Nsc(self, Type):

        # checked
        carrier = self.return_carrer(Type, opposite=True)

        # print  (self.return_dopant('hole')) * self.Z('hole')

        return (self.return_dopant('electron') * self.Z('electron')) + (
            self.return_dopant('hole') * self.Z('hole') +
            carrier)

    def Nsceff(self, Type):

        # checked
        carrier = self.return_carrer(Type, opposite=True)

        if Type == 'electron':
            Na = self.G(Type)
            Na *= self.return_dopant('hole') * self.Z('hole')
            Nd = self.return_dopant('electron') * self.Z('electron')
        elif Type == 'hole':
            Nd = self.G(Type)
            Nd *= self.return_dopant('electron') * self.Z('electron')
            Na = self.return_dopant('hole') * self.Z('hole')

        return Na + Nd + carrier / self.F(Type)

    def Z(self, Type):
        """
        accounts for high doping effects - clustering
        """
        # Done
        i = self.type_dic[Type]
        return 1. + 1. / (self.c[i] +
                          (self.Nref2[i] / self.return_dopant(Type))**2.)

    def uc(self, Type):
        """
        excess carrier scattering
        """
        # Done

        i = self.type_dic[Type]

        return self.umin[i] * self.umax[i] / (self.umax[i] - self.umin[i])

    def return_carrer(self, Type, opposite=True):

        if opposite:
            if Type == 'hole':
                Type = 'electron'
            else:
                Type = 'hole'

        if Type == 'hole':
            carrier = self.p
        elif Type == 'electron':
            carrier = self.n
        return carrier

    def return_dopant(self, Type, opposite=True):

        if Type == 'hole':
            dopant = self.N_a
        elif Type == 'electron':
            dopant = self.N_d
        return dopant

    def carrier_sum(self):

        return self.p + self.n

    def G(self, Type):
        """
        Accounts for minority impurity scattering
        """

        i = self.type_dic[Type]
        a = 1.
        b = - self.s1 / \
            (self.s2 + (self.T / 300. / self.mr[i])
             ** self.s4 * self.P(Type))**self.s3
        c = self.s5 / \
            ((300. / self.T / self.mr[i])**self.s7 * self.P(Type))**self.s6
        return a + b + c

    def P(self, Type):
        # Done
        return 1. / (self.fCW / self.PCW(Type) + self.fBH / self.PBH(Type))

    def PCW(self, Type):
        # Done
        return 3.97e13 * (
            1. / (self.Nsc(Type)) * ((self.T / 300.)**(3.)))**(2. / 3.)

    def PBH(self, Type):
        # Done
        i = self.type_dic[Type]
        return 1.36e20 / self.carrier_sum() * (
            self.mr[i] * (self.T / 300.0)**2.0)

    def F(self, Type):
        """
        Accounts for electron-hole scattering
        """
        # done
        i = self.type_dic[Type]
        # uses Since True == 1 and False == 0 in python
        j = (not i) * 1

        return (self.r1 * self.P(Type)**self.r6
                + self.r2 + self.r3 * self.mr[i] / self.mr[j]
                ) / (
            self.P(Type)**(self.r6) + self.r4 +
            self.r5 * self.mr[i] / self.mr[j])

if __name__ == "__main__":

    a = Mobility('Si')
    print a.Models.sections()
    a.available_models()

    dn = np.logspace(10, 20)
    Nd = 1e16
    mob_e, mob_h = a.dorkel(a.vals, Nd, dn, dn + Nd, 300.)
    # print a.Radiative.ni, a.Auger.ni

    plt.plot(dn, mob_e)
    plt.plot(dn, mob_h)
    plt.ylim(bottom=0)

    plt.semilogx()
    plt.show()
