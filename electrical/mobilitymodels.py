
import numpy as np
import matplotlib.pylab as plt
import sys
import os
import ConfigParser


def add_mobilities(self, mobility_list):
    imobility = 0

    for i in mobility_list:
        imobility += 1. / i

    return 1. / imobility



def CaugheyThomas( vals,  impurity, min_car_den, **kwargs):
    '''
    emperical form for one temperature taken from:
    D. M. Caughey and R. E. Thomas, Proc. U.E.E., pp. 2192,
    2193 (Dec. 1977).

    inputs:
        impurty: the number of impurities (cm^-3)
        min_carr_den: the number of minoirty carrier densities (cm^-3)
        maj_car:  the majority carrier type
        temp: temperature (K)
    output:
        mobility (cm^2 V^-1 s^-1)

    '''

    mu = vals['mu_min'] + (vals['mu_max'] - vals['mu_min']
                           ) / (1. + (impurity / vals['nr'])**vals['alpha'])
    return mu 

def dorkel(vals, impurity, min_car_den, maj_car_den, temp, carrier):
    '''
    not consistent with PVlihthouse at high injection

    inputs:
        impurty: the number of impurities (cm^-3)
        min_carr_den: the number of minoirty carrier densities (cm^-3)
        maj_car_den: the number of majority carrier densities (cm^-3)
        temp: temperature (K)
    output:
         electron mobility (cm^2 V^-1 s^-1)
         hole mobility (cm^2 V^-1 s^-1)
    '''

    # this relatves the carrier to the extension in the variable name
    if carrier == 'electron':
        carrier = 'e'
    elif carrier == 'hole':
        carrier = 'h'

    # determine hole dependent carrier partial mobilities
    mu_L = lattice_mobility(vals, temp, carrier)
    mu_i = impurity_mobility(vals, impurity, temp, carrier)

    # determine both carrier scattering mobilities
    mu_css = carrier_scattering_mobility(
        vals, min_car_den, maj_car_den, temp)

    # determine sudo function
    X = np.sqrt(6. * mu_L * (mu_i + mu_css) / (mu_i * mu_css))

    # combine partial moblities into total
    mu = mu_L * (1.025 / (1. + (X / 1.68)**(1.43)) - 0.025)

    return mu

def lattice_mobility(vals, temp, carrier):
    ''' 
    due to scattering of acoustic phonons
    '''
    mu_L = vals['mul0' + carrier] * \
        (temp / vals['temp0'])**(-vals['alpha' + carrier])
    return mu_L

def impurity_mobility( vals, impurity, temp, carrier):
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

def impurity_neutral():
    pass

def carrier_scattering_mobility( vals, min_car_den, maj_car_den, temp):
    '''
    The coefficient in B is 2e17 (equation 3) and not 2e7 (Equation 7) as presented in the paper

    '''

    A = np.log(1. + 8.28e8 * temp**2 / (min_car_den * maj_car_den)**(1. / 3))
    B = 2e17 * temp**(3. / 2) / np.sqrt(min_car_den * maj_car_den)
    mu_css = B / A
    return mu_css



class unified_mobility():

    """
    Thaken from: 

    [1] D. B. M. Klaassen, 
    "A unified mobility model for device simulation-I. Model equations and concentration dependence"
     Solid. State. Electron., vol. 35, no. 7, pp. 953-959, Jul. 1992. 

    [2] D. B. M. Klaassen,
    "A unified mobility model for device simulation-II. Temperature dependence of carrier mobility and lifetime,"
    Solid. State. Electron., vol. 35, no. 7, pp. 961-967, Jul. 1992.

    This is the Klaassen's mobility model, for which the calculations  with two exceptions: 
        (i) r5 is set to -0.8552 rather than -0.01552 (see Table 2 of [1]), 
        (ii) Eq. A3 of [1] is adjusted such that PCWe is determined with Ne,sc rather than (Z^3 Ni) 
         and PCWh is determined with Nh,sc rather than (Z^3 Ni);

    these changes give a better fit to the solid calculated lines in Figures 6 and 7 of [1], which better fits the experimental data. 
    These modifications are also contained in Sentaurus's version of Klaassen's model [5].
    Klaassen's mobility model fits reasonably with experimental data over an estimated temperature range of 100 - 450 K.
    Its accuracy is greatest at 300 K (see [1,2]).
    """

    # these are the values for phosphorous and boron respectively.

    # Original value
    # r5 = -0.01552, changing this means changing 2 equations as well

    # a switch used for different types
    # change to hle and electron for clarity
    type_dic = {'hole': 1, 'electron': 0}

    def __init__(self, vals, impurity, dn, N_a, N_d, temp, carrier):
        '''
        Both these need to be fixed
        impurity isn't used
        Na: acceptor ataoms, (creates p-type matterial)
        Nd: donor ataoms, (creates n-type matterial)
        '''

        self.Na = N_a
        self.Nd = N_d
        self.update_carriers(dn, N_d, N_a)
        return self.mobility(carrier)


    def update_carriers(self, deltan, p0, n0):

        # finding the majority carriers
        p0 = self.return_dopant('hole') - self.return_dopant('electron')
        if p0.all() > 0:
            n0 = self.ni**2 / p0
            # print 'p-type'
        else:
            # print 'n-type'
            n0 = -p0
            p0 = self.ni**2 / p0

        self.p = deltan + p0
        self.n = deltan + n0


    def mobility(self, carrier):

        self.update_carriers(deltan)

        return 1. / (1. / self.uDCS('carrier') + 1. / self.uLS('carrier'))


    def uLS(self, Type):
        i = self.type_dic[Type]
        return self.vals['umax'][i] * (300. / self.temp)**self.vals['theta'][i]

    def uDCS(self, Type):

        i = self.type_dic[Type]
        print 'here'
        # print  self.Nsc(Type) ,'\n', self.Nsceff(Type) ,'\n', (
        #     self.Nref[i] / self.Nsc(Type))**(self.alpha[i]), (
        #     self.uc(Type) * self.carrier_sum() / self.Nsceff(Type))

        return self.un(Type) * self.Nsc(Type) / self.Nsceff(Type) * (
            self.vals['Nref'][i] / self.Nsc(Type))**(self.vals['alpha'][i]) + (
            self.uc(Type) * self.carrier_sum() / self.Nsceff(Type))

    def un(self, Type):
        """
        majority dopant scattering (with screening)
        """
        # Done
        i = self.type_dic[Type]

        return self.vals['umax'][i] * self.vals['umax'][i] / (self.vals['umax'][i] - self.vals['umin'][i])

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

        return self.vals['umin'][i] * self.vals['umax'][i] / (self.vals['umax'][i] - self.vals['umin'][i])

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
            dopant = self.Na
        elif Type == 'electron':
            dopant = self.Nd
        return dopant

    def carrier_sum(self):

        return self.p + self.n

    def G(self, Type):
        """
        Accounts for minority impurity scattering
        """

        i = self.type_dic[Type]
        a = 1.
        b = - self.vals['s1'] / \
            (self.vals['s2'] + (self.temp / 300. / self.vals['mr'][i])
             ** self.vals['s4'] * self.P(Type))**self.vals['s3']
        c = self.vals['s5'] / \
            ((300. / self.temp / self.vals['mr'][i])**self.vals['s7'] * self.P(Type))**self.vals['s6']
        return a + b + c

    def P(self, Type):
        # Done
        return 1. / (self.vals['fCW'] / self.PCW(Type) + self.vals['fBH'] / self.PBH(Type))

    def PCW(self, Type):
        # Done
        return 3.97e13 * (
            1. / (self.Nsc(Type)) * ((self.temp / 300.)**(3.)))**(2. / 3.)

    def PBH(self, Type):
        # Done
        i = self.type_dic[Type]
        return 1.36e20 / self.carrier_sum() * (
            self.vals['mr'][i] * (self.temp / 300.0)**2.0)

    def F(self, Type):
        """
        Accounts for electron-hole scattering
        """
        # done
        i = self.type_dic[Type]
        # uses Since True == 1 and False == 0 in python
        j = (not i) * 1

        return (self.vals['r1'] * self.P(Type)**self.vals['r6']
                + self.vals['r2'] + self.vals['r3'] * self.vals['mr'][i] / self.vals['mr'][j]
                ) / (
            self.P(Type)**(self.vals['r6']) + self.vals['r4'] +
            self.vals['r5'] * self.vals['mr'][i] / self.vals['mr'][j])


