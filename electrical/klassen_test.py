'''
Created on Aug 19, 2012

@author: mattias

This is about mobility, here type (Type) doesn't not stand
for the wafer type but rather the carrier type.

I can't remember if this is the sentarous model or not
'''

# import tkFileDialog, sys, re, numpy,scipy
import numpy as np
import time
import matplotlib.pylab as plt
import os
from scipy.optimize import minimize
# import glob
# import os
# import scipy.optimize


###
# This is for testing new things, like the mobility model im about to right
###

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
    These modifications are also contained in Sentaurus's version of Klaassen's model .
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

    Temp = 300.
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

    # check that G doesn't get large
    min_G_check = {'temp': -1, 'P': 300, 'val': 1000}

    def update_carriers(self, deltan):

        # finding the majority carriers
        self.p0 = self.return_dopant('hole') - self.return_dopant('electron')
        # print self.return_dopant('hole')- self.return_dopant('electron')

        if np.all(self.p0 > 0):
            self.n0 = self.ni**2 / self.p0
            # print 'p-type'
        elif np.all(self.p0 < 0):
            # print 'n-type'
            self.n0 = -self.p0
            self.p0 = self.ni**2 / self.p0
        else:
            print 'how can there be no net dopants?'
            print self.p0

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
        return self.umax[i] * (300. / self.Temp)**self.theta[i]

    def uDCS(self, Type):
        """
        I think this is mu i"""
        i = self.type_dic[Type]

        return self.un(Type) * self.Nsc(Type) / self.Nsceff(Type) * (
            self.Nref[i] / self.Nsc(Type))**(self.alpha[i]) + (
            self.uc(Type) * self.carrier_sum() / self.Nsceff(Type))

    def un(self, Type):
        """
        majority dopant scattering (with screening)
        """
        # Done
        i = self.type_dic[Type]

        return self.umax[i] * self.umax[i] / (self.umax[i] - self.umin[i]) * (self.Temp / 300.)**(3. * self.alpha[i] - 1.5)

    def Nsc(self, Type):
        """
        This is Z**2 N_i 
        """
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

        return self.umin[i] * self.umax[i] / (self.umax[i] - self.umin[i]) * (300. / self.Temp)**0.5

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

    def return_dopant(self, Type):
        # print Type
        if Type == 'hole':
            dopant = self.N_a
        elif Type == 'electron':
            dopant = self.N_d
        return dopant

    def carrier_sum(self):

        return self.p + self.n

    def find_min_G(self):
        '''
        As mentioned in 1992:
        UNIFIED MOBILITY MODEL FOR DEVICE SIMULATION--II. TEMPERATURE DEPENDENCE OF
        CARRIER MOBILITY AND LIFETIME

        The function of G should not increase again at low values of P. This prevents that
        '''

        def test_G(P, Type):

            a = 1.
            b = - self.s1 / \
                (self.s2 + (self.Temp / 300. / self.mr)
                 ** self.s4 * P)**self.s3
            c = self.s5 / \
                ((300. / self.Temp / self.mr)**self.s7 * P)**self.s6

            b = b[self.type_dic[Type]]
            c = c[self.type_dic[Type]]

            return (a + b + c)

        # print type(test_G(1))
        # print minimize(test_G,0.5)
        # print minimize(test_G, 0.5, args=('electron',))
        temp = minimize(test_G, 0.5, args=('electron',))['x'][0]
        temp2 = minimize(test_G, 0.5, args=('hole',))['x'][0]

        self.min_G_check['P'] = [temp, temp2]

        self.min_G_check['G'] = [
                            test_G(self.min_G_check['P'], 'electron'),
                            test_G(self.min_G_check['P'], 'hole')
                                ]

        self.min_G_check['temp'] = self.Temp


    def G(self, Type):
        """
        Accounts for minority impurity scattering
        """

        if self.min_G_check['temp'] != self.Temp:
            self.find_min_G()

        P = self.P(Type)
        i = self.type_dic[Type]


        P[self.min_G_check['P'][i] > P] = self.min_G_check['P'][i]

        a = 1.
        b = - self.s1 / \
            (self.s2 + (self.Temp / 300. / self.mr[i])
             ** self.s4 * P)**self.s3
        c = self.s5 / \
            ((300. / self.Temp / self.mr[i])**self.s7 * P)**self.s6
        return a + b + c

    def P(self, Type):
        # Done
        return 1. / (self.fCW / self.PCW(Type) + self.fBH / self.PBH(Type))

    def PCW(self, Type):
        """
        Eq. A3 of [1] is adjusted such that PCWe is determined with Ne,sc rather than (Z^3 Ni) 
        and PCWh is determined with Nh,sc rather than (Z^3 Ni)
        """

        return 3.97e13 * (
            1. / (self.Nsc(Type)) * ((self.Temp / 300.)**(3.)))**(2. / 3.)

    def PBH(self, Type):
        # Done
        i = self.type_dic[Type]
        return 1.36e20 / self.carrier_sum() * (
            self.mr[i] * (self.Temp / 300.0)**2.0)

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


class Mobility_Si():

    N_a = 0
    N_d = 1e16


def Check_Doping():
    plt.figure('Doping Check')
    Mob = Mobility_Klassen()
    Mob.N_a = np.array([0])
    Mob.N_d = np.logspace(12, 20)

    Mob.N_d = np.array([0])
    Mob.N_a = np.logspace(12, 20)

    deltan = np.array([1e14])

    plt.plot(Mob.N_a, Mob.mobility_electron(deltan), label='electron type')
    plt.plot(Mob.N_a, Mob.mobility_hole(deltan), label='hole type')

    # Compare against data from PVlighthouse

    folder = os.path.join(os.getcwd(), 'Si')
    fname = 'Klassen_1e14_carriers.dat'

    data = np.genfromtxt(os.path.join(folder, fname), names=True)

    plt.plot(data['Ndop'], data['ue'], '--',
             label='electron - PV-lighthouse')
    plt.plot(data['Ndop'], data['uh'], '--',
             label='hole - PV-lighthouse')

    plt.semilogx()
    plt.legend(loc=0)
    plt.grid(True)


def Check_Carriers():
    plt.figure('Carrier Check')
    Mob = Mobility_Klassen()

    Mob.N_d = np.array([0])
    Mob.N_a = np.array([1e14])

    deltan = np.logspace(10, 20)

    plt.plot(deltan, Mob.mobility_electron(deltan), label='electron mobility')
    plt.plot(deltan, Mob.mobility_hole(deltan), label='hole mobility')

    plt.semilogx()
    plt.legend(loc=0)
    plt.grid(True)

    folder = os.path.join(os.getcwd(), 'Si')
    fname = 'Klassen_1e14_dopants.dat'
    data = np.genfromtxt(os.path.join(folder, fname), names=True)
    # print data.dtype.names
    plt.plot(data['deltan'], data['ue'], '--',
             label='electron - PV-lighthouse')
    plt.plot(data['deltan'], data['uh'], '--', label='hole - PV-lighthouse')
    plt.legend(loc=0)
    # plt.show()

    plt.figure('components')
    print Mob.uLS('electron'), Mob.uLS('hole')
    print Mob.un('electron'), Mob.un('hole')
    print Mob.uc('electron'), Mob.uc('hole')
    plt.plot(deltan, Mob.Nsc('electron'))
    plt.plot(deltan, Mob.Nsc('hole'))
    plt.plot(deltan, Mob.Nsceff('electron'))
    plt.plot(deltan, Mob.Nsceff('hole'))

    plt.loglog()

def check_G():
    plt.figure('Carrier Check')
    Mob = Mobility_Klassen()

    for temp in [300,400,500,600]:
        Mob.Temp = temp
        Mob.N_d = np.array([0])
        Mob.N_a = np.array([1e14])

        deltan = np.logspace(10, 20)
        Mob.update_carriers(np.logspace(12, 20))

        plt.plot(Mob.P('electron'), Mob.G('electron'), label='electron mobility')
        plt.plot(Mob.P('hole'), Mob.G('hole'), label='hole mobility')
        plt.plot(Mob.min_G_check['P'], Mob.min_G_check[
                 'G'], 'ro', label='hole mobility')

        plt.loglog()
        plt.legend(loc=0)
        plt.grid(True)

    # folder = os.path.join(os.getcwd(), 'Si')
    # fname = 'Klassen_1e14_dopants.dat'
    # data = np.genfromtxt(os.path.join(folder, fname), names=True)
    # print data.dtype.names
    # plt.plot(data['deltan'], data['ue'], '--',
             # label='electron - PV-lighthouse')
    # plt.plot(data['deltan'], data['uh'], '--', label='hole - PV-lighthouse')
    plt.legend(loc=0)


def check_temp():
    plt.figure('Temp Check')
    Mob = Mobility_Klassen()

    for Temp in [450.]:
        Mob.Temp = Temp
        Mob.N_d = np.array([0.])
        Mob.N_a = np.array([1e14])

        deltan = np.logspace(10, 20)

        plt.plot(deltan, Mob.mobility_electron(
            deltan), label='electron mobility')
        plt.plot(deltan, Mob.mobility_hole(deltan), label='hole mobility')

        plt.semilogx()
        plt.legend(loc=0)
        plt.grid(True)

        folder = os.path.join(os.getcwd(), 'Si')
        fname = 'Klassen_1e14_temp-450.dat'
        data = np.genfromtxt(os.path.join(folder, fname), names=True)

        plt.plot(data['deltan'], data['ue'], '--',
                 label='electron - PV-lighthouse')
        plt.plot(
            data['deltan'], data['uh'], '--', label='hole - PV-lighthouse')
        plt.legend(loc=0)

        # Mob.mobility_electron(deltan)
        # plt.plot(Mob.P('hole'),Mob.G('hole'))
        # plt.loglog()


if __name__ == '__main__':
    # check_G()
    # Check_Doping()
    # Check_Carriers()
    check_temp()
    plt.show()
