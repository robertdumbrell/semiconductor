#!/usr/local/bin/python
# UTF-8

import numpy as np
import semiconductor.general_functions.carrierfunctions as GF


def add_mobilities(self, mobility_list):
    imobility = 0

    for i in mobility_list:
        imobility += 1. / i

    return 1. / imobility


def CaugheyThomas(vals, Na, Nd, nxc, **kwargs):
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
    impurity = Na, Nd
    mu = vals['mu_min'] + (vals['mu_max'] - vals['mu_min']
                           ) / (1. + (impurity / vals['nr'])**vals['alpha'])
    return mu


def dorkel(vals, Na, Nd, nxc, temp, carrier, **kwargs):
    '''
    inputs:
        impurty: the number of impurities (cm^-3)
        min_carr_den: the number of minoirty carrier densities (cm^-3)
        maj_car_den: the number of majority carrier densities (cm^-3)
        temp: temperature (K)
    output:
         electron mobility (cm^2 V^-1 s^-1)
         hole mobility (cm^2 V^-1 s^-1)
    '''

    impurity = Na + Nd

    ne, nh = GF.get_carriers(Na,
                             Nd,
                             nxc,
                             temp=temp)
    print Na, Nd, nxc, temp

    if np.all(nh < ne):
        nxc = nh
        maj_car_den = ne
    elif np.all(nh >= ne):
        maj_car_den = nh
        nxc = ne
    else:
        sys.exit("Mixed doping types not permitted")

    # this relates the carrier to the extension in the variable name
    if carrier == 'electron':
        carrier = 'e'
    elif carrier == 'hole':
        carrier = 'h'
    else:
        sys.exit("inappropriate values for carrier passed")

    # determine hole dependent carrier partial mobilities
    mu_L = lattice_mobility(vals, temp, carrier)
    mu_i = impurity_mobility(vals, impurity, temp, carrier)

    # determine both carrier scattering mobilities
    mu_css = carrier_scattering_mobility(
        vals, nxc, maj_car_den, temp)

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


def impurity_mobility(vals, impurity, temp, carrier):
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


def carrier_scattering_mobility(vals, nxc, maj_car_den, temp):
    '''
    The coefficient in B is 2e17 (equation 3) and not 2e7 (Equation 7) as presented in the paper

    '''

    A = np.log(1. + 8.28e8 * temp**2 / (nxc * maj_car_den)**(1. / 3))
    B = 2e17 * temp**(3. / 2) / np.sqrt(nxc * maj_car_den)
    mu_css = B / A
    return mu_css


# below this are the functions for klaassen's model


def unified_mobility(vals, Na, Nd, nxc, temp, carrier):
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
    type_dic = {'hole': 'h', 'electron': 'e'}

    if carrier in type_dic:
        carrier = type_dic[carrier]
    else:
        print 'incorrect input for carrier input'

    # Things to fix up
    # ni = ni

    # the only thing ni is used for, this can be factored out so these values
    # are passed to this function

    ne, nh = GF.get_carriers(Na,
                             Nd,
                             nxc,
                             temp=temp)

    return 1. / (
        1. / uDCS(carrier, vals, nh, ne, Na, Nd, temp) +
        1. / uLS(carrier, vals, temp))


def uLS(carrier, vals, temp):

    return vals['umax_' + carrier] * (300. / temp)**vals['theta_' + carrier]


def uDCS(carrier, vals, nh, ne, Na, Nd, temp):
    carrier_sum = nh+ne

    return un(carrier, vals, temp) * Nsc(carrier, vals, nh, ne, Na, Nd) / Nsceff(carrier, vals, nh, ne, Na, Nd, temp) * (
        vals['nref_' + carrier] / Nsc(carrier, vals, nh, ne, Na, Nd))**(vals['alpha_' + carrier]) + (
        uc(carrier, vals, temp) * carrier_sum / Nsceff(carrier, vals, nh, ne, Na, Nd, temp))


def un(carrier, vals, temp):
    """
    majority dopant scattering (with screening)
    """
    # Done
    return vals['umax_' + carrier] * vals['umax_' + carrier] / (vals['umax_' + carrier] - vals['umin_' + carrier]) *\
        (temp / 300.)**(3. * vals['alpha_' + carrier] - 1.5)


def Nsc(carrier, vals, nh, ne, Na, Nd):

    # checked
    car_den = return_carrer(carrier, nh, ne, opposite=True)

    return return_dopant('e', Na, Nd) * Z('e', vals, Na, Nd) + (
        return_dopant('h', Na, Nd) * Z('h', vals, Na, Nd) +
        car_den)


def Nsceff(carrier, vals, nh, ne, Na, Nd, temp):

    # checked

    car_den = return_carrer(carrier, nh, ne, opposite=True)

    if carrier == 'e':
        N_a = G(carrier, vals, nh, ne, Na, Nd, temp)
        N_a *= return_dopant('h', Na, Nd) * Z('h', vals, Na, Nd)
        N_d = return_dopant('e', Na, Nd) * Z('e', vals, Na, Nd)
    elif carrier == 'h':
        N_d = G(carrier, vals, nh, ne, Na, Nd, temp)
        N_d *= return_dopant('e', Na, Nd) * Z('e', vals, Na, Nd)
        N_a = return_dopant('h', Na, Nd) * Z('h', vals, Na, Nd)

    else:
        print 'Something has gone wrong in the code'

    # plt.figure('test')
    # print ' starting:'
    # print N_d, N_a, car_den, carrier
    # print F(carrier, vals, p,n, Na, Nd, temp)
    # plt.loglog()
    # plt.show()
    return N_a + N_d + car_den / F(carrier, vals, nh, ne, Na, Nd, temp)


def Z(carrier, vals, Na, Nd):
    """
    accounts for high doping effects - clustering
    """

    return 1. + 1. / (vals['c_' + carrier] +
                      (vals['nref2_' + carrier] / return_dopant(carrier, Na, Nd))**2.)


def G(carrier, vals, nh, ne, Na, Nd, temp):
    """
    Accounts for minority impurity scattering
    """

    P_value = P(carrier, vals, nh, ne, Na, Nd, temp)

    a = 1.
    b = - vals['s1'] / \
        (vals['s2'] + (temp / 300. / vals['mr_' + carrier])
         ** vals['s4'] * P_value)**vals['s3']
    c = vals['s5'] / \
        ((300. / temp /
          vals['mr_' + carrier])**vals['s7'] * P_value)**vals['s6']
    return a + b + c


def P(carrier, vals, nh, ne, Na, Nd, temp):
    return 1. / (vals['fcw'] / PCW(carrier, vals, nh, ne, Na, Nd, temp)
                 + vals['fbh'] / PBH(carrier, vals, temp, ne + nh))


def PCW(carrier, vals, nh, ne, Na, Nd, temp):
    # Done
    return 3.97e13 * (
        1. / (Nsc(carrier, vals, nh, ne, Na, Nd)) * ((temp / 300.)**(3.)))**(2. / 3.)


def PBH(carrier, vals, temp, carrier_sum):
    # Done
    return 1.36e20 / carrier_sum * (
        vals['mr_' + carrier] * (temp / 300.0)**2.0)


def F(carrier, vals, nh, ne, Na, Nd, temp):
    """
    Accounts for electron-hole scattering
    """
    # done
    # uses Since True == 1 and False == 0 in python

    switch = {'e': 'h', 'h': 'e'}

    return (vals['r1'] * P(carrier, vals, nh, ne, Na, Nd, temp)**vals['r6']
            + vals['r2'] + vals['r3'] *
            vals['mr_' + carrier] / vals['mr_' + switch[carrier]]
            ) / (
        P(carrier, vals, nh, ne, Na, Nd, temp)**(vals['r6']) + vals['r4'] +
        vals['r5'] * vals['mr_' + carrier] / vals['mr_' + switch[carrier]])


def uc(carrier, vals, temp):
    """
    excess carrier scattering
    """
    # Done

    return vals['umin_' + carrier] * vals['umax_' + carrier] / (
        vals['umax_' + carrier] - vals['umin_' + carrier]) * (300. / temp)**0.5


def return_carrer(carrier, nh, ne, opposite=True):

    if opposite:
        switch = {'e': 'h', 'h': 'e'}
        carrier = switch[carrier]

    if carrier == 'h':
        car_den = nh
    elif carrier == 'e':
        car_den = ne
    return car_den


def return_dopant(carrier, Na, Nd):

    if carrier == 'h':
        dopant = np.array([Na]).flatten()
    elif carrier == 'e':
        dopant = np.array([Nd]).flatten()

    return dopant
