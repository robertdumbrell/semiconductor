
import numpy as np
import os


def auger_dopants(vals, min_car_den, ne0, nh0):
    '''This is the classic auger model that only depends on doping

    The model requires two cosntants
    cn and cp

    It requires the  the dark carrier concentrations of the
    sample to be known
    '''
    Ce = vals['cn']
    Ch = vals['cp']

    nh = nh0 + min_car_den
    ne = ne0 + min_car_den

    R = Ce * ne**2 * nh + Ch * ne * nh**2

    return min_car_den / R


def auger(vals, min_car_den, ne0, nh0):
    '''This is the classic auger model that includes 
    the impact of excess carriers

    it requires 3 constants. One for holes, 
    one for electrons and one for carrier to carrier interaction

    It requires the  the dark carrier concentrations of the
    sample to be known
    '''

    nh = nh0 + min_car_den
    ne = ne0 + min_car_den

    Ce = vals['ced'] * ne0 / \
        (ne0 + nh) + vals['ccc'] / 2 * nh / (nh + ne0)
    Ch = vals['chd'] * nh0 / \
        (nh0 + ne) + vals['ccc'] / 2 * ne / (ne + nh0)

    R = (Ce * ne + Ch * nh) * (ne * nh - nh0 * ne0)

    return min_car_den / R


def coulomb_enhanced_auger(vals, min_car_den, ne0, nh0):
    '''The coulomb enhanced auger model, as proposed by Kerr2002

    The model requires the inputs
    K_n, K_p, K_Delta, delta, L_eeh, N_eeh, K_eeh, P_ehh, L_ehh, K_ehh

    It requires the  the dark carrier concentrations of the
    sample to be known

    There may be an older form from Schmit. 
    '''
    # then need to get doping and delta n
    nh = nh0 + min_car_den
    ne = ne0 + min_car_den

    g_eeh = 1. + vals['l_eeh'] * (1. - np.tanh(
        (ne0 / vals['k_eeh'])**vals['n_eeh']))
    g_ehh = 1. + vals['l_ehh'] * (1. - np.tanh(
        (nh0 / vals['k_ehh'])**vals['p_ehh']))

    # Then the lifetime is provided by this
    return min_car_den / (ne * nh - ne0 * nh0) /\
        (vals['k_n'] * g_eeh * ne0 +
         vals['k_p'] * g_ehh *
         nh0 +
         vals['k_delta'] *
         min_car_den**vals['delta']
         )


def auger_notimplimented(*args):
    """ A dummy class for models with incomplete parameters"""
    pass


def Richter2012_check(vals, func, fig, ax):

    folder = os.path.join(os.path.dirname(__file__), 'Si', 'check_data')
    files = ['Richter_dop_1e15.csv',
             'Richter_dop_1e17.csv',
             'Richter_dop_1e19.csv']

    min_car_den = np.logspace(10, 20)

    for fname, doping in zip(files, [1e15, 1e17, 1e19, 1e19]):
        ne0 = 1e20 / doping
        nh0 = doping
        taus = func(vals, min_car_den, ne0, nh0)

        data = np.genfromtxt(
            os.path.join(folder, fname), delimiter=',', names=True)

        ax.plot(min_car_den, taus)

        ax.plot(data['min_car_den'], data['tau'], '--')
    ax.set_xlabel('Excess carrer density (cm^-3)')
    ax.set_ylabel('$\tau$ (cm$^{-3}$)')
    ax.loglog()
