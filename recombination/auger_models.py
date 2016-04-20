
import numpy as np
import os


def auger_dopants(vals, nxc, ne0, nh0):
    '''This is the classic auger model that only depends on doping

    The model requires two cosntants
    cn and cp

    It requires the  the dark carrier concentrations of the
    sample to be known
    '''
    Ce = vals['cn']
    Ch = vals['cp']

    nh = nh0 + nxc
    ne = ne0 + nxc

    R = Ce * ne**2 * nh + Ch * ne * nh**2

    return nxc / R


def auger(vals, nxc, ne0, nh0):
    '''This is the classic auger model that includes 
    the impact of excess carriers

    it requires 3 constants. One for holes, 
    one for electrons and one for carrier to carrier interaction

    It requires the  the dark carrier concentrations of the
    sample to be known
    '''

    nh = nh0 + nxc
    ne = ne0 + nxc

    Ce = vals['ced'] * ne0 / \
        (ne0 + nh) + vals['ccc'] / 2 * nh / (nh + ne0)
    Ch = vals['chd'] * nh0 / \
        (nh0 + ne) + vals['ccc'] / 2 * ne / (ne + nh0)

    R = (Ce * ne + Ch * nh) * (ne * nh - nh0 * ne0)

    return nxc / R


def coulomb_enhanced_auger(vals, nxc, ne0, nh0):
    '''The coulomb enhanced auger model, as proposed by Kerr2002

    The model requires the inputs
    K_n, K_p, K_Delta, delta, L_eeh, N_eeh, K_eeh, P_ehh, L_ehh, K_ehh

    It requires the  the dark carrier concentrations of the
    sample to be known

    There may be an older form from Schmit. 
    '''
    # then need to get doping and delta n
    nh = nh0 + nxc
    ne = ne0 + nxc

    g_eeh = 1. + vals['l_eeh'] * (1. - np.tanh(
        (ne0 / vals['k_eeh'])**vals['n_eeh']))
    g_ehh = 1. + vals['l_ehh'] * (1. - np.tanh(
        (nh0 / vals['k_ehh'])**vals['p_ehh']))

    # Then the lifetime is provided by this
    return nxc / (ne * nh - ne0 * nh0) /\
        (vals['k_n'] * g_eeh * ne0 +
         vals['k_p'] * g_ehh *
         nh0 +
         vals['k_delta'] *
         nxc**vals['delta']
         )


def auger_notimplimented(*args):
    """ A dummy class for models with incomplete parameters"""
    pass


def Richter2012_check(vals, func, fig, ax):

    folder = os.path.join(os.path.dirname(__file__), 'Si', 'check_data')
    files = ['Richter_dop_1e15.csv',
             'Richter_dop_1e17.csv',
             'Richter_dop_1e19.csv']

    nxc = np.logspace(10, 20)

    for fname, doping in zip(files, [1e15, 1e17, 1e19, 1e19]):
        ne0 = 1e20 / doping
        nh0 = doping
        taus = func(vals, nxc, ne0, nh0)

        data = np.genfromtxt(
            os.path.join(folder, fname), delimiter=',', names=True)

        ax.plot(nxc, taus)

        ax.plot(data['nxc'], data['tau'], '--')
    ax.set_xlabel('Excess carrer density (cm^-3)')
    ax.set_ylabel('$\tau$ (cm$^{-3}$)')
    ax.loglog()
