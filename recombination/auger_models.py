
import numpy as np
import os
from semiconductor.helper.helper import Webplotdig_JSONreader


def auger_dopants(vals, nxc, ne0, nh0, **args):
    '''
    This is the classic auger model that only depends on doping

    The model requires two constants
    Cn and Cp

    It requires the  the dark carrier concentrations of the
    sample to be known
    '''
    Ce = vals['cn']
    Ch = vals['cp']

    nh = nh0 + nxc
    ne = ne0 + nxc

    R = Ce * ne**2 * nh + Ch * ne * nh**2

    return nxc / R


def auger(vals, nxc, ne0, nh0, **args):
    '''
    This is the classic auger model that includes
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


def coulomb_enhanced_auger_Altermatt(vals, nxc, ne0, nh0, temp):
    '''
    The coulomb enhanced auger model, as proposed by Altermatt
    This uses equations 2, 3, 5, 6 from 10.1063/1.370784

    This model include the temperature dependence of Cn and Cp
    '''

    # then need to get doping and delta n
    nh = nh0 + nxc
    ne = ne0 + nxc

    gmaxn = vals['k_gmaxn'] * temp**vals['p_gmaxn']
    gmaxh = vals['k_gmaxh'] * temp**vals['p_gmaxh']

    # enhancement factors
    g_eeh = 1. + (gmaxn - 1) * (1. - np.tanh((
        (ne / vals['k_eeh'])**vals['n_eeh'])))
    g_ehh = 1. + (gmaxh - 1) * (1. - np.tanh(
        ((nh / vals['k_ehh'])**vals['p_ehh'])))

    # Cn can be considered temp independent for 70 - 400K
    Cn = vals['k_n']
    Cn *= g_eeh

    # Ch can not
    Ch = vals['k_h0'] + vals['k_h1'] * temp + vals['k_h2'] * temp**2
    Ch *= g_ehh

    # the auger recombination rate is given by
    R = Cn * (ne**2 * nh - ne0**2 * nh0) + \
        Ch * (ne * nh**2 - ne0 * nh0**2)

    # Then the lifetime is provided by this
    return nxc / R


def coulomb_enhanced_auger_Glunz(vals, nxc, ne0, nh0, **args):
    '''
    The coulomb enhanced auger model, as proposed by Glunz
    This uses equations 3 - 5 from 10.1063/1.370784.

    This work extend the work of Altermatt to include the
    impact of the excess carrier density (high injection effects)
    '''

    # then need to get doping and delta n
    nh = nh0 + nxc
    ne = ne0 + nxc

    gmaxn = vals['k_gmaxn'] * 300.**vals['p_gmaxn']
    gmaxh = vals['k_gmaxh'] * 300.**vals['p_gmaxh']

    # enhancement factors
    g_eeh = 1. + (gmaxn - 1) * (1. - np.tanh((
        (ne / vals['k_eeh'])**vals['n_eeh'])))
    g_ehh = 1. + (gmaxh - 1) * (1. - np.tanh(
        ((nh / vals['k_ehh'])**vals['p_ehh'])))

    Cn = g_eeh * vals['k_n'] * ne0 / (ne0 + nxc) + \
        vals['k_a'] / 2 * nxc / (ne0 + nxc)

    Ch = g_ehh * vals['k_h'] * nh0 / (nh0 + nxc) + \
        vals['k_a'] / 2 * nxc / (nh0 + nxc)

    # the auger recombination rate is given by
    R = Cn * (ne**2 * nh - ne0**2 * nh0) + \
        Ch * (ne * nh**2 - ne0 * nh0**2)

    # Then the lifetime is provided by this
    return nxc / R


def coulomb_enhanced_auger_Kerrsimple(vals, nxc, ne0, nh0, **args):
    '''
    The coulomb enhanced auger model, as proposed by Kerr2002
    This uses equation 20 from 10.1063/1.1432476

    This model was a demonstration of the model of suggested by Glunz
    This model underestimates Auger compared to earlier models
    '''
    nh = nh0 + nxc
    ne = ne0 + nxc

    R = (ne * nh - ne0 * nh0) *\
        (vals['k_n'] * ne0**vals['p_n'] +
         vals['k_h'] * nh0**vals['p_h'] +
         vals['k_nxc'] * nxc**vals['p_nxc']
         )

    # Then the lifetime is provided by this
    return nxc / R


def coulomb_enhanced_auger_Kerr(vals, nxc, ne0, nh0, **args):
    '''
    The coulomb enhanced auger model, as proposed by Kerr2002
    This uses equation 23 from 10.1063/1.1432476

    While more complicated it has better agreement with previous Auger models
    '''

    # then need to get doping and delta n
    nh = nh0 + nxc
    ne = ne0 + nxc

    g_eeh = 1. + vals['l_eeh'] * (1. - np.tanh(
        (ne0 / vals['k_eeh'])**vals['n_eeh']))
    g_ehh = 1. + vals['l_ehh'] * (1. - np.tanh(
        (nh0 / vals['k_ehh'])**vals['p_ehh']))

    # the auger recombination rate is given by
    R = (ne * nh - ne0 * nh0) *\
        (vals['k_n'] * g_eeh * ne0 +
         vals['k_h'] * g_ehh * nh0 +
         vals['k_nxc'] *
         nxc * (nxc + vals['n_nxc']) / (nxc + vals['d_nxc'])
         )

    # Then the lifetime is provided by this
    return nxc / R


def coulomb_enhanced_auger_Richter(vals, nxc, ne0, nh0, **args):
    '''
    The coulomb enhanced auger model, as proposed by Kerr2002
    '''

    # then need to get doping and delta n
    nh = nh0 + nxc
    ne = ne0 + nxc

    g_eeh = 1. + vals['l_eeh'] * (1. - np.tanh(
        (ne0 / vals['k_eeh'])**vals['n_eeh']))
    g_ehh = 1. + vals['l_ehh'] * (1. - np.tanh(
        (nh0 / vals['k_ehh'])**vals['p_ehh']))

    # the auger recombination rate is given by
    R = (ne * nh - ne0 * nh0) *\
        (vals['k_n'] * g_eeh * ne0 +
         vals['k_p'] * g_ehh *
         nh0 +
         vals['k_delta'] *
         nxc**vals['delta']
         )

    # Then the lifetime is provided by this
    return nxc / R


def auger_notimplimented(*args):
    """
    A dummy class for models with incomplete parameters
    """
    pass


def Altermatt1997_check(vals, func, fig, ax):
    """
    A function that checks the Altermatt's model
    This is incomplete
    """

    print 'Have not placed in check data for Altermatt\'s model yet'
    nxc = np.logspace(13, 20)
    for Nd in [1e15, 1e18, 1e20]:

        # get digitized data

        taus = func(vals, nxc, float(Nd), 0, 300)

        ax.plot(nxc, taus, '--', label=Nd)

    ax.set_xlabel('Excess carrer density (cm$^{-3}$)')
    ax.set_ylabel(r'$\tau$ (cm$^{-3}$)')
    ax.loglog()


def Glunz1999_check(vals, func, fig, ax):
    """
    There is no published data for Glunz mode,
    This function thus just plots Glunz model
    """

    nxc = np.logspace(13, 20)
    for Nd in [1e15, 1e18, 1e20]:

        # get digitized data

        taus = func(vals, nxc, float(Nd), 0)

        ax.plot(nxc, taus, '--', label=Nd)

    ax.set_xlabel('Excess carrer density (cm$^{-3}$)')
    ax.set_ylabel(r'$\tau$ (cm$^{-3}$)')
    ax.loglog()


def Kerr2002_simple_check(vals, func, fig, ax):
    """
    A function that checks the Kerr 2002 simple against published data
    """
    folder = os.path.join(os.path.dirname(__file__), 'Si', 'check_data')
    file = 'Kerr2002_simple.json'
    fpath = os.path.join(folder, file)
    # with open(fpath, 'r') as data_file:
    data = Webplotdig_JSONreader(fpath)

    for Nd in data.getDatasetNames():

        # get digitized data
        dataset = data.getDatasetByName(Nd)
        nxc, tau = data.getDatasetValues(dataset).T

        taus = func(vals, nxc, float(Nd), 0)

        p1, = ax.plot(nxc, tau, '.', label=Nd)
        color = p1.get_color()
        ax.plot(nxc, taus, '--', c=color, label=Nd)

    ax.set_xlabel('Excess carrer density (cm$^{-3}$)')
    ax.set_ylabel(r'$\tau$ (cm$^{-3}$)')
    ax.loglog()


def Kerr2002_complex_check(vals, func, fig, ax):
    """
    A function that checks the Kerr2002_complex model against published data
    """
    folder = os.path.join(os.path.dirname(__file__), 'Si', 'check_data')
    file = 'Kerr2002.json'
    fpath = os.path.join(folder, file)
    # with open(fpath, 'r') as data_file:
    data = Webplotdig_JSONreader(fpath)

    for Nd in data.getDatasetNames():

        # get digitized data
        dataset = data.getDatasetByName(Nd)
        nxc, tau = data.getDatasetValues(dataset).T

        taus = func(vals, nxc, float(Nd), 0)

        p1, = ax.plot(nxc, tau, '.', label=Nd)
        color = p1.get_color()
        ax.plot(nxc, taus, '--', c=color, label=Nd)

    ax.set_xlabel('Excess carrer density (cm$^{-3}$)')
    ax.set_ylabel(r'$\tau$ (cm$^{-3}$)')
    ax.loglog()


def Richter2012_check(vals, func, fig, ax):
    """
    A function that checks theRichter model against published data
    """

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

        p1, = ax.plot(nxc, taus)
        color = p1.get_color()

        ax.plot(data['nxc'], data['tau'], '--', c=color)
    ax.set_xlabel('Excess carrer density (cm$^{-3}$)')
    ax.set_ylabel(r'$\tau$ (cm$^{-3}$)')
    ax.loglog()
