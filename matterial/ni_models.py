
import numpy as np
import scipy.constants as Const


def ni_temp(vals, temp, **kargs):
    """
     This form comes from Bludau, Onton, and
     Heinke3 and Macfarlane et a1.31 as cited by Green,3 is
     given by
    """
    ni = vals['a'] * (temp)**vals['power'] * \
        np.exp(- vals['eg'] / temp)
    # print self.ni, - vals['eg'], 2., Const.k, self.temp

    return ni


def ni_temp_eg(vals, temp, doping, Eg, *args):
    """
     This form comes from Bludau, Onton, and
     Heinke3 and Macfarlane et a1.31 as cited by Green,3 
    """

    ni = vals['a'] * temp**vals['power'] * \
        np.exp(- Eg / 2. * Const.e / Const.k / temp)
    # print vals, Eg, ni

    return ni
