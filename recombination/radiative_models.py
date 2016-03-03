
import numpy as np




def Roosbroeck(vals, min_car_den, Nh_0, Ne_0,  temp = None, B=None):

    # B = Roosbroeck_B(vals, B)

    Nh = Nh_0 + min_car_den
    Ne = Ne_0 + min_car_den

    R = B * (Ne * Nh - Ne_0 * Nh_0)

    return min_car_den / R


def Roosbroeck_with_screening_B(vals, min_car_den, doping, temp, Blow):
    """
    This is the roosbroeck model that accounts for many things
    It needs temperature, min_car_den, doping and blow to be defined
    """

    bmin = vals['rmax'] + (vals['rmin'] - vals['rmax']) / (
        1. + (temp / vals['r1'])**vals['r2'])
    b1 = (vals['smax'] + (vals['smin'] - vals['smax']) / (
        1. + (temp / vals['s1'])**vals['s2'])) * 2
    b3 = (vals['wmax'] + (vals['wmin'] - vals['wmax']
                          ) / (
        1. + (temp / vals['w1'])**vals['w2'])) * 2

    # print bmin

    B = Blow * (bmin + (vals['bmax'] - bmin) / (
        1. + ((2. * min_car_den + doping) / b1
              )**vals['b2']
        + ((2. * min_car_den + doping) / b3)**vals['b4']))

    return B


def Roosbroeck_with_screening(vals, min_car_den, Nh_0, Ne_0, temp, B):
    """
    This is the roosbroeck model that accounts for many things
    It needs temperature, min_car_den, doping and blow to be defined
    """

    Roosbroeck_with_screening_B(vals, min_car_den, np.amax(Nh_0, Ne_0), temp, B)

    return Roosbroeck(min_car_den, Nh_0, Ne_0, B)


def cubic_loglog_parm(vals, temp):

    if np.any(temp > vals['temp_limit']):
        # print temp
        print 'Out of Blow model temperature range'
        print 'the ', vals['temp_limit'], 'is provided'

    B = vals['a'] * np.log10(temp)**3 +\
        vals['b'] * np.log10(temp)**2 \
        + vals['c'] * np.log10(temp)**1 \
        + vals['d']
    return 10**B


