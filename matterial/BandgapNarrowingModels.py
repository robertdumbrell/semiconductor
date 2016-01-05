
import numpy as np

def apparent_BGN(vals, doping, *args):
    '''
    It returns the 'apparent BGN'. This estimates the real bandgap narrowing, 
    but uses the accurate boltzman statastics. N is the net dopant concentration. 
        This model incorporates degeneracy as well as band gap narrowing,
    which is why it determines an 'apparent BGN' abd not an actual BGN.
    The apparent occurs as Boltzmann statistics are used rather than Fermi-Dirac satastic.
    '''

    if type(doping).__module__ != np.__name__:
        doping = np.array(doping)

    BGN = np.zeros(np.array(doping).shape)

    index = doping > vals['n_onset']
    BGN[index] = (
        vals['de_slope'] * np.log(doping[index] / vals['n_onset']))

    return BGN

def not_implimented(vals, doping, *args):
    print 'model not implimented, returning 0 values'
    return np.zeros(doping.shape)

def BGN(vals, doping, *args):
    '''
    Returns the BGN when applied for carriers with fermi distribution.
    This estimates the real bandgap narrowing, 
    but uses boltzman stats
    where N is the net dopant concentration. 
    '''
    # BGN = np.zeros(doping.shape)

    bgn = vals['de_slope']\
        * np.power(np.log(doping / vals['n_onset']), vals['b'])\
        + vals['de_offset']

    # ensures no negitive values
    if bgn.size > 1 :
        bgn[bgn < 0] = 0

    return bgn