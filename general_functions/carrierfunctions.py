#!/usr/local/bin/python
# UTF-8

import numpy as np
from semiconductor.matterial.ni import IntrinsicCarrierDensity as NI


def get_carriers(Na, Nd, min_car_den,
                 temp=300,matterial='Si', ni_author=None,  ni=None):
    '''
    returns the carrier density given the doping and ni
    and the excess carriers

    input:
    Na
    Nd 
    min_car_den
    temp
    ni: (optional)
        provide  a values so this function doesn't caculate ni

    returns ne, nh

    '''
    if ni is None:
        ni = NI(matterial=matterial).update_ni(author=ni_author, temp=temp)
    ne, nh = Nd - Na, Na - Nd
    print ni
    if np.all(Na < Nd):
        ne += min_car_den
        nh = min_car_den + ni**2/ne
    elif np.all(Na > Nd):
        nh += min_car_den
        ne = min_car_den + ni**2/nh
    else:
        print 'determination of total carrier connc didn\'t work'

    return ne, nh
