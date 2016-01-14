#!/usr/local/bin/python
# UTF-8

import numpy as np
from semiconductor.matterial.ni import IntrinsicCarrierDensity as NI


def get_carriers(Na, Nd, dn, temp=300, ni_model=None, matterial='Si', ni=None):
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
        ni = NI(matterial=matterial).update_ni(model=ni_model, temp=temp)
    ne, nh = Nd - Na, Na - Nd

    if np.all(Na < Nd):
        ne += dn
        nh = dn + ni
    elif np.all(Na > Nd):
        nh += dn
        ne = dn + ni

    return ne, nh
