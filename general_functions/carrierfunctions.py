#!/usr/local/bin/python
# UTF-8

import numpy as np
from semiconductor.matterial.ni import IntrinsicCarrierDensity as NI


def get_carriers(Na, Nd, nxc,
                 temp=300,matterial='Si', ni_author=None,  ni=None):
    '''
    returns the carrier densities given the number of ionised dopants and ni
    and the excess carriers

    input:
    Na = the number of ionised acceptors
    Nd = the number of ionised donors
    nxc = the excess carrier density. In this function assume 
                  deltap = deltan
    temp = temperature
    ni: (optional)
        provide  a values so this function doesn't caculate ni

    returns ne, nh

    '''
    if ni is None:
        ni = NI(matterial=matterial).update(author=ni_author, temp=temp)
    
    # Calculated on the assumption that at thermal equilibrium in the 
    # dark n0p0 = ni**2, and that charge neutrality holds. Usually 
    # simplified to saying the majority carrier density ~ the doping and min
    # carrier denisty is the number of excess carriers. The below version 
    # more accurately incorporates ni though, which is particularly important 
    # for temperature dependent measurements.
    maj_car_den = (0.5 * (np.abs(Nd - Na) + np.sqrt((Nd - Na)**2 + 4 * ni**2
                  ))) + nxc 
    nxc = (ni**2 / maj_car_den) + nxc
    
    if np.all(Na < Nd):
        ne = maj_car_den        
        nh = nxc        
    elif np.all(Na >= Nd):
        nh = maj_car_den       
        ne = nxc
    else:
        print 'determination of total carrier connc didn\'t work'
    
    return ne, nh
        