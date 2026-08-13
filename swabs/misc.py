#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from astropy import units as un, constants as const
import numpy as np
import matplotlib.pyplot as plt
import os

e = const.e.gauss.value * un.cm**1.5 * un.g**0.5/un.s
G = (un.erg/un.cm**3)**0.5

plot_props = {}

def speed_to_rel_corr_energy(v):
    """
    Converts speed to relativistc-corrected energy in keV
    :param v: DESCRIPTION
    :type v: TYPE
    :raises Warning: DESCRIPTION
    :return: DESCRIPTION
    :rtype: TYPE

    """
    if type(v) != un.quantity.Quantity:
        print("Assuming velocity is in units of cm/s")
        v *= un.cm/un.s
    elif type(v) == un.quantity.Quantity:
        assert('speed' in v.unit.physical_type )
    
    if v > const.c.cgs:
        E = np.nan
        raise Warning("speed of electron higher than speed of light")
    else:
        E = (const.m_e * const.c**2 * ((1/(1 - (v/const.c)**2))**0.5 - 1)).to('keV')
    return E


def calc_thermal_electron_speed(temp):
    """
    Calculates thermal electron speed
    :param temp: thermal temperature
    :type temp: float or astropy quantity with units temperature
    :return: velocity in cm/s
    :rtype: astropy quantity

    """
    if type(temp) != un.quantity.Quantity:
        print("Assuming temperature is in Kelvin")
        temp *= un.K
    elif type(temp) == un.quantity.Quantity:
        assert(temp.unit.physical_type == 'temperature')
    v = np.sqrt(2 * temp * const.k_B/const.m_e)
    return v.to('cm/s')


def density_to_frequency(density):
    """
    Converts electron number density to fundamental plasma frequency
    :param density: electron number density
    :type density: float or astropy quantity with units inverse volume
    :return: frequency in MHz
    :rtype: astropy quantity

    """
    if type(density) != un.quantity.Quantity:
        print("No units provided for density; assuming in cubic centimeters")
        density *= un.cm**-3
    freq = 8.978 * un.kHz.to('MHz') * np.sqrt(density.to('cm**-3').value) * un.MHz
    return freq


def frequency_to_density(frequency):
    """
    Converts fundamental plasma frequency to electron number density
    :param frequency: plasma frequency
    :type frequency: float or astropy quantity with units frequency
    :return: plasma density in units of inverse cubic centimetre
    :rtype: astropy quantity

    """
    if type(frequency) != un.quantity.Quantity:
        print("No units provided, assuming in MHz")
        frequency *= un.MHz
    dens = ((frequency/(8.978*un.kHz)).to(''))**2
    return dens * un.cm**-3


def check_units(kwarg_dict, default_dict):
    """
    Checks that units of keyword arguments are correct
    :param kwarg_dict: a dictionary of keyword arguments
    :type kwarg_dict: dictionary
    :param default_dict: the dictionary to check the units against
    :type default_dict: dictionary
    :return: unit-corrected dictionary
    :rtype: dictionary

    """
    for k in kwarg_dict.keys():
        if type(kwarg_dict[k]) == un.quantity.Quantity:
            assert(kwarg_dict[k].unit.is_equivalent(default_dict[k].unit)), f"{kwarg_dict[k]} should be in units of {default_dict[k].unit.physical_type}"
            
        elif type(kwarg_dict[k]) != un.quantity.Quantity:
            if type(default_dict[k]) == un.quantity.Quantity:
                kwarg_dict[k] = kwarg_dict[k] * default_dict[k].unit
                
    return kwarg_dict
