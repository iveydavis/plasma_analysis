#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 17 11:23:12 2026

@author: idavis
"""
from astropy import units as un, constants as const
import numpy as np
import matplotlib.pyplot as plt

e = const.e.gauss.value * un.cm**1.5 * un.g**0.5/un.s
G = (un.erg/un.cm**3)**0.5

plot_props = {}

def speed_to_rel_corr_energy(v):
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
    if type(temp) != un.quantity.Quantity:
        print("Assuming temperature is in Kelvin")
        temp *= un.K
    elif type(temp) == un.quantity.Quantity:
        assert(temp.unit.physical_type == 'temperature')
    v = np.sqrt(2 * temp * const.k_B/const.m_e)
    return v.to('cm/s')


def density_to_frequency(density):
    if type(density) != un.quantity.Quantity:
        print("No units provided for density; assuming in cubic centimeters")
        density *= un.cm**-3
    freq = 9.8 * un.kHz.to('MHz') * np.sqrt(density.to('cm**-3').value) * un.MHz
    return freq

def frequency_to_density(frequency):
    if type(frequency) != un.quantity.Quantity:
        print("No units provided, assuming in MHz")
        frequency *= un.MHz
    dens = ((frequency/(9.8*un.kHz)).to(''))**2
    return dens * un.cm**-3
