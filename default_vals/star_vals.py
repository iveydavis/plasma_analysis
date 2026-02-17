#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 17 11:18:33 2026

@author: idavis
"""
from astropy import units as un, constants as const

default_star_vals = {
    "R_star":0.94*const.R_sun.cgs,
    "M_star":0.95*const.M_sun.cgs,
    "T_phot": 5600*un.K,
    "rotation_period":10*un.day}
