#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 17 11:21:07 2026

@author: idavis
"""
from misc import un, const, G
from star import Star

default_corona_vals = {
    "n0":4e10/un.cm**3,
    "temp":1e7*un.K,
    "star": Star(),
    "R_star":0.94*const.R_sun.cgs,
    "M_star":0.95*const.M_sun.cgs,
    "B0":1000*G,
    "r_res":300,
    "mass_fraction":1,
    "r_max":40*const.R_sun.cgs,
    }   
