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
    "mean_mol_weight":1,
    "r_max":40*const.R_sun.cgs,
    "mass_fraction":1
    }

default_polytrope_vals = { 
    "n0":4e10/un.cm**3,
     "T0":1e7*un.K,
     "R":0.94*const.R_sun.cgs,
     "M":0.95*const.M_sun.cgs,
     "r_res":300,
     "r_max":30*const.R_sun.cgs,
     "mean_mol_weight":1,
     "poly_idx": 1.1,
     "mass_fraction":1
     }   
