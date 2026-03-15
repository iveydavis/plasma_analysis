#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from misc import un, const, G
from star import Star

default_corona_vals = {
    "n0":4e10/un.cm**3,
    "temp":1e7*un.K,
    "star": Star(),
    "B0":1000*G,
    "r_res":300,
    "mean_mol_weight":1,
    "mass_fraction":1,
    "r_max":40*const.R_sun.cgs
    }

default_polytrope_vals = { 
    "n0":4e10/un.cm**3,
     "T0":1e7*un.K,
     "R_star":0.94*un.R_sun,
     "M_star":0.95*un.M_sun,
     "r_res":300,
     "r_max":30*un.R_sun,
     "mean_mol_weight":1,
     "mass_fraction":1,
     "poly_idx": 1.1
     }   
