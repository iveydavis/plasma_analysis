#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from misc import un, const, G
from star import Star
import copy

default_isothermal_vals = {
    "n0":4e10/un.cm**3,
    "T0":1e7*un.K,
    "star": Star(),
    "B0":1000*G,
    "r_res":300,
    "mean_mol_weight":1,
    "mass_fraction":1,
    "r_max":40*const.R_sun.cgs
    }

default_polytrope_vals = copy.deepcopy(default_isothermal_vals)
default_polytrope_vals.update({"poly_idx": 1.1})

