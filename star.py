#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 17 11:17:20 2026

@author: idavis
"""
from astropy import units as un, constants as const
import numpy as np
from default_vals.star_vals import default_star_vals

class Star:
    def __init__(self, **kwargs):
        for k in list(default_star_vals.keys()):
            if k not in list(kwargs.keys()):
                kwargs.update({k:default_star_vals[k]})
        
        self.R_star = kwargs['R_star']
        self.M_star = kwargs['M_star']
        self.rotation_period = kwargs['rotation_period']
        self.T_phot = kwargs['T_phot']
        return