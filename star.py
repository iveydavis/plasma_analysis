#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from misc import check_units
from default_vals.star_vals import default_star_vals


class Star:
    def __init__(self, **kwargs):
        """
        Defines characteristics of the star (not the corona). Defaults to values for EK Draconis
        :Keyword arguments:
            ** R_star: radius of the star, assumed in units of R_sun if not defined
            ** M_star: mass of the star, assumed in units of M_sun if not defined
            ** T_phot: photospheric temperature of the star, assumed in units Kelvin if not defined
        """
        for k in list(default_star_vals.keys()):
            if k not in list(kwargs.keys()):
                kwargs.update({k:default_star_vals[k]})
        for k in kwargs.keys():
            if k not in list(default_star_vals.keys()):
                raise Warning(f"{k} not recognised keyword argument")
        
        kwargs = check_units(kwargs, default_star_vals)
        
        self.R_star = kwargs['R_star'].cgs
        self.M_star = kwargs['M_star'].cgs
        self.T_phot = kwargs['T_phot'].cgs
        return
    