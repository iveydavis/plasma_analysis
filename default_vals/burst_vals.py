#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 17 11:25:48 2026

@author: idavis
"""
from misc import const, un
from corona import Corona

default_iii_vals = {'corona': Corona(process=True),
                    'vb': 0.5*const.c.cgs,
                    'starting_height_factor': 0.1,
                    'spec_idx': 3}

    
default_ii_vals = {
    'corona': Corona(process=True),
    'vb': 500*un.km/un.s,
    'starting_height_factor': 0.4,
    'starting_width_factor': 0.25,
    'width_growth_factor': 1}