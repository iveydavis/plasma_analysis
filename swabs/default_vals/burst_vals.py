#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from swabs.misc import const, un
from swabs.corona import Corona

default_iii_vals = {
    'vb': 0.5*const.c,
    'starting_height_factor': 0.1,
    'starting_width_factor': 0.1,
    'width_growth_factor': 0.5}

    
default_ii_vals = {
    'vb': 500*un.km/un.s,
    'starting_height_factor': 0.1,
    'starting_width_factor': 0.1,
    'width_growth_factor': 0.5}
