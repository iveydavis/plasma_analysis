#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pytest
from astropy import units as un

from swabs.misc import density_to_frequency, frequency_to_density
from swabs.star import Star


def test_star_uses_default_values():
    """A Star built with no arguments falls back to the EK Draconis defaults."""
    star = Star()

    assert star.R_star == (0.94 * un.R_sun).cgs
    assert star.M_star == (0.95 * un.M_sun).cgs
    assert star.T_phot == (5600 * un.K).cgs


def test_star_accepts_an_override():
    """A keyword argument should replace the default for that field only."""
    star = Star(R_star=2.0 * un.R_sun)

    assert star.R_star == (2.0 * un.R_sun).cgs
    # The fields that were not overridden still hold their defaults.
    assert star.M_star == (0.95 * un.M_sun).cgs


def test_density_and_frequency_round_trip():
    """density -> plasma frequency -> density should return the original value."""
    density = 1e8 * un.cm ** -3

    frequency = density_to_frequency(density)
    round_tripped = frequency_to_density(frequency)

    assert round_tripped.to("cm**-3").value == pytest.approx(density.value)
