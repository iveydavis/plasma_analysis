#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 12 11:10:40 2026

@author: idavis
"""
from astropy import units as un, constants as const
import numpy as np
import matplotlib.pyplot as plt

e = const.e.gauss.value * un.cm**1.5 * un.g**0.5/un.s
G = (un.erg/un.cm**3)**0.5


#####################################################################


def speed_to_rel_corr_energy(v):
    if type(v) != un.quantity.Quantity:
        print("Assuming velocity is in units of cm/s")
        v *= un.cm/un.s
    elif type(v) == un.quantity.Quantity:
        assert('speed' in v.unit.physical_type )
    
    if v > const.c.cgs:
        E = np.nan
        raise Warning("speed of electron higher than speed of light")
    else:
        E = (const.m_e * const.c**2 * ((1/(1 - (v/const.c)**2))**0.5 - 1)).to('keV')
    return E


def calc_thermal_electron_speed(temp):
    if type(temp) != un.quantity.Quantity:
        print("Assuming temperature is in Kelvin")
        temp *= un.K
    elif type(temp) == un.quantity.Quantity:
        assert(temp.unit.physical_type == 'temperature')
    v = np.sqrt(2 * temp * const.k_B/const.m_e)
    return v.to('cm/s')


def density_to_frequency(density):
    if type(density) != un.quantity.Quantity:
        print("No units provided for density; assuming in cubic centimeters")
        density *= un.cm**-3
    freq = 9.8 * un.kHz.to('MHz') * np.sqrt(density.to('cm**-3').value) * un.MHz
    return freq


#####################################################################


default_star_vals = {
    "R_star":0.94*const.R_sun.cgs,
    "M_star":0.95*const.M_sun.cgs,
    "T_phot": 5600*un.K,
    "rotation_period":10*un.day}

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


#####################################################################


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

class Corona:
    def __init__(self, process=False, **kwargs):
        for k in list(default_corona_vals.keys()):
            if k not in list(kwargs.keys()):
                kwargs.update({k:default_corona_vals[k]})
        self.n0 = kwargs['n0']
        self.temp = kwargs['temp']
        self.star = kwargs['star']
        self.B0 = kwargs['B0']
        self.mass_fraction = kwargs['mass_fraction']
        self.r_res = kwargs['r_res']
        self.r_max = kwargs['r_max']
        
        self.r_vec = np.linspace(self.star.R_star + 0.0001*self.star.R_star, self.r_max, self.r_res)
        self.wind_profile = np.zeros(len(self.r_vec))
        self.dterms = np.zeros(len(self.r_vec))
        self.right_terms = np.zeros(len(self.r_vec))
        self.mass_density_profile = np.zeros(len(self.r_vec))
        self.number_density_profile = np.zeros(len(self.r_vec))
        
        self.c0 = np.sqrt(const.k_B * self.temp/ (const.m_p * self.mass_fraction)).to('km/s')
        self.rc = (const.G * self.star.M_star/2/self.c0**2).to('R_sun')
        self.vthermal = np.sqrt(2 *const.k_B * self.temp/(const.m_p * self.mass_fraction) )
        if process:
            self.get_parker_solutions()
        return
    
    
    def get_parker_solutions(self, starting_resolution = 5000, max_iter=500, prec=1e-3):
        velocities = np.zeros(len(self.r_vec)) * un.km/un.s
        
        for i,r in enumerate(self.r_vec):
            right_term = 4 * np.log(r/self.rc) + 4 * self.rc/r - 3
            
            if i == 0 or np.isnan(velocities[i-1]):
                v_vec = np.linspace(0, self.c0.to('km/s').value, 5000) * un.km/un.s
            else:
                v_vec = np.linspace(velocities[i-1].to("km/s").value, velocities[i-1].to('km/s').value + 100,5000) * un.km/un.s
                
            left_term = (v_vec/self.c0)**2 - np.log(v_vec**2/self.c0**2)    
            
            count = 0
            res = starting_resolution
            while np.nanmax(left_term) < right_term or np.nanmin(left_term) > right_term: 
                v_vec = np.linspace(v_vec[0].to("km/s").value, v_vec[-1].to("km/s").value + 1, starting_resolution) * un.km/un.s
                if v_vec[0] <= 0:
                    v_vec[0] = 0
                left_term = (v_vec/self.c0)**2 - np.log(v_vec**2/self.c0**2)
            
            while np.nanmin(np.abs(left_term - right_term)) > prec and count <= max_iter:
                res += 1000
                if i == 0:
                    vmin = v_vec[0].to('km/s').value
                elif i != 0:
                    vmin = velocities[i-1].to('km/s').value
                    
                v_vec = np.linspace(vmin, v_vec[-1].to("km/s").value + 1, res) * un.km/un.s    
                left_term = (v_vec/self.c0)**2 - np.log(v_vec**2/self.c0**2)
                count += 1
                
            self.right_terms[i] = right_term
            dterm = np.abs(right_term - left_term)
            
            diff_idxs = np.where(dterm <= 5*prec)[0]
            if len(diff_idxs) > 1:
                vels = v_vec[diff_idxs]
                if i == 0:
                    min_idx = diff_idxs[0]
                elif i != 0:
                    vdiff = vels - velocities[-1]
                    bad_idxs = np.where(vdiff < 0)[0]
                    vdiff[bad_idxs] = np.nan
                    diff_idx = np.where(vdiff >= 0)[0]
                    min_idx = diff_idxs[diff_idx[0]]
                    
                v = v_vec[min_idx]
            elif len(diff_idxs) == 1:
                min_idx = diff_idxs
                v = v_vec[min_idx]
            elif len(diff_idxs) == 0:
                min_idx = np.nanmin(dterm)
                v = np.nan * un.km/un.s
                
            self.dterms[i] = dterm[min_idx]
            velocities[i] = v
        
        diff = velocities - np.roll(velocities, 1)
        failed_idx = np.where(diff == 0)[0]
        if len(failed_idx) > 0:
            i0 = failed_idx[0]
            ie = failed_idx[-1]
            buff = int(self.r_res/1000) + 2
            
            dx = ie  - i0 + 2*buff
            dy = velocities[ie + buff] - velocities[i0 - buff]
            slope = dy/dx
            velocities[i0 - buff: ie+buff] = slope * np.linspace(1, dx, dx) + velocities[i0 - buff-1]
            
        density = self.n0 * self.mass_fraction*const.m_p*(self.r_vec[0]/self.r_vec)**2 * (velocities[0]/velocities)
        B_field = self.B0/((self.r_vec/self.star.R_star).to(''))**3
        alf_speed = (B_field/np.sqrt(density * 4 * np.pi)).to('km/s')        
        try:
            alf_idx = np.where(alf_speed < velocities)[0][0]
            density[alf_idx:] = density[alf_idx] *  (self.r_vec[alf_idx]/self.r_vec[alf_idx:])**2
            velocities[alf_idx:] = velocities[alf_idx]
            B_field[alf_idx:] = B_field[alf_idx] *  (self.r_vec[alf_idx]/self.r_vec[alf_idx:])**2
        except:
            print("No Alfven radius found")
            
        self.mag_field_profile = B_field
        self.mass_density_profile = density
        self.number_density_profile = self.mass_density_profile/(const.m_p * self.mass_fraction)
        self.velocity_profile = velocities
        self.alfven_speed = (self.mag_field_profile/np.sqrt(4 * np.pi * self.mass_density_profile)).to('km/s')
        
        
    def get_density(self, dist_from_surface, return_idxs=False):
        r = self.r_vec - self.star.R_star
        if dist_from_surface > self.r_vec[-1]:
            densities = np.nan /un.cm**3
            idx = []
        else:
            diff = np.abs(r - dist_from_surface)
            idx = np.where(diff == np.nanmin(diff))[0][0]
            densities = self.number_density_profile[idx]

        if not return_idxs:
            return np.nanmean(densities)
        elif return_idxs:
            return np.nanmean(densities), idx

    
    def calc_debye_lengths(self):
        freqs = 9800 * un.Hz * 2 * np.pi * np.sqrt(self.number_density_profile.to('cm**-3').value)
        vTe = calc_thermal_electron_speed(self.temp)
        self.debye_lengths = (vTe/freqs).to('cm')
        return
    
    
    def plot(self, prop='velocity'):
        properties = ['velocity', 'mass_density', 'number_density', 'magnetic_field', 'alfven_speed']
        if type(prop) != list:
            prop = [prop]
                
        if type(prop) == list:
            for p in prop:
                p = p.lower()
                assert(p in properties), f"{p} not recognized property: {properties}"
        
        fig, ax = plt.subplots(len(prop), sharex=True)
        if len(prop) == 1:
            ax = [ax]
        for i, p in enumerate(prop):
            if p.lower() == 'velocity':
                dat = self.velocity_profile.to('km/s')
                ylabel = "Velocity [km/s]"
            elif p.lower() == 'mass_density':
                dat = self.mass_density_profile.to('g/cm**3')
                ylabel = r'Mass density [g/cm$^3$]'
            elif p.lower() == 'number_density':
                dat = self.number_density_profile
                ylabel = r'Number density [cm$^-$$^3$]'
            elif p.lower() == 'magnetic_field':
                dat = self.mag_field_profile 
                ylabel = 'Magnetic field [G]'
            elif p.lower() == 'alfven_speed':
                dat = self.alfven_speed
                ylabel = 'Alfven velocity [km/s]'
                
            ax[i].semilogy((self.r_vec/self.star.R_star).to(''), dat)
            
            ax[i].set_ylabel(ylabel)
        ax[-1].set_xlabel(r"Distance [R$_\star$]")
   
    
#####################################################################


default_iii_vals = {'corona': Corona(process=True),
                    'vb': 0.5*const.c.cgs,
                    'starting_height_factor': 0.1,
                    'spec_idx': 3}

class Type_III_Burst:
    def __init__(self, **kwargs):
        for k in list(default_iii_vals.keys()):
            if k not in list(kwargs.keys()):
                kwargs.update({k:default_iii_vals[k]}) 
        self.corona = kwargs['corona']
        self.vb = kwargs['vb']
        self.spec_idx = kwargs['spec_idx']
        self.starting_height = kwargs['starting_height_factor'] * self.corona.star.R_star
        return
    
    
    def get_frequency_profile(self, duration=10*un.min, t_res:int = 1000, propagation='radial', do_return=False):
        propagation = propagation.lower()
        propagations = ['radial', 'spiral']
        assert(propagation in propagations), f"{propagation} not a recognized propagation profile ({propagations})"
        
        locs = locals()
        for lk in locs.keys():
            self.__dict__.update({lk:locs[lk]})
            
        tvec = np.linspace(0*un.s, duration, t_res)
        freqs = np.zeros(t_res)*un.MHz
        dists = np.zeros(t_res) * un.cm
        if propagation == 'radial' or 'spiral':
            for i, t in enumerate(tvec):
                dist = t * self.vb + self.starting_height
                density = self.corona.get_density(dist)
                freqs[i] = density_to_frequency(density)
                dists[i] = dist
                
        self.frequencies = freqs
        self.times = tvec
        self.distances = dists
        if do_return:
            return tvec, freqs
        
        
    def make_dynamic_spectrum(self, duration=10*un.min, t_res = 1*un.s, nu_min=30*un.MHz, nu_max=168*un.MHz, nu_res=1*un.MHz, duration_factor=1, duration_exp=0.95):
        nu_min = nu_min.to('MHz').value
        nu_max = nu_max.to('MHz').value
        tvec = np.arange(0,duration.to('s').value + 1, t_res.to('s').value) * un.s
        nuvec = np.arange(nu_min, nu_max+1, nu_res.to('MHz').value)
        
        subband_extents = [[nuvec[0], nuvec[1]]]
        i = 1
        while i < len(nuvec) - 1:
            subband_extents.append([nuvec[i], nuvec[i] + 1])
            i += 1
        sb = np.array(subband_extents)
        
        dyn_spec = np.zeros((len(tvec), len(nuvec))) 
        for i, t in enumerate(tvec):
            dist = self.starting_height + self.vb * t
            density, idx = self.corona.get_density(dist, return_idxs=True)
            nu = density_to_frequency(density).to('MHz').value
            
            if nu_min <= nu <= nu_max:
                diff = np.transpose(sb - nu)
                idx = np.where((diff[0] <= 0) & (diff[1] >= 0))[0][0]
                
                t_dec= 10**7.7 * (nu * 1e6)**(-duration_exp) * duration_factor
                dt = np.round(t_dec/t_res.to('s').value).astype(int)
                if t + dt*un.s > tvec[-1]:
                    idx_t = -1
                else:
                    idx_t = i + dt
                dyn_spec[i:idx_t,idx] = 1
                if np.nansum(dyn_spec[i-1:i, idx+1]) == 0:
                    try:
                        freq_min_idx = np.where(dyn_spec[i-1,:] == 1)[0][0]
                    except:
                        freq_min_idx = -1
                    dyn_spec[i:idx_t,idx:freq_min_idx] = 1
        idxs = np.where(dyn_spec == 0)
        dyn_spec[idxs] = np.nan
        self.dyn_spec = dyn_spec.transpose()
        self.extents = [0, duration.to('min').value, nu_min, nu_max]
        return
    
    def plot_dyn_spec(self, cmap='magma', interpolation='nearest'):
        assert('dyn_spec' in self.__dict__.keys()), "No dyn spec has been made to plot"
        
        fig, ax = plt.subplots()
        ax.imshow(self.dyn_spec, aspect='auto', origin='lower', extent=self.extents, interpolation=interpolation, cmap=cmap)
        ax.tick_params(axis = 'both', which='both', labelsize=15)
        ax.set_ylim(self.extents[2], self.extents[3])
        ax.set_xlim(self.extents[0], self.extents[1])
        ax.set_ylabel("Frequency [MHz]", fontsize = 18, multialignment='center')
        ax.set_xlabel("Time since flare onset [min]", fontsize = 18)
        fig.tight_layout()
        return


#####################################################################

    
default_ii_vals = {
    'corona': Corona(process=True),
    'vb': 500*un.km/un.s,
    'starting_height_factor': 0.4,
    'starting_width_factor': 0.25,
    'width_growth_factor': 1}

class Type_II_Burst:
    def __init__(self, **kwargs):
        for k in list(default_ii_vals.keys()):
            if k not in list(kwargs.keys()):
                kwargs.update({k:default_ii_vals[k]}) 
        
        assert(kwargs['starting_height_factor'] > kwargs['starting_width_factor']/2)
        self.corona = kwargs['corona']
        self.vb = kwargs['vb']
        self.starting_height = kwargs['starting_height_factor'] * self.corona.star.R_star
        self.starting_width = kwargs['starting_width_factor']  * self.corona.star.R_star
        self.width_growth_factor = kwargs['width_growth_factor']
        return
    
    def get_frequency_profile(self, duration=30*un.min, t_res:int = 1000, propagation='radial', do_return=False):
        propagation = propagation.lower()
        propagations = ['radial', 'spiral']
        assert(propagation in propagations), f"{propagation} not a recognized propagation profile ({propagations})"
        
        locs = locals()
        for lk in locs.keys():
            self.__dict__.update({lk:locs[lk]})
        
        tvec = np.linspace(0*un.s, duration, t_res)
        freqs = np.zeros(t_res)*un.MHz
        dists = np.zeros(t_res) * un.cm
        if propagation == 'radial' or 'spiral':
            for i, t in enumerate(tvec):
                dist = t * self.vb + self.starting_height
                density = self.corona.get_density(dist)
                freqs[i] = density_to_frequency(density)
                dists[i] = dist
                
        self.frequencies = freqs
        self.times = tvec
        self.distances = dists
        if do_return:
            return tvec, freqs
        

    def make_dynamic_spectrum(self, t_res=1*un.s, nu_res=1*un.MHz, duration=30*un.min, nu_min = 30*un.MHz, nu_max = 168*un.MHz):
        locs = locals()
        for lk in locs.keys():
            self.__dict__.update({lk:locs[lk]})
        
        nu_min = nu_min.to('MHz').value
        nu_max = nu_max.to('MHz').value
        tvec = np.arange(0,duration.to('s').value + 1, t_res.to('s').value) * un.s
        nuvec = np.arange(nu_min, nu_max+1, nu_res.to('MHz').value)
        
        subband_extents = [[nuvec[0], nuvec[1]]]
        i = 1
        while i < len(nuvec) - 1:
            subband_extents.append([nuvec[i], nuvec[i] + 1])
            i += 1
        sb = np.array(subband_extents)
        
        dist_vec = self.starting_height + self.vb * tvec
        width_vec = self.starting_width + self.width_growth_factor * (self.starting_width * (dist_vec/dist_vec[0] -1))
        dyn_spec = np.zeros((len(tvec), len(nuvec)))*np.nan
        
        velocities = self.corona.velocity_profile
        v_A = self.corona.alfven_speed
        for i, t in enumerate(tvec):
            d = dist_vec[i]
            w = width_vec[i]
            d_min = d - w/2
            d_max = d + w/2
            
            dens_min, d_max_idx = self.corona.get_density(d_max, return_idxs=True)
            dens_max, d_min_idx = self.corona.get_density(d_min, return_idxs=True)
            
            freq_min = density_to_frequency(dens_min).to('MHz').value
            freq_max = density_to_frequency(dens_max).to('MHz').value
            
            if (nu_min <= freq_max <= nu_max) or (nu_min <= freq_min <= nu_max):
                if freq_max <= nu_max:
                    diff = np.abs(sb.transpose()[1] - freq_max)
                    sb_idx1 = np.where(diff == np.nanmin(diff))[0][0]
                elif freq_max > nu_max:
                    sb_idx1 = len(subband_extents)
                
                if freq_min >= nu_min:
                    diff = np.abs(freq_min - sb.transpose()[0])
                    sb_idx0 = np.where(diff == np.nanmin(diff))[0][0]
                elif freq_min < nu_min:
                    sb_idx0 = 0
                
                if sb_idx1 == sb_idx0:
                    dyn_spec[i,sb_idx0] = -1
                else:
                    dyn_spec[i,sb_idx0:sb_idx1] = -1
                    
            if freq_max >= nu_max and freq_min <= nu_min:
                dyn_spec[i,0:-1] = -1

            super_sonic_idxs = np.where((velocities+v_A)[d_min_idx:d_max_idx] < self.vb)[0]
            if len(super_sonic_idxs) != 0:
                dist_subset = self.corona.r_vec[d_min_idx:d_max_idx][super_sonic_idxs] - self.corona.star.R_star
  
                dens_max, d_min_idx = self.corona.get_density(dist_subset[0], return_idxs=True)
                dens_min, d_max_idx = self.corona.get_density(dist_subset[-1], return_idxs=True)
                
                freq_max = density_to_frequency(dens_max).to('MHz').value
                freq_min = density_to_frequency(dens_min).to('MHz').value
                
                if (nu_min <= freq_max <= nu_max) or (nu_min <= freq_min <= nu_max):
                    if freq_max <= nu_max:
                        diff = np.abs(sb.transpose()[1] - freq_max)
                        sb_idx1 = np.where(diff == np.nanmin(diff))[0][0]
                    elif freq_max > nu_max:
                        sb_idx1 = len(subband_extents)
                    
                    if freq_min >= nu_min:
                        diff = np.abs(freq_min - sb.transpose()[0])
                        sb_idx0 = np.where(diff == np.nanmin(diff))[0][0]
                    elif freq_min < nu_min:
                        sb_idx0 = 0
                    
                    if sb_idx1 == sb_idx0:
                        dyn_spec[i,sb_idx0] = 1
                    else:
                        dyn_spec[i,sb_idx0:sb_idx1] = 1
                        
                if freq_max >= nu_max and freq_min <= nu_min:
                    dyn_spec[i,0:-1] = 1
        
        self.dyn_spec = dyn_spec.transpose()
        self.extents = [0, duration.to('min').value, nu_min, nu_max]
        return
    
    
    def plot_dyn_spec(self, interpolation='nearest', cmap='magma'):
        assert('dyn_spec' in self.__dict__.keys()), "No dyn spec has been made to plot"
        
        fig, ax = plt.subplots()
        ax.imshow(self.dyn_spec, aspect='auto', origin='lower', extent=self.extents, interpolation=interpolation, cmap=cmap)
        ax.tick_params(axis = 'both', which='both', labelsize=15)
        ax.set_ylim(self.extents[2], self.extents[3])
        ax.set_xlim(self.extents[0], self.extents[1])
        ax.set_ylabel("Frequency [MHz]", fontsize = 18, multialignment='center')
        ax.set_xlabel("Time since flare onset [min]", fontsize = 18)
        fig.tight_layout()
        return fig, ax