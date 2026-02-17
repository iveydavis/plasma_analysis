#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 12 11:10:40 2026

@author: idavis
"""

import matplotlib.pyplot as plt
from misc import un, const, G, density_to_frequency, np
from default_vals.burst_vals import default_ii_vals, default_iii_vals


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