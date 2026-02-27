#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 12 11:10:40 2026

@author: idavis
"""

from misc import un, density_to_frequency, np, plt, frequency_to_density
from default_vals.burst_vals import default_ii_vals, default_iii_vals

C_sun = 4.22e36 * un.Hz
class Burst:
    def __init__(self, **kwargs):
        self.corona = kwargs['corona']
        self.vb = kwargs['vb']
        self.starting_height = kwargs['starting_height_factor'] * self.corona.star.R_star
        return
    
    def get_frequency_profile(self):
        return
    
    def make_dynamic_spectrum(self):
        return
    
    def plot_dyn_spec(self, cmap='magma', interpolation='nearest'):
        
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
    
    
class Type_III_Burst(Burst):
    def __init__(self, **kwargs):
        for k in list(default_iii_vals.keys()):
            if k not in list(kwargs.keys()):
                kwargs.update({k:default_iii_vals[k]}) 
        super().__init__(corona = kwargs['corona'], vb = kwargs['vb'], starting_height_factor=kwargs['starting_height_factor'])
        return 
     
        
    def get_frequency_profile(self, duration=10*un.min, t_res:int = 1*un.s, do_return=False):
        locs = locals()
        for lk in locs.keys():
            self.__dict__.update({lk:locs[lk]})
            
        tvec =  np.arange(0,duration.to('s').value + 1, t_res.to('s').value) * un.s
        freqs = np.zeros(len(tvec))*un.MHz
        dists = np.zeros(len(tvec)) * un.cm
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
            
        dyn_spec = np.zeros((len(tvec), len(nuvec)))
        densities = self.corona.number_density_profile
        distances = self.corona.r_vec
        velocities = self.corona.velocity_profile
        
        for i, nu in enumerate(nuvec):
            density = frequency_to_density(nu * un.MHz)
            diff = np.abs(densities - density)
            idx = np.where(diff == np.nanmin(diff))[0]

            t = self.corona.r_vec[idx]/self.vb
            t_diff = np.abs(tvec - t)
            t_idx = np.where(t_diff == np.nanmin(t_diff))[0][0]
            
            if t_idx == len(tvec) -1:
                continue
            
            if idx < len(densities) - 1:
                nu_c = (9.8 * un.kHz * np.sqrt((C_sun/(velocities[idx] * distances[idx]**2)).to('cm**-3').value)).to('Hz').value
                t_dec = (10**7.71 /nu_c**0.95)[0]#*un.s
                dt = np.round(t_dec/t_res.to('s').value).astype(int)
                te_idx = t_idx + dt
        
                if te_idx > len(tvec):
                    te_idx = len(tvec)
                dyn_spec[t_idx:te_idx, i] = 1

        idxs = np.where(dyn_spec == 0)
        dyn_spec[idxs] = np.nan
        self.dyn_spec = dyn_spec.transpose()
        self.extents = [0, duration.to('min').value, nu_min, nu_max]
        return
    
    
    def plot_dyn_spec(self, cmap='magma', interpolation='nearest'):
        assert('dyn_spec' in self.__dict__.keys()), "No dyn spec has been made to plot"
        super().plot_dyn_spec(cmap=cmap, interpolation=interpolation)
        return
#####################################################################


class Type_II_Burst(Burst):
    def __init__(self, **kwargs):
        for k in list(default_ii_vals.keys()):
            if k not in list(kwargs.keys()):
                kwargs.update({k:default_ii_vals[k]}) 
        
        assert(kwargs['starting_height_factor'] > kwargs['starting_width_factor']/2)
        super().__init__(corona = kwargs['corona'], vb = kwargs['vb'], starting_height_factor=kwargs['starting_height_factor'])
        self.starting_width = kwargs['starting_width_factor']  * self.corona.star.R_star
        self.width_growth_factor = kwargs['width_growth_factor']
        return
    
       
    def get_frequency_profile(self, duration=30*un.min, t_res=1*un.s, propagation='radial', do_return=False):
        propagation = propagation.lower()
        propagations = ['radial', 'spiral']
        assert(propagation in propagations), f"{propagation} not a recognized propagation profile ({propagations})"
        
        locs = locals()
        for lk in locs.keys():
            self.__dict__.update({lk:locs[lk]})
        
        tvec =  np.arange(0,duration.to('s').value + 1, t_res.to('s').value) * un.s
        freqs = np.zeros(len(tvec))*un.MHz
        dists = np.zeros(len(tvec)) * un.cm
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
        
    
    def get_subband_idxs(self, freq_min, freq_max=None):
        nu_min = self.subbands[0][0]
        nu_max = self.subbands[-1][-1]
        if (nu_min <= freq_max <= nu_max) or (nu_min <= freq_min <= nu_max):
            if freq_max <= nu_max:
                diff = np.abs(self.subbands.transpose()[1] - freq_max)
                sb_idx1 = np.where(diff == np.nanmin(diff))[0][0]
            elif freq_max > nu_max:
                sb_idx1 = len(self.subbands)
            
            if freq_min >= nu_min:
                diff = np.abs(freq_min - self.subbands.transpose()[0])
                sb_idx0 = np.where(diff == np.nanmin(diff))[0][0]
            elif freq_min < nu_min:
                sb_idx0 = 0
                
        if freq_max >= nu_max and freq_min <= nu_min:
            sb_idx0 = 0
            sb_idx1 = -1
            
        if freq_min > nu_max or freq_max < nu_min:
            sb_idx0 = None
            sb_idx1 = None
        
        return sb_idx0, sb_idx1
    
    
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
        while i < len(nuvec)-1:
            subband_extents.append([nuvec[i], nuvec[i + 1]])
            i += 1
        sb = np.array(subband_extents)
        self.subbands = sb
        
        dist_vec = self.starting_height + self.vb * tvec
        width_vec = self.starting_width*(1  + self.width_growth_factor * (dist_vec/dist_vec[0] -1))
        v_tot = self.corona.velocity_profile + self.corona.alfven_speed
        
        if dist_vec[-1] + width_vec[-1]/2 > self.corona.r_vec[-1].cgs - self.corona.star.R_star.cgs:
            max_dist = self.corona.star.R_star + dist_vec[-1] + width_vec[-1]/2 
            self.corona.update_max_distance(max_dist)
        
        dyn_spec = np.zeros((len(tvec), len(sb)))#*np.nan   
        
        numins = []
        numaxs = []
        for i, t in enumerate(tvec):
            d = dist_vec[i]
            w = width_vec[i]
            d_min = d - w/2
            d_max = d + w/2
            
            dens_min, d_max_idx = self.corona.get_density(d_max, return_idxs=True)
            dens_max, d_min_idx = self.corona.get_density(d_min, return_idxs=True)
            
            try:
                d_max_idx = int(d_max_idx)
            except:
                d_max_idx = -1
            
            freq_min = density_to_frequency(dens_min).to('MHz').value
            freq_max = density_to_frequency(dens_max).to('MHz').value
            
            sb_idx0, sb_idx1 = self.get_subband_idxs(freq_min, freq_max)
            
            if sb_idx0 != None:
                if sb_idx1 == sb_idx0:
                    dyn_spec[i,sb_idx0] = -1
                else:
                    dyn_spec[i,sb_idx0:sb_idx1+1] = -1
            
            numins.append(freq_min)
            numaxs.append(freq_max)
            super_sonic_idxs = np.where(v_tot[d_min_idx:d_max_idx] < self.vb)[0]
            if len(super_sonic_idxs) != 0:
                dist_subset = self.corona.r_vec[d_min_idx:d_max_idx][super_sonic_idxs] - self.corona.star.R_star
  
                dens_max = self.corona.get_density(dist_subset[0], return_idxs=False)
                dens_min = self.corona.get_density(dist_subset[-1], return_idxs=False)
                
                freq_max = density_to_frequency(dens_max).to('MHz').value
                freq_min = density_to_frequency(dens_min).to('MHz').value                
                idx0, idx1 = self.get_subband_idxs(freq_min, freq_max)
                
                if idx0 != None:
                    if idx0 == idx1:
                        dyn_spec[i, idx0] = 1
                    else:
                        dyn_spec[i, idx0:idx1+1] = 1
                    
        self.min_freqs = np.array(numins) * un.MHz
        self.max_freqs = np.array(numaxs) * un.MHz        
        self.dyn_spec = dyn_spec.transpose()
        self.extents = [0, duration.to('min').value, nu_min, nu_max]
        self.tvec = tvec
        self.nuvec = nuvec
        return
    
    def plot_dyn_spec(self, cmap='magma', interpolation='nearest'):
        assert('dyn_spec' in self.__dict__.keys()), "No dyn spec has been made to plot"
        super().plot_dyn_spec(cmap=cmap, interpolation=interpolation)
        return
