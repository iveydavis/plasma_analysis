#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 17 11:16:18 2026

@author: idavis
"""
from misc import un, const, np, plt, calc_thermal_electron_speed
from default_vals.corona_vals import default_corona_vals

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
        self.vthermal = np.sqrt(2 *const.k_B * self.temp/(const.m_e * self.mass_fraction) )
        if process:
            self.get_parker_solutions()
        return
    
    
    def get_parker_solutions(self, starting_resolution = 5000, max_iter=500, prec=1e-3, find_open=True, alf_dist=None):
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
            buff = int(self.r_res/1000) + 2
            i0 = failed_idx[0] - buff
            ie = failed_idx[-1] + buff
            if i0 < 1:
                i0 = 1
            
            dx = ie  - i0 
            dy = velocities[ie] - velocities[i0]
            slope = dy/dx
            velocities[i0:ie] = slope * np.linspace(1, dx, dx) + velocities[i0-1]
            
        density = self.n0 * self.mass_fraction*const.m_p*(self.r_vec[0]/self.r_vec)**2 * (velocities[0]/velocities)
        B_field = self.B0/((self.r_vec/self.star.R_star).to(''))**3
        alf_speed = (B_field/np.sqrt(density * 4 * np.pi)).to('km/s')   
        if find_open:
            try:
                if alf_dist is None:
                    alf_idx = np.where(alf_speed < velocities)[0][0]
                elif alf_dist is not None:
                    alf_idx = np.where(self.r_vec > alf_dist)[0][0]
                density[alf_idx:] = density[alf_idx] *  (self.r_vec[alf_idx]/self.r_vec[alf_idx:])**2
                velocities[alf_idx:] = velocities[alf_idx]
                B_field[alf_idx:] = B_field[alf_idx] *  (self.r_vec[alf_idx]/self.r_vec[alf_idx:])**2
                self.alf_dist = self.r_vec[alf_idx]
            except:
                print("No Alfven radius found")
            
        self.mag_field_profile = B_field
        self.mass_density_profile = density
        self.number_density_profile = self.mass_density_profile/(const.m_p * self.mass_fraction)
        self.velocity_profile = velocities
        self.alfven_speed = (self.mag_field_profile/np.sqrt(4 * np.pi * self.mass_density_profile)).to('km/s')
        
        
    def get_density(self, dist_from_surface, return_idxs=False, predict=True):
        r = self.r_vec - self.star.R_star
        if dist_from_surface > r[-1]:
            idx = []
            if not predict:
                densities = np.nan /un.cm**3    
            elif predict:
                dens0 = self.number_density_profile[-1]
                d0 = r[-1]
                densities = (dens0 * (d0/dist_from_surface)**2).cgs
                
        else:
            diff = np.abs(r - dist_from_surface)
            idx = np.where(diff == np.nanmin(diff))[0][0]
            densities = self.number_density_profile[idx]

        if not return_idxs:
            return np.nanmean(densities)
        elif return_idxs:
            return np.nanmean(densities), idx

    def update_max_distance(self, new_max_dist):
        dr = int((self.r_vec[1] - self.r_vec[0]).cgs.value)
        new_r_vec = np.arange(self.r_vec[0].cgs.value, new_max_dist.cgs.value + dr, dr)*un.cm
        
        new_mag_prof = np.zeros(len(new_r_vec)) * self.mag_field_profile.unit
        new_mass_prof = np.zeros(len(new_r_vec)) * self.mass_density_profile.unit
        new_num_prof = np.zeros(len(new_r_vec)) * self.number_density_profile.unit
        new_velocity_prof = np.zeros(len(new_r_vec)) * self.velocity_profile.unit
        
        new_mag_prof[:self.r_res] = self.mag_field_profile 
        new_mass_prof[:self.r_res] = self.mass_density_profile
        new_num_prof[:self.r_res] = self.number_density_profile
        new_velocity_prof[:self.r_res] = self.velocity_profile
        
        d0 = self.r_vec[-1].cgs
        new_mag_prof[self.r_res:] = new_mag_prof[self.r_res-1] * (d0/new_r_vec[self.r_res:])**2
        new_mass_prof[self.r_res:] = new_mass_prof[self.r_res-1] * (d0/new_r_vec[self.r_res:])**2
        new_num_prof[self.r_res:] = new_num_prof[self.r_res-1] * (d0/new_r_vec[self.r_res:])**2
        new_velocity_prof[self.r_res:] = self.velocity_profile[-1]
        new_alf_prof = (new_mag_prof/np.sqrt(4 * np.pi * new_mass_prof)).to('km/s')
        
        self.r_res = len(new_r_vec)
        self.r_max = new_max_dist.cgs
        self.r_vec = new_r_vec
        self.mag_field_profile = new_mag_prof
        self.velocity_profile = new_velocity_prof
        self.alfven_speed = new_alf_prof
        self.mass_density_profile = new_mass_prof
        self.number_density_profile = new_num_prof
        return
    
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