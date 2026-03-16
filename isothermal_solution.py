#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from corona import Corona
from default_vals.corona_vals import default_isothermal_vals
from misc import const, un, check_units, np

class Isothermal_Corona(Corona):
    def __init__(self, process=False, **kwargs):
        """
        Class for solving for winds and estimating coronal properties
        :Keyword arguments:
            ** star: Star class instance with properties of the star to model the corona for
            ** temp: base coronal temperature, assumed in units Kelvin if not defined
            ** n0: base coronal number density, assumed in units per cubic cm if not defined
            ** r_res (int): number of of elements in the spatial array
            ** r_max: maximum physical distance in the spatial array, assumed in units of R_sun if not defined
            ** mean_mol_weight: mean molecular weight of gas, default is 0.6
            ** mass_fraction: isotopic mass fraction of the gas, default is 1 (hydrogen)
        """
        for k in list(default_isothermal_vals.keys()):
            if k not in list(kwargs.keys()):
                kwargs.update({k:default_isothermal_vals[k]})
        for k in list(kwargs.keys()):
            if k not in list(default_isothermal_vals.keys()):
                raise Warning(f"{k} not recognized keyword argument")
                
        for k in kwargs:
            self.__dict__.update({k:kwargs[k]})
        
        check_units(kwargs, default_isothermal_vals)
        
        super().__init__(kwargs)
        
        self.rc = (const.G * self.star.M_star/2/self.c0**2).to('R_sun')
        if process:
            self.calc_wind_solutions(sol_type='parker')
        return


    def calc_c0(self):
        return np.sqrt(const.k_B * self.T0/ (const.m_p * self.mean_mol_weight)).to('km/s')
    
    
    def calc_wind_solution(self, starting_resolution:int=5000, max_iter:int=500, prec=1e-3, find_open:bool=True, alf_dist=None):
        """
        Gets the solution for an isothermal/Parker wind
        :param starting_resolution: starting velocity-space resolution to search for solution at given distance, defaults to 5000
        :type starting_resolution: int, optional
        :param max_iter: Maximum number of iterations to try increasing resolution to reach precision, defaults to 500
        :type max_iter: int, optional
        :param prec: precision of difference between left and right hand side of the equation, defaults to 1e-3
        :type prec: int, float, optional
        :param find_open: Describes whether to try to calculate the Alfven/open-field distance,, defaults to True
        :type find_open: bool, optional
        :param alf_dist: A value to force the Alfven/open-field distance to be, defaults to None
        :type alf_dist: un.quantity.Quantity, optional
        """
        
        self.poly_idx = 1
        velocities = np.zeros(len(self.r_vec)) * un.km/un.s
        dterms = np.zeros(len(self.r_vec))
        right_terms = np.zeros(len(self.r_vec))
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
                
            right_terms[i] = right_term
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
                
            dterms[i] = dterm[min_idx]
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
        
        super()._get_other_wind_properties(velocities, find_open=find_open, alf_dist=alf_dist)
        self.temperature_profile  = np.zeros(len(self.r_vec)) + self.T0
        return
        
    
    def get_density(self, dist_from_surface, return_idxs:bool=False, predict:bool=True):
        super().get_density(dist_from_surface, return_idxs, predict)
        return

    
    def update_max_distance(self, new_max_dist):
        super().update_max_distance(new_max_dist)
        return
    
    
    def calc_debye_lengths(self):
        super().calc_debye_lengths()
        return
    
    
    def plot(self, prop='velocity'):
        super().plot(prop)
        return
    
        
    def save(self, outpath):
        super().save(outpath)