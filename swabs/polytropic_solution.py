#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from swabs.misc import un, const, np, plt, check_units
from swabs.default_vals.corona_vals import default_polytrope_vals 
from swabs.corona import Corona

class Polytropic_Corona(Corona):
    def __init__(self, **kwargs):
        """
        Class for solving polytropic wind equation
        :Keyword arguments:
            ** star: instance of Star class with relevant information
            ** B0: base magnetic field strength
            ** T0: base coronal temperature, assumed in units Kelvin if not defined
            ** n0: base coronal number density, assumed in units per cubic cm if not defined
            ** r_res (int): number of of elements in the spatial array
            ** r_max: maximum physical distance in the spatial array, assumed in units of R_sun if not defined
            ** mean_mol_weight: mean molecular weight of gas, default is 0.6
            ** mass_fraction: isotopic mass fraction of the gas, default is 1 (hydrogen)
            ** poly_idx: polytropic index of the gas
        """
        
        for k in default_polytrope_vals:
            if k not in kwargs.keys():
                kwargs.update({k:default_polytrope_vals[k]})
                
        for k in list(kwargs.keys()):
            if k not in list(default_polytrope_vals.keys()):
                raise Warning(f"{k} not recognized keyword argument")
                
        for k in kwargs:
            self.__dict__.update({k:kwargs[k]})
        kwargs = check_units(kwargs, default_polytrope_vals)
        super().__init__(kwargs)
        self.calc_cg()
        assert((self.c0/self.cg)**2 > (self.poly_idx-1)/(2*self.poly_idx)), "No accelerating solution exists"

        return
    
    def calc_c0(self):
        return ((self.poly_idx * const.k_B * self.T0/(const.m_p*self.mean_mol_weight)).cgs)**0.5
        
    
    def calc_cg(self):
        """
        Calculates gravitational sound speed; updates property cg
        """
        self.cg =np.sqrt(const.G * self.star.M_star /2/self.star.R_star).cgs
        return
    
    
    def calc_sc(self, starting_res:int=5000, prec_thres=1e-3):
        """
        Calculates the normalised critical distance and updates property sc
        :param starting_res: the initial resolution in sc space to solve within, defaults to 5000
        :type starting_res: int, optional
        :param prec_thres: Precision threshold for the difference between the left and right hand sides of the equation, defaults to 1e-3
        :type prec_thres: float, int, optional
        """
        
        gam = self.poly_idx
        res = starting_res
        end = 10
        diff = np.array([-1])
        
        Cg = self.cg
        C0 = self.c0
        while all(diff < 0) or all(diff > 0):
            sc_vals = np.linspace(prec_thres, end, res)
            
            logl11 = 4/(gam -1) * np.log10(Cg/C0) + 2*(2*gam -3)/(gam -1)*np.log10(sc_vals)
            linl11 = 10**logl11
            l11 = 0.5 * (linl11 -1)
            
            l2 = 1/(gam - 1) * (C0**2/Cg**2 * sc_vals -1)
            r = 2 * (sc_vals-1)
            
            diff = (r - (l11 + l2)).value
            end += 10
            inf_idxs = np.where(np.isinf(diff))[0]
            diff[inf_idxs] = np.nan
        
        while np.nanmin(np.abs(diff)) > prec_thres:
            res += 100
            sc_vals = np.linspace(sc_vals[0], sc_vals[-1], res)
            logl11 = 4/(gam -1) * np.log10(Cg/C0) + 2*(2*gam -3)/(gam -1)*np.log10(sc_vals)
            linl11 = 10**logl11
            l1 = 0.5 * (linl11 -1)
            l2 = 1/(gam - 1) * (C0**2/Cg**2 * sc_vals -1)
            l = l1 + l2
            r = 2 * (sc_vals-1)
            
            diff = (r - l).value
            inf_idxs = np.where(np.isinf(diff))[0]
            diff[inf_idxs] = np.nan
        
        idx = np.where(np.abs(diff) == np.nanmin(np.abs(diff)))[0]
        self.sc = sc_vals[idx][0]
        if self.sc < 1:
            self.sc = 1
        return

        
    def calc_v0(self):
        """
        Calculates the velocity of the wind at the base of the corona, updates property v0
        """
        
        if 'sc' not in self.__dict__.keys():
            self.calc_sc()
            
        Cg = self.cg
        C0 = self.c0
        gam = self.poly_idx
        sc = self.sc
        
        V02 = Cg**2 * (Cg/C0)**(4/(gam - 1)) * sc**((3*gam - 5)/(gam - 1))
        self.v0 = V02**0.5
        return 
    
    
    def calc_vc(self):
        """
        Calculates the velocity of the wind at the critical distance sc, updates property vc
        """
        
        if 'sc' not in self.__dict__.keys():
            self.calc_sc()
        
        self.vc = (self.cg**2/self.sc)**0.5
        return
    
    
    def calc_wind_solution(self, starting_res:int=5000, prec_vel:un.quantity.Quantity=10*un.km/un.s, max_iter:int=50, find_open=True, alf_dist=None):
        """
        Calculates the wind speed profile, updates property v_prof
        :param starting_res: initial resolution in velocity space to solve for v at a given distance s
        :param prec_vel: precision threshold for velocity difference between left and right side of the equation, defaults to 10*un.km/un.s
        :type prec_vel: un.quantity.Quantity, optional
        :param max_iter: Number of times to increase resolution to solve within velocity precision, defaults to 50
        :type max_iter: int, optionals
        """
        
        self.calc_v0()
        self.calc_vc()
            
        s_vec = self.r_vec/self.star.R_star
        vc = self.vc
        sc = self.sc
        if sc < 1:
            sc = 1
        gam = self.poly_idx
        self.velocity_profile = np.zeros(len(s_vec))*vc.unit
        vmax = 100
        vmin = 0
        for i,s in enumerate(s_vec):
            v_vec = np.linspace(vmin,vmax,5000) * un.km/un.s
            
            l = 0.5 * (v_vec**2 - vc**2)
            
            r1 = -vc**2/(gam -1) * ((vc*sc**2/(v_vec*s**2))**(gam-1) -1)
            r2 = 2 * self.cg**2 * (1/s - 1/sc)
            
            diff = l - (r1 + r2)
            inf_idxs = np.where(np.isinf(diff))[0]
            diff[inf_idxs] = np.nan
            length = starting_res
            diff0 = np.nanmin(np.abs(diff))
            niter = 0
            
            while np.nanmin(np.abs(diff))**0.5 > prec_vel and niter < max_iter:
                v_vec = np.linspace(vmin, vmax, length) * un.km/un.s
                l = 0.5 * (v_vec**2 - vc**2)
                
                r1 = -vc**2/(gam -1) * ((vc*sc**2/(v_vec*s**2))**(gam-1) -1)
                r2 = 2 * self.cg**2 * (1/s - 1/sc)
                
                diff = l - (r1 + r2)
                inf_idxs = np.where(np.isinf(diff))[0]
                diff[inf_idxs] = np.nan
                
                if np.nanmin(np.abs(diff)) >= 0.1*diff0:                    
                    vmax += prec_vel.to('km/s').value
                    diff0 = np.nanmin(np.abs(diff))
                    length += 100
                else:
                    length += 100
                niter += 1

            idx = np.where(np.abs(diff) == np.nanmin(np.abs(diff)))[0]
            if len(idx) > 1:
                diffv = np.abs(v_vec - self.velocity_profile[i-1])
                idx = np.where(diffv == np.nanmin(diffv))[0][0]
            
            self.velocity_profile[i] = v_vec[idx]
            dv = (self.velocity_profile[i] - self.velocity_profile[i-1]).to('km/s').value
            if dv == 0 or dv > prec_vel.to('km/s').value:
                dv = prec_vel.to('km/s').value
            vmax = self.velocity_profile[i].to('km/s').value + 0.1*dv
            vmin = self.velocity_profile[i].to('km/s').value     
        super()._get_other_wind_properties(self.velocity_profile, find_open=find_open, alf_dist=alf_dist)
        self.temperature_profile = self.T0 * (self.n0/self.number_density_profile)**(1-self.poly_idx)
        return
    
    
    def get_density(self, dist_from_surface, return_idxs:bool=False, predict:bool=True):
        sol = super().get_density(dist_from_surface, return_idxs, predict)
        return sol

    
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

