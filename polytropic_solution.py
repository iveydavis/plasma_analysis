#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 27 14:57:02 2026

@author: idavis
"""

from misc import un, const, np, plt
from default_vals.corona_vals import default_polytrope_vals 
# https://arxiv.org/pdf/2209.03508


class polytropic_wind:
    def __init__(self, **kwargs):
        for k in default_polytrope_vals:
            if k not in kwargs.keys():
                kwargs.update({k:default_polytrope_vals[k]})
        
        self.M = kwargs['M']
        self.R = kwargs['R']
        self.T0 = kwargs['T0']
        self.n0 = kwargs['n0']
        self.mass_fraction = kwargs['mass_fraction']
        self.rho0 = self.n0 * const.m_p * self.mass_fraction
        self.poly_idx = kwargs['poly_idx']
        self.calc_cg()
        self.calc_c0()
        assert((self.C0/self.Cg)**2 > (self.poly_idx-1)/(2*self.poly_idx)), "No accelerating solution exists"
        self.r_vec = np.linspace(self.R, kwargs['r_max'], kwargs['r_res'])
        

    def calc_cg(self):
        self.Cg =np.sqrt(const.G * self.M /2/self.R).cgs
        return
    
    
    def calc_c0(self):
        self.C0 = ((self.poly_idx * const.k_B * self.T0/const.m_p).cgs)**0.5
        return


    def calc_sc(self, starting_res = 5000, prec_thres=1e-3):
        gam = self.poly_idx
        res = starting_res
        end = 10
        diff = np.array([-1])
        
        Cg = self.Cg
        C0 = self.C0
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
        return

        
    def calc_v0(self):
        Cg = self.Cg
        C0 = self.C0
        gam = self.poly_idx
        sc = self.sc
        
        V02 = Cg**2 * (Cg/C0)**(4/(gam - 1)) * sc**((3*gam - 5)/(gam - 1))
        self.v0 = V02**0.5
        return 
    
    
    def calc_vc(self):
        self.vc = (self.Cg**2/self.sc)**0.5
    
    def calc_wind_prof(self, prec_thres=1e-3, prec_vel=10*un.km/un.s, max_iter=50):
        if 'sc' not in self.__dict__.keys():
            self.calc_sc()
        if 'vc' not in self.__dict__.keys():
            self.calc_vc()
            
        s_vec = self.r_vec/self.R
        vc = self.vc
        sc = self.sc
        gam = self.poly_idx
        self.v_prof = np.zeros(len(s_vec))*vc.unit
        vmax = 100
        vmin = 0
        for i,s in enumerate(s_vec):
            v_vec = np.linspace(vmin,vmax,5000) * un.km/un.s
            
            l = 0.5 * (v_vec**2 - vc**2)
            
            r1 = -vc**2/(gam -1) * ((vc*sc**2/(v_vec*s**2))**(gam-1) -1)
            r2 = 2 * self.Cg**2 * (1/s - 1/sc)
            
            diff = l - (r1 + r2)
            inf_idxs = np.where(np.isinf(diff))[0]
            diff[inf_idxs] = np.nan
            length = 5000
            diff0 = np.nanmin(np.abs(diff))
            niter = 0
            # print(i)
            while np.nanmin(np.abs(diff))**0.5 > prec_vel and niter < max_iter:
                # print("iter")
                v_vec = np.linspace(vmin, vmax, length) * un.km/un.s
                l = 0.5 * (v_vec**2 - vc**2)
                
                r1 = -vc**2/(gam -1) * ((vc*sc**2/(v_vec*s**2))**(gam-1) -1)
                r2 = 2 * self.Cg**2 * (1/s - 1/sc)
                
                diff = l - (r1 + r2)
                inf_idxs = np.where(np.isinf(diff))[0]
                diff[inf_idxs] = np.nan
                
                if np.nanmin(np.abs(diff)) >= 0.1*diff0:                    
                    vmax += 2*prec_vel.to('km/s').value
                    diff0 = np.nanmin(np.abs(diff))
                    length += 100
                else:
                    length += 100
                niter += 1

            idx = np.where(np.abs(diff) == np.nanmin(np.abs(diff)))[0]
            
            self.v_prof[i] = v_vec[idx]
            dv = (self.v_prof[i] - self.v_prof[i-1]).to('km/s').value
            if dv == 0:
                dv = prec_vel.to('km/s').value
            vmax = self.v_prof[i].to('km/s').value + 1.05*dv
            vmin = self.v_prof[i].to('km/s').value
            
        return
    
    def get_density_profile(self):
        if 'v_prof' not in self.__dict__.keys():
            self.calc_wind_prof()
        C = self.n0 * self.v_prof[0] * self.r_vec[0]**2
        n = C/(self.v_prof * self.r_vec**2)
        self.number_density_profile = n
        self.mass_density_profile = n * const.m_p * self.mass_fraction


