#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from misc import un, const, np, plt, check_units
from default_vals.corona_vals import default_polytrope_vals 

class Polytropic_Wind:
    def __init__(self, **kwargs):
        """
        Class for solving polytropic wind equation
        :Keyword arguments:
            ** R_star: radius of the star, assumed in units of R_sun if not defined
            ** M_star: mass of the star, assumed in units of M_sun if not defined
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
        
        assert(type(self.r_res) == int)
        kwargs = check_units(kwargs, default_polytrope_vals)
        
        self.rho0 = self.n0 * const.m_p * self.mass_fraction
        self.calc_cg()
        self.calc_c0()
        assert((self.C0/self.Cg)**2 > (self.poly_idx-1)/(2*self.poly_idx)), "No accelerating solution exists"
        self.r_vec = np.linspace(self.R_star, self.r_max, self.r_res)
        
        
        for k in self.__dict__.keys():
            if type(self.__dict__[k]) == un.quantity.Quantity:
                self.__dict__[k] = self.__dict__[k].cgs
        return
        

    def calc_cg(self):
        """
        Calculates gravitational sound speed; updates property Cg
        """
        
        self.Cg =np.sqrt(const.G * self.M_star /2/self.R_star).cgs
        return
    
    
    def calc_c0(self):
        """
        Calculates initial thermal sound speed; updates property C0
        """
        
        self.C0 = ((self.poly_idx * const.k_B * self.T0/(const.m_p*self.mean_mol_weight)).cgs)**0.5
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
        """
        Calculates the velocity of the wind at the base of the corona, updates property v0
        """
        
        if 'sc' not in self.__dict__.keys():
            self.calc_sc()
            
        Cg = self.Cg
        C0 = self.C0
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
        self.vc = (self.Cg**2/self.sc)**0.5
        return
    
    
    def calc_wind_prof(self, starting_res:int=5000, prec_vel:un.quantity.Quantity=10*un.km/un.s, max_iter:int=50):
        """
        Calculates the wind speed profile, updates property v_prof
        :param starting_res: initial resolution in velocity space to solve for v at a given distance s
        :param prec_vel: precision threshold for velocity difference between left and right side of the equation, defaults to 10*un.km/un.s
        :type prec_vel: un.quantity.Quantity, optional
        :param max_iter: Number of times to increase resolution to solve within velocity precision, defaults to 50
        :type max_iter: int, optional
        """
        
        self.calc_v0()
        self.calc_vc()
            
        s_vec = self.r_vec/self.R_star
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
            length = starting_res
            diff0 = np.nanmin(np.abs(diff))
            niter = 0
            
            while np.nanmin(np.abs(diff))**0.5 > prec_vel and niter < max_iter:
                v_vec = np.linspace(vmin, vmax, length) * un.km/un.s
                l = 0.5 * (v_vec**2 - vc**2)
                
                r1 = -vc**2/(gam -1) * ((vc*sc**2/(v_vec*s**2))**(gam-1) -1)
                r2 = 2 * self.Cg**2 * (1/s - 1/sc)
                
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
                diffv = np.abs(v_vec - self.v_prof[i-1])
                idx = np.where(diffv == np.nanmin(diffv))[0][0]
            
            self.v_prof[i] = v_vec[idx]
            dv = (self.v_prof[i] - self.v_prof[i-1]).to('km/s').value
            if dv == 0 or dv > prec_vel.to('km/s').value:
                dv = prec_vel.to('km/s').value
            vmax = self.v_prof[i].to('km/s').value + 0.1*dv
            vmin = self.v_prof[i].to('km/s').value            
        return
    
    
    def get_density_profile(self):
        """
        Calculates the density profile based on the wind. Updates properties 
        number_density_profile and mass_density_profile
        """
        
        if 'v_prof' not in self.__dict__.keys():
            self.calc_wind_prof()
        C = self.n0 * self.v_prof[0] * self.r_vec[0]**2
        n = C/(self.v_prof * self.r_vec**2)
        self.number_density_profile = n
        self.mass_density_profile = n * const.m_p * self.mass_fraction
        
        return


    def get_temperature_profile(self):
        """
        Calculates the temperature profile. Updates property T
        """
        
        if "number_density_profile" not in self.__dict__.keys():
            self.get_density_profile()
            
        self.T = self.T0 * (self.n0/self.number_density_profile)**(1-self.poly_idx)
        
        return
    
    
    def plot(self, prop:str, list = 'velocity'):
        """
        Plots a given property as function of distance
        :param prop: Property or properties to be plotted, defaults to 'velocity'
        :type prop: str or list of str
        :return: figure and axis objects
        :rtype: Figure, AxesSubplot or array of AxesSubplot 
        """
        
        props = ['velocity', 'density', 'temperature']
        if type(prop) == str:
            prop = [prop]
        for p in prop:
            assert(p in props)
            
        fig, axs = plt.subplots(len(prop), sharex=True)
        if len(prop) == 1:
            axs = [axs]
            
        r = self.r_vec.cgs/self.R_star.cgs
        for i,p in enumerate(prop):
            if p == 'velocity':
                data = self.v_prof.to('km/s')
                ylabel = 'Veloicty [km/s]'
            if p == 'density':
                data = self.number_density_profile.cgs
                ylabel = r'Density [cm$^-$$^-3$]'
            if p == 'temperature':
                data = self.T
                ylabel = 'Temperature [K]'
            
            
            axs[i].plot(r, data)
            axs[i].set_ylabel(ylabel)
        axs[i].set_xlabel(r"Distance [R$_\star$]")

        return fig, axs
