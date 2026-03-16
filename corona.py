#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from misc import un, const, np, plt, calc_thermal_electron_speed
import os

class Corona:
    def __init__(self, process=False, **kwargs):
        self.r_vec = np.linspace(self.star.R_star + 0.0001*self.star.R_star, self.r_max, self.r_res)
        self.velocity_profile = np.zeros(len(self.r_vec))
        
        self.mass_density_profile = np.zeros(len(self.r_vec))
        self.number_density_profile = np.zeros(len(self.r_vec))
        self.temperature_profile = np.zeros(len(self.r_vec))
        self.c0 = self.calc_c0()
        return
    
    def calc_c0(self):
        """
        Calculates the base isothermal speed and updates the c0 property
        """
        return
    
    
    def calc_wind_solution(self):
        """
        Calculates the wind speed profile and then also determines the density, temperature, and magnetic
        field profiles
        """
        return
    
    def _get_other_wind_properties(self, velocities, find_open=True, alf_dist=None):
        """
        Calculates the properties of the wind besides the velocity, and updates
        the velocity profile based on the opened field condition
        :param velocities: velocity profile that informs the solutions
        :type velocities: ndarray of un.quantity.Quantity
        :param find_open: , defaults to True
        :param find_open: Describes whether to try to calculate the Alfven/open-field distance,, defaults to True
        :type find_open: bool, optional
        :param alf_dist: A value to force the Alfven/open-field distance to be, defaults to None
        :type alf_dist: un.quantity.Quantity, optional

        """
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
        return
    
    def update_max_distance(self, new_max_dist):
        """
        Increases the maximum distance that the coronal properties
        are calculated for assuming that the last quantity was in the open-field
        regime. Maintains the same spatial resolution as the original r_vec property
        :param new_max_dist: The new maximum distance to extend the wind properties 
        out to
        :type new_max_dist: un.quantity.Quantity
        """
        assert(new_max_dist > self.r_vec[-1]), "New maximum distance should be further than the original maximum distance"
        
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
        new_temp_prof = self.temperature_profile[0] * (new_num_prof[0]/new_num_prof)**(1-self.poly_idx)
        
        self.r_res = len(new_r_vec)
        self.r_max = new_max_dist.cgs
        self.r_vec = new_r_vec
        self.mag_field_profile = new_mag_prof
        self.velocity_profile = new_velocity_prof
        self.alfven_speed = new_alf_prof
        self.mass_density_profile = new_mass_prof
        self.number_density_profile = new_num_prof
        self.temperature_profile = new_temp_prof
        return
    
    
    def get_density(self, dist_from_surface, return_idxs:bool=False, predict:bool=True):
        """
        Returns the number density based on a distance from the stellar surface provided
        :param dist_from_surface: The distance from surface to retrieve the number density from
        :type dist_from_surface: un.quantity.Quantity
        :param return_idxs: describes whether to return the index of the vector that the density was retrieved from, defaults to False
        :type return_idxs: bool, optional
        :param predict: If the distance is further than the profiles have been
        calculated for, then prediction will try to estimat the density assuming
        the field is opened. Otherwise, it returns NaN, defaults to True
        :type predict: bool, optional
        :return: density, (idx)
        :rtype: un.quantity.Quantity, (int)
        """
        
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
        
    
    def calc_debye_lengths(self):
        """
        Calculates debye lengths as function of distance. Updates debye_lengths property
        """
        
        freqs = 9800 * un.Hz * 2 * np.pi * np.sqrt(self.number_density_profile.to('cm**-3').value)
        vTe = calc_thermal_electron_speed(self.temperature_profile)
        self.debye_lengths = (vTe/freqs).to('cm')
        return
    
    
    def plot(self, prop='velocity'):
        """
        Plots a given property as function of distance
        :param prop: Property or properties to be plotted, defaults to 'velocity'
        :type prop: str or list of str
        :return: figure and axis objects
        :rtype: Figure, AxesSubplot or array of AxesSubplot 
        """

        properties = ['velocity', 'mass_density', 'number_density', 'magnetic_field', 'alfven_speed', 'temperature']
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
            elif p.lower() == 'temperature':
                dat = self.temperature_profile
                ylabel = 'Temperature [K]'
                
            ax[i].semilogy((self.r_vec/self.star.R_star).to(''), dat)
            
            ax[i].set_ylabel(ylabel)
        ax[-1].set_xlabel(r"Distance [R$_\star$]")
        return fig, ax
    
    
    def save(self, outpath:str=None):
        """
        Saves the wind solution
        :param outpath: name of the path and file name to save data to, defaults to None
        :type outpath: str, optional
        """
        if outpath is None:
            dens = f"{int(self.n0.cgs.value/1e9)}e9"
            temp = f"{int(self.T0.cgs.value/1e6)}MK"
            field = f"{int(self.B0.cgs.value)}G"
            outpath = f"{os.getcwd()}/corona_n0{dens}_temp{temp}_field{field}.npz"
        np.savez(outpath, self.__dict__)    
        return

    
def load(file:str):
    """
    Loads corona file
    :param file: pathname of the file
    :type file: str
    :return: the corona data
    :rtype: Corona

    """
    dat = np.load(file, allow_pickle = True)
    d = dat['arr_0'].flatten()[0]
    if 'process' in d.keys():
        d['process'] = False
    elif 'process' not in d.keys():
        d.update({'process':False})
        
    cor = Corona(d)
    return cor