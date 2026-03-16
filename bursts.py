#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from misc import un, const, density_to_frequency, np, plt, os
from corona import Corona
from default_vals.burst_vals import default_ii_vals, default_iii_vals
import copy
from scipy.optimize import curve_fit


class Burst:
    def __init__(self, corona:Corona, burst_type:str, **kwargs):
        """
        Class for estimating type II or III burst structure in dynamic spectra 
        and estimating related properties
        :param burst_type: type of burst, either 'ii' or 'iii'
        :type burst_type: str
        :Keyword Arguments:
            corona: Corona class instance that contains solutions for wind properties
            vb: speed of the beam. If burst_type is 'ii', the default is 500km/s;
                if it is 'iii', the default is 0.5c
            starting_height_factor: fraction of the stellar radius above the 
                stellar surface that the burst starts at
            starting_width_factor: fraction of the stellar radius that the beam
                length starts as
            width_growth_factor: the rate of growth of the beam.

        """
        
        btype = burst_type.lower()
        assert(btype in ['ii', 'iii'] )
        if btype == 'ii':
            default_vals = default_ii_vals
        elif btype == 'iii':
            default_vals = default_iii_vals
            
        for k in list(default_vals.keys()):
            if k not in list(kwargs.keys()):
                kwargs.update({k:default_vals[k]}) 
        for k in list(kwargs.keys()):
            if k not in list(default_vals.keys()):
                raise Warning(f"{k} not recognized keyword argument")
                
        for k in kwargs:
            self.__dict__.update({k:kwargs[k]})
        self.corona = corona
        self.starting_height = self.starting_height_factor * self.corona.star.R_star
        self.starting_width = self.starting_width_factor  * self.corona.star.R_star
        self.burst_type = btype
        self.set_vmin()
        return
    
    
    def set_vmin(self, vmin = None):
        """
        Defines the minimum velocity that the beam must exceed to generate a burst.
        For a type II burst, this is the wind-frame Alfven velocity. For type III
        burst it is 3 times the background thermal velocity.
        Updates the vmin property.
        :param vmin: Allows for manual definition of the minimum velocity. If None,
            the minimum velocity is calculated based on the type of burst, defaults to None
        :type vmin: array of un.quantity.Quantity, optional
        """
        if vmin is None:
            if self.burst_type == 'ii':
                self.vmin = self.corona.velocity_profile + self.corona.alfven_speed
            elif self.burst_type == 'iii':
                self.vmin = np.sqrt(2 * const.k_B * self.corona.temperature_profile/ const.m_e) * 3
                
        elif vmin is not None:
            assert(len(vmin) == len(self.corona.velocity_profile))
            self.vmin = vmin
        return
    
    
    def make_frequency_profile(self, duration:un.quantity.Quantity=30*un.min, 
                               t_res:un.quantity.Quantity=1*un.s, 
                               do_return:bool=False):
        """
        Makes a profile of frequency excited by the center of the beam as a function of time.
        Updates porperties profile_freqs, profile_times, and profile_distances
        :param duration: Length of time to calculate the profile for, defaults to 30*un.min
        :type duration: un.quantity.Quantity, optional
        :param t_res: time resolution, defaults to 1*un.s
        :type t_res: un.quantity.Quantity, optional
        :param do_return: option to return the values instead of just updating properties, defaults to False
        :type do_return: bool, optional
        :return: if do_return, returns the time and frequency arrays for the profile
        :rtype: (ndarray, ndarray)
        """
        
        tvec =  np.arange(0,duration.to('s').value + 1, t_res.to('s').value) * un.s
        freqs = np.zeros(len(tvec))*un.MHz
        dists = np.zeros(len(tvec)) * un.cm
        for i, t in enumerate(tvec):
            dist = t * self.vb + self.starting_height
            density = self.corona.get_density(dist)
            freqs[i] = density_to_frequency(density)
            dists[i] = dist
                
        self.profile_freqs = freqs
        self.profile_times = tvec
        self.profile_distances = dists
        if do_return:
            return tvec, freqs
    
    
    def get_subband_idxs(self, freq_min, freq_max):
        """
        Gets the frequency indices in a dynamic spectrum for a range of frequencies
        :param freq_min: minimum frequency of the frequency range to identify indices for
        :type freq_min: un.quantity.Quantity
        :param freq_max: maximum frequency of the range to identify indices for
        :type freq_max: un.quantity.Quantity, optional
        :return: starting and ending indices of the dynamic spectrum that satisfy
            the provided frequency range
        :rtype: (int, int)
        """
        
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
    
    
    def make_dynamic_spectrum(self, t_res:un.quantity.Quantity=1*un.s, 
                              nu_res:un.quantity.Quantity=1*un.MHz, 
                              duration:un.quantity.Quantity=30*un.min, 
                              nu_min:un.quantity.Quantity=30*un.MHz, 
                              nu_max:un.quantity.Quantity=168*un.MHz):
        """
        Creates a dynamic spectrum of the burst. Updates properties min_freqs,
        max_freqs, dyn_spec, extents, dynspec_times, dynspec_freqs
        :param t_res: Time resolution of the dynamic spectrum, defaults to 1*un.s
        :type t_res: un.quantity.Quantity, optional
        :param nu_res: frequency resolution of the dynamic spectrum, defaults to 1*un.MHz
        :type nu_res: un.quantity.Quantity, optional
        :param duration: duration of the dynamic spectrum, defaults to 30*un.min
        :type duration: un.quantity.Quantity, optional
        :param nu_min: Minimum frequency of the dynamic spectrum, defaults to 30*un.MHz
        :type nu_min: un.quantity.Quantity, optional
        :param nu_max: Maximum frequency of the dynamic spectrum, defaults to 168*un.MHz
        :type nu_max: un.quantity.Quantity, optional

        """
        
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
        
        if dist_vec[-1] + width_vec[-1]/2 > self.corona.r_vec[-1].cgs - self.corona.star.R_star.cgs:
            max_dist = self.corona.star.R_star + dist_vec[-1] + width_vec[-1]/2 
            self.corona.update_max_distance(max_dist)
            self.set_vmin()
        
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
            super_sonic_idxs = np.where(self.vmin[d_min_idx:d_max_idx] < self.vb)[0]
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
        self.extents = [0*un.min, duration.to('min'), nu_min*un.MHz, nu_max*un.MHz]
        self.dynspec_times = tvec
        self.dynspec_freqs = nuvec*un.MHz
        return
    
    
    def fit_drift_rate(self, freq_range:tuple=None, return_profile:bool=False):
        """
        Fits the drift rate to estimate de-dispersion parameters. Updates
        drift_rate_params
        :param freq_range: The frequency range to do the fit over. If None, then 
            it does it over the entire frequency range from make_frequency_profile, 
            defaults to None
        :type freq_range: tuple, optional
        :param return_profile: If True, then it returns the frequencies and the 
            profile produced from the fit, defaults to False
        :type return_profile: bool, optional
        :return: If return_profile is True, returns the frequencies and the fit profile
        :rtype: tuple of ndarray
        """
        
        assert('profile_freqs' in self.__dict__.keys()), "Need to run make_frequency_profile() first"
        
        freqs = self.profile_freqs
        times = self.profile_times
        
        dnudt = (freqs - np.roll(freqs, 1))/(times[1] - times[0])
        dnudt = dnudt[1:]
        freqs = freqs[1:]
        
        if freq_range is not None:
            try:
                idxs = np.where((freqs <= np.nanmax(freq_range)) & (freqs <= np.nanmin(freq_range)))[0]
                freqs = freqs[idxs]
                dnudt = dnudt[idxs]
            except Exception as e:
                raise Exception(e)
                
        fit, __ = curve_fit(drift_v_freq, freqs, dnudt)
        self.drift_fit_params = fit
        
        if return_profile:
            return freqs, drift_v_freq(freqs, fit[0], fit[1])
        return
        
        
    def fit_duration(self, freq_range:tuple=None, return_profile:bool=False):
        """
        Fits the burst duration as a functino of frequency. Updates
        duration_fit_params
        :param freq_range: The frequency range to do the fit over. If None, then 
            it does it over the entire frequency range from make_frequency_profile, 
            defaults to None
        :type freq_range: tuple, optional
        :param return_profile: If True, then it returns the frequencies and the 
            profile produced from the fit, defaults to False, defaults to False
        :type return_profile: bool, optional
        :return: If return_profile is True, returns the frequencies and the fit profile
        :rtype: tuple of ndarray
        """
        
        assert('dyn_spec' in self.__dict__.keys()), "Need to run make_frequency_profile() first"
        ds = copy.deepcopy(self.dyn_spec)
        dt = self.dynspec_times[1] - self.dynspec_times[0]
        if freq_range is not None:
            try:
                idx_min, idx_max = self.get_subband_idxs(np.nanmin(freq_range), np.nanmax(freq_range))
            except Exception as e:
                raise Exception(e)
        elif freq_range is None:
            idx_min = 0
            idx_max = -1
        
        ds[ds < 0] = 1
        duration = np.nansum(ds, axis = 1) * dt.to('s').value
        
        duration = duration[idx_min:idx_max]
        freqs = self.dynspec_freqs[idx_min:idx_max]
        if len(freqs) > len(duration):
            freqs = freqs[:-1]
        fit, __ = curve_fit(duration_v_freq, freqs, duration)
        self.duration_fit_params = fit
        
        if return_profile:
            return freqs, duration_v_freq(freqs, fit[0], fit[1])
        return
    
    
    def dedisperse(self, a=None, alpha=None):
        """
        Dedisperses the dynamic spectrum based on the equation:
            delta_t = 1/(a * (alpha -1 ))* (nu^(1-alpha) - nu0^(1-alpha))
        :param a: defaults to None
        :type a: float, int
        :param alpha: defaults to None
        :type alpha: float, int
        :return: the de-dispersed dynamic spectrum
        :rtype: 2-dimensional ndarray.

        """
        if a is None or alpha is None:
            self.fit_dispersion_properties()
        if a is None:
            a = self.drift_fit_params[0]
        if alpha is None:
            alpha = self.drift_fit_params[1]
        
        freqs_flip = np.flip(self.dynspec_freqs[:-1])
        ds = np.flip(copy.deepcopy(self.dyn_spec), axis = 0)
        
        dt = (self.dynspec_times[1] - self.dynspec_times[0]).to('s').value
        exp = 1 - alpha
        nu0 = freqs_flip[0]
        deltat_vals = []
        for j, nu in enumerate(freqs_flip):
            deltat = 1/(a*(alpha -1)) * (nu**exp - nu0**exp)
            dn = int(deltat/dt)
            ds[j] = np.roll(ds[j], -dn)
            if dn != 0:
                ds[j,-dn:] = 0
            deltat_vals.append(deltat)
            
        ds = np.flip(ds, axis = 0)
        return ds
    
    
    def plot_dyn_spec(self, cmap:str='magma', interpolation:str='nearest', time_unit:str='min', freq_unit:str='MHz'):
        """
        Plots the dynamic spectrum saved in property dyn_spec
        :param cmap: Colomap for pyplot.imshow, defaults to 'magma'
        :type cmap: str, optional
        :param interpolation: Interpolation method for pyplot.imshow, defaults to 'nearest'
        :type interpolation: str, optional
        :param time_unit: Unit of the time axis, defaults to 'min'
        :type time_unit: str, optional
        :param freq_unit: Unit of the frequency axis, defaults to 'MHz'
        :type freq_unit: str, optional
        :return: figure and axis objects of the dynamic spectrum
        :rtype: Figure, AxesSubplot or array of AxesSubplot 

        """
        fig, ax = plt.subplots()
        extents = copy.deepcopy(self.extents)
        extents[0] = extents[0].to(time_unit).value
        extents[1] = extents[1].to(time_unit).value
        
        extents[2] = extents[2].to(freq_unit).value
        extents[3] = extents[3].to(freq_unit).value
        
        ax.imshow(self.dyn_spec, aspect='auto', origin='lower', extent=extents, interpolation=interpolation, cmap=cmap)
        ax.tick_params(axis = 'both', which='both', labelsize=15)
        ax.set_ylim(extents[2], extents[3])
        ax.set_xlim(extents[0], extents[1])
        ax.set_ylabel(f"Frequency [{freq_unit}]", fontsize = 18, multialignment='center')
        ax.set_xlabel(f"Time since flare onset [{time_unit}]", fontsize = 18)
        fig.tight_layout()
        return fig, ax
    
    
    def save(self, outpath):
        """
        Saves the burst solution
        :param outpath: name of the path and file name to save data to, defaults to None
        :type outpath: str, optional
        """
        if outpath is None:
            if self.burst_type == 'iii':
                vb = f"{np.round((self.vb.cgs/const.c.cgs),2)}c"
            elif self.burst_type == 'ii':
                vb = f"{int(self.vb.to('km/s')).value}"
            outpath = f"{os.getcwd()}/type{self.burst_type}_burst_vb_{vb}.npz"
        np.savez(outpath, self.__dict__)
        return 
    

def load(file):
    """
    Loads burst file
    :param file: pathname of the file
    :type file: str
    :return: the burst data
    :rtype: Burst

    """
    dat = np.load(file, allow_pickle = True)
    d = dat['arr_0'].flatten()[0]
    burst = Burst(d)
    return burst


def duration_v_freq(nu, a, alpha):
    return a * nu**-alpha


def drift_v_freq(nu, a, alpha):
    return -a * nu**alpha
