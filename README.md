### Overview
This is simple 1-D code for evaluating how type II and III bursts look in various coronae/winds. The work assumes:

* The wind velocity follows a polytropic profile ($\gamma=1$ returns the E. Parker 1965 solution)
* The magnetic field follows dipolar scaling and is "opened" at the Alfv\'en radius
* The emitting region of a burst's accelerator is the entire beam (_e.g_., CME or electron beam) which increases linearly in extent with distance. The velocity of the beam is constant

There are three classes:

`Star`: Stellar information including radius (`R_star`), mass (`M_star`), effective photospheric temperature (`T_eff`). The default values are estimates for the G-dwarf EK Draconis

`Corona`: Base information about the corona, including the base number density density (`n0`), an isothermal temperature (`temp`), base magnetic field strength (`B0`), and the mass fraction (`mass_fraction`). Also requires an instance of the `Star` class. Properties `r_res` and `r_max` can also be defined to set the number of spatial array elements and maximum distance, respectively, to derive wind values.  The main function of this class is `get_wind_solution`, which derives the wind speed, density, temperature, and magnetic field as a function of distance.

`Burst`: For estimating the shape of a type II and III bursts in a dynamic spectrum. The required value is the burst type ('ii' or 'iii').  Takes values for the speed of the CME or electron beam (`vb`), the height that the beam is initially accelerated in units of stellar radii (`starting_height_factor`), the initial width of the beam in units of stellar radii (`starting_width_factor`) and the growth rate of the beam with distance (`width_growth_factor`). It also requires an instance of the `Corona` class. The main function is the `make_dynamic_spectrum`, which calculates both the frequencies and whether the material exceeds the minimum velocity needed to excite plasma emission. It can be easily plotted with `plot_dyn_spec`
