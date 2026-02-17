### Overview
This is simple 1-D code for evaluating how type II and III bursts look in various coronae/winds. The work assumes:

* The wind velocity follows a Parker spiral and the accelerators for the type II and III bursts follow this spiral
* The magnetic field follows dipolar scaling and is "opened" at the Alfv\'en radius
* The emitting region of a type II burst is associated with a CME which increases linearly with distance
* The duration of a type III burst at a given frequency follows the power-law from Alvarez & Haddock 1973 (some flexibility is provided to alter this)

There are currently four classes:

`Star`: Stellar information including radius (`R_star`), mass (`M_star`), effective photospheric temperature (`T_eff`). The default values are estimates for the G-dwarf EK Draconis

`Corona`: Base information about the corona, including the base number density density (`n0`), an isothermal temperature (`temp`), base magnetic field strength (`B0`), and the mass fraction (`mass_fraction`). Also requires an instance of the `Star` class. Properties `r_res` and `r_max` can also be defined to set the resolution and maximum distance, respectively, to derive wind values.  The main function of this class is `get_parker_solutions`, which derives the the solution to the parker wind model (equation 17 in [Parker, E. 1965](https://ui.adsabs.harvard.edu/abs/1965SSRv....4..666P/abstract)) given some precision criterion.

`Type_II_Burst`: For estimating the shape of a type II burst in a dynamic spectrum. Takes values for the speed of the CME (`vb`), the height that the CME is initially accelerated in units of stellar radii (`starting_height_factor`), the initial width of the CME in units of stellar radii (`starting_width_factor`) and the growth rate of the CME with distance (`width_growth_factor`). It also requires an instance of the `Corona` class. The main function is the `make_dynamic_spectrum`, which calculates both the frequencies and whether the material is super-Alfvenic. It can be easily plotted with `plot_dyn_spec`

`Type_III_Burst` For estimating the shape of a type III burst in a dynamic spectrum. Takes values for the speed of the electron beam (`vb`) and the height that the beam is initially accelerated in units of stellar radii (`starting_height_factor`) It also requires an instance of the `Corona` class. The main function is the `make_dynamic_spectrum`, which be easily plotted with `plot_dyn_spec`
