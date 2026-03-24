### Overview
This is simple 1-D code for evaluating how type II and III bursts look in various coronae/winds. The work currently assumes:

* The wind velocity follows a polytropic profile ($\gamma=1$ returns the E. Parker 1965 solution) that does not depend on stellar rotation
* The magnetic field follows dipolar scaling and is "opened" at the Alfv\'en radius
* The emitting region of a burst's accelerator is the entire beam (_e.g_., CME or electron beam) which increases linearly in extent with distance. 
* The velocity of the beam is constant

### To install and use package:

* Run `git clone https://github.com/iveydavis/plasma_analysis.git ` or directly download code
* In the `plasma_analysis` path, run `pip install .`
* From python: `from swabs import star, isothermal_solution, polytropic_solution, bursts `
