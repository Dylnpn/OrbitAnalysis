r"""
constants.py

These are physical constants that are used throughout the orbitalAnalysis package. 

All of the following values are given in SI units, and are obtained from current acceptable astrodynamics references. 
Gravitational discrepancies like J2 oblateness are neglected in order to keep the program simple. Additionally 
drag, and perturbations due to a third body are neglected as well.
"""

# Earth's gravitational parameter (mu) in m^3/s^2 also (mu = G * M_earth)
# This is used in orbital velocity and energy calculations
MU_EARTH = 3.986004418e14  # [m^3/s^2]

# Earth's radius in meters
# This is used for altitude calculations and conversions
R_EARTH = 6.371e6  # [m]

# Graviation acceleration at Earth's surface
# Used for the rocket equation and propellant mass calculations
G0 = 9.80665  # [m/s^2]
