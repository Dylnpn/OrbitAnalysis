"""
ep.py

This is the low-thrust (EP) orbit transger mode l for the orbitalAnalysis package.

Purpose
-------
This scripts provides a simplified and numerically-integrated EP transfer model that is used to be
compared against the impulsive chemical baseline with the Hohmann + plane chnage transfer. 

Here, the RK4 numerical method is used as the low thrust model is treated as a coupled ODE system.

Model
-----
A near-circular approximation is used and assumed, which is propagated using:

    y(t) = [a9t), i(t), m(t)]

Here:
- a(t) = semimajor axis [m], for circular orbits a=r
- i(t) = inclination [rad]
- m(t) = mass of spacecraft [kg]

For EP, thrust is assumed to be constant, with magnitude T. The thrust is decomposed in two directions:
- tangential (prograde), creating acceleration to raise or lower the semimajor axis
- normal (out of plane), creating an acceleration that changes the iniclination 

Orbital Dynamics
-----------------
Assume circular approximations.

mu is the gravitational parameter of the Earth.

1). The Semi-major axis chage due to tangential acceleration a_t is:
    da/dt = (2 a^(3/2) / sqrt(mu)) * a_t

2). The Inclination chnage due to the norma acceleration a_n si:
    di/dt = sqrt(a / mu) * a_n

3). The Mass flow rate is obtained from the rocket equation with constatn thrust and Isp

Flow
----
At each time t:
- calculate remaining inclination chnage di_left
- calculate time left t_left as per the user's "desired time"
- choose the normal accleration a_n needed to achieve inclination time within desired time
- limit normal acceleration a_n to available based on thrust and mass (T/m)
- the rest of the acceleration is used in the tangential direction a_t, prograde

Might not be optimal, but reasonable for the scope of the project.

Feasibility
-----------
The model has the capability of finding if the maneuver is not feasible based on two common facts:
1). Not enough acceleration to perform inclination change within the limited time
2). not enough time or thrust to reach the target altitude within the capability of the integration

Units
-----
- Distance: meters [m]
- Time: seconds [s]
- Angle: [degree] input and [radian] calculations
- Mass: [kg]
- Thrust: [N]
- Isp: [s]
"""

import math 
# Avoid having to define __init__, __repr__, etc. for simple classes
from dataclasses import dataclass

from .constants import MU_EARTH, R_EARTH, G0
from .rungeKutta4 import rk4Integrator

@dataclass
class LowThrustResult:
    """
    This class contains the results obtained for a low-thrust EP transfer simulation.

    Contains
    --------
    feasible: bool
        TRUE:If the target altitude is reached within simulaiton time limit and inclination feasibility checks passed
    message: str
        Summary of the succes or failure with warnings
    TOF: float
        Time of Flight until end of intergration or success [s].
    DeltV: float
        Delta V obtained from integrating (T/m) over time [m/s]
        This is used for comparison with chemical propulsion.
    m0: float
        Initial mass [kg]
    mf: float
        Final mass [kg]
    propSpent: float
        Propellant consumed = m0 - mf [kg]
    history: dict
        Dictionary containing results for plotting.
        - t: time [s]
        - a: semimajor axis [m]
        - inc: inclination [rad]
        - m: mass [kg]
        - a_t: tangential accel [m/s^2]
        - a_n: normal accel [m/s^2]
    """

    # initialize variable types
    feasible: bool
    message: str
    TOF: float
    DeltV: float
    m0: float
    mf: float
    propSpent: float
    history: dict


def lowThrustTransferSim(h1, h2, inc1, inc2, tlim, m0, thrust, isp, dt, mu=MU_EARTH, rBody=R_EARTH):
    """
    Parameters
    ----------
    h1, h2: float
        initial and final orbital altitudes [m].
    inc1, inc2: float
        Initial and final inclinations [deg].
    tlim: float
        Transfer time desired by user [s]. Schedules inclination chang following the Flow described.
    m0: float
        Intial spacercaft mass [kg].
    thrust: float
        Thrust with constant maginute [N].
    isp: float
        Specific Impulse [s].
    dt: float
        Step size for integratio [s].
    mu: float
        Gravitational parameter [m^3/s^2].
    rBody: float
        Earth radius [m]. 

    RETURNS
    --------


    
    """

