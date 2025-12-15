r"""
chem.py

This is the chemical impulsive orbital transfer analysys module for the orbitalAnalysis package.

A simple and idealized chemical-propulsion system is assumed. The results from the chemical propulsion model will be comapred 
to the electric propulsion model (EP) in order to highlight the differences between the two systems. 

For chemical propulsion, the delta V maneuvers are assumed to be instantaneous (impulse)

Assumptions
-----------
- Earth is central body 
- Initial and final orbits are circular
- Transfer is a standard two-impulse Hohmann transfer
- Inclination changes are instantaneous and combined with the delta V maneuvers
- No perurbations are considered
- Minimum energy estimates

Units
-----
- Distance: meters [m]
- Time: seconds [s
- Velocity: meters per second [m/s]
- Angle: degrees for input and radians for calculations
"""

import math
# Avoid having to define __init__, __repr__, etc. for simple classes
from dataclasses import dataclass
from .constants import MU_EARTH, R_EARTH

def hohmann_transfer(r1, r2, mu=MU_EARTH):
    """
    This function computes the delta V and time of flight for coplanar Hohmann transfer.

    Parameters
    ----------
    r1: float
        Initial orbital radius from Earth's center [m]
    r2: float
        Final orbital radius
    mu = float
        Gravitational parameter of the central body [m^3/s^2]

    Returns
    -------
    dv1: float
        Delta V for the first maneuver [m/s]
    dev2: float
        Delta V for the second maneuver [m/s]
    TOF: float
        Time of flight for the transfer [s] 
    """
    # make sure orbital radii are not negative
    if r1 <= 0 or r2 <= 0:
        raise ValueError("Orbital radii must be positive values.")
    
    # Calculate circular orbital velocities
    v1 = math.sqrt(mu/r1)
    v2 = math.sqrt(mu/r2)

    # Semi-major axis of transfer orbit
    at = (r1 + r2) / 2

    # Transfer ellipse: velocities at apogee and perigee
    va = math.sqrt(mu*(2/r1 - 1/at)) # velocity at perigee
    vp = math.sqrt(mu*(2/r2 - 1/at)) # velocity at apogee

    dv1 = abs(vp - v1) # Delta V for first maneuver
    dv2 = abs(v2 - va) # Delta V for second maneuver

    # Time of Flight (TOF) for half the ellipse period 
    TOF = math.pi*math.sqrt(at**3/mu)

    return dv1, dv2, TOF

def planeChangedV(v, delta_i):
    """
    This function computes the delta-V required to change inclination impulsively.

    Parameters
    ----------
    v: float
        Orbital velocity where the plane change takes place [m/s]
    delta_i: float 
        Change in inclination [degrees]
    
    Returns
    -------
    dvPlane: float
        Delta V needed for the plane change [m/s]

    Note that for an instantaneous plane change at a constant speed V, 

        .. math::
            \Delta V = 2 V \sin{\left(\frac{\Delta i}{2}\right)}
    """

    # Make sure velocity is positive
    if v < 0: 
        raise ValueError("velocity must be non-negative")
    
    # Turn inclination degrees to radians
    delta_iRad = math.radians(delta_i)

    # return delta V for plane change
    return 2*v*math.sin(abs(delta_iRad/2))

@dataclass
class ChemicalTransferResult:
    """
    This class stores the chemical transfer results. 
    
    Contains
    ----------
    dvTotal: float
        Total delta V needed for transfer [m/s].
    TOF: float
        Time of flight for the transfer [s] for Hohmann transfer.
    location: str
        Location of where the plane change occurs.
    dvBook: dict
        Contains a dictionary with all individual delta V maneuvers [m/s].
    """

    dvTotal: float
    TOF: float
    location: str
    dvBook: dict

def HohmannWithPlaneChange(h1, h2, inc1, inc2, mu=MU_EARTH, rbody=R_EARTH):
    """
    This function computes the chemical Hohmann transfer with the inclination change.

    The following are the different options for when to perform inclination change:
        1). Plane change in the initial circular orbit
        2). Plane change in the final circular orbit
        3). Plane change at slowed point of the transfer ellipse (apoapsis for orbit raising, periapsis for orbit lowering)
    The option with the minimum total delta V is selected.

    PARAMETERS
    ----------
    h1: float
        Intial orbital altitude [m]
    h2: float
        Final orbital altitude [m]
    inc1: float
        Initial orbital inclination [degrees]
    inc2: float
        Final orbital inclination [degrees]
    mu: float
        Gravitational parameter of the central body [m^3/s^2]
    rbody: float
        Radius of the central body [m]

    RETURNS
    -------
    ChemicalTrasferResult
        A ChemicalTransferResult object containing the results of the transfer, wwhich are
        the delta V, TOF, and location of the plane change.

    """

    # Compute orbital radii
    r1 = rbody + h1
    r2 = rbody + h2

    # Inclination change
    delta_i = inc2 - inc1

    # Perform coplanar Hohmann transfer 
    dv1, dv2, TOF = hohmann_transfer(r1, r2, mu)

    # Calculate circular orbital velocities
    vcirc1 = math.sqrt(mu/r1)
    vcirc2 = math.sqrt(mu/r2)

    # Trasnfer ellipse velocities
    at = 0.5 * (r1+r2)
    vp = math.sqrt(mu*(2/r1 - 1/at))
    va = math.sqrt(mu*(2/r2 - 1/at))

    # Option 1: Plane chnage @ initial orbit
    DeltVPlane1 = planeChangedV(vcirc1, delta_i)
    DeltVTotal1 = dv1 + dv2 + DeltVPlane1

    # Option 2: Plane change @ final orbit
    DeltVPlane2 = planeChangedV(vcirc2, delta_i)
    DeltVTotal2 = dv1 + dv2 + DeltVPlane2

    # Option 3: Plane chnage @ slowest point of transfer ellipse 
    slowestV = min(vp, va)
    DeltVPlane3 = planeChangedV(slowestV, delta_i)
    DeltVTotal3 = dv1 + dv2 + DeltVPlane3

    options = [("plane change at INITIAL orbit", DeltVTotal1, DeltVPlane1),
               ("plane change at FINAL orbit", DeltVTotal2, DeltVPlane2),
               ("plane change at slowest point of transfer", DeltVTotal3, DeltVPlane3)]

    locations, DeltVTotal, DeltVPlane = min(options, key=lambda x: x[1])

    # Dictionary with the delta V break down
    dVBook = {"dv1Hohmann":dv1,"dv2Hohmann":dv2,"dvPlane": DeltVPlane}
    
    return ChemicalTransferResult(dvTotal = DeltVTotal, TOF = TOF, location = locations, dvBook = dVBook)