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
    LowThrustResult
        Result of integration that can be plotted.
    """
    # validating user input
    if tlim <= 0:
        raise ValueError("tlim must be > 0.")
    if m0 <= 0:
        raise ValueError("m0 must be > 0.")
    if thrust <= 0:
        raise ValueError("thrust must be > 0.")
    if isp <= 0: 
        raise ValueError("Isp must be > 0.")
    if dt <= 0:
        raise ValueError("dt must be > 0.")
    
    # compute radii from altitude
    a0 = rBody + h1
    atarget = rBody + h2

    # inclination in radias
    i0 = math.radians(inc1)
    itarget = math.radians(inc2)

    # Inclination change required 
    diTotal = itarget - i0

    # Mass flow rate with constant thrust assumption 
    # dm/dt = -T /(g0 * Isp)
    mdot = thrust / (G0 * isp)

    # -------------------------------------------------------------------------------
    # The following is a check to make sure the inclination scheduling is feasible
    #
    # Using di/dt = sqrt(a/mu)*a_n
    # To complete |diTotal| within time tlim, then
    #   |a_n_needed| ~ |diTotal| / tlim / sqrt(a/mu)
    #
    # The max available acceleration magnitude at the start is:
    #   aTotalMax = T/m0
    #
    # If the required normal acceleration exceeds the total acceleration, NOT FEASIBLE
    # --------------------------------------------------------------------------------
    aTotal0 = thrust/m0
    a_n_needed0 = abs(diTotal) / tlim / math.sqrt(a0 /mu)

    if a_n_needed0 > aTotal0:
        return LowThrustResult(feasible = False, message=(
            "NOT FEASIBLE (Inclination): The required normal acceleration exceeds" \
            "the vailable thrust/mass at the start. To avoid this do either of the following:" \
            "   1). Increase Thrust" \
            "   2). Increase allowed time" \
            "   3). Reduce inclination change" \
            "   4). Reduce mass." 
            ),
            TOF = 0.0,
            DeltV = 0,
            m0 = 0,
            mf = 0,
            propSpent= 0,
            history= {},
        )

    # ---------------------------------------------------------------------------------------
    # Flow Logic:
    # Choose a_n for each step in order to close the remaining inclination by time tlim.
    # Make aTotal either (+) or (-) so that the acceleration availble is not exceeded. 
    # The acceleration leftover becomes a_t (tangential) and always prograde. 
    # ---------------------------------------------------------------------------------------
    def dynamics(t,y):
        a, inc, m = y

        # When mass hits zero because no prop, stop dynamic computation to avoid dividing by zero
        if m <= 0.0:
            return [0.0, 0.0, 0.0]
        
        # Total acceleration mag. available
        aTotal = thrust/m

        # Time left based on desired flow 
        tRemain = max(tlim - t, 1e-6) # avoids zero

        # Inclination change still left to perform
        diRemain = itarget - inc

        # Normal acceleration that is required to finish the inclination chnage by the time limit tlim
        a_n_need = diRemain / tRemain / math.sqrt(a /mu)

        # Use available thrust
        a_n = max(-aTotal, min(aTotal, a_n_need))

        # The tangential acceleration uses leftover magnitude now
        a_t = math.sqrt(max(aTotal * aTotal - a_n * a_n, 0.0))

        # ODEs for near circular motion 
        dadt = (2 * (a ** 1.5) / math.sqrt(mu)) * a_t
        didt = math.sqrt(a / mu)* a_n
        dmdt = -mdot

        return [dadt, didt, dmdt]
    
    # End Event: This stopes the integrator when the desired semi major axis is reached
    def event(t,y):
        return y[0] - atarget # always negative UNTIL target reaches 
    
    # Initial state 
    y0  = [a0, i0, m0]

    # Integration time limiter:
    # For EP transfers, they usually take longer than the desired time if the thrust is low,
    # to combat that we allow twice the desired time as the integrator time limiter. 
    tf = 2 * tlim

    # Call runge-Kutta 4 to solve ODE
    result = rk4Integrator(dynamics, y0, 0, tf, dt, event=event)

    # Sort solution histories 
    tHist = result.t
    aHist = [state[0] for state in result.y]
    iHist = [state[1] for state in result.y]
    mHist = [max(state[2], 0) for state in result.y] # avoids negative mass

    # CHECK: Did we actually reach the target
    reachTarget = (aHist[-1] >= atarget - 1e-3) and result.EndEvent

    # Delta V is approximated by integrating aTotal = T / m over time spent:
    # dv = integral(T/m)dt
    dv = 0
    atHist = [0] # store for plotting later
    anHist = [0]

    # loop to start integration 
    for k in range(1, len(tHist)):
        dt = tHist[k] - tHist[k-1]

        # Midpoint mass for the aTotal estimate
        # The integration is discrete, so we're doing apprixmations at the midpoint of the time step!
        mMid = 0.5 * (mHist[k] + mHist[k-1])
        aMid = 0.5 * (aHist[k] + aHist[k-1])
        iMid = 0.5 * (iHist[k] + iHist[k-1])
        tMid = 0.5 * (tHist[k] + tHist[k-1])

        if mMid <= 0:
            atHist.append(0)
            anHist.append(0)
            continue 

        aTotal = thrust /mMid
        dv += aTotal * dt

        # Acceleration split is recomputed at the midpoint 
        tRemain = max(tlim - tMid, 1e-6)
        diRemain = itarget - iMid
        a_n_need = diRemain / tRemain / math.sqrt(aMid / mu)
        a_n = max(-aTotal, min(aTotal, a_n_need))
        a_t = math.sqrt(max(aTotal * aTotal - a_n * a_n, 0))

        atHist.append(a_t)
        anHist.append(a_n)
    
    TOF = tHist[-1]
    mf = mHist[-1]
    propSpent = m0 - mf

    # Display message for user
    if reachTarget:
        msg = "YES FEASIBLE: reached target altitude within the integration time limit."
    else: 
        msg = (" NOT FEASIBLE: Could not reach target altitude within the integration time limit. Do one of the following:"\
            "   1). Increase Thrust" \
            "   2). Increase allowed time" \
            "   3). Reduce inclination change" \
            "   4). Reduce mass.")
        
    # Prompt warning for the inclination accuracy
    incError = abs(math.degrees(iHist[-1] - itarget))
    if reachTarget and incError > 0.25:
        msg += f"Warning: final inclination error is {incError:.2f} degrees (Flow Logic)."

    return LowThrustResult(
        feasible= reachTarget,
        message = msg,
        TOF = TOF,
        DeltV = dv,
        m0 = m0,
        mf = mf,
        propSpent = propSpent,
        history={"t": tHist, "a": aHist, "inc":iHist, "m":mHist,"a_t":atHist,"a_n":anHist},
    )

