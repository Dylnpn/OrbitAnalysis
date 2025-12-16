"""
trade.py

This script compares chem. and EP propulsion to logic through a recommendation. 
Part of the orbitalAnalysis package.

Purpose
-------
This script compares the two transfer methods:

1). Chemical Propulsion (Impulse):
    - Hohmann  transfer and plane chnage (idealized and simplified)
    - Computes total deltaV and transfer time

2). Electric Propulsion (Low-Thrust) model:
    - Numerically intergates ODEs using Runge Kutta 4th Order.
    - Computes propellant used, time and feasibility given user's needs.

A recommendatin is provided: a propulsion system is recommended, which minimizes the 
propellant mass used and is subject to feasibility with respect to the user's desired trasnfer time.
"""
# Avoid having to define __init__, __repr__, etc. for simple classes
from dataclasses import dataclass

from .chem import HohmannWithPlaneChange
from .ep import lowThrustTransferSim

@dataclass
class TradeResults:
    """
    The propulsion trade study results are contained here.
    """
    recommendation: str
    reason: str

    # Chemical results
    chemDeltVtotal : float
    chemTOF: float
    chemStrat: str
    chemDeltVBreakdown: dict

    # EP Result
    epFeasible: bool
    epMessage: str
    epTOF: float
    epPropSpent: float
    epDeltV: float

def chemicalVsEP(h1,h2,inc1,inc2,tlim,m0,epThrust,epIsp,ep_dt = 10,):
    """
    Comaparing chemical and EP transfers, ultimately recommending a system. 

    Parameters
    ----------
    h1, h2; float
        Initial and final altidues [m]
    inc1, inc2: float
        intial and final inclination [deg]
    tlim: float
        Desired transfer time for EP [s]
    m0: float
        inital mass of spacecraft [kg]
    epThrust: float
        EP Thrust [N]
    epIsp: float
        EP Isp [s]
    ep_dt: flat
        RK4 step size for EP integration [s]
    
    RETURNS
    -------
    TradeResults
        This contains the results and a final recommendation message for the user.
    """

    # ---------------------------
    # Chemical Case (Impulsive)
    # ---------------------------
    chem = HohmannWithPlaneChange(h1, h2, inc1, inc2)

    # ---------------------------
    # EP Case (Numerical)
    # ---------------------------
    ep = lowThrustTransferSim(h1 = h1, h2=h2, inc1 = inc1, inc2 = inc2, tlim = tlim, m0 = m0, thrust = epThrust, isp = epIsp, dt = ep_dt,)

    # ---------------------------
    # Recommendation Logic
    # ---------------------------
    # If EP is NOT feasible, chemical is recommended by default
    if not ep.feasible:
        recommendation = "Chemical Propulsion"
        reason = ("EP case is NOT feasible given the thrust and time constraints."
        "Recommendation: Impulsive chemical transfer as a baseline.")
    else:
        # EP YES Feasible: compare propellant uused vs the equivalent of chemical.
        #
        # Note that computing propellant mass for a chemical rocket requires engine design analysis
        # as well as Isp. For this reason, we focus mostly on the fact that EP consumes much less 
        # propellant as a general rule of thumb, but this comes at the cost of time. EP is recommends if YES
        # feasible and prop is low. 
        recommendation = "Electric Propulsion"
        reason = ("EP case YES feasible and minimizes the propellant mass compared to a "
        "chemical propulsion system. However, this comes at the cost of longer transfer time.")

        return TradeResults(recommendation = recommendation, reason = reason,
                            
                            chemDeltVtotal=chem.dvTotal,
                            chemTOF = chem.TOF,
                            chemStrat = chem.location,
                            chemDeltVBreakdown = chem.dvBook,

                            epFeasible = ep.feasible,
                            epMessage = ep.message,
                            epTOF = ep.TOF,
                            epPropSpent = ep.propSpent,
                            epDeltV = ep.DeltV
        )

