"""
cli.py

This is the Command-Line Interface (CLI) for the orbitalAnalysis package.

Purpose
-------
This script allows a user to run the propulsion trade study directly from a terminal:

    orbitalalysis --h1-km 400 --h2-km 800 --inc1-deg 28.5 --inc2-deg 28.5\
        --t-days 120 --m0-kg 500 --ep-thrust 2.0 --ep-isp 16000 --dt 20

This avoids complexity by being able to interact with the program without editing every
Python script. 

"""

import argparse 

from .trade import chemicalVsEP

def build_parser():
    """
    The argument parser is created and returned.
    The CLI is easier to test and extend when put in a different function like this. 
    """
    p = argparse.ArgumentParser(
        prog = "orbitalanalysis",
        description = "Compare chemical (impulsive) vs EP (low-thrust) orbit transfers."
    )

    # Geometry inputs for orbit
    p.add_argument("--h1-km", type=float, required=True, help = "Initial altitude [km].")
    p.add_argument("--h2-km", type=float, required=True, help = "Final altitude [km].")
    p.add_argument("--inc1-deg", type=float, required=True, help = "Initial inclination [deg].")
    p.add_argument("--inc2-deg", type=float, required=True, help = "Final inclination [deg].")

    # Constrain inputs for transfer
    p.add_argument("--t-days", type=float, required=True, help ="Desired transfer time [days].")

    # Spacecraft / EP inputs 
    p.add_argument("--m0-kg", type=float, required=True, help = "Initial spacecraft mass [kg].")
    p.add_argument("--ep-thrust", type=float, required=True, help = "EP Thrust [N].")
    p.add_argument("--ep-isp", type=float, required=True, help = "EP specific impulse [s].")

    # Numerical control
    p.add_argument("--dt", type=float, default =10, help = "RK4 step size for EP integration [s].")

    return p

def main(argv=None):
    """
    This is the CLI entry point.

    Parameters
    ----------
    argv: list[str] or None
        Optional list of arguments for testing. 
    """

    parser = build_parser()
    args = parser.parse_args(argv)

    # convert to SI units
    h1 = args.h1_km * 1000
    h2 = args.h2_km * 1000
    tlim = args.t_days * 86400

    # Run Trade
    res = chemicalVsEP(h1=h1, h2=h2, inc1=args.inc1_deg, inc2=args.inc2_deg, tlim = tlim, m0 = args.m0_kg, epThrust = args.ep_thrust, epIsp = args.ep_isp, ep_dt = args.dt,)

    # Print a report
    print("\n=== orbitalAnalysis Trade Study ===")
    print(f"Recommendation: {res.recommendation.upper()}")
    print(f"Reason: {res.reason}\n")

    print("Chemical (Impulsive):")
    print(f"    Total Delta V: {res.chemDeltVtotal:.2f} m/s")
    print(f"    Hohmann TOF: {res.chemTOF/3600:.2f} hr")
    print(f"    Plane Change: {res.chemStrat}")
    print(f"    Breakdown: {res.chemDeltVBreakdown}\n")

    print("Electric Propulsiion (Low-Thrust):")
    print(f"    Feasible: {res.epFeasible}")
    print(f"    Message: {res.epMessage}")
    print(f"    TOF: {res.epTOF/3600:.3f} hr")
    print(f"    Propellant Spent: {res.epPropSpent:.3f} kg")
    print(f"    Total Delta V: {res.epDeltV:.2f} m/s\n")

if __name__ == "__main__":
    main()