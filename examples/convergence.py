"""
convergence.py

This is a step size convergence study for the EP integrator.

This script essentially runs the same EP transfer with different fixed time steps (dt)
and compares each result to a refernece solution computed with the smallest dt in the list. 

The error vs dt is plotted on a log-log scale to show numerical convergence. 
"""

import numpy as np
import matplotlib.pyplot as plt

from orbitalanalysis.ep import lowThrustTransferSim

def run(dt):
    """
    A singular EP simulaiton is ran for a fized dt, returning key scalars.
    """

    res = lowThrustTransferSim(h1 = 400e3, h2 = 800e3, inc1 = 28.5, inc2 = 28.5, tlim = 120*86400, m0 = 500, thrust = 2, isp = 1600, dt = dt)

    if not res.feasible:
        raise RuntimeError(f"Run NOT feasible at dt={dt}. Message:{res.message}")
    
    # These convergence metrics depend on integrgated history:
    # - time to reach target (event time)
    # - prop used (mass history)
    # - integrated found Delta V (post integral)
    return{"TOF":res.TOF, "propSpent":res.propSpent, "DeltV": res.DeltV}

def main():
    # dt should go from coarse to fine
    dtList = np.array([1000,900,800,700,600,500,400,300,200,100,50,25,10,5,2,1,0.5,0.25], dtype=float)

    results = []
    for dt in dtList:
        out = run(float(dt))
        results.append(out)
        print(f"dt={dt:>6.1f} s | TOF={out['TOF']/86400:>9.4f} days | "
            f"prop={out['propSpent']:.6f} kg | dv={out['DeltV']:.2f} m/s")
        
    # Reference = smallest dt which is the last
    ref = results[-1]

    tof_err = np.array([abs(r["TOF"] - ref["TOF"]) for r in results])
    prop_err = np.array([abs(r["propSpent"] - ref["propSpent"]) for r in results])
    dv_err = np.array([abs(r["DeltV"] - ref["DeltV"]) for r in results])

    # Avoid log(0) issues for very small errors
    eps = 1e-16
    tof_err = np.maximum(tof_err, eps)
    prop_err = np.maximum(prop_err, eps)
    dv_err = np.maximum(dv_err, eps)

    # Common plot settings
    FIGSIZE = (7.0, 5.0)
    LINEWIDTH = 2.0
    MARKERSIZE = 6
    TITLE_FONTSIZE = 14
    LABEL_FONTSIZE = 12
    TICK_FONTSIZE = 10
    GRID_ALPHA_MAJOR = 0.6
    GRID_ALPHA_MINOR = 0.3

    # -------------------------------
    # PLOT 1: TOF Error vs dt
    # -------------------------------
    plt.figure(figsize=FIGSIZE)

    plt.loglog(
        dtList,
        tof_err,
        marker="o",
        linewidth=LINEWIDTH,
        markersize=MARKERSIZE
    )

    plt.gca().invert_xaxis()

    plt.xlabel("Time step, Δt [s]", fontsize=LABEL_FONTSIZE)
    plt.ylabel(r"$|\,\mathrm{TOF}(\Delta t) - \mathrm{TOF}_{\mathrm{ref}}\,|$ [s]",
            fontsize=LABEL_FONTSIZE)

    plt.title("Electric Propulsion Convergence:\nTime-of-Flight Error vs Time Step",
            fontsize=TITLE_FONTSIZE)

    plt.tick_params(axis="both", which="both", labelsize=TICK_FONTSIZE)

    plt.grid(True, which="major", linestyle="-", alpha=GRID_ALPHA_MAJOR)
    plt.grid(True, which="minor", linestyle="--", alpha=GRID_ALPHA_MINOR)

    plt.tight_layout()
    plt.savefig("ep_convergence_tof.png", dpi=300)
    plt.show()


    # -------------------------------
    # PLOT 2: Propellant Error vs dt
    # -------------------------------
    plt.figure(figsize=FIGSIZE)

    plt.loglog(
        dtList,
        prop_err,
        marker="o",
        linewidth=LINEWIDTH,
        markersize=MARKERSIZE
    )

    plt.gca().invert_xaxis()

    plt.xlabel("Time step, Δt [s]", fontsize=LABEL_FONTSIZE)
    plt.ylabel(r"$|\,m_{\mathrm{prop}}(\Delta t) - m_{\mathrm{prop,ref}}\,|$ [kg]",
            fontsize=LABEL_FONTSIZE)

    plt.title("Electric Propulsion Convergence:\nPropellant Mass Error vs Time Step",
            fontsize=TITLE_FONTSIZE)

    plt.tick_params(axis="both", which="both", labelsize=TICK_FONTSIZE)

    plt.grid(True, which="major", linestyle="-", alpha=GRID_ALPHA_MAJOR)
    plt.grid(True, which="minor", linestyle="--", alpha=GRID_ALPHA_MINOR)

    plt.tight_layout()
    plt.savefig("ep_convergence_prop.png", dpi=300)
    plt.show()

    print("\nSaved figures:")
    print("  ep_convergence_tof.png")
    print("  ep_convergence_prop.png")


if __name__ == "__main__":
    main()


