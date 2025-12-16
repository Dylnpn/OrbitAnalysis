"""
trade_sweep.py

Parametric sweep for EP trade study: propellant used vs allowed transfer time.

This script is meant to generate an “engineering insight” plot for your report:
- EP is mass-efficient but can be infeasible if the allowed time is too short,
  especially for large inclination changes.

We sweep t_days and run the EP simulation each time.
"""

import numpy as np
import matplotlib.pyplot as plt

from orbitalanalysis.ep import lowThrustTransferSim


def run_ep_case(t_days, h1_m, h2_m, inc1_deg, inc2_deg, m0_kg, thrust_n, isp_s, dt_s):
    """Run EP sim for a single time value and return (feasible, prop_used_kg, tof_days)."""
    res = lowThrustTransferSim(
        h1=h1_m,
        h2=h2_m,
        inc1=inc1_deg,
        inc2=inc2_deg,
        tlim=t_days * 86400.0,
        m0=m0_kg,
        thrust=thrust_n,
        isp=isp_s,
        dt=dt_s,
    )
    return res.feasible, res.propSpent, res.TOF / 86400.0, res.message


def sweep_time(t_days_list, case_label, h1_m, h2_m, inc1_deg, inc2_deg, m0_kg, thrust_n, isp_s, dt_s):
    """
    Sweep over transfer time and return arrays for plotting.
    Infeasible cases are stored as NaN for propellant so the curve breaks cleanly.
    """
    prop_used = np.full_like(t_days_list, np.nan, dtype=float)
    tof_used = np.full_like(t_days_list, np.nan, dtype=float)
    feasible_mask = np.zeros_like(t_days_list, dtype=bool)

    for k, t_days in enumerate(t_days_list):
        feasible, prop, tof_days, msg = run_ep_case(
            t_days=float(t_days),
            h1_m=h1_m,
            h2_m=h2_m,
            inc1_deg=inc1_deg,
            inc2_deg=inc2_deg,
            m0_kg=m0_kg,
            thrust_n=thrust_n,
            isp_s=isp_s,
            dt_s=dt_s,
        )

        feasible_mask[k] = feasible
        if feasible:
            prop_used[k] = prop
            tof_used[k] = tof_days

        # Optional: print status for debugging / transparency
        print(f"[{case_label}] t={t_days:6.1f} days | feasible={feasible} | "
              f"prop={prop:10.4f} kg | tof_used={tof_days:8.2f} days | {msg}")

    return prop_used, tof_used, feasible_mask


def main():
    # -----------------------
    # Fixed scenario settings
    # -----------------------
    h1_m = 400e3
    h2_m = 800e3

    m0_kg = 500.0
    thrust_n = 2.0
    isp_s = 1600.0

    # Integration step size (keep fixed for the sweep)
    dt_s = 20.0

    # -----------------------
    # Sweep definition (days)
    # -----------------------
    # Use a range that will show both infeasible and feasible regions.
    t_days_list = np.array([1,5, 7, 10, 14, 20, 30, 45, 60, 90, 120, 150, 200], dtype=float)

    # Two cases: no plane change vs significant plane change
    cases = [
        {"label": "Δi = 0° (28.5 → 28.5)", "inc1": 28.5, "inc2": 28.5},
        {"label": "Δi = 23.1° (28.5 → 51.6)", "inc1": 28.5, "inc2": 51.6},
    ]

    # Run sweeps
    outputs = []
    for c in cases:
        prop_used, tof_used, feasible_mask = sweep_time(
            t_days_list=t_days_list,
            case_label=c["label"],
            h1_m=h1_m,
            h2_m=h2_m,
            inc1_deg=c["inc1"],
            inc2_deg=c["inc2"],
            m0_kg=m0_kg,
            thrust_n=thrust_n,
            isp_s=isp_s,
            dt_s=dt_s,
        )
        outputs.append((c["label"], prop_used, feasible_mask))

    # -----------------------
    # Plot: propellant vs time
    # -----------------------
    # ---- Style constants (tweak once, applies everywhere)
    FIGSIZE = (7.5, 5.0)
    LINEWIDTH = 2.0
    MARKERSIZE = 5.5
    TITLE_FONTSIZE = 14
    LABEL_FONTSIZE = 12
    TICK_FONTSIZE = 10
    LEGEND_FONTSIZE = 10
    GRID_ALPHA_MAJOR = 0.6
    GRID_ALPHA_MINOR = 0.3

    # Safety: ensure arrays are numpy arrays for masking/indexing
    t_days = np.asarray(t_days_list, dtype=float)

    plt.figure(figsize=FIGSIZE)

    # Track global feasible min to place infeasible markers consistently
    global_feasible_min = np.inf
    for (_, prop_used, feasible_mask) in outputs:
        prop = np.asarray(prop_used, dtype=float)
        feas = np.asarray(feasible_mask, dtype=bool)
        if np.any(feas):
            global_feasible_min = min(global_feasible_min, np.nanmin(prop[feas]))

    # Fallback if nothing is feasible anywhere
    if not np.isfinite(global_feasible_min):
        global_feasible_min = 0.0

    # Place infeasible markers slightly below the smallest feasible value
    y_mark = 0.9 * global_feasible_min if global_feasible_min > 0 else 0.0

    for (label, prop_used, feasible_mask) in outputs:
        prop = np.asarray(prop_used, dtype=float)
        feas = np.asarray(feasible_mask, dtype=bool)

        # Plot feasible portion as a continuous line; NaNs break the line
        prop_plot = prop.copy()
        prop_plot[~feas] = np.nan

        plt.plot(
            t_days,
            prop_plot,
            marker="o",
            markersize=MARKERSIZE,
            linewidth=LINEWIDTH,
            label=label
        )

        # Mark infeasible points clearly at a consistent baseline
        if np.any(~feas):
            plt.plot(
                t_days[~feas],
                np.full(np.sum(~feas), y_mark),
                marker="x",
                linestyle="none",
                markersize=MARKERSIZE + 1,
                label=None
            )

    # Labels and title
    plt.xlabel("Allowed transfer time [days]", fontsize=LABEL_FONTSIZE)
    plt.ylabel("EP propellant used [kg]", fontsize=LABEL_FONTSIZE)
    plt.title("EP Trade Sweep: Propellant vs Allowed Time (400 km → 800 km)",
            fontsize=TITLE_FONTSIZE)

    # Ticks/grid
    plt.tick_params(axis="both", which="both", labelsize=TICK_FONTSIZE)
    plt.minorticks_on()
    plt.grid(True, which="major", linestyle="-", alpha=GRID_ALPHA_MAJOR)
    plt.grid(True, which="minor", linestyle="--", alpha=GRID_ALPHA_MINOR)

    # Legend: put it outside so it never covers data
    plt.legend(
        frameon=True,
        fontsize=LEGEND_FONTSIZE,
        loc="center left",
        bbox_to_anchor=(1.02, 0.5),
        borderaxespad=0.0
    )

    plt.tight_layout()
    plt.savefig("trade_sweep_ep_prop_vs_time.png", dpi=300, bbox_inches="tight")
    plt.show()
    print("\nSaved figure: trade_sweep_ep_prop_vs_time.png")


if __name__ == "__main__":
    main()
