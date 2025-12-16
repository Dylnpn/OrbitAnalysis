# orbitalAnalysis

**orbitalAnalysis** is a Python package that compares **chemical (impulsive)** and **electrical propulsion (low-thrust)** for orbital transfers. This package demonstrates the use of a 4th order Runge Kutta solver, feasibility analysis, and a propulsion trade study. 

---

## Overview
The software allows the user to enter an initial and final orbit, inclination change, and desired transfer time. This software does the folloing:
- Computes an analytical chemical baseline that includes a Hohmann + plane change transfer.
- Simulates low-thrust electric propulsion (EP) transfer numerically.
- Evaluates feasibility for the EP case.
- Makes recommendation for a propulsion system to be used based on propellant efficiency. 

---

## Features

- Chemical Hohmann transfer with inclination change.
- Low-thrust EP dynamics using a fixed-step RK4 integrator.
- Event stopping simulation at target for EP.
- Trade study and recommendation logic for propulsion system.
- Command-line interface (CLI).
- Step-size convergence verification.
- Parametric trade sweep that plots propellant over time.

---

## Governing Orbital Mechanics

### Two-Body Dynamics with Thrust

Spacecraft motion is governed by the two-body equation with an applied thrust acceleration:

$$
\ddot{\mathbf{r}} = -\frac{\mu}{r^3}\mathbf{r} + \mathbf{a}_{\text{thrust}}
$$

where:
- $\mathbf{r}$ is the inertial position vector  
- $\mu$ is Earth’s gravitational parameter  
- $\mathbf{a}_{\text{thrust}}$ is the thrust acceleration  

---

## Impulsive Transfer Model (Chemical)

For a circular-to-circular Hohmann transfer:

$$
\Delta v_1 = \sqrt{\frac{\mu}{r_1}}
\left(
\sqrt{\frac{2r_2}{r_1+r_2}} - 1
\right)
$$

$$
\Delta v_2 = \sqrt{\frac{\mu}{r_2}}
\left(
1 - \sqrt{\frac{2r_1}{r_1+r_2}}
\right)
$$

$$
\Delta v_{\text{total}} = \Delta v_1 + \Delta v_2
$$

Propellant mass is computed using the Tsiolkovsky rocket equation:

$$
\Delta v = g_0 I_{sp} \ln\left(\frac{m_0}{m_f}\right)
$$

---

## Low-Thrust Propulsion Model (Electric)

For continuous tangential thrust:

$$
\mathbf{a}_{\text{thrust}} = \frac{T}{m}\hat{\mathbf{v}}
$$

Mass flow rate:

$$
\dot{m} = -\frac{T}{g_0 I_{sp}}
$$

Total propellant consumed over transfer duration $t_f$:

$$
m_{\text{prop}} = \int_0^{t_f} |\dot{m}| \, dt
$$

---

## Numerical Method

### State Vector

The propagated spacecraft state is:

$$
\mathbf{x} =
\begin{bmatrix}
\mathbf{r} \\
\mathbf{v} \\
m
\end{bmatrix}
$$

where:
- $\mathbf{r}$ is position  
- $\mathbf{v}$ is velocity  
- $m$ is spacecraft mass  

---

### Governing Ordinary Differential Equations

$$
\dot{\mathbf{r}} = \mathbf{v}
$$

$$
\dot{\mathbf{v}} = -\frac{\mu}{r^3}\mathbf{r} + \mathbf{a}_{\text{thrust}}
$$

$$
\dot{m} = -\frac{T}{g_0 I_{sp}}
$$

---

### Runge–Kutta Time Integration

The low-thrust equations of motion are propagated using a **fourth-order Runge–Kutta (RK4)** integration scheme with fixed time step $\Delta t$.

The state update is:

$$
\mathbf{x}_{k+1} =
\mathbf{x}_k +
\frac{\Delta t}{6}
\left(
\mathbf{k}_1 + 2\mathbf{k}_2 + 2\mathbf{k}_3 + \mathbf{k}_4
\right)
$$

with intermediate slopes:

$$
\mathbf{k}_1 = f(\mathbf{x}_k)
$$

$$
\mathbf{k}_2 = f\left(\mathbf{x}_k + \frac{\Delta t}{2}\mathbf{k}_1\right)
$$

$$
\mathbf{k}_3 = f\left(\mathbf{x}_k + \frac{\Delta t}{2}\mathbf{k}_2\right)
$$

$$
\mathbf{k}_4 = f\left(\mathbf{x}_k + \Delta t\,\mathbf{k}_3\right)
$$

where $f(\mathbf{x})$ represents the state derivative function.

---
### Convergence Verification

A numerical convergence is verified using **time-step refinement**, which reduces the $\Delta t$ until changes in the followin are negligible:
- Final orbital radius
- Total propellant spent

This ensures the numerical error doesn't affect the propulsion trade. Two plots are produced to ensure convergence:
    1. TOF error vs time step
    2. Prop used error vs time step
    
These plots show error convergence towards zero as the time step decreases (finer calculations).


---

### Sweep Trade Analysis

A Trade sweep is performed, comparing propellant used and allowed time for specified orbital altitude and inclination changes. The parameters to be analyzed must be manually chnaged in the "trade_sweep.py" file. 

---

 ## Installation
 Download the source code or clone the repo locally.
 In the project root directory, open a terminal and create/activate a fresh conda environment (or reuse an existing one):
 ```
 bash
 conda create-n orbitTool python=3.14
 conda activate orbitTool
 pip install -e .
 ```

 ---

 ## Quickstart
To run the trade study:
Replace all "<#>" with preferred input:
 ```
orbitAnalysis --h1-km <#> --h2-km <#> --inc1-deg <#> --inc2-deg <#> --t-days <#> --m0-kg <#> --ep-thrust <#> --ep-isp <#> --dt <#>

-----------------------------------------------------------------
Example 1:
orbitAnalysis --h1-km 200 --h2-km 300 --inc1-deg 30 --inc2-deg 32 --t-days 40 --m0-kg 1000 --ep-thrust 40 --ep-isp 1600 --dt 1
--> Should return EP YES Feasible
-----------------------------------------------------------------
Example 2: 
orbitAnalysis --h1-km 200 --h2-km 700 --inc1-deg 30 --inc2-deg 40 --t-days 4 --m0-kg 1000 --ep-thrust 2 --ep-isp 160 --dt 1
--> Should return chem NOT EP
-----------------------------------------------------------------
 ```
To run the convergence plot use the following. Note that you must go into "convergence.py" if you wiish to change the step size. The output should be two plots showing convergence. 
 ```
 python examples/convergence.py
 ```
To run the sweep trade analysis plot, run the following line. Note that in order to change the analysis trade parameters, you must go into "trade_sweep.py". The output should show EP prop used vs transfer time. 
 ```
 python examples/trade_sweep.py
 ``` 

---

## Project Structure
orbitAnalysis/\
|\
|---chem.py         # chemical (impulsive) simulation\
|---cli.py          # Facilitated user interface\
|---constants.py    # constants used\
|---ep.py           # EP numerical simulation\
|---rungeKutta4.py  # RK4 integration utilities\
|---trade.py        # compares prop systems and makes rec.\
|\
|---examples/\
| |---convergence.py # plots convergence \
| |---trade_sweep.py # plots trade analysis\
|\
|---pyproject.toml\
|---README.md\


