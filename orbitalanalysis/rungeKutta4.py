"""
rungeKutta4.py

This scripts containes the most important part of the package: 4th Order Runge-Kutta numerical solver.
The 4th-Order Runge-Kutta (RK4) integrator is a classex fixed-step numerical method used to solve 
ordinary differential equations (ODEs).

This package containes the use of RK4 because of the use of eelctric propulsion as one of our propulsion systems.
EP systems are low thrust, so they do not impose an instantanous delta V on a vehicle in orbit. Low thrust orbital
transfers are oftentimes modeled as ODEs due to the continues state chnage with time. These types of ODEs do not 
have a closed-form solution most of the time when thrust and mass play a key role. 

Note that structure was inspired by MATLAB's ODE113.

Assumptions
-----------
- Fixed time step (dt) for easy convergence
- "Event" function to stop integration after orbital transfer is performed

Units
-----
User is responsible for inputting adequate units.
"""

# Avoid having to define __init__, __repr__, etc. for simple classes
from dataclasses import dataclass

@dataclass 

class RK4Output:
    """
    This class stores the RK4 numerical solutions.

    Contains
    ---------
    t: list of floats list[float]
        Time stamps for the solution 
    y: list of list of floats list[list[float]]
        y[k] is a full state at each time t[k[]]
    EndEvent: bool
        If integration stops because condition is met, then TRUE
    """

    # Initialize lists and bool
    t: list
    y: list
    EndEvent: bool

def rk4Integrator(f, y0, t0, tf, dt, event=None):
    """
    The 4th-Order Runge Kutta numerical solver lives here. It integrates an ODE system.

    Parameters
    ----------
    f: is callable
    f(t, y) with dy_dt. This is a dynamic function with two variables
        - t is time (float)
        - y is the state (list[float])
        - dy_dt is the time derivative of the state (list[float])
    y0: list[float]
        Intial condition of state vector
    t0: float
        Initial time. 
    tf: float
        Final time in the integrator
    dt: float
        Integration step size (fixed)
    event:
        The integration stops when g crosses from negative to more than or qual to 0. 
        This is ideal when used for checking if altitude has been reached. 
        g(t, y) to float. 

    RETURNS
    -------
    RK4Output
        Time (history), state (history), bool of state (did it stop?)

    Numerical Method
    ----------------
    RK4 updates with:
    y_{n+1} = y_n + 9dt/6) * (k1 + 2*k2 + 2*k3 + k4)
    with:
        k1 = f(t, y)
        k2 = f(t + dt/2, y + dt*k1/2)
        k3 = f(t + dt/2, y + dt*k2/2)
        k4 = f(t + dt, y + dt*k3)

    Note that each ki uses the result from k{i-1}
    No adaptive step sizing for simplicity of project!
    """

    # make sure change in time is logical and not negative
    if dt <= 0:
        raise ValueError("dt is negative... dt MUST BE > 0")
    
    # to avoid rounding errors and error propagation, convert states to floats
    y = [float(v) for v in y0]
    t = float(t0)

    # Store arrays of times and states
    tHist = [t]
    yHist = [y.copy()]

    # EndEvent to stop integration (keeps track) 
    past_g = event(t,y) if event is not None else None # ternary operator (if event is not none, evaluate to the value of event, otherwise evaluate it to None)
    EndEvent = False
    # if event isnt None, return current value in event
    # if event IS none, return None
    # Detects if sign is negative
   
    # create function to allow vector addition since python cannot do it on its own
    def vecAddition(a, b, scale =1.0):
        """
        Return the vector: a + scale*b
        """
        return [ai + scale * bi for ai, bi in zip(a,b)]
    # zip(a,b) pairs corresponding elements in the same branched arrays 
    # so if you have:
    #   vecAddition([1, 2], [3, 4]) --> [4,6]

    # Ensure the while loop doesn't bug out due to small tolerance
    timeBuffer = 1e-11
    while t < tf - timeBuffer:

        # The time interval [t0, tf] isnt usually an exact multiple of dt,
        # so what this lin tells us is that we stop at a time lower than tf, never above it; this 
        # then allows h = tf - t, which stops at tf. It serves as a check to avoid bugs. 
        h = min(dt, tf - t)

        k1 = f(t,y)
        k2 = f(t+0.5*h, vecAddition(y, k1, 0.5*h))
        k3 = f(t+0.5*h, vecAddition(y, k2, 0.5*h))
        k4 = f(t + h, vecAddition(y, k3, h))

        # Average slope weighted
        slope = [(k1i + 2*k2i + 2*k3i + k4i)/6 for k1i, k2i, k3i, k4i in zip(k1, k2, k3, k4)]

        # Move forward in state and time
        y = vecAddition(y, slope, h)
        t = t + h
        # Add these values to the histories
        tHist.append(t)
        yHist.append(y.copy())

        # Check if the EndEvent condition has been met: stop if we go from g > 0 to  g >= 0
        if event is not None:
            g = event(t, y)
            # Stop event!!!
            if past_g is not None and (past_g < 0.0) and (g >= 0.0):
                EndEvent = True
                break
            # Don't stop, proceed with past event state
            past_g = g

    return RK4Output(t = tHist, y=yHist, EndEvent = EndEvent)