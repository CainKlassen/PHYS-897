import numpy as np
import matplotlib.pyplot as plt
from CourseWork.Assignments.common import generic_plots

def rk1(T0, t0 = 0, tf = 50, dt = 1, C = 4*10**8, alpha = 0.32, tau = 0.57):
    """

    :param T0: Initial surface temperature
    :param t0: Starting time in years (defaults to 0)
    :param tf: Ending time in years (defaults to 50)
    :param dt: Timestep in years (defaults to 1)
    :param C: Heat capacity (defaults to 4*10**8)
    :param alpha: Albedo (defaults to 0.32)
    :param tau: Emissivity (defaults to 0.57)
    :return:
    """

    N = np.round((tf-t0)/dt).astype(int) + 1 # Find the number of iterations
    t = np.linspace(t0, tf, N)  # Create time array for plotting

    # Convert time array into years for calculations
    dt = dt * 31556926

    Ts = np.zeros(N) # Create an empty array to hold values
    Ts[0] = T0 # Set first value to initial condition

    # Define constants
    sigma = 5.67*10**(-8)
    Q = 341.2 # W/m2

    # Loop thorugh time steps
    for i in range(0,N-1):
        # Calculate radiation both incoming and outgoing
        asr = (1-alpha)*Q
        olr = tau*sigma*Ts[i]**4

        # Update the surface temperature according to model
        Ts[i+1] = Ts[i] + (asr - olr)*dt/C

    return t, Ts

# ---------------------- Plot 1 -------------------

# Create a figure to plot on
fig = plt.figure()
fig.suptitle("")

# Calculate state from the first C - value
t, Ts1 = rk1(288,
            0,
            50,
            1,
            4*10**8,
            0.32,
            0.57)

# Plot using a custom function of mine
generic_plots.plot_Generic_Scatter(t,
                         Ts1,
                         xlabel = "Time [Years]",
                         ylabel = "Temperature [K]",
                         title = "Equilibrium In The Energy Balance Model",
                         DataLabel = "C = 4e8 [J/K]",
                         PlotID = fig,
                         ReadytoShow = False,
                         Markerparam = (2, None, None),
                         xlim = None,
                         ylim = None,
                         Line = (None, 2, None) )

# Caluclate state from the second C-value
t, Ts2 = rk1(288,
            0,
            50,
            1,
            8*10**8,
            0.32,
            0.57)

# Plot both on the same figure
generic_plots.plot_Generic_Scatter(t,
                         Ts2,
                         xlabel = "Time [Years]",
                         ylabel = "Temperature [K]",
                         title = "Equilibrium In The 0-Dim Energy Balance Model",
                         DataLabel = "C = 8e8 [J/K]",
                         PlotID = fig,
                         ReadytoShow = False,
                         Markerparam = (2, None, None),
                         xlim = None,
                         ylim = None,
                         Line = (None, 2, None) )

Fig2 = plt.figure()
generic_plots.plot_Generic_Scatter(t,
                         Ts1-Ts2,
                         xlabel = "Time [Years]",
                         ylabel = "Difference in Temperature [K]",
                         title = "Temperature Difference Between Two Heat Capacities",
                         DataLabel = ".",
                         PlotID = Fig2,
                         ReadytoShow = True,
                         Markerparam = (2, None, None),
                         xlim = None,
                         ylim = None,
                         Line = (None, 2, None) )

