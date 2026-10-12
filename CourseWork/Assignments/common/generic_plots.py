import matplotlib.pyplot as plt
import numpy as np

# ------------------------------------------------------------------------------------------------------------------- #
# -------------------------------------------- FUNCTION GENERIC_HIST() ---------------------------------------------- #
# ------------------------------------------------------------------------------------------------------------------- #
def plot_Generic_Hist(data, XlimThreshold, N_bins = 10, units = "", xlabel = "", ylabel = "",
                      title = "", DataLabel = "", ProbabilityDensity = False, PlotRange = None,
                      PlotID = None, Plot_as_Scatter = False, color = None, ReadytoShow = True):
    """
    This function takes an arbitrary set of data and creates a histogram that fits my general preference
    without needing to re-do all the plotting code every time I want to make a new histogram

    Parameters
    ----------
    data - This is a (1xn) dimensional array that contains the data to make a histogram out of

    XlimThreshold - This parameter is useful when the input data has a wide spread with tails that quickly
                     drop off but don't go exactly to zero. For example, in a Gaussian distribution that peaks
                     at 10,000 over a range of 5,000 bins, the user could opt to choose only to include the
                     bins which have more than 100 counts in them. This might reduce the x limits to 500 bins
                     which makes the output easier to analyze as it is no longer dominated by trailing zeros

    N_bins - This variable allows the user to choose how many bins exist within the histogram

    units - This is a string that will be appended to the x-axis in order to indicate the bin resolution
                 of the histogram

    xlabel - This parameter is a string that will be added to the x label. It is left empty by default

    ylabel - This parameter is a string that will be added to the y label. It is left empty by default

    title - This parameter is a string that will be used as the title for the histogram plot. It is empty
            by default

    ProbabilityDensity - Setting this parameter to "True" results in a distribution that is normalized such
                         that the area under the curve is 1. If unset, it defaults to "False"

    PlotRange - This parameter sets the range for the x-axis. If it is unset by the user, it defaults to
                the range of which the data spans (max - min)/#ofbins. The expected data type is a tuple
                with form (datamin, datamax)

    PlotID - If the user wishes to plot this histogram on an existing figure they can do so by passing the
             plot ID here

    Plot_as_scatter - Setting this to True results in a histogram containing points located at the middle
                      of each bin rather than a bar

    color - This parameter is used for setting the plot color as dictated by plt.histogram()

    ReadyToShow - Setting this parameter to True causes the plot to show. It can be set to False to
                  suppress intermediate outputs if a single figure is to contain multiple plots

    Returns
    -------
    The output of this function is a histogram of the input data over a user specified range. This is intended to
    make histogram plotting simpler without needing to recall all the plotting and formatting code every time I
    want to make a plot.

    """

    # If the user doesn't want to specify a range then default to the datas max and min
    if PlotRange is None:
        PlotRange = (float(np.min(data)), float(np.max(data)))

    # Create a mask to remove trailing zeros from the DN values array

    hist_data, bin_edges = np.histogram(data, bins = N_bins, range=PlotRange, density=ProbabilityDensity)
    nonzero_idx = np.where(hist_data > XlimThreshold)[0]

    if len(nonzero_idx) == 0:
        raise ValueError("\n" "\n" "NO HISTOGRAM BINS HAVE A GREATER VALUE THAN XlimThreshold!!!" 
              rf" THE MAXIMUM IS {np.max(hist_data):.0f}")

    start, end = nonzero_idx[0], nonzero_idx[-1]

    binwidth = (PlotRange[-1] - PlotRange[0])/N_bins

    # Create a new figure for the plot if specified
    if PlotID is None:
        PlotID = plt.figure()
    else:
        plt.figure(PlotID)

    # Plot as a scatter plot if requested
    if Plot_as_Scatter:
        bin_centers = (bin_edges[:-1] + bin_edges[1:])/2
        plt.scatter(bin_centers[start:end+1],hist_data[start:end+1], s=1, label=DataLabel)
    # Else Plot as a bar histogram
    else:
        plt.hist(data, bins=N_bins, density = ProbabilityDensity, range = PlotRange, label=DataLabel, color=color)

    # Trim the Axes and apply the labels
    plt.title(title)
    plt.xlabel(xlabel + " " + rf"({binwidth:.2f}" + units + "/bin)")
    plt.ylabel(ylabel)
    plt.xlim((bin_edges[start], bin_edges[end+1]))
    plt.grid(True, alpha=0.3)

    plt.legend()

    # Show plot if ready
    if ReadytoShow:
        plt.tight_layout()
        plt.show()
        plt.close()

    return PlotID

# ------------------------------------------------------------------------------------------------------------------- #
# ------------------------------------------------- END OF FUNCTION ------------------------------------------------- #
# ------------------------------------------------------------------------------------------------------------------- #



# ------------------------------------------------------------------------------------------------------------------- #
# ----------------------------------------- FUNCTION PLOT_GENERIC_HEATMAP() ----------------------------------------- #
# ------------------------------------------------------------------------------------------------------------------- #
def plot_Generic_Heatmap(data,
                         xlabel = "",
                         ylabel = "",
                         title = "",
                         ColorbarLabel = "",
                         XYlim = None,
                         PlotID = None,
                         ReadytoShow = True,
                         minpercentile = 0,
                         maxpercentile = 100,
                         swapy = False):
    """
    This function takes an arbitrary set of data and creates a histogram that fits my general preference
    without needing to re-do all the plotting code every time I want to make a new histogram

    Parameters
    ----------
    data - This is a (1xn) dimensional array that contains the data to make a histogram out of

    xlabel - This parameter is a string that will be added to the x label. It is left empty by default

    ylabel - This parameter is a string that will be added to the y label. It is left empty by default

    title - This parameter is a string that will be used as the title for the histogram plot. It is empty
            by default

    ColorbarLabel - This input expects a string and is used to set the label for the colorbar that is shown
                    alongside the heatmap

    XYlim - This parameter allows the user to set the units for the X and Y axis if they don't want to
            default to integer indexes. The expected data type is a list of form [xmin, xmax, ymin, ymax]

    PlotID - If the user wishes to plot this histogram on an existing figure they can do so by passing the
             plot ID here

    ReadyToShow - Setting this parameter to True causes the plot to show. It can be set to False to
                  suppress intermediate outputs if a single figure is to contain multiple plots

    Minpercentile - This parameter is a float and is used to determine the lower percentile of data to
                    include in the heatmap. Increasing this removes low valued outliers

    Maxpercentile - This parameter is a float and is used to determine the upper percentile of data to
                    include in the heatmap. Increasing this removes high valued outliers

    swapy - Changes the y axis from ascending to descending to match pixel indexing convention

    Returns
    -------
    The output of this function is a heat map of the input data over a user specified range. This is intended to
    make 2D plotting simpler without needing to recall all the plotting and formatting code every time I
    want to make a plot.

    """

    # If the user doesn't want to specify a range then default to all values
    varmin, varmax = np.nanpercentile(data, [minpercentile, maxpercentile])

    # Create a new figure for the plot if specified
    if PlotID is None:
        PlotID = plt.figure()
    else:
        plt.figure(PlotID)

    # Check if y axis should be swapped
    if swapy:
        swapy = "upper"
        XYlim[2], XYlim[3] = XYlim[3], XYlim[2]
    else:
        swapy = "lower"

    # Plot the heatmap
    plt.imshow(data, extent = XYlim, vmin = varmin, vmax = varmax, origin = swapy)

    # Trim the Axes and apply the labels
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.colorbar(label=ColorbarLabel, fraction=0.046, pad=0.04)

    # Show plot if ready
    if ReadytoShow:
        plt.tight_layout()
        plt.show()
        plt.close()

    return PlotID

# ------------------------------------------------------------------------------------------------------------------- #
# ------------------------------------------------- END OF FUNCTION ------------------------------------------------- #
# ------------------------------------------------------------------------------------------------------------------- #



# ------------------------------------------------------------------------------------------------------------------- #
# ----------------------------------------- FUNCTION plot_Generic_Scatter() ----------------------------------------- #
# ------------------------------------------------------------------------------------------------------------------- #
def plot_Generic_Scatter(xdata, ydata, xlabel = "", ylabel = "", title = "", DataLabel = "",
                         PlotID = None, ReadytoShow = True, Markerparam = (None, None, None), xlim = None,
                         ylim = None, Line = (None, None, None) ):
    """

    Parameters
    ----------
    xdata (np.array) - This contains an array of x axis data to be plotted. Must be the same length as ydata

    ydata (np.array) - This contains an array of y axis data to be plotted. Must be the same length as xdata

    xlabel (string) - This parameter is a string that will be added to the x label. It is left empty
                      by default

    ylabel (string) - This parameter is a string that will be added to the y label. It is left empty by default

    title (string) - This parameter is a string that will be used as the title for the histogram plot. It
                     is empty by default

    DataLabel (string) - This input expects a string and is used to set the label for the data that was
                         just plotted

    PlotID (IDK) - If the user wishes to plot this histogram on an existing figure they can do so by
                   passing the plot ID here

    ReadyToShow (Bool) - Setting this parameter to True causes the plot to show. It can be set to False to
                         suppress intermediate outputs if a single figure is to contain multiple plots

    Markerparam (float, string, string) - This parameter sets the size of each plot point as dictated by
                                          plt.scatter() along with the marker size and style in the form

                                               e.g (markersize, color, markerstyle)

    xlim (float) - This parameter expects a tuple in the form (min,max) to set the x limit of the
                   scatter plot

    ylim (float) - This parameter expects a tuple in the form (min,max) to set the y limit of the
                   scatter plot

    Line (float, float, string) - Indicate whether the user wants a line connecting the points. This is a
                                 tuple where the first parameter is a Bool saying if to make the line,
                                 the second is a float showing line width, and the third is a string
                                 dictating the linestyle in the form

                                       e.g (color, linewidth, linestyle)

    Returns
    -------
    The output is a generic scatter plot with a figure produced using standard plt.scatter()

    """

    # If needing to plot on an existing figure, use the ID given by the user
    if PlotID is None:
        PlotID = plt.figure()
    else:
        plt.figure(PlotID)

    # Create the scatter plot and apply labels
    plt.scatter(xdata, ydata, label = DataLabel, s=Markerparam[0], c=Markerparam[1], marker=Markerparam[2])
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)

    if any(x is not None for x in Line):
        plt.plot(xdata, ydata, color = Line[0], linewidth = Line[1], linestyle = Line[2])

    # Apply the limits on the plot
    plt.xlim(xlim)
    plt.ylim(ylim)
    plt.grid(True, alpha=0.3)

    # Apply the legend
    plt.legend()

    # Show plot if ready
    if ReadytoShow:
        plt.tight_layout()
        plt.show()
        plt.close()

    return PlotID


# ------------------------------------------------------------------------------------------------------------------- #
# ------------------------------------------------- END OF FUNCTION ------------------------------------------------- #
# ------------------------------------------------------------------------------------------------------------------- #








# ------------------------------------------------------------------------------------------------------------------- #
# -------------------------------------------- FUNCTION () -------------------------------------------- #
# ------------------------------------------------------------------------------------------------------------------- #
# FUNCTION DECLARATION HERE
#
# Description:
#
#
# Inputs:
#
#
# Output:
#
#
# ---------------------------------------------------------------------------------------------------------------------

# CODE HERE



# ------------------------------------------------------------------------------------------------------------------- #
# ------------------------------------------------- END OF FUNCTION ------------------------------------------------- #
# ------------------------------------------------------------------------------------------------------------------- #