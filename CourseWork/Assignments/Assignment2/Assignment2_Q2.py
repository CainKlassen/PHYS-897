import numpy as np
import xarray as xr
import matplotlib.pyplot as plt

import CourseWork.Assignments.common.generic_plots

# ------------------ NCEP DATA --------------------------------
url = "http://www.esrl.noaa.gov/psd/thredds/dodsC/Datasets/ncep.reanalysis.derived/surface_gauss/air.2m.mon.ltm.nc"
ncep_air2m = xr.open_dataset(url, decode_times=False)

Ts = ncep_air2m['air']

# Find the mean weight by dividing cos(lat) by mean of cos(lat)
coslat = np.cos(np.deg2rad(Ts.lat))
mean_weight = coslat / coslat.mean(dim='lat')

Ts_avg = (Ts * mean_weight).mean(dim=('lat', 'lon')).values

# --------------------- CESM PREDICTIONS ---------------------
# Data path for the pre-industrial control CESM data
cesm_data_path = "http://thredds.atmos.albany.edu:8080/thredds/dodsC/CESMA/"
atm_control = xr.open_dataset(cesm_data_path + "cpl_1850_f19/concatenated/cpl_1850_f19.cam.h0.nc")

# Find the mean weight by dividing cos(lat) by mean of cos(lat)
coslat = np.cos(np.deg2rad(atm_control.lat))
mean_weight = coslat / coslat.mean(dim='lat')

Ts_avg_CESM = (atm_control['TS']*mean_weight).mean(dim=('lat','lon')).isel(time=slice(0,12)).values
# ----------------------------------------

months = [
    "Jan", "Feb", "Mar", "Apr", "May", "June",

    "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"
]

# ------------------------ Part a)
fig = plt.figure()

CourseWork.Assignments.common.generic_plots.plot_Generic_Scatter(xdata = months,
                                                                 ydata = Ts_avg - 273.15,
                                                                 xlabel = 'Time (Months)',
                                                                 ylabel= "Average Surface Temperature " + f"($^o$C)",
                                                                 title= "Global Average Surface Temperature Dynamics",
                                                                 DataLabel="NCEP",
                                                                 PlotID=fig,
                                                                 ReadytoShow=False,
                                                                 Markerparam=(4,None,None),
                                                                 xlim = None,
                                                                 ylim = None,
                                                                 Line = (None,3,None),)

CourseWork.Assignments.common.generic_plots.plot_Generic_Scatter(xdata = months,
                                                                 ydata = Ts_avg_CESM - 273.15,
                                                                 xlabel = 'Time (Months)',
                                                                 ylabel= "Average Surface Temperature " + f"($^o$C)",
                                                                 title= "Global Average Surface Temperature Dynamics",
                                                                 DataLabel="CESM",
                                                                 PlotID=fig,
                                                                 ReadytoShow=True,
                                                                 Markerparam=(4,None,None),
                                                                 xlim = None,
                                                                 ylim = (12,17),
                                                                 Line = (None,3,None),)