# %% Cell 1: Imports
import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
import climlab
import importlib
from CourseWork.Assignments.common import generic_plots
importlib.reload(generic_plots)

# %% Cell 2: Find server data
#  This code is used just to create the skew-T plot of global, annual mean air temperature
time_coder = xr.coders.CFDatetimeCoder(use_cftime=True)

## We usually default to the data served from NOAA ESRL
ncep_url = "http://www.esrl.noaa.gov/psd/thredds/dodsC/Datasets/ncep.reanalysis.derived/"
ncep_air = xr.open_dataset( ncep_url + "pressure/air.mon.1981-2010.ltm.nc", decode_times=time_coder)

## Alternate data source if the NOAA ESRL server is down
# url = 'http://apdrc.soest.hawaii.edu:80/dods/public_data/Reanalysis_Data/NCEP/NCEP/clima/pressure/air'
# air = xr.open_dataset(url)
# The name of the vertical axis is different than the NOAA ESRL version..
#ncep_air = air.rename({'lev': 'level'})

#  Take global, annual average
coslat = np.cos(np.deg2rad(ncep_air.lat))
weight = coslat / coslat.mean(dim='lat')
Tglobal = (ncep_air.air * weight).mean(dim=('lat','lon','time'))

# Get the water vapor data from CESM output
cesm_data_path = "http://thredds.atmos.albany.edu:8080/thredds/dodsC/CESMA/"
atm_control = xr.open_dataset(cesm_data_path + "cpl_1850_f19/concatenated/cpl_1850_f19.cam.h0.nc")
# Take global, annual average of the specific humidity
weight_factor = atm_control.gw / atm_control.gw.mean(dim='lat')
Qglobal = (atm_control.Q * weight_factor).mean(dim=('lat','lon','time'))

# %% Cell 3: Create a ClimLab Model
#  Make a model on same vertical domain as the GCM
mystate = climlab.column_state(lev=Qglobal.lev, water_depth=2.5)
#  Build the radiation model -- just like we already did
rad = climlab.radiation.RRTMG(name='Radiation',
                              state=mystate,
                              specific_humidity=Qglobal.values,
                              timestep = climlab.constants.seconds_per_day,
                              albedo = 0.25,  # surface albedo, tuned to give reasonable ASR for cloud-free model
                             )
#  Now create the convection model
conv = climlab.convection.ConvectiveAdjustment(name='Convection',
                                               state=mystate,
                                               adj_lapse_rate=6.5,
                                               timestep=rad.timestep,
                                              )
#  Here is where we build the model by coupling together the two components
rcm = climlab.couple([rad, conv], name='Radiative-Convective Model')
rcm.integrate_years(5)

# %% Cell 4: 2xCO2 Model
rcm_2xCO2 = climlab.process_like(rcm)
rcm_2xCO2.name = 'Radiative-Convective Model (2xCO2)'

print(rcm_2xCO2.subprocess['Radiation'].absorber_vmr['CO2'])

rcm_2xCO2.subprocess['Radiation'].absorber_vmr['CO2'] *= 2

print(rcm_2xCO2.subprocess['Radiation'].absorber_vmr['CO2'])


# %% Cell 5: Atmosphere Adjusted Radiative Forcing

rcm_2xCO2_adjusted = climlab.process_like(rcm_2xCO2)
rcm_2xCO2_adjusted.name = 'Radiative-Convective Model (2xCO2 atmosphere-adjusted)'
for n in range(1000):
    rcm_2xCO2_adjusted.step_forward()
    # hold surface temperature fixed
    rcm_2xCO2_adjusted.Ts[:] = rcm.Ts[:]

DeltaR_adjusted = (rcm_2xCO2_adjusted.ASR - rcm_2xCO2_adjusted.OLR) - (rcm.ASR - rcm.OLR)
print(" ")
print("DeltaR (Adjusted)")
print(DeltaR_adjusted)

# %% Cell 6:
rcm_2xCO2_eq = climlab.process_like(rcm_2xCO2)
rcm_2xCO2_eq.name = 'Radiative-Convective Model (2xCO2 equilibrium)'
rcm_2xCO2_eq.integrate_years(5)

DeltaT_0 = rcm_2xCO2_eq.Ts - rcm.Ts

print(" ")
print("Equilibrium Climate Sensitivity")
print(DeltaT_0)
print("-" * 20)
print("Climate Response Parameter (No feedback)")
print(DeltaR_adjusted/DeltaT_0)

# %% Scratch Cell
