import numpy as np
import xarray as xr

# Data path for the pre-industrial control CESM data
cesm_data_path = "http://thredds.atmos.albany.edu:8080/thredds/dodsC/CESMA/"

# Isolate the atmospheric data specifically
atm_control = xr.open_dataset(cesm_data_path + "cpl_1850_f19/concatenated/cpl_1850_f19.cam.h0.nc")

# Find the mean weight by dividing cos(lat) by mean of cos(lat)
coslat = np.cos(np.deg2rad(atm_control.lat))
mean_weight = coslat / coslat.mean(dim='lat')

def take_global_average(dataname:str, xarray_data):
    """
    :param dataname: String containing the DataArray dictionary name
    :return: Prints the dataname as well as it's value
    """

    # Print data variable name
    print(dataname + ":")
    # Print global average of data variable
    print((xarray_data * mean_weight).mean(dim=('lat','lon','time')).values)
    print("-"*20)

# Solar Radiation Budget
solar_insolation = atm_control['SOLIN']
toa_reflected = atm_control['FSUTOA']
surface_reflected = atm_control['FSDS'] - atm_control['FSNS']
surface_absorbed = atm_control['FSNS']
clouds_reflected = toa_reflected - surface_reflected
toa_absorbed = solar_insolation - toa_reflected
total_absorbed = toa_absorbed - surface_absorbed

take_global_average("Solar Insolation", solar_insolation)
take_global_average("Reflected Insolation", toa_reflected)
take_global_average("Reflected at Surface", surface_reflected)
take_global_average("Absorbed at Surface", surface_absorbed)
take_global_average("Reflected by Clouds", clouds_reflected)
take_global_average("Absorbed at Top Of Atmosphere", toa_absorbed)
take_global_average("Absorbed by Atmosphere", total_absorbed)

# Longwave Radiation Budget
OLR = atm_control['FLUT']
surface_olr_up = atm_control['FLDS'] + atm_control['FLNS']
surface_olr_down = atm_control['FLDS']

take_global_average("Outgoing Longwave Radiation (OLR)", OLR)
take_global_average("Upward Longwave Emission From Surface", surface_olr_up)
take_global_average("Downward Longwave Emission From Surface", surface_olr_down)

# Other surface fluxes
thermals = atm_control['SHFLX']
evapotranspiration = atm_control['LHFLX']

take_global_average("Thermals", thermals)
take_global_average("Evapotranspiration", evapotranspiration)

