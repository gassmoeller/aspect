# Python script to load the solution data into a
# numpy array, interpolate the data to a new
# uniform grid, and plot the results.

# Load modules
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
from mpl_toolkits import mplot3d
import scipy.interpolate as interpolate
import pyshtools as pysh

# Relative (from current directory) path to the
# solution data folder.
mypath = os.path.dirname(__file__)
output_directory = os.path.join(mypath, 'output-spherical_harmonic_expansion')
solution_directory = os.path.join(output_directory, 'solution')

# Time step to analyze.
output_number = '00000'

text_filename = os.path.join(output_directory, 'heat_flux.' + output_number)

#heat_flux = np.genfromtxt(text_filename, delimiter=' ', skip_header=1)
heat_flux = pd.read_csv(text_filename, delimiter=' ', skiprows=1, header=0, names=['x','y','z','heat_flux'])

heat_flux['r'] = heat_flux.apply(lambda row: np.sqrt(row.x*row.x+row.y*row.y+row.z*row.z), axis = 1)
heat_flux['latitude'] = heat_flux.apply(lambda row: 90 - np.arccos(row.z/row.r) * 180 / np.pi, axis = 1)
heat_flux['longitude'] = heat_flux.apply(lambda row: 180 + np.arctan2(row.y,row.x) * 180 / np.pi, axis = 1)

bottom_heat_flux = heat_flux[heat_flux.r < 4e6].to_numpy()

# Gauss Legendre sample points
lmax = 32
latglq, longlq = pysh.shtools.GLQGridCoord (lmax)

# Equally spaced sample points
#lat_i = np.linspace(-90, 90, 51)
#long_i = np.linspace(-180, 180, 51)

# Generate data interpolation grid. 
Lon, Lat = np.meshgrid(longlq,latglq)

# create interpolation function
#function = interpolate.interp2d(bottom_heat_flux[:,5], bottom_heat_flux[:,6], bottom_heat_flux[:,3], kind='linear')
#function = interpolate.interp2d(bottom_heat_flux[:,5], bottom_heat_flux[:,6], bottom_heat_flux[:,3], kind='linear')
#bottom_heat_flux_interpolated = function(lat_i, long_i)

# Plot heat flux field
fig, ax = plt.subplots(2,2)

# Cartesian scatter
# ax0 = fig.gca(projection='3d')
# c = ax0.scatter(bottom_heat_flux[:,0],bottom_heat_flux[:,1],bottom_heat_flux[:,2],c=bottom_heat_flux[:,3])

# Geographic scatter
c = ax[0,0].scatter(bottom_heat_flux[:,6],bottom_heat_flux[:,5],c=bottom_heat_flux[:,3],s=1)
ax[0,0].set_aspect('equal', 'box')
ax[0,0].set_xlabel('Longitude')
ax[0,0].set_ylabel('Latitude')
ax[0,0].set_title('Heat flux on cell centers', fontsize=10)
fig.colorbar(c, ax=ax[0,0])

#ax1 = fig.add_subplot(212)
#d = ax1.scatter(Lon, Lat,c=bottom_heat_flux_interpolated)

# Geographic contour
# Interpolate heat flux onto new grid
heat_flux_structured = interpolate.griddata((bottom_heat_flux[:,5:]), bottom_heat_flux[:,3], (Lat, Lon), method='cubic', fill_value=0.0)
d = ax[1,0].pcolormesh(Lon,Lat,heat_flux_structured, shading='auto')
ax[1,0].set_aspect('equal', 'box')
ax[1,0].set_xlabel('Longitude')
ax[1,0].set_ylabel('Latitude')
ax[1,0].set_title('Heat flux on GLQ points', fontsize=10)
fig.colorbar(d, ax=ax[1,0])


# Spherical harmonic expansion
zero, w = pysh.expand.SHGLQ (lmax)
cilm = pysh.expand.SHExpandGLQ (heat_flux_structured, w, zero, norm=2, csphase=-1)
power_per_l = pysh.spectralanalysis.spectrum(cilm,normalization='schmidt')
#power_per_l = np.arange(cilm.shape[1])
degrees = np.arange(cilm.shape[1])
e = ax[0,1].plot(degrees, power_per_l)
ax[0,1].set_xlim(0, lmax)
ax[0,1].set_ylim(3e-7, 3e-3)
ax[0,1].set_yscale('log')
ax[0,1].set_xlabel('SPH degree')
ax[0,1].set_ylabel('Power')
ax[0,1].set_title('Heat flux spectrum', fontsize=10)

# Expand spherical harmonic expansion onto new grid
heat_flux_grid = pysh.expand.MakeGridGLQ(cilm, zero, lmax, norm=2, csphase=-1)
d = ax[1,1].pcolormesh(Lon,Lat,heat_flux_grid, shading='auto')
ax[1,1].set_aspect('equal', 'box')
ax[1,1].set_xlabel('Longitude')
ax[1,1].set_ylabel('Latitude')
ax[1,1].set_title('Heat flux from SPH on GLQ points', fontsize=10)
fig.colorbar(d, ax=ax[1,1])

# Output spherical harmonic expansion
coeffs = pysh.SHCoeffs.from_array(cilm, normalization='schmidt', csphase=-1, lmax=lmax, name='heatflux')
coeffs.info()
coeffs.to_file(os.path.join(mypath, 'heat_flux_coeffs.txt'))
coeffs.to_netcdf(os.path.join(mypath, 'heat_flux_coeffs.nc'))
#plt.show()
plt.savefig("heat_flux.png", dpi=300)
plt.close()