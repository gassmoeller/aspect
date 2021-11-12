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

import vtk as vtk; from vtk.util import numpy_support
import scipy.interpolate as interpolate


def load_vtk(filename):
  # Load vtu data (pvtu directs to vtu files)
  reader = vtk.vtkXMLPUnstructuredGridReader()
  reader.SetFileName(filename)
  reader.Update()

  # Get the coordinates of nodes in the mesh
  nodes_vtk_array= reader.GetOutput().GetPoints().GetData()

  # Convert nodal vtk data to a numpy array
  nodes_numpy_array = vtk.util.numpy_support.vtk_to_numpy(nodes_vtk_array)

  # Extract x, y and z coordinates from numpy array 
  x,y = nodes_numpy_array[:,0] , nodes_numpy_array[:,1]

  # Determine the number of scalar fields contained in the .pvtu file
  number_of_fields = reader.GetOutput().GetPointData().GetNumberOfArrays()

  # Determine the name of each field and place it in an array.
  field_names = []
  for i in range(number_of_fields):
    field_names.append(reader.GetOutput().GetPointData().GetArrayName(i))

  # Determine the index of the field strain_rate
  idx = field_names.index("heat_flux_map")

  # Extract values of strain_rate
  field_vtk_array = reader.GetOutput().GetPointData().GetArray(idx)
  strain_rate     = numpy_support.vtk_to_numpy(field_vtk_array)
  return x, y, strain_rate

# Relative (from current directory) path to the
# solution data folder.
mypath = os.path.dirname(__file__)
output_directory = os.path.join(mypath, 'output-spherical_harmonic_expansion')
solution_directory = os.path.join(output_directory, 'solution')

# Time step to analyze. This number corresponds to 
# the numbering of .pvtu files. In this instance
# '00000.0213' is the last nonlinear iteration
# of the first time step, but note that the 
# number of nonlinear iterations may vary slightly
# between machines, or significantly when additional
# parameters are altered.
number = '00000'

vtu_filename = os.path.join(solution_directory, 'solution-' + number + '.pvtu')
text_filename = os.path.join(output_directory, 'heat_flux.' + number)

#heat_flux = np.genfromtxt(text_filename, delimiter=' ', skip_header=1)
heat_flux = pd.read_csv(text_filename, delimiter=' ', skiprows=1, header=0, names=['x','y','z','heat_flux'])

heat_flux['r'] = heat_flux.apply(lambda row: np.sqrt(row.x*row.x+row.y*row.y+row.z*row.z), axis = 1)
heat_flux['latitude'] = heat_flux.apply(lambda row: 90 - np.arccos(row.z/row.r) * 180 / np.pi, axis = 1)
heat_flux['longitude'] = heat_flux.apply(lambda row: np.arctan2(row.y,row.x) * 180 / np.pi, axis = 1)

bottom_heat_flux = heat_flux[heat_flux.r < 4e6].to_numpy()


lat_i = np.linspace(-90, 90, 51)
long_i = np.linspace(-180, 180, 51)
# Generate data interpolation grid. 
Lat, Lon = np.meshgrid(lat_i,long_i)

# create interpolation function
#function = interpolate.interp2d(bottom_heat_flux[:,5], bottom_heat_flux[:,6], bottom_heat_flux[:,3], kind='linear')
#function = interpolate.interp2d(bottom_heat_flux[:,5], bottom_heat_flux[:,6], bottom_heat_flux[:,3], kind='linear')
#bottom_heat_flux_interpolated = function(lat_i, long_i)

# Plot heat flux field
fig = plt.figure()

# Cartesian scatter
# ax0 = fig.gca(projection='3d')
# c = ax0.scatter(bottom_heat_flux[:,0],bottom_heat_flux[:,1],bottom_heat_flux[:,2],c=bottom_heat_flux[:,3])

# Geographic scatter
ax0 = fig.add_subplot(211)
c = ax0.scatter(bottom_heat_flux[:,6],bottom_heat_flux[:,5],c=bottom_heat_flux[:,3])

#ax1 = fig.add_subplot(212)
#d = ax1.scatter(Lon, Lat,c=bottom_heat_flux_interpolated)

# Geographic contour
# Interpolate heat flux onto new grid
heat_flux_structured = interpolate.griddata((bottom_heat_flux[:,5:]), bottom_heat_flux[:,3], (Lat,Lon), method='cubic')
ax1 = fig.add_subplot(212)
d = ax1.pcolormesh(Lon,Lat,heat_flux_structured)


ax0.set_aspect('equal', 'box')
ax0.set_xlabel('Longitude')
ax0.set_ylabel('Latitude')
ax0.set_title('Heat flux on cell centers', fontsize=10)

ax1.set_aspect('equal', 'box')
ax1.set_xlabel('Longitude')
ax1.set_ylabel('Latitude')
ax1.set_title('Heat flux on GLQ points', fontsize=10)

fig.colorbar(c, ax=ax0)
fig.colorbar(d, ax=ax1)

plt.show()
plt.savefig("heat_flux.png", dpi=300)
plt.close()