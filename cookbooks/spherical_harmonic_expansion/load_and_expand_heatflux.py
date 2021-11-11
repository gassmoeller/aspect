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
from scipy.interpolate import griddata


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
bottom_heat_flux = heat_flux[heat_flux.r < 4e6].to_numpy()


xi = np.linspace(np.min(bottom_heat_flux[:,0]), np.max(bottom_heat_flux[:,0]), 101)
yi = np.linspace(np.min(bottom_heat_flux[:,1]), np.max(bottom_heat_flux[:,1]), 101)
zi = np.linspace(np.min(bottom_heat_flux[:,2]), np.max(bottom_heat_flux[:,2]), 2)

# Generate data interpolation grid. 
X, Y, Z = np.meshgrid(xi,yi,zi)

# Interpolate heat flux onto new grid
heat_flux_structured = griddata((bottom_heat_flux[:,0:3]), bottom_heat_flux[:,3], (X,Y,Z), method='nearest')

# Plot strain-rate field
fig = plt.figure()
ax0 = fig.gca(projection='3d')

c = ax0.scatter(bottom_heat_flux[:,0],bottom_heat_flux[:,1],bottom_heat_flux[:,2],c=bottom_heat_flux[:,3])
ax0.set_aspect('equal', 'box')
ax0.set_xlabel('Horizontal Position (m)')
ax0.set_ylabel('Vertical Position (m)')
ax0.set_title('Strain Rate Second Invariant (1/s)', fontsize=10)
fig.colorbar(c, ax=ax0)
plt.show()
#plt.savefig("strain_rate_field.png", dpi=300)
#plt.close()