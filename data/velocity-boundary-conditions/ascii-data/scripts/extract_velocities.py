import numpy as np
import sys
import glob

from paraview import simple
from paraview import servermanager

def generate_sphere_octant_coordinates(resolution_radius, resolution_lateral, boundary_name, filename="coords"):
    
    radius_bottom = 3491000
    radius_surface = 6171000
    
    coords = []
    scoords = []
    radiuss = np.linspace(radius_bottom,radius_surface,resolution_radius)
    lats = np.linspace(0,np.pi/2,resolution_lateral)
    lons = np.linspace(0,np.pi/2,resolution_lateral)

    if boundary_name == "west":
        lons = 0,
    elif boundary_name == "east":
        lons = np.pi/2,
    elif boundary_name == "south":
        lats = np.pi/2,
    elif boundary_name == "inner":
        radiuss = radius_bottom,
    elif boundary_name == "outer":
        radiuss = radius_surface,
                                
    for theta in lats:
        for phi in lons:
            for r in radiuss:
                coords.append([r*np.cos(phi)*np.sin(theta),
                               r*np.sin(phi)*np.sin(theta),
                               r*np.cos(theta)])
                scoords.append([r,phi,theta])
     
    header = "x,y,z"   
    np_coords = np.array(coords)
    np.savetxt(filename  + "_" + boundary_name + ".csv", np_coords, delimiter=',', header=header, footer="", comments='')
    return filename, np.array(scoords)



""" Usage: This script expects 5 input arguments:
Input model file (*.pvd)
Output data name (is extended with the boundary's name and saved as csv and vtp)
Radial resolution (number of points)
Lateral resolution (number of points)
Property (velocity, temperature or the name of another existing property in the input model
"""

input_data = str(sys.argv[1])
output_data = str(sys.argv[2])

resolution_radius = int(sys.argv[3])
resolution_lateral = int(sys.argv[4])

properties = ["Points0","Points1","Points2"]

if sys.argv[5] == "velocity":
    properties += ["velocity0","velocity1","velocity2"] 
elif sys.argv[5] == "temperature":
    properties += ["T",] 
else:
    properties += [sys.argv[5],] 

resolutions = resolution_radius,resolution_lateral,resolution_lateral

names = "west","east","south","inner","outer"

model = simple.OpenDataFile(input_data)
model.UpdatePipeline()

transformed_model = simple.Transform()
transformed_model.Input = model
transformed_model.Transform.Rotate = [0,0,0]

for name in names:
    input_coord_filename,scoords = generate_sphere_octant_coordinates(resolution_radius, resolution_lateral, name)
    input_coord_filename += "_"
    reader = simple.OpenDataFile(input_coord_filename + name + ".csv")

    tabletopoints = servermanager.filters.TableToPoints()
    tabletopoints.XColumn = 'x'
    tabletopoints.YColumn = 'y'
    tabletopoints.ZColumn = 'z'
    tabletopoints.Input = reader
    tabletopoints.UpdatePipeline()
    
    resample = simple.ResampleWithDataset()
    resample.Input = transformed_model
    resample.Source = tabletopoints
    
    writer_vtp = simple.XMLPolyDataWriter()
    writer_vtp.Input = resample
    writer_vtp.Writealltimestepsasfileseries = 1
    writer_vtp.FileName = output_data + "_" + name + ".vtp"
    writer_vtp.UpdatePipeline()
    
    writer = simple.DataSetCSVWriter()
    writer.Input = resample
    writer.WriteAllTimeSteps = 1
    writer.FileName = output_data "_" + name + ".csv"
        
    writer.UpdatePipeline()
    
    files = glob.glob(output_data "_" + name + ".*.csv")
    
    for file in files:
        reformat_data = np.genfromtxt(file, delimiter=',',names=True)
        reformatted_data = reformat_data[properties]
        reformatted_data = reformatted_data.view((float, len(reformatted_data.dtype.names)))

        if name == "west":
            np.savetxt(file,np.hstack((scoords[:,[0,2]],reformatted_data[:,3:])), fmt='%.10g', delimiter=' ', header="# POINTS: " + str(resolutions[0]) + " " + str(resolutions[2]), footer="", comments='')
        if name == "east":
            np.savetxt(file,np.hstack((scoords[:,[0,2]],reformatted_data[:,3:])), fmt='%.10g', delimiter=' ', header="# POINTS: " + str(resolutions[0]) + " " + str(resolutions[2]), footer="", comments='')           
        if name == "south":
            np.savetxt(file,np.hstack((scoords[:,[0,1]],reformatted_data[:,3:])), fmt='%.10g', delimiter=' ', header="# POINTS: " + str(resolutions[0]) + " " + str(resolutions[1]), footer="", comments='')        
        if name == "inner":
            np.savetxt(file,np.hstack((scoords[:,[1,2]],reformatted_data[:,3:])), fmt='%.10g', delimiter=' ', header="# POINTS: " + str(resolutions[1]) + " " + str(resolutions[2]), footer="", comments='')
        if name == "outer":
            np.savetxt(file,np.hstack((scoords[:,[1,2]],reformatted_data[:,3:])), fmt='%.10g', delimiter=' ', header="# POINTS: " + str(resolutions[1]) + " " + str(resolutions[2]), footer="", comments='')            
