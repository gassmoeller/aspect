import numpy as np
import sys
import glob

from paraview import simple
from paraview import servermanager

input_data = str(sys.argv[1])
input_coords = str(sys.argv[2])
output_data = str(sys.argv[3])

names = "west","east","south","inner","outer"

model = simple.OpenDataFile(input_data)
model.UpdatePipeline()

for name in names:
    reader = simple.OpenDataFile(input_coords + name + ".csv")
    
    tabletopoints = servermanager.filters.TableToPoints()
    tabletopoints.XColumn = 'x'
    tabletopoints.YColumn = 'y'
    tabletopoints.ZColumn = 'z'
    tabletopoints.Input = reader
    tabletopoints.UpdatePipeline()
    
    resample = simple.ResampleWithDataset()
    resample.Input = model
    resample.Source = tabletopoints
    
    writer_vtp = simple.XMLPolyDataWriter()
    writer_vtp.Input = resample
    writer_vtp.Writealltimestepsasfileseries = 1
    writer_vtp.FileName = output_data + name + ".vtp"
    writer_vtp.UpdatePipeline()
    
    writer = simple.DataSetCSVWriter()
    writer.Input = resample
    writer.WriteAllTimeSteps = 1
    writer.FileName = output_data + name + ".csv"
        
    writer.UpdatePipeline()
    
    files = glob.glob(output_data + name + ".*.csv")
    
    for file in files:
        reformat_data = np.genfromtxt(file, delimiter=',', names=True)    
        reformatted_data = reformat_data[["Points0","Points1","Points2","velocity0","velocity1","velocity2"]]
        np.savetxt(file,reformatted_data, delimiter=' ', header="", footer="", comments='')        
        
