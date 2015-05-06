import numpy as np
import sys

resolution_radius = int(sys.argv[1])
resolution_lateral = int(sys.argv[2])

radius_bottom = 3391000
radius_surface = 6371000
model_depth = radius_surface - radius_bottom

# western boundary
header = "x,y,z"

coords = []
radiuss = np.arange(radius_bottom,radius_surface,model_depth/resolution_radius)
lats = np.arange(0,np.pi/2,np.pi/(2*resolution_lateral))
phi = 0

for theta in lats:
    for r in radiuss:
        coords.append([r*np.cos(phi)*np.sin(theta),
                       r*np.sin(phi)*np.sin(theta),
                       r*np.cos(theta)])
        
np_coords = np.array(coords)
np.savetxt("coords_west.csv", np_coords, delimiter=',', header=header, footer="", comments='')    

    
# eastern boundary
header = "x,y,z"
coords = np.zeros((resolution_radius,resolution_lateral,3))

coords = []
radiuss = np.arange(radius_bottom,radius_surface,model_depth/resolution_radius)
lats = np.arange(0,np.pi/2,np.pi/(2*resolution_lateral))
phi = np.pi/2

for theta in lats:
    for r in radiuss:
        coords.append([r*np.cos(phi)*np.sin(theta),
                       r*np.sin(phi)*np.sin(theta),
                       r*np.cos(theta)])
        
np_coords = np.array(coords)
np.savetxt("coords_east.csv", np_coords, delimiter=',', header=header, footer="", comments='')        

# southern boundary
header = "x,y,z"

coords = []
theta = np.pi / 2
radiuss = np.arange(radius_bottom,radius_surface,model_depth/resolution_radius)
lons = np.arange(0,np.pi/2,np.pi/(2*resolution_lateral))

for phi in lons:
    for r in radiuss:
        coords.append([r*np.cos(phi)*np.sin(theta),
                       r*np.sin(phi)*np.sin(theta),
                       r*np.cos(theta)])
        
np_coords = np.array(coords)
np.savetxt("coords_south.csv", np_coords, delimiter=',', header=header, footer="", comments='')        

# inner boundary
header = "x,y,z"

coords = []
theta = np.arange(0,np.pi/2,np.pi/(2*resolution_lateral))
radiuss = radius_bottom
lons = np.arange(0,np.pi/2,np.pi/(2*resolution_lateral))

for theta in lats:
    for phi in lons:
        coords.append([r*np.cos(phi)*np.sin(theta),
                       r*np.sin(phi)*np.sin(theta),
                       r*np.cos(theta)])

np_coords = np.array(coords)
np.savetxt("coords_inner.csv", np_coords, delimiter=',', header=header, footer="", comments='')        

# outer boundary
header = "x,y,z"

coords = []
theta = np.arange(0,np.pi/2,np.pi/(2*resolution_lateral))
radiuss = radius_surface
lons = np.arange(0,np.pi/2,np.pi/(2*resolution_lateral))

for theta in lats:
    for phi in lons:
        coords.append([r*np.cos(phi)*np.sin(theta),
                       r*np.sin(phi)*np.sin(theta),
                       r*np.cos(theta)])
        
np_coords = np.array(coords)
np.savetxt("coords_outer.csv", np_coords, delimiter=',', header=header, footer="", comments='')        


