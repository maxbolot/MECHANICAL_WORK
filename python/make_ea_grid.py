import numpy as np
from pathlib import Path

# Grid dimensions matching the 0.25-degree setup
nx = 1440
ny = 720

# Longitudes: standard 0.25 degree spacing
xvals = np.linspace(0.125, 360 - 0.125, nx)
xbounds = np.zeros((nx, 2))
xbounds[:, 0] = np.linspace(0, 360 - 0.25, nx)
xbounds[:, 1] = np.linspace(0.25, 360, nx)

# Latitudes: Equal-Area spacing (uniform in sin(lat))
sin_lat_bounds = np.linspace(-1.0, 1.0, ny + 1)
# Cell centers are the midpoints in sine space to ensure exact area conservation
sin_lat_vals = (sin_lat_bounds[:-1] + sin_lat_bounds[1:]) / 2.0

# Convert back to degrees for CDO
yvals = np.degrees(np.arcsin(sin_lat_vals))
ybounds = np.zeros((ny, 2))
ybounds[:, 0] = np.degrees(np.arcsin(sin_lat_bounds[:-1]))
ybounds[:, 1] = np.degrees(np.arcsin(sin_lat_bounds[1:]))

# Write the CDO grid description file under output/grids.
project_root = Path(__file__).resolve().parent.parent
output_dir = project_root / 'output' / 'grids'
output_dir.mkdir(parents=True, exist_ok=True)
output_path = output_dir / 'ea_25km_grid.txt'

with output_path.open('w') as f:
    f.write("gridtype  = lonlat\n")
    f.write(f"gridsize  = {nx * ny}\n")
    f.write(f"xsize     = {nx}\n")
    f.write(f"ysize     = {ny}\n")
    f.write("xname     = lon\n")
    f.write("xunits    = degrees_east\n")
    f.write("yname     = lat\n")
    f.write("yunits    = degrees_north\n")

    f.write("xvals     = " + " ".join([f"{x:.5f}" for x in xvals]) + "\n")
    
    f.write("xbounds   = ")
    for i in range(nx):
        f.write(f"{xbounds[i,0]:.5f} {xbounds[i,1]:.5f} ")
    f.write("\n")

    f.write("yvals     = " + " ".join([f"{y:.5f}" for y in yvals]) + "\n")
    
    f.write("ybounds   = ")
    for i in range(ny):
        f.write(f"{ybounds[i,0]:.5f} {ybounds[i,1]:.5f} ")
    f.write("\n")

print(f"Equal-area grid description saved to '{output_path}'")