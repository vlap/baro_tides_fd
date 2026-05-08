import numpy as np
import netCDF4 as nc
import os

os.makedirs('data/ETOPO', exist_ok=True)
fn = 'data/ETOPO/ETOPO1_Ice_g_gmt4.grd'

# Create a small dummy ETOPO file
# 1 minute resolution would be 21601 x 10801, too big.
# I'll use 1 degree for the test file but name it the same.
nlons = 360
nlats = 180

with nc.Dataset(fn, 'w', format='NETCDF4') as ds:
    ds.title = 'Dummy ETOPO for Testing'
    ds.createDimension('x', nlons)
    ds.createDimension('y', nlats)
    
    x = ds.createVariable('x', 'f8', ('x',))
    y = ds.createVariable('y', 'f8', ('y',))
    z = ds.createVariable('z', 'i2', ('y', 'x'))
    
    x[:] = np.linspace(-180, 180, nlons)
    y[:] = np.linspace(-90, 90, nlats)
    
    x.long_name = 'Longitude'
    x.actual_range = np.array([-180.0, 180.0])
    x.units = 'degrees_east'
    
    y.long_name = 'Latitude'
    y.actual_range = np.array([-90.0, 90.0])
    y.units = 'degrees_north'
    
    z.long_name = 'Height'
    z.actual_range = np.array([-10000.0, 8000.0])
    z.units = 'meters'
    
    # Flat ocean with some land
    depth = np.full((nlats, nlons), -4000, dtype='i2')
    depth[70:110, 160:200] = 500 # "Island"
    z[:] = depth

print(f"Created dummy ETOPO file at {fn}")
