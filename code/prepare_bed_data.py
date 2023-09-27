"""
Prepare the bed elevation data (Adrien's bed, GPR obs) to be in the right format and same vertical reference as Pleiades DEMs, i.e ellipsoid.
- Convert the bed provided by Adrien in matlab format into a Gtiff in UTM32.
- The vertical datum is shifted from IGN RAF09 to ellipsoid.
- The same is done for the GPR observations

NB: There is a vertical shift of 6.6 m and horizontal east/north shift of (4.04 m, -8.88 m ) between our 20190825 Pleiades DEM and the one provided by Adrien, meaning that all Pleiades DEM are probably not exactly aligned with all other datasets...
See old script code/archived/prepare_bed_data.py.

Author: Amaury Dehecq
"""

import os

import geopandas as gpd
import geoutils as gu
import matplotlib.pyplot as plt
import numpy as np
import rasterio as rio
import xdem
from scipy.io import loadmat

# - Load Adrien's bed - #
bed_elmer = loadmat("data/ice_thx/orig/Bed_Elmer.mat")
bed_x = bed_elmer["X"]
bed_y = bed_elmer["Y"]
bed_z = bed_elmer["Zbed"]

# Calculate image transform and convert to Raster
dx = bed_x[0, 1] - bed_x[0, 0]
dy = bed_y[1, 0] - bed_y[0, 0]
transform = rio.transform.from_origin(np.min(bed_x), np.max(bed_y), dx, dy)
bed_rst = gu.Raster.from_array(np.flipud(bed_z), transform=transform, crs="EPSG:27572")

# Reproject to UTM32 and save - only used for comparison with raw obs
# bed_rst_utm = bed_rst.reproject(dst_crs="EPSG:32632")
# os.makedirs("data/ice_thx/processed", exist_ok=True)
# print("Saving file data/ice_thx/processed/dem_bed_adrien_utm32_unshifted.tif")
# bed_rst_utm.save("data/ice_thx/processed/dem_bed_adrien_utm32_unshifted.tif")

# - Change vertical datum, from IGN90 to ellipsoid (presumably used for Pleiades) - #
# Does not work with original CRS, first need to convert to UTM32...
bed_rst_utm = bed_rst.reproject(dst_crs="EPSG:32632", dst_res=bed_rst.res)
bed_rst_vref = xdem.DEM(bed_rst_utm, vcrs="fr_ign_RAF09.tif")
bed_rst_vref.to_vcrs("Ellipsoid")
(bed_rst_vref - bed_rst_utm).show()
plt.title("Difference RAF09 - ellipsoid")
plt.show()

# Save
os.makedirs("data/ice_thx/processed", exist_ok=True)
print("Saving file data/ice_thx/processed/dem_bed_adrien_utm32.tif")
bed_rst_utm.save("data/ice_thx/processed/dem_bed_adrien_utm32.tif")

# - Calculate thickness - #
# Calculate thickness for fictive date
dem_med_date = xdem.DEM("output/dem_2017-02-15/dem_2017-02-15_interp_filtered.tif")
thick = dem_med_date.reproject(bed_rst_vref) - bed_rst_vref

# Mask few pixels with negative values and re-interpolate between valid pixels and edges set at 0
thick_fixed = thick.copy()
thick_fixed.set_mask(thick_fixed.data < 0)
thick_fixed.data[thick.data.mask] = 0
thick_fixed.data.mask[thick.data.mask] = False
thick_fixed.data = rio.fill.fillnodata(thick_fixed.data)
thick_fixed.data.mask = thick.data.mask

# Plot results
fig = plt.figure(figsize=(18, 6))
plt.suptitle("Calculated thickness")
ax1 = plt.subplot(131)
thick.show(ax=ax1, title="Full range")
ax2 = plt.subplot(132, sharex=ax1, sharey=ax1)
thick.show(ax=ax2, vmin=-30, vmax=30, cmap="RdYlBu", title="Zoom around 0")
ax3 = plt.subplot(133, sharex=ax1, sharey=ax1)
thick_fixed.show(ax=ax3, vmin=-30, vmax=30, cmap="RdYlBu", title="After fixing negative values")
plt.show()

# Save
print("Saving file data/ice_thx/processed/thickness_adrien_2017-02-15_utm32.tif")
thick_fixed.save("data/ice_thx/processed/thickness_adrien_2017-02-15_utm32.tif")


# -- Change vertical datum for the GPR observations -- #

# Load GPR data
bed_obs_ds = gpd.read_file("data/ice_thx/orig/zbed_arg_measured_UTM32N.shp")
zsurf = bed_obs_ds.PixelValue
bed_obs_z = bed_obs_ds.Field_3
bed_obs_x = bed_obs_ds.geometry.x
bed_obs_y = bed_obs_ds.geometry.y

# Create vertical coordinate reference systems, for datum change (same as for Adrien's bed)
# RAF09
raf09 = xdem.vcrs._vcrs_from_user_input("fr_ign_RAF09.tif")
raf09_ccrs = xdem.vcrs._build_ccrs_from_crs_and_vcrs(bed_obs_ds.crs, raf09)

# Ellipsoid
ellipsoid = xdem.vcrs._vcrs_from_user_input("Ellipsoid")
ellipsoid_ccrs = xdem.vcrs._build_ccrs_from_crs_and_vcrs(bed_obs_ds.crs, ellipsoid)

# Apply vertical transformation
bed_obs_z_shifted = xdem.vcrs._transform_zz(raf09_ccrs, ellipsoid_ccrs, bed_obs_z, bed_obs_y, bed_obs_z)

# Update pandas' dataset and save
bed_obs_ds_shifted = bed_obs_ds.copy()
bed_obs_ds_shifted.Field_3 = bed_obs_z_shifted
print("Saving file data/ice_thx/processed/zbed_arg_measured_UTM32N_shifted.shp")
bed_obs_ds_shifted.to_file("data/ice_thx/processed/zbed_arg_measured_UTM32N_shifted.shp")
