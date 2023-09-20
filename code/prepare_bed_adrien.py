"""
Convert the bed provided by Adrien in .dat format into a Gtiff in UTM32.
The DEM is coregistered to the 2019-08-25 Pleiades DEM outside Argentiere extent to be coherent with all other Pleiades DEMs used in this study.
There is a ~60 m vertical shift between both. Most of it is due to the fact that Adrien's data is expressed w.r.t. to the geoid (most likely IGN RAF09) while Luc's Pleiades DEMs are w.r.t. the ellipsoid, which causes a ~53 m difference. The rest may be due to the fact that Luc's Pleiades DEMs are only coregistered relative to each other, but most likely free floating.
This means the Pleiades elevations are not true elevation!

Author: Amaury Dehecq
"""

import os

import geoutils as gu
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import rasterio as rio
import xdem

# Read input data
dem_bed_adrien = pd.read_csv("data/IceTHX/DEMbed_argentiere.dat", delimiter=" ", names=["X", "Y", "Z"])
bed_x = dem_bed_adrien.X.values
bed_y = dem_bed_adrien.Y.values
bed_z = dem_bed_adrien.Z.values

# Add missing 2 in front of all y values (due to surfer)
bed_y = bed_y + 2e6

# Reshape to 2D
ncols = len(np.unique(bed_x))
nrows = len(np.unique(bed_y))
bed_x = bed_x.reshape((nrows, ncols))
bed_y = bed_y.reshape((nrows, ncols))
bed_z = bed_z.reshape((nrows, ncols))

# Calculate image transform and convert to Raster
dx = bed_x[0, 1] - bed_x[0, 0]
dy = bed_y[1, 0] - bed_y[0, 0]
transform = rio.transform.from_origin(np.min(bed_x), np.max(bed_y), dx, dy)
bed_rst = gu.Raster.from_array(np.flipud(bed_z), transform=transform, crs="EPSG:27572")

# Reproject to UTM32 and save - only used for comparison with raw obs
# bed_rst_utm = bed_rst.reproject(dst_crs="EPSG:32632")
# os.makedirs("outputs/IceTHX", exist_ok=True)
# bed_rst_utm.save("outputs/IceTHX/dem_bed_adrien_utm32_unshifted.tif")

# optionally change vertical datum, from IGN90 to ellipsoid (presumably used for Pleiades)
# Does not work with original CRS, first need to convert to UTM32...
# Not used here as coregistration seem to work well without taking the vertical datum into account
# bed_rst_utm = bed_rst.reproject(dst_crs="EPSG:32632")
# bed_rst_test = xdem.DEM(bed_rst_utm)
# bed_rst_test.to_vcrs(src_vcrs="Ellipsoid", dst_vcrs="fr_ign_RAF09.tif")
# (bed_rst_test - bed_rst_utm).show()
# plt.title("Difference RAF09 - ellipsoid")
# plt.show()

# Coregister Adrien's bed to Pleiades 2019-08-25, which is the DEM used outside of Argentiere
# Use Nuth & Kaab + vertical shift
print("\nCoregistering...")
pleiades_dem = xdem.DEM("data/DEMs/20190825_DEM_4m_nuth_x-1.60_y+6.74_z+9.61-0.02_align_corr.tif")
bed_rst_coreg, coreg_method, stats, inlier_mask = xdem.coreg.dem_coregistration(
    bed_rst, pleiades_dem, verbose=True, plot=True, shp_list=["data/GIS/Argentiere_outline_201809.shp"]
)

# Plot
ddem_coreg = pleiades_dem - bed_rst_coreg.reproject(pleiades_dem)
ddem_coreg.show(vmin=-5, vmax=5, cmap="coolwarm_r")
plt.title("Pleiades 2019-08-25 - Adrien's bed")
plt.show()

# Reproject to UTM32 and save
bed_rst_utm = bed_rst_coreg.reproject(dst_crs="EPSG:32632")
bed_rst_utm.save("outputs/IceTHX/dem_bed_adrien_utm32.tif")

# Calculate thickness for fictive date and save
dem_med_date = xdem.DEM("outputs/dh/dem_interp_2017_02_15.tif")
thick = dem_med_date - bed_rst_coreg.reproject(dem_med_date)
thick.data[thick.data < 0] = 0  # Filter a few negative thickness values
thick.save("outputs/IceTHX/thickness_adrien_2017-02-15_utm32.tif")
