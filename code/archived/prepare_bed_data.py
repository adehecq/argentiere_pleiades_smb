"""
Prepare the bed elevation data (Adrien's bed, GPR obs) to be in the right format and coregistered to the Pleiades DEMs.
- Convert the bed provided by Adrien in .dat format into a Gtiff in UTM32.
- The DEM is coregistered to the 2019-08-25 Pleiades DEM outside Argentiere extent to be coherent with all other Pleiades DEMs used in this study.
There is a ~60 m vertical shift between both. Most of it is due to the fact that Adrien's data is expressed w.r.t. to the geoid (most likely IGN RAF09) while Luc's Pleiades DEMs are w.r.t. the ellipsoid, which causes a ~53 m difference. The rest may be due to the fact that Luc's Pleiades DEMs are only coregistered relative to each other, but most likely free floating.
This means the Pleiades elevations are not true elevation!
- Apply the same vertical shift to the GPR observations.

Author: Amaury Dehecq
"""

import os

import geopandas as gpd
import geoutils as gu
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import rasterio as rio
import xdem
from scipy.io import loadmat
from shapely.geometry import Polygon


# --- Read input data --- #

# - Load Adrien's bed - #
dem_bed_adrien = pd.read_csv("data/IceTHX/orig/DEMbed_argentiere.dat", delimiter=" ", names=["X", "Y", "Z"])
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

# - Load Argentiere contour, to crop outside Elmer boundary - #
arg_contour = loadmat("data/IceTHX/orig/contour_arg2003.mat")
arg_contour_x, arg_contour_y = arg_contour["A"].T[0], arg_contour["A"].T[1]

# Add missing 2 in front of all y values (due to surfer)
arg_contour_y = arg_contour_y + 2e6

# Convert to Vector
geom = Polygon(np.transpose((arg_contour_x, arg_contour_y)))
arg_contour_vect = gu.Vector(gpd.GeoDataFrame(geometry=[geom], crs="EPSG:27572"))

# Plot
ax = plt.subplot(111)
bed_rst.show(ax=ax)
arg_contour_vect.show(ec="k", fc="none", ax=ax)
plt.show()

# Convert to UTM32
arg_contour_vect = arg_contour_vect.reproject(dst_crs="EPSG:32632")

# Reproject to UTM32 and save - only used for comparison with raw obs
# bed_rst_utm = bed_rst.reproject(dst_crs="EPSG:32632")
# os.makedirs("output/IceTHX", exist_ok=True)
# bed_rst_utm.save("output/IceTHX/dem_bed_adrien_utm32_unshifted.tif")

# optionally change vertical datum, from IGN90 to ellipsoid (presumably used for Pleiades)
# Does not work with original CRS, first need to convert to UTM32...
bed_rst_utm = bed_rst.reproject(dst_crs="EPSG:32632")
bed_rst_vref = xdem.DEM(bed_rst_utm, vcrs="fr_ign_RAF09.tif")
bed_rst_vref.to_vcrs("Ellipsoid")
(bed_rst_vref - bed_rst_utm).show()
plt.title("Difference RAF09 - ellipsoid")
plt.show()

# Coregister Adrien's bed to Pleiades 2019-08-25, which is the DEM used outside of Argentiere
# Use Nuth & Kaab + vertical shift
print("\nCoregistering...")
pleiades_dem = xdem.DEM("data/dh/DEMs/20190825_DEM_4m_shift_H-V_clip.tif")
bed_rst_coreg, coreg_method, stats, inlier_mask = xdem.coreg.dem_coregistration(
    bed_rst_vref, pleiades_dem, verbose=True, plot=True, shp_list=["data/GIS/Argentiere_outline_201809.shp"]
)

# Plot
ddem_coreg = pleiades_dem - bed_rst_coreg.reproject(pleiades_dem)
ax = plt.subplot(111)
ddem_coreg.show(vmin=-5, vmax=5, cmap="coolwarm_r", ax=ax)
arg_contour_vect.show(fc="none", ec="k", ax=ax)
plt.title("Pleiades 2019-08-25 - Adrien's bed")
plt.show()

# Reproject to UTM32 and save
bed_rst_utm = bed_rst_coreg.reproject(dst_crs="EPSG:32632", silent=True)
os.makedirs("data/IceTHX/processed", exist_ok=True)
bed_rst_utm.save("data/IceTHX/processed/dem_bed_adrien_utm32.tif")

# Calculate thickness for fictive date and save
# TODO: check issue with interpolated DEM (issue #2) and re-run
dem_med_date = xdem.DEM("output/dem_2017-02-15/dem_2017-02-15_interp_filtered.tif")
thick = dem_med_date - bed_rst_coreg.reproject(dem_med_date)
thick.data[thick.data < 0] = 0  # Filter a few negative thickness values
thick.save("data/IceTHX/processed/thickness_adrien_2017-02-15_utm32.tif")

# -- Apply same transformation to the GPR observations -- #

# Load GPR data
bed_obs_ds = gpd.read_file("data/IceTHX/orig/zbed_arg_measured_UTM32N.shp")
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

bed_obs_z_shifted = xdem.vcrs._transform_zz(raf09_ccrs, ellipsoid_ccrs, bed_obs_z, bed_obs_y, bed_obs_z)

# Apply also vertical shift estimated during coreg
vshift = coreg_method.pipeline[0]._meta["vshift"] + coreg_method.pipeline[1]._meta["vshift"]
bed_obs_z_shifted += vshift

# Update pandas' dataset and save
bed_obs_ds_shifted = bed_obs_ds.copy()
bed_obs_ds_shifted.Field_3 = bed_obs_z_shifted
bed_obs_ds_shifted.to_file("data/IceTHX/processed/zbed_arg_measured_UTM32N_shifted.shp")
