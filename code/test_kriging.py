"""
Attempt to extrapolate the bed obsrvations to the whole glacier using ordinary kriging with a deterministic thickness model, based on slope and velocity. The steps are the following:
- load thickness observations, slope and velocity
- bin the thickness obs by slope and velocity
- use this empirical relationship to calculate a modeled thickness over the whole glacier
- remove the model to the observations
- fit a variogram model to teh residuals
- run ordinary kriging to teh residuals
 
TODO:
- smooth slope
- fit a polynomial relationship between thickness, slope and velocity

Author: Amaury Dehecq
"""

import os
import time

import geopandas as gpd
import geoutils as gu
import gstools as gs
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xdem
from scipy.signal import fftconvolve
from sklearn import linear_model


def mean_filter_nan(array_in, kernel_width):
    """
    Apply a mean filter to an array array_in, on a square kernel of size kernel_width.
    The method account for NaN values in array_in and return the count of valid pixels used in the convolution.
    Inputs:
    - array_in - 2D array, the array on which the mean filter is to be applied
    - kernel_width - int, the size of the kernel
    Outputs:
    - array of same shape as array_in, containing the result of the convolution
    - array of same shape as array_in, containing the count of valid pixels averaged
    """
    # Copy of array_in where nan are replaced by 0
    array_zero = array_in.copy()
    array_zero[np.isnan(array_in)] = 0

    # Convolution with square kernel
    kernel = np.ones((kernel_width, kernel_width), dtype="f")
    array_conv = fftconvolve(array_zero, kernel, mode="same")

    # Array containing 0 where array_in is NaN, 1 otherwise
    array_mask = 0 * array_in.copy() + 1
    array_mask[np.isnan(array_in)] = 0

    # Convolution -> count of valid pixels in the kernel
    array_count = fftconvolve(array_mask, kernel, mode="same")

    # Final convolution result with NaN accounted for
    array_out = array_conv / array_count

    return array_out, array_count


# -- parameters -- #
outdir = os.path.join("output", "thickness_kriging")


# --- Load input data --- #

# Load Argentiere glacier outline
arg_outline = gu.Vector("data/gis/outline_pleiades_2020-09-08.shp")

# Load thickness obs
zbed_ds = gpd.read_file("data/ice_thx/processed/zbed_arg_measured_UTM32N_shifted.shp")
zsurf = zbed_ds.PixelValue
zbed = zbed_ds.Field_3
zbed_x = zbed_ds.geometry.x
zbed_y = zbed_ds.geometry.y
Transect = zbed_ds.Transect.astype(float)

# remove some observations 

# Calculate thickness and remove bad values
zsurf = zsurf.replace(0, np.nan)
H_obs = zsurf - zbed
H_obs[H_obs < 0] = np.nan

valid_obs = np.where(np.isfinite(H_obs))
H_obs = H_obs.values[valid_obs]
zbed_x = zbed_x.values[valid_obs]
zbed_y = zbed_y.values[valid_obs]

# Load velocity data
velocity = gu.Raster("output/velocity/v_mean_after_30_15°.tif")

# Load Pleiades fictive DEM
dem = gu.Raster("output/dh_results/meanDEM-2017_02_15.tif")

# Crop DEM and velocity to Argentiere extent, reproject on same grid as velocity
left, bottom, right, top = list(arg_outline.bounds)
velocity.crop([left - 200, bottom - 200, right + 200, top + 200])
# velocity.crop(arg_outline)
dem = dem.reproject(velocity, resampling="average")

# Calculate slope
slope = xdem.terrain.slope(dem)

# Mask pixels outside glaciers and smooth
sm_length = 400
gl_mask = arg_outline.create_mask(slope)
slope_arg = np.where(gl_mask.data & ~slope.data.mask, slope.data.data, np.nan)
slope_sm, count = mean_filter_nan(slope_arg, int(sm_length / slope.res[0]))
slope_sm[~gl_mask.data] = np.nan
slope_sm[count <= 1] = np.nan
slope_sm = slope.copy(new_array=slope_sm)

# Calculate tangent of slope, for the V vs tau relationship
slope_tan = np.tan(slope_sm * np.pi / 180.)

# Plot input data
fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(18, 6))

arg_outline.show(fc="none", ec="k")
sc = axes[0].scatter(zbed_x, zbed_y, c=H_obs)
cb = plt.colorbar(sc)
cb.set_label("Thickness (m)")

velocity.show(ax=axes[1], cbar_title="Velocity (m/yr)", vmax=150)
arg_outline.show(fc="none", ec="k", ax=axes[1])

slope_sm.show(ax=axes[2], cbar_title="Slope (degree)")
arg_outline.show(fc="none", ec="k", ax=axes[2])

plt.tight_layout()
plt.show()

# --- Prepare data --- #

# Extract slope at location of thickness obs
slope_obs = slope_tan.value_at_coords(zbed_x, zbed_y)
slope_obs[slope_obs < 0] = np.nan  # remove values outside glacier outlines
plt.hist(slope_obs)

# Extract velocity at location of thickness obs
vel_obs = velocity.value_at_coords(zbed_x, zbed_y)
plt.hist(vel_obs)

# --- Model H as a function of slope and velocity --- #

# Calculate mean H as a function of slope and velocity
# slope_bins = np.tan(np.array([0, 2, 4, 6, 8, 10, 20, 40, np.max(slope)]) * np.pi / 180)
# velocity_bins = [0, 10, 20, 30, 40, 50, np.max(vel_obs)]
slope_bins = np.percentile(slope_obs[np.isfinite(slope_obs)], [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100])
velocity_bins = np.percentile(vel_obs[np.isfinite(vel_obs)], [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100])

df = xdem.spatialstats.nd_binning(
    values=H_obs,
    list_var=[slope_obs, vel_obs],
    list_var_names=["slope", "velocity"],
    statistics=["count", "mean", "median"],
    list_var_bins=[slope_bins, velocity_bins],
)

# Minimum number of samples in each bin to be considered valid
min_count = 3

# Plot relationship against each variable
fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(12, 6))

xdem.spatialstats.plot_1d_binning(
    df,
    var_name="slope",
    statistic_name="mean",
    label_var="Slope tangent",
    label_statistic="Mean thickness (m)",
    ax=axes[0],
    min_count=min_count,
)

xdem.spatialstats.plot_1d_binning(
    df,
    var_name="velocity",
    statistic_name="mean",
    label_var="Velocity (m/yr)",
    label_statistic="Mean thickness (m)",
    ax=axes[1],
    min_count=min_count,
)

plt.tight_layout()
plt.show()

# Plot 2D relationship -> some data gaps
xdem.spatialstats.plot_2d_binning(
    df,
    var_name_1="slope",
    var_name_2="velocity",
    statistic_name="mean",
    label_var_name_1="Slope tangent",
    label_var_name_2="Velocity (m/yr)",
    label_statistic="Mean thickness (m)",
    min_count=min_count,
)
plt.tight_layout()
plt.show()

# Option 1 - Interpolate model for H
# H_model = xdem.spatialstats.interp_nd_binning(
#     df, list_var_names=["slope", "velocity"], statistic="mean", min_count=min_count
# )

# - Option 2 - Polynomial fit of relationship - #

# Test 1D, but not used
# Vs slope
df_sub = df[np.logical_and(df.nd == 1, np.isfinite(pd.IntervalIndex(df["slope"]).mid))].copy()
xx1 = np.log10(pd.IntervalIndex(df_sub["slope"][df_sub["count"].values > 5]).mid)
yy1 = np.log10(df_sub["median"][df_sub["count"].values > 5])
poly_slope = np.polyfit(xx1, yy1, deg=1)

# Vs velocity
df_sub = df[np.logical_and(df.nd == 1, np.isfinite(pd.IntervalIndex(df["velocity"]).mid))].copy()
xx2 = np.log10(pd.IntervalIndex(df_sub["velocity"][df_sub["count"].values > 5]).mid)
yy2 = np.log10(df_sub["median"][df_sub["count"].values > 5])
poly_vel = np.polyfit(xx2, yy2, deg=1)

plt.figure(figsize=(12, 6))

plt.subplot(121)
x_reg = [np.min(xx1), np.max(xx1)]
y_reg = np.polyval(poly_slope, x_reg)

plt.plot(xx1, yy1, "r+")
plt.plot(x_reg, y_reg, "k-")
plt.xlabel("Log(tan(slope))")
plt.ylabel("Log(thickness)")

plt.subplot(122)
x_reg = [np.min(xx2), np.max(xx2)]
y_reg = np.polyval(poly_vel, x_reg)

plt.plot(xx2, yy2, "r+")
plt.plot(x_reg, y_reg, "k-")
plt.xlabel("Log(velocity)")
plt.ylabel("Log(thickness)")

plt.show()

# Multivariate regression
df_sub = df[np.logical_and(df.nd == 2, np.isfinite(pd.IntervalIndex(df["slope"]).mid), np.isfinite(pd.IntervalIndex(df["velocity"]).mid))].copy()
xx1 = np.log10(pd.IntervalIndex(df_sub["slope"][df_sub["count"].values > 5]).mid)
xx2 = np.log10(pd.IntervalIndex(df_sub["velocity"][df_sub["count"].values > 5]).mid)
yy = np.log10(df_sub["median"][df_sub["count"].values > 5])

lin_reg = linear_model.LinearRegression()
lin_reg.fit(np.transpose([xx1, xx2]), yy)

# Predicition at obs
valid_obs = np.where(np.isfinite(slope_obs) & np.isfinite(vel_obs))
H_modeled_obs = np.nan * np.zeros_like(H_obs)
H_modeled_obs[valid_obs] = 10**lin_reg.predict(np.transpose((np.log10(slope_obs[valid_obs]), np.log10(vel_obs[valid_obs]))))

# Generate full map of modeled H
x1 = np.log10(slope_tan.data)
x2 = np.log10(velocity.data)
valid_x = np.where(np.isfinite(x1) & np.isfinite(x2))
H_modeled_tmp = 10**lin_reg.predict(np.transpose((x1[valid_x], x2[valid_x])))
H_modeled = np.nan * np.zeros_like(slope_tan.data)
H_modeled[valid_x] = H_modeled_tmp
H_modeled = gu.Raster.from_array(H_modeled, transform=velocity.transform, crs=velocity.crs, nodata=-9999)

# Calculate Pearson correlation coeff
r = np.corrcoef(H_obs[np.isfinite(H_modeled_obs)], H_modeled_obs[np.isfinite(H_modeled_obs)])[0, 1]

# - Plot -

# Plot one vs the other
plt.plot(H_obs, H_modeled_obs, ls="", marker="+", c="C0")
plt.xlabel("Observed thickness (m)")
plt.ylabel("Modeled thickness (m)")
plt.axline((0, 0), slope=1, color="k", ls="--")  # Plot 1:1 line
plt.title(f"R2 = {r**2:.2f}")
plt.show()

# Calculate model residuals
res = H_obs - H_modeled_obs

# Plot map along with map of residuals
fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(18, 6))
arg_outline.show(fc="none", ec="k", ax=axes[0])
sc = axes[0].scatter(zbed_x, zbed_y, c=res)
cb = plt.colorbar(sc)
cb.set_label("Residual obs - model (m)")
axes[1].hist(res, bins=40)
axes[1].set_xlabel("Obs - model residual (m)")
H_modeled.show(vmin=0, vmax=500, cbar_title="Modeled thickness (m)", ax=axes[2])
arg_outline.show(fc="none", ec="k", ax=axes[2])
plt.show()

# - Calculate variogram on H residuals using xdem -

# Sample each 50 m distance until 1 km, then every 200 m
vg_bins = np.append(np.arange(0, 1500, 50), np.arange(1500, 5000, 200))
# vg_bins = 30*np.array([np.sqrt(2)**i for i in range(1, 16)])  # xdem's default

vg_output = xdem.spatialstats.sample_empirical_variogram(
    res, coords=np.array((zbed_x, zbed_y)), subsample_method="cdist_point", bin_func=vg_bins
)

func_vgm, params_vgm = xdem.spatialstats.fit_sum_model_variogram(
    list_models=["Spherical"], empirical_variogram=vg_output
)
print("\n## Variogram model with xdem")
print(params_vgm)

xdem.spatialstats.plot_variogram(
    vg_output,
    list_fit_fun=[
        func_vgm,
    ],
    list_fit_fun_label=[
        "Variogram model",
    ],
)
plt.tight_layout()
plt.show()

# - Test variogram and kriging with gstools - #

# experimental variogram
bin_center, gamma = gs.vario_estimate((zbed_x, zbed_y), res)

# Fit spherical model, with nugget
vg_model = gs.Spherical(dim=2)
vg_model.fit_variogram(bin_center, gamma, nugget=True)
print("\n## Variogram model with gstools")
print(vg_model)

# plot
ax = vg_model.plot(x_max=max(bin_center))
ax.scatter(bin_center, gamma)
plt.xlabel("Lag Distance")
plt.ylabel("Variance")
plt.title(f"Sill: {vg_model.var:g} m2, range: {vg_model.len_scale:.0f} m, nugget: {vg_model.nugget} m")
plt.show()

# - Ordinary kriging - #

# Add conditiosn to the obs positions
conditions = np.array([zbed_x, zbed_y, res])

# create glacier contour
mask1 = arg_outline.create_mask(velocity)
mask2 = arg_outline.create_mask(velocity, buffer=velocity.res[0])
final_mask = mask2 & ~mask1

# Add "fake" observations at 0 on glacier edges
idx_contour = np.where(final_mask.data)
res_contour = np.ones(np.size(idx_contour, 1)) * 0.0
grid_x, grid_y = velocity.coords()

cond_x = np.append(conditions[0], grid_x[idx_contour][::2])  # downsampling for tests
cond_y = np.append(conditions[1], grid_y[idx_contour][::2])
cond_z = np.append(conditions[2], res_contour[::2])
# res = np.concatenate((res, res_contour))
# zbed_x = np.concatenate((zbed_x, idx_contour[0]))
# zbed_y = np.concatenate((zbed_y, idx_contour[1]))

# prepare kriging
krig = gs.krige.Ordinary(vg_model, cond_pos=[cond_x, cond_y], cond_val=cond_z)

print(f"Starting kriging at {time.strftime('%H:%M:%S', time.localtime())}")
t0 = time.time()

# downsampling factor
dw = 1

# with structured mesh, on the whole scene -> slower
# x_coords = grid_x[0, :]
# y_coords = grid_y[:, 0]
# krig_res, krig_var  = krig.structured([x_coords[::dw], y_coords[::dw]])
# krig_final = krig_res.T  # x and y are reversed so need to transpose output

# over glacier only
mask = arg_outline.create_mask(velocity)
mask_dw = mask.data[::dw, ::dw]
grid_x_arg = grid_x[::dw, ::dw][mask_dw]
grid_y_arg = grid_y[::dw, ::dw][mask_dw]

krig_res, krig_var = krig.unstructured([grid_x_arg, grid_y_arg])
krig_final, krig_var_final = np.zeros(mask_dw.shape), np.zeros(mask_dw.shape)
krig_final[mask_dw] = krig_res
krig_var_final[mask_dw] = krig_var
print(f"Took {time.time() - t0} s")

# Plot
plt.figure(figsize=(12, 6))
plt.subplot(121)
plt.imshow(krig_final, vmin=-150, vmax=150)
cb = plt.colorbar()
cb.set_label("Kriging of the residuals (m)")

plt.subplot(122)
plt.imshow(np.sqrt(krig_var_final), vmin=0, vmax=70)
cb = plt.colorbar()
cb.set_label("Kriging standard error (m)")

plt.tight_layout()
plt.show()

# Calculate final thickness and save #

# crude oversampling, for tests only
if dw > 1:
    # duplicate values
    krig_final = np.repeat(np.repeat(krig_final, dw, axis=0), dw, axis=1)
    krig_var_final = np.repeat(np.repeat(krig_var_final, dw, axis=0), dw, axis=1)

    # remove extra rows/columns that were added
    krig_final = krig_final[: H_modeled.shape[0], : H_modeled.shape[1]]
    krig_var_final = krig_var_final[: H_modeled.shape[0], : H_modeled.shape[1]]

# Add back modeled H
H_final = H_modeled + krig_final
err_final = H_modeled * 0 + np.sqrt(krig_var_final)

# Mask outside glacier
H_final.set_mask(~mask)
err_final.set_mask(~mask)

# Plot
plt.figure(figsize=(12, 6))
ax1 = plt.subplot(121)
H_final.show(vmin=0, vmax=500, cbar_title="Krigged thickness (m)", ax=ax1)
arg_outline.show(fc="none", ec="k", ax=ax1)

ax2 = plt.subplot(122)
err_final.show(vmin=0, cbar_title="Kriging error (m)", ax=ax2)
arg_outline.show(fc="none", ec="k", ax=ax2)

plt.tight_layout()
plt.savefig(os.path.join(outdir, "thickness_from_kriging.png"), dpi=200)
plt.show()

# Save
os.makedirs(outdir, exist_ok=True)
H_final.save(os.path.join(outdir, "Argentiere_thickness_from_kriging.tif"))
err_final.save(os.path.join(outdir, "Argentiere_thickness_err_from_kriging.tif"))
