#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script used to simulate 100 realistic bed topographies, compatible with the observations, using Sequential Gaussian Simulations (SGS).
All tunable arguments are at the beginning of main.

Authors: Auguste Basset, Amaury Dehecq
"""
import os
import time

import geopandas as gpd
import geoutils as gu
import gstools as gs
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xdem
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.interpolate import griddata
from scipy.signal import fftconvolve
from sklearn import linear_model
from tqdm import tqdm


def run_srf(ref_raster, vg_model, conditions, downsampling, n_simu, mask=None):
    """
    Simulate conditioned random fields using a reference raster grid, a variogram model and conditions.

    ref_raster: the reference raster to use for the output grid
    vg_model: the variogram model as estimated with gstools
    conditions: a 2D array of shape (3, N) containing the x, y and obs at each known point
    n_simu: number of simulations to run
    mask: if not None, simulations will be calculated only where mask is True. Must be of same shape as ref_raster.
    """
    # Run kriging on input conditions and variogram model and create CondSRF instance
    krig = gs.krige.Ordinary(vg_model, cond_pos=conditions[:2], cond_val=conditions[2])
    cond_srf = gs.CondSRF(krig)

    # Set position of pixels where to run simulations
    x_model, y_model = ref_raster.coords()
    if mask is None:
        x_model = x_model[0]
        y_model = y_model[:, 0]
        cond_srf.set_pos([x_model[::downsampling], y_model[::downsampling]], "structured")
    else:
        x_model = x_model[::downsampling, ::downsampling]
        y_model = y_model[::downsampling, ::downsampling]
        mask_dw = mask[::downsampling, ::downsampling]
        cond_srf.set_pos([x_model[mask_dw], y_model[mask_dw]], "unstructured")

    # seeded ensemble generation
    seed = gs.random.MasterRNG(20170519)
    for i in tqdm(range(n_simu)):
        #print(float(i)/n_simu)
        cond_srf(seed=seed(), store=[f"fld{i}", False, False])

    return cond_srf

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

# ----- Main ----- #

# - A list of arguments that can be changed -

# Whether or not to apply a log transformation to the residuals
use_log = False

# Which bed to use (only one can be true)
bed_sia = True 
bed_elmer = False
bed_IGM = False
bed_farinotti = False

# output directory
if use_log:
    if bed_sia:
        outdir = os.path.join("output", "thickness_SIA", "simulated_beds_log")
    if bed_elmer:
        outdir = os.path.join("output", "simulated_beds_elmer_log")
    if bed_IGM:
        outdir = os.path.join("output", "thickness_IGM", "simulated_beds_log")
    if bed_farinotti:
        outdir = os.path.join("output", "thickness_Farinotti", "simulated_beds_log")
else:
    if bed_sia:
        outdir = os.path.join("output", "thickness_SIA", "simulated_beds")
    if bed_elmer:
        outdir = os.path.join("output", "simulated_beds_elmer")
    if bed_IGM:
        outdir = os.path.join("output", "thickness_IGM", "simulated_beds")
    if bed_farinotti:
        outdir = os.path.join("output", "thickness_Farinotti", "simulated_beds")
        
os.makedirs(outdir, exist_ok=True)

# Resolution and number of simulations
# Running at 1/2 resolution for 100 simulation takes ~5h on 16 cores
downsampling = 4
n_simu = 100

# --- Load input data --- #

# Load Argentiere glacier outline
arg_outline = gu.Vector(r"data/gis/outline_pleiades_2020-09-08.shp")

# Load thickness obs
zbed_ds = gpd.read_file(r"data/ice_thx/processed/zbed_arg_measured_UTM32N_shifted.shp")
zbed = zbed_ds.Field_3
zbed_x = zbed_ds.geometry.x
zbed_y = zbed_ds.geometry.y

# load data from fictive DEM to be coherent with methodology
mnt_fict = gu.Raster("output/dh_results/meanDEM-2017_02_15.tif")
zsurf = mnt_fict.value_at_coords(zbed_x, zbed_y)

# Calculate thickness and remove bad values
zsurf[zsurf == 0] = np.nan
zsurf = pd.Series(zsurf.flatten())
H_obs = zsurf - zbed
H_obs[H_obs < 0] = np.nan

valid_obs = np.where(np.isfinite(H_obs))
H_obs = H_obs.values[valid_obs]
zbed_x = zbed_x.values[valid_obs]
zbed_y = zbed_y.values[valid_obs]
zbed = zbed.values[valid_obs]

# Load or calculate modeled thickness and crop to glacier extent with 200 m buffer
if bed_elmer:
    H_model = gu.Raster(r"data/ice_thx/processed/thickness_adrien_2017-02-15_utm32.tif")
if bed_sia:
    H_model = gu.Raster(r"output/thickness_SIA/Argentiere_thickness_SIA.tif")
if bed_IGM:
    H_model = gu.Raster(r"C:/Users/kneibm/Documents/CAIRN/Accu_Argentiere/IGM/ThxInversionGJ/2023_10_24_MK/thkresult.tif")
if bed_farinotti:
    H_model = gu.Raster(r"data/ice_thx/orig/RGI60-11.03638_thickness_Faconsensus.tif")

    
left, bottom, right, top = list(arg_outline.bounds)
H_model.crop([left - 200, bottom - 200, right + 200, top + 200])
H_model.show()
# H_model.crop(arg_outline)

# --- Prepare data --- #

# Extract modeled H at obs
H_model_obs = H_model.value_at_coords(zbed_x, zbed_y).squeeze()
H_model_obs[H_model_obs == H_model.nodata] = np.nan

# Calculate Pearson correlation coeff
r = np.corrcoef(H_obs[np.isfinite(H_model_obs)], H_model_obs[np.isfinite(H_model_obs)])[0, 1]

# Calculate model residuals
res = H_obs - H_model_obs
res_log = np.log(H_obs) - np.log(H_model_obs)
res_log[np.isinf(res_log)] = np.nan  # when H = 0, res are infinite

# - Plot -

fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(12, 12))

# Modeled H
H_model.show(ax=axes[0, 0], cbar_title="Modeled thickness (m)", vmin=0, vmax=400)
arg_outline.show(fc="none", ec="k", ax=axes[0, 0])

# Residuals map
cmap = "coolwarm"
vmax = 150
arg_outline.show(fc="none", ec="k", ax=axes[0, 1])
sc = axes[0, 1].scatter(zbed_x, zbed_y, c=res, vmin=-vmax, vmax=vmax, cmap=cmap)
# make sure cmap is same height as subplot
divider = make_axes_locatable(axes[0, 1])
cax = divider.append_axes("right", size="5%", pad=0.05)
norm = matplotlib.colors.Normalize(vmin=-vmax, vmax=150)
cb = matplotlib.colorbar.ColorbarBase(cax, cmap=cmap, norm=norm)
# cb = plt.colorbar(sc)
cb.set_label("Residual obs - model (m)")

# Residuals scatterplot
plt.sca(axes[1, 0])
plt.plot(H_obs, H_model_obs, ls="", marker="+", c="C0")
plt.xlabel("Observed thickness (m)")
plt.ylabel("Modeled thickness (m)")
plt.axline((0, 0), slope=1, color="k", ls="--")  # Plot 1:1 line
plt.title(f"R2 = {r**2:.2f}")

# Residuals histogram
plt.sca(axes[1, 1])
plt.hist(res, bins=40)
plt.gca().set_xlabel("Obs - model residual (m)")

plt.tight_layout()
plt.savefig(os.path.join(outdir, "resiudals_plot.png"), dpi=200)
plt.show()

# --- Sequential Gaussian Simulations --- #

# Compute a variogram and fit that variogram to an exponential model
if use_log:
    bin_center, gamma = gs.vario_estimate((zbed_x, zbed_y), res_log)
else:
    bin_center, gamma = gs.vario_estimate((zbed_x, zbed_y), res)

# -- Tests to find the best covariance model --
# Define a set of models to test
# models = {
#     "Gaussian": gs.Gaussian,
#     "Exponential": gs.Exponential,
#     "Matern": gs.Matern,
#     "Stable": gs.Stable,
#     "Rational": gs.Rational,
#     "Circular": gs.Circular,
#     "Spherical": gs.Spherical,
#     "SuperSpherical": gs.SuperSpherical,
#     "JBessel": gs.JBessel,
# }
models = {
    "Gaussian": gs.Gaussian,
    "Exponential": gs.Exponential,
    "Matern": gs.Matern,
}

scores = {}

# Iterate over all models, fit their variogram and calculate the r2 score.
# plot the estimated variogram
plt.figure()
plt.scatter(bin_center, gamma, color="k", label="data")
ax = plt.gca()
# fit all models to the estimated variogram
for model in models:
    fit_model = models[model](dim=2)
    para, pcov, r2 = fit_model.fit_variogram(bin_center, gamma, nugget=True, return_r2=True)
    fit_model.plot(x_max=3000, ax=ax)
    scores[model] = r2

# Create a ranking based on the score and determine the best models
ranking = sorted(scores.items(), key=lambda item: item[1], reverse=True)
print("RANKING by Pseudo-r2 score")

for i, (model, score) in enumerate(ranking, 1):
    print(f"{i:>6}. {model:>15}: {score:.5}")
plt.show()

# -- Select "best" model and fit again --
fit_model = gs.Exponential(dim=2)
fit_model.fit_variogram(bin_center, gamma, nugget=True)

ax = fit_model.plot(x_max=max(bin_center))
ax.scatter(bin_center, gamma)
plt.xlabel("Lag Distance")
plt.ylabel("Variance")
plt.title(f"Sill: {fit_model.var:g} m2, range: {fit_model.len_scale:.0f} m")
plt.savefig(os.path.join(outdir, "variogram_model.png"), dpi=200)
plt.show()

# create glacier contour
mask1 = arg_outline.create_mask(H_model)
mask2 = arg_outline.create_mask(H_model, buffer=H_model.res[0])
final_mask = mask2 & ~mask1

# Add "fake" observations at 0 on glacier edges
idx_contour = np.where(final_mask.data)
res_contour = np.ones(np.size(idx_contour, 1)) * 0.0001
res = np.concatenate((res, res_contour))
res_log = np.concatenate((res_log, np.log(res_contour)))
grid_x, grid_y = H_model.coords()
zbed_x = np.concatenate((zbed_x, grid_x[idx_contour]))
zbed_y = np.concatenate((zbed_y, grid_y[idx_contour]))

# -- Run simulations --

conditions = np.array([zbed_x, zbed_y, res_log])

print(f"Starting SGS at {time.strftime('%H:%M:%S', time.localtime())}")
t0 = time.time()
cond_srf = run_srf(H_model, fit_model, conditions, downsampling=downsampling, n_simu=n_simu)
print(f"Took {time.time() - t0} s")

# Sum (downsampled) simulated bed + reference bed
H_simu = np.zeros((n_simu, *H_model.shape))
for k in range(n_simu):
    if use_log:
        H_simu[k][::downsampling, ::downsampling] = cond_srf.all_fields[k].T + np.log(
            H_model.data[::downsampling, ::downsampling]
        )
    else:
        H_simu[k][::downsampling, ::downsampling] = (
            cond_srf.all_fields[k].T + H_model.data[::downsampling, ::downsampling]
        )

# -- Interpolate missing values (due to downsampling) with linear interpolation -- #

# Indexes where simulations where run
idx_run = np.where(H_simu[0] != 0)

# Grid (over glaciers) where to interpolate
rowmin, rowmax, colmin, colmax = gu.raster.get_valid_extent(np.where(mask1.data == 0, np.nan, 1))
row_grid_glacier, col_grid_glacier = np.meshgrid(
    np.arange(rowmin, rowmax + 1), np.arange(colmin, colmax + 1), indexing="ij"
)

# Interpolate for each simulation
H_simu = np.copy(H_simu)
for i in range(n_simu):
    z = H_simu[i][H_simu[i] != 0]
    interp = griddata(np.transpose(idx_run), z, (row_grid_glacier, col_grid_glacier), method="linear")
    if use_log:
        H_simu[i][row_grid_glacier, col_grid_glacier] = np.exp(interp)
    else:
        H_simu[i][row_grid_glacier, col_grid_glacier] = interp

    # masking off glacier
    H_simu[i][~mask1.data] = np.nan

# Save to raster
print(f"Saving output files in {outdir}")
for k in range(n_simu):
    raster = gu.Raster.from_array(H_simu[k], H_model.transform, H_model.crs, nodata=H_model.nodata)
    raster.save(os.path.join(outdir, f"simulated_thickness_{k}.tif"))

# plotting first 4 simulations
fig, ax = plt.subplots(2, 2, sharex=True, sharey=True)
ax = ax.flatten()
for i in range(4):
    im = ax[i].imshow(H_simu[i], vmin=0, vmax=500)
    if i == 0:
        cbar = fig.colorbar(im, ax=ax, orientation="vertical", fraction=0.05, pad=0.05)
plt.savefig(os.path.join(outdir, "simu_samples.png"), dpi=200)
plt.show()

# -- Calculating mean, median and std --

thx_std = np.std(H_simu, axis=0)
thx_mean = np.mean(H_simu, axis=0)
thx_median = np.median(H_simu, axis=0)

# Saving to file
thx_std_rst = gu.Raster.from_array(thx_std, H_model.transform, H_model.crs, nodata=H_model.nodata)
thx_std_rst.save(os.path.join(outdir, "simulated_thickness_std.tif"))

thx_mean_rst = gu.Raster.from_array(thx_mean, H_model.transform, H_model.crs, nodata=H_model.nodata)
thx_mean_rst.save(os.path.join(outdir, "simulated_thickness_mean.tif"))

thx_median_rst = gu.Raster.from_array(thx_median, H_model.transform, H_model.crs, nodata=H_model.nodata)
thx_median_rst.save(os.path.join(outdir, "simulated_thickness_median.tif"))

# Plotting std
plt.imshow(thx_std)
cb = plt.colorbar()
cb.set_label("Std (m)")
plt.show()
