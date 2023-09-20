#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script used to simulate 100 realistic bed topographies, compatible with the observations, using Seuqnetial Gaussian Simulations (SGS).

Authors: AUguste Basset, Amaury Dehecq
"""
import geopandas as gpd
import geoutils as gu
import gstools as gs
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas
import pyvista as pv
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.interpolate import griddata


def insert_zeros(array):

    original_shape = array.shape
    new_shape = (original_shape[0], original_shape[1] + 1, original_shape[2])
    new_array = np.zeros(new_shape)
    new_array[:, :-1, :] = array

    return new_array


def run_srf(H_model, vg_model, downsampling, ens_no):

    global H_simu
    # Create the kriging grid relative to the considered zone
    nrows, ncols = H_model.res

    grid = pv.UniformGrid()
    # Grid origin
    grid.origin = (H_model.bounds.left, H_model.bounds.top, 0)
    # Cell sizes
    grid.spacing = (ncols // downsampling, nrows // downsampling, 1)
    # Number of cells in each direction
    grid.dimensions = (H_model.width // downsampling, H_model.height // downsampling, 1)

    # conditions
    cond_pos = list([list(zbed_x), list(zbed_y)])
    cond_val = list(res)

    # Kriging and creating a condSrf class object
    krig = gs.krige.Ordinary(vg_model, cond_pos=cond_pos, cond_val=cond_val)
    krig.mesh(grid, name="Residuals (obs-mod)")
    # grid now has the residual scalar field and kriging variance as data arrays.
    cond_srf = gs.CondSRF(krig)

    # same output positions for all ensemble members
    mask = arg_outline.create_mask(H_model)
    x_model, y_model = H_model.coords()
    x_model = x_model[0]
    y_model = y_model[:, 0]

    cond_srf.set_pos([x_model[::downsampling], y_model[::downsampling]], "structured")

    # seeded ensemble generation
    seed = gs.random.MasterRNG(20170519)
    for i in range(ens_no):
        cond_srf(seed=seed(), store=[f"fld{i}", False, False])

    # fig, ax = plt.subplots(4, 5, sharex=True, sharey=True)
    # ax = ax.flatten()
    # for i in range(ens_no):
    #     # on peut afficher les champs simulés seuls ou bien replacés sur le glacier
    #     im = ax[i].imshow(cond_srf[i].T, origin="lower")
    #     if i==1:
    #         cbar = fig.colorbar(im, ax=ax, orientation = 'vertical', fraction = 0.05, pad = 0.05)
    # plt.show()

    return cond_srf


# --------------------------------------------------------------- Main ---------------------------------------------------------
# --- Load input data --- #

# Load Argentiere glacier outline
arg_outline = gu.Vector(r"C:\Users\Aug\Desktop\Stage_IGE\Data\GIS\outline_bed_v4.shp")

# Load thickness obs
zbed_ds = gpd.read_file(r"C:\Users\Aug\Desktop\Stage_IGE\Data\GIS\Zbed_Arg_Measured_UTM32N.shp")
zbed = zbed_ds.Field_3
zbed_x = zbed_ds.geometry.x
zbed_y = zbed_ds.geometry.y

# load data from fictive DEM to be coherent with methodology
mnt_fict = gu.Raster(r"C:\Users\Aug\Desktop\Stage_IGE\Data\dH\MNT\MNT_fict.tif")
zsurf = mnt_fict.value_at_coords(zbed_x, zbed_y)

# Calculate thickness and remove bad values
zsurf[zsurf == 0] = np.nan
zsurf = pandas.Series(zsurf.flatten())
H_obs = zsurf - zbed
H_obs[H_obs < 0] = np.nan

valid_obs = np.where(np.isfinite(H_obs))
H_obs = H_obs.values[valid_obs]
zbed_x = zbed_x.values[valid_obs]
zbed_y = zbed_y.values[valid_obs]

# Load modeled bed
H_model = gu.Raster(r"/home/bassetau/Documents/Auguste/Accu_Argentiere/IceTHX/Bed_Adrien.tif")
H_model.crop(arg_outline)

# --- Prepare data --- #

# Extract modeled H at obs
H_model_obs = H_model.value_at_coords(zbed_x, zbed_y)
H_model_obs[H_model_obs == H_model.nodata] = np.nan

# Calculate Pearson correlation coeff
r = np.corrcoef(H_obs[np.isfinite(H_model_obs)], H_model_obs[np.isfinite(H_model_obs)])[0, 1]

# Calculate model residuals
res = np.log(H_obs) - np.log(H_model_obs)


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
plt.show()

## -------------- Gaussian simulation-----------------##

# Compute a variogram and fit that variogram to an exponential model
bin_center, gamma = gs.vario_estimate((zbed_x, zbed_y), res)

# ------------------Test to find the best covariance model-----------
# Define a set of models to test
models = {
    "Gaussian": gs.Gaussian,
    "Exponential": gs.Exponential,
    "Matern": gs.Matern,
    "Stable": gs.Stable,
    "Rational": gs.Rational,
    "Circular": gs.Circular,
    "Spherical": gs.Spherical,
    "SuperSpherical": gs.SuperSpherical,
    "JBessel": gs.JBessel,
}
scores = {}

# Iterate over all models, fit their variogram and calculate the r2 score.
# plot the estimated variogram
plt.scatter(bin_center, gamma, color="k", label="data")
ax = plt.gca()
# fit all models to the estimated variogram
for model in models:
    fit_model = models[model](dim=2)
    para, pcov, r2 = fit_model.fit_variogram(bin_center, gamma, return_r2=True)
    fit_model.plot(x_max=3000, ax=ax)
    scores[model] = r2

# Create a ranking based on the score and determine the best models
ranking = sorted(scores.items(), key=lambda item: item[1], reverse=True)
print("RANKING by Pseudo-r2 score")

for i, (model, score) in enumerate(ranking, 1):
    print(f"{i:>6}. {model:>15}: {score:.5}")
plt.show()
# --------------------------------------------------------------------

# fit_model = gs.Exponential(dim = 2, var = 2658.4120995073145, len_scale = 440.7028787096736, nugget = 356.8240249489487)
fit_model = gs.Exponential(dim=2)
fit_model.fit_variogram(bin_center, gamma, nugget=True)

ax = fit_model.plot(x_max=max(bin_center))
ax.scatter(bin_center, gamma)
gs.covmodel.plot.plot_variogram(fit_model, ax=ax)
plt.xlabel("Lag Distance")
plt.ylabel("Variogram")
plt.show()


# file_zone = r'C:\Users\Aug\Desktop\Stage IGE\Data\IceTHX\Argentiere_part2.gpkg'

# fields = zonal_srf_bed(zbed_ds, H_model, file_zone, fit_model).all_fields
# x_model_zone, y_model_zone = zonal_coord(H_model)
# H_cond.value_at_coords(x_model_zone, y_model_zone) = fields

# create glacier contour
mask1 = arg_outline.create_mask(H_model)
mask2 = arg_outline.create_mask(H_model, buffer=10)
final_mask = mask2 & ~mask1
idx_contour = np.where(final_mask.data)
res_contour = np.ones(np.size(idx_contour, 1)) * 0.001
res = np.concatenate((res, res_contour))
# res now contains residuals on measurements sites and on glacier contour
# modify zbed_x and zbed_y accordingly
zbed_x = np.concatenate((zbed_x, idx_contour[0]))
zbed_y = np.concatenate((zbed_y, idx_contour[1]))

# modify idx_contour and mask1 for what follows
idx_contour = np.asarray(idx_contour)
idx_contour = np.delete(idx_contour, -1, axis=1)
idx_contour = np.delete(idx_contour, -1, axis=1)
idx_contour = np.delete(idx_contour, -1, axis=1)
final_mask_array = np.delete(final_mask.data, -1, axis=0)

# run simulations
downsampling = 2
ens_no = 100
cond_srf = run_srf(H_model, fit_model, downsampling=downsampling, ens_no=ens_no)

# replace simulations on the glacier with reference value of bed topography
# first, find (x,y) coordinates of on-glacier pixels based on every pixels
x_model, y_model = H_model.coords()
x_model = x_model[0]
y_model = y_model[:, 0]
y_model = np.delete(y_model, H_model.shape[0] - 1, 0)
on_glacier_coords = []
for i in x_model:
    for j in y_model:
        if H_model.value_at_coords(i, j) != H_model.nodata:
            on_glacier_coords.append((i, j))

# Retrieve (x,y) on-glacier coordinates
x_on_glacier = []
y_on_glacier = []
for elem in on_glacier_coords:
    x_on_glacier.append(elem[0])
    y_on_glacier.append(elem[1])
i_on_glacier, j_on_glacier = H_model.xy2ij(x_on_glacier, y_on_glacier)

# we ran simulation for :
x_run = x_model[::downsampling]
x_run = x_run.tolist()
y_run = y_model[::downsampling]
y_run = y_run.tolist()

# retrieve coords of simulation AND on-glacier
x_run_on_glacier = []
y_run_on_glacier = []
for x in x_run:
    for y in y_run:
        if (x, y) in on_glacier_coords:
            x_run_on_glacier.append(x)
            y_run_on_glacier.append(y)

# convert it to indexes
i_run_on_glacier, j_run_on_glacier = H_model.xy2ij(x_run_on_glacier, y_run_on_glacier)
i_run_on_glacier = np.unique(i_run_on_glacier)
j_run_on_glacier = np.unique(j_run_on_glacier)

i_run = []
j_run = []
for i in x_run:
    for j in y_run:
        k, l = H_model.xy2ij(i, j)
        i_run.append(k)
        j_run.append(l)

# i_run = np.unique(i_run)
# j_run = np.unique(j_run)

# initialize simulated beds arrays
H_simu = np.zeros((ens_no, np.size(y_model), np.size(x_model)))
H_model_temp = H_model.data
H_model_temp = np.delete(H_model_temp, -1, axis=0)

for k in range(ens_no):
    compt_i = 0
    compt_j = 0
    H_simu[k][~final_mask_array] = np.nan
    for i in i_run_on_glacier:
        for j in j_run_on_glacier:
            H_simu[k][i][j] = cond_srf.all_fields[k][compt_j][compt_i] + np.log(H_model_temp[i][j])
            compt_j += 1
        compt_j = 0
        compt_i += 1

# Now we have to interpolate simulated residuals over the all glacier
# xx, yy = np.meshgrid(i_run, j_run, indexing='ij')

idx_run = np.where(~np.isnan(H_simu[0]))
idx_run = np.asarray(idx_run)

grid_i, grid_j = np.meshgrid(np.unique(i_on_glacier), np.unique(j_on_glacier))
for i in range(ens_no):
    z = H_simu[i][~np.isnan(H_simu[i])]
    func = griddata(np.transpose(idx_run), z, (grid_i, grid_j), method="linear")
    H_simu[i] = np.exp(np.transpose(func))
H_simu = insert_zeros(H_simu)

# on redécoupe selon les limites du glacier
for k in range(ens_no):
    raster = gu.Raster.from_array(H_simu[k], H_model.transform, H_model.crs, nodata=H_model.nodata)
    raster.set_mask(~mask1)
    var = k
    raster.save(f"simulated_bed_{var}.tif")
    H_simu[k] = raster.to_xarray()
    H_simu[k][H_simu[k] < 0] = np.nan

# # Now we can go back to physical values the values and add reference bed
# for k in range(ens_no):
#     H_simu[k] = np.exp(H_simu[k])
#     H_simu[k][H_simu[k]<0] = np.nan

# plotting
fig, ax = plt.subplots(2, 2, sharex=True, sharey=True)
ax = ax.flatten()
for i in range(4):
    im = ax[i].imshow(H_simu[i + 12], vmin=0, vmax=500)
    if i == 0:
        cbar = fig.colorbar(im, ax=ax, orientation="vertical", fraction=0.05, pad=0.05)
plt.show()
