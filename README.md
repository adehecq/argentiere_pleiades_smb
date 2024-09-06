# argentiere_pleiades_smb
Repository to calculate the specific mass balance of Argentiere glacier in the French Alps using Pleiades satellite data.

Main reference: Kneib, M., Dehecq, A., Gilbert, A., Basset, A., Miles, E. S., Jouvet, G., Jourdain, B., Ducasse, E., Beraud, L., Rabatel, A., Mouginot, J., Carcanade, G., Laarman, O., Brun, F., and Six, D.: Distributed surface mass balance of an avalanche-fed glacier, EGUsphere [preprint], https://doi.org/10.5194/egusphere-2024-1733, 2024.

# Dependencies

Data to run the codes in available on the GLACIOCLIM website: https://glacioclim.osug.fr/-Acces-a-des-donnees-elaborees-

The Python scripts use the xdem and geoutils libraries: https://xdem.readthedocs.io/ and https://geoutils.readthedocs.io/

A number of these MATLAB scripts were derived from the work by Miles et al. (2021) and are available on Evan Miles' github: https://github.com/miles916/grid_continuity_SMB. 

These scripts call some functions developed in ImGRAFT (Messerli & Grinsted, 2015): https://github.com/grinsted/ImGRAFT

We used the IGM model for one of the ice thickness modelling approaches (Jouvet et al., 2023): https://github.com/jouvetg/igm 

We used the Elmer/Ice model for the forward modelling (Gagliardini et al., 2013): https://github.com/elmercsc/elmerfem. 

We have also made the SnowSlide snow redistribution scheme available online (Bernhardt & Schulz, 2010): https://github.com/OGGM/Snowslide. 

# References
Bernhardt, M., & Schulz, K. (2010). SnowSlide: A simple routine for calculating gravitational snow transport. Geophysical Research Letters, 37(11). https://doi.org/10.1029/2010GL043086
Gagliardini, O., Zwinger, T., Gillet-Chaulet, F., Durand, G., Favier, L., de Fleurian, B., Greve, R., Malinen, M., Martín, C., Råback, P., Ruokolainen, J., Sacchettini, M., Schäfer, M., Seddik, H., & Thies, J. (2013). Capabilities and performance of Elmer/Ice, a new-generation ice sheet model. Geoscientific Model Development, 6(4), 1299–1318. https://doi.org/10.5194/gmd-6-1299-2013
Jouvet, G., & Cordonnier, G. (2023). Ice-flow model emulator based on physics-informed deep learning. Journal of Glaciology, 1–15. https://doi.org/10.1017/jog.2023.73
Miles, E., McCarthy, M., Dehecq, A., Kneib, M., Fugger, S., & Pellicciotti, F. (2021). Health and sustainability of glaciers in High Mountain Asia. Nature Communications 2021 12:1, 12(1), 1–10. https://doi.org/10.1038/s41467-021-23073-4
Messerli, A., & Grinsted, A. (2015). Image georectification and feature tracking toolbox: ImGRAFT. Geoscientific Instrumentation, Methods and Data Systems, 4(1), 23–34. https://doi.org/10.5194/gi-4-23-2015

# Contacts 
Marin Kneib, marin.kneib@univ-grenoble-alpes.fr
Amaury Dehcq, amaury.dehecq@univ-grenoble-alpes.fr

# Content
- IGM_inversion:
		- igm_run.py: classic igm_run script
		- params.json: parameters and modules used for the flux inversion.
		- sensitivity_igm.py: script to run iteratively the inversion, changing the weights at each iteration.
- environment_xdem.yml: environment used to run Sequential Gaussian Simulations.
- test_kriging.py: attempt to extrapolate the bed observations to the whole glacier using ordinary kriging with a deterministic thickness model, based on slope and velocity.
- Thickness_SIA.py: script to calculate the bed from surface velocity observations using the Shallow Ice Approximation (SIA) calibrated with GPR observations.
- prepare_bed_data.py: pre-processing script for bed elevation data.
- simulate_bed_with_sgs.py: script to run Sequential Gaussian Simulations.
- UsefulCodes:
		- geotiffcrop_shp.m: MATLAB function to crop raster with shapefile.
- FluxCalcsSimple.m: Function to estimate surface mass balance distribution for a glacier from inputs of ice thickness, thinning, and velocity, based on the continuity equation. Function adapted from Miles et al. (2021): 10.1038/s41467-021-23073-4
- FluxCalcsUncertainty.m: Function to estimate surface mass balance distribution for a glacier from inputs of ice thickness, thinning, and velocity, based on the continuity equation. Function adapted from Miles et al. (2021): 10.1038/s41467-021-23073-4. In addition to FluxCalcsSimple.m, this function allows to conduct Monte-Carlo simulations to account for uncertainties.
- FluxIGM_SMB_Argentiere.m: Master script to estimate surface mass balance distribution for a glacier from inputs of IGM flux inversion and thinning (inspired from the scripts by Miles et al, 2021).
- master_ContinuitySMB_Argentiere.m - Master script to estimate surface mass balance distribution for a glacier from inputs of ice thickness, thinning, and velocity, based on the continuity equation (adapted from Miles et al, 2021)
- SMB_GLACIOCLIM.m: script to pre-process the SMB GLACIOCLIM measurements.
- tendancy_dh.m: script to derive mean and filtered elevation change from a set of DEMs.
