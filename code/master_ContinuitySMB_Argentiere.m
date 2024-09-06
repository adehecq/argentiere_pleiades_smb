%master_ContinuitySMB_Argentiere.m - Master script to estimate surface mass balance
%distribution for a glacier from inputs of ice thickness, thinning, and
%velocity, based on the continuity equation (adapted from Miles et al, 2021)
%
% New glaciers require this template to be coded with paths for the full set of
% inputs, and any additional needed preprocessing. 
%
% Other m-files required: C2xyz.m, index_nanfill.m, remove_small_zones.m,
% segment_Gmask_slope2.m, subset_geo.m, through_fluxes.m, zonal_aggregate.m
%
% Author: Evan Miles
% Work address: Swiss Federal Research Institute WSL
% Email: evan.miles@wsl.ch
% Updated: Marin Kneib
% Email: marin.kneib@univ-grenoble-alpes.fr

%% ------------- BEGIN CODE --------------
clear
close all

%set home directory
homedir = 'C:\Users\kneibm\Documents\CAIRN\Accu_Argentiere\argentiere_pleiades_smb\code\';
cd(homedir)

% path to useful codes
addpath('../UsefulCodes/ImGRAFT')
addpath('../UsefulCodes/')
addpath('MilesSMB/')
addpath(homedir)
addpath('C:\Users\kneibm\Documents\CAIRN\Accu_Argentiere\argentiere_pleiades_smb\SMB_VanTricht2021')

%glacier ID
P.Glacier = 'Argentiere';

%title of output directory
resultsfolder = '..\output\smb\';

% which configuration (ice thickness, velocity, dh/dt...)
thx_elmer = 0;
thx_farinotti = 0;
thx_sia = 1;
thx_igm = 0;

vel_pleiades = 1;
vel_millan = 0;
vel_igm = 0;

dh_pleiades_20182023 = 1;
dh_pleiades = 1;
dh_hugonnet = 0;

%% settings
P.plotouts=0; %output plots or not
P.exports=0; %save geotiffs or not
P.exports_indiv = 0; % Save every single SMB geotiff or not
P.DX = 20; %resolution to run calculations at 
P.segdist=300; %effective linear distance between flowbands
P.Vfilter=0; %swtich to smooth velocity data or not
P.dhfilter=0; %switch to peform 2x 3-sigma outlier removal (from overall dataset - only if erroneous pixels are common)
P.THXfilter=0; %switch to apply simply Gaussian filter to thickness data (useful for field thickness measurements)
P.umult=0; %switch for column-average velocity [0.8-1] is physical range. '0' estimates based on THX dist. '2' estimates for each pixel. 
P.Vreproj=0; %0 if velocities are provided oriented in the correct coordinate system, 1 if they are in the source data projection, [2 if they are true north], [3 to determine from slope]
P.fdivfilt=1; %use VanTricht gradient filters (2), just flux filter (1) or not at al (0)
P.uncertainty=100; %0 for 'simple run without uncertainty, otherwise N runs; note that you will then need to setup perturbations within the N structure
P.Valongslope=0; %1 for velocity along slope, 0 for 'normal' (i.e. as in Van Tricht, 2021; Miles 2021) 

%% inputs     

% velocity
if vel_millan
    V.pathx = fullfile(homedir,'..\data\velocity\Millan2022\VX_RGI-11_2021July01_clip.tif');
    V.pathy = fullfile(homedir,'..\data\velocity\Millan2022\VY_RGI-11_2021July01_clip.tif');
    V.mult=1; %scale to convert input units to m/a
    v_suff = '_vMillan';
end
if vel_pleiades
    V.pathx = fullfile(homedir,'..\output\velocity\vx_mean_after_30_15°.tif');
    V.pathy = fullfile(homedir,'..\output\velocity\vy_mean_after_30_15°.tif');
    V.mult=1; %scale to convert input units to m/a
    v_suff = '_vPleiades';
end
if vel_igm
    V.pathx = 'C:\Users\kneibm\Documents\CAIRN\Accu_Argentiere\IGM\ThxInversionGJ\2024_01_30_MK\full_50m_2024_02_13\ubar.tif';
    V.pathy = 'C:\Users\kneibm\Documents\CAIRN\Accu_Argentiere\IGM\ThxInversionGJ\2024_01_30_MK\full_50m_2024_02_13\vbar.tif';
    V.mult=1; %scale to convert input units to m/a
    v_suff = '_vIGM';
end

% dh
if dh_pleiades_20182023
    DH.path = fullfile(homedir,'..\output\dh_results\tendancy_filt_2018-2023.tif');
    DH.mult=1; %scale to convert input units to m/a
    dh_suff = '_dhPleiades_2018-2023';
end
if dh_pleiades
    DH.path = fullfile(homedir,'..\output\dh_results\tendancy_filt.tif');
    DH.mult=1; %scale to convert input units to m/a
    dh_suff = '_dhPleiades';
end
if dh_hugonnet
    DH.path = fullfile(homedir,'..\data\dh\Hugonnet_2010-2020_utm32N.tif');
    DH.mult=1; %scale to convert input units to m/a
    dh_suff = '_dhHugonnet';
    % TODO
end

% thickness
if thx_elmer
    if P.uncertainty==0
        THX.path = 'C:\Users\kneibm\Documents\CAIRN\Accu_Argentiere\Elmer\BedInversion_2024_01_18\processed\thickness_adrien_2017-02-15_utm32_NewMaskAccu.tif';
    else
        THX.path = fullfile(homedir,'..\output\thickness_Elmer1\simulated_beds\');
    end
    thx_suff = '_tElmer';
end
if thx_farinotti
    if P.uncertainty==0
        THX.path = fullfile(homedir,'..\data\ice_thx\orig\RGI60-11.03638_thickness_Faconsensus.tif');
    else
        THX.path = fullfile(homedir,'..\output\thickness_Farinotti\simulated_beds\');
    end
    thx_suff = '_tFarinotti';
end
if thx_sia
    if P.uncertainty==0
        THX.path = fullfile(homedir,'..\output\thickness_SIA\Argentiere_thickness_SIA.tif');
    else
        THX.path = fullfile(homedir,'..\output\thickness_SIA\simulated_beds\');
    end
    thx_suff = '_tSIA';
end
if thx_igm
    if P.uncertainty==0
        THX.path = 'C:\Users\kneibm\Documents\CAIRN\Accu_Argentiere\IGM\ThxInversionGJ\2024_01_30_MK\full_50m_bed_constrained\thkresult.tif';
    else
        THX.path = fullfile(homedir,'..\output\thickness_igm_upd\simulated_beds\');
    end
    thx_suff = '_tIGM';
end

% DEM (always the same)
DEM.path = fullfile(homedir,'..\output\dh_results\meanDEM-2017_02_15.tif');

%% initialization
addpath(genpath([homedir,'MilesSMB']))
addpath(genpath([homedir,'MilesSMB\cbrewer']))

if P.fdivfilt==0
    filt_suff = '_fNone';
end
if P.fdivfilt==1
    filt_suff = '_fFluxFilter';
end
if P.fdivfilt==2
    filt_suff = '_fVanTricht';
end
if P.Valongslope==1
    AS_suff = '_vAlongSlope';
else
    AS_suff = '';
end
if P.uncertainty > 0
    str = string(P.uncertainty);
    UNC_suff = ['_uncertainty_N=',str{:}];
else
    UNC_suff = '';
end

date = string(datetime("today","Format","uuuu-MM-dd"));
outdir = [resultsfolder,date{:},v_suff,dh_suff,thx_suff,filt_suff,AS_suff,UNC_suff];
mkdir(outdir)

outtitle1 = [num2str(P.DX) 'mgrid'];
    
%% Load all input data, resample to common grid, etc
if P.uncertainty==0
    [THX,DEM,V,DH]=load_subset_all(P,THX,DEM,V,DH);
else
%     [THX,DEM,V,DH]=load_subset_all(P,THX,DEM,V,DH);
    [THX,DEM,V,DH]=load_subset_all_SGS(P,THX,DEM,V,DH);
end
% clip to glacier outlines
if P.uncertainty == 0
    outlines = shaperead('C:\Users\kneibm\Documents\CAIRN\Accu_Argentiere\argentiere_pleiades_smb\data\gis\outline_pleiades_2020-09-08.shp');
    dgeot = geotiffinfo(THX.path);
    [~,THX.data] = geotiffcrop_shp(outlines,THX.data,dgeot);
end

THX.data(THX.data>700 | THX.data<-200) = NaN;
%% Cleaning raw data if needed

%identify likely errors
V.iERR=(abs(V.Uraw)>400)|(abs(V.Vraw)>400);
if P.Vfilter==1 %gaussian low-pass filter removing extreme variations
    V.Uraw=imgaussfilt(V.Uraw,5);
    V.Vraw=imgaussfilt(V.Vraw,5);
end

if P.THXfilter==1% filter thickness data
    THX.dataR=THX.data;
    MASK = (THX.data==0);
    THX.data=imgaussfilt(THX.dataR,4); %gaussian low-pass filter. Important for thickness maps derived from field data
    THX.data(MASK)=0;
end

if P.dhfilter==1 %if noisy or gappy, apply some preprocessing to the dH
    DH.data2=DH.data;
    DH.errthresh1=3.*nanstd(DH.data(:));
    DH.data2(abs(DH.data2)>DH.errthresh1)=NaN;
    DH.errthresh2=3.*nanstd(DH.data2(:));
    DH.data2(abs(DH.data2)>DH.errthresh2)=NaN;
    DH.dH3=imgaussfilt(DH.data2); %smooths dH slightly, expands NaN around bad DH data
else
    DH.dH3=DH.data;
end

%% REPROJECT AND RESAMPLE ALL INPUTS (to THX coordinate system), derive MASK, etc
if P.uncertainty == 0
    N = resample_inputs_mk(P.DX,THX,V,DEM,DH);
else
%     N = resample_inputs_mk(P.DX,THX,V,DEM,DH);
    N = resample_inputs_mk_SGS(P.DX,THX,V,DEM,DH);
end
N.THX(N.THX<-200) = NaN; % remove results from interpolation on the edges
N.fdivfilt=P.fdivfilt;
N.uncertainty=P.uncertainty;
N.Valongslope = P.Valongslope;
N.exports_indiv = P.exports_indiv;
N.Glacier = P.Glacier;
N.info.GeoTIFFTags.GeoKeyDirectoryTag = THX.info.GeoTIFFTags.GeoKeyDirectoryTag;


%% fill gaps in reprojected inputs if needed

%gap-fill dH based on elevation
N.DH0=N.DH; %unfilled values
N.DH(N.MASK)= index_nanfill(N.DH(N.MASK),N.DEM(N.MASK));
N.DH((N.MASK==0))=0;

%% zone segmentation

N.zones = uint16(segment_Gmask_slope2(N.DEM,N.MASK,P.DX,P.segdist)); %segments glacier mask and DEM into approximately uniformly-spaced elevation bands

%% determine column-average velocity correction
if vel_igm
    N.umult = 1;
    N.sigUmult = 0;
else
    N.umult = 0.9;
    N.sigUmult = 0.1;
end

%% Uncertainties on velocities and dh
if vel_millan
    N.sigU = 10;
    N.sigV = 10;
else
    N.sigU = 2.4;
    N.sigV = 2.4;
end

if dh_hugonnet
    N.sigdH = 1;
else
    N.sigdH = 0.07;
end
N.sigdensity_mixedzone = 0.15;
N.density_mixedzone = 0.75;

%% SMB calculations
cd(outdir)

diary('output.txt')
if N.uncertainty>1
    N=FluxCalcsUncertainty(N); %calculates fluxes 
else
    N=FluxCalcsSimple(N); %calculates fluxes 
end    
diary off;

%% Get misfit values
% Read the saved output from the file
fid = fopen('output.txt', 'r');
output = textscan(fid, '%s', 'Delimiter', '\n');
fclose(fid);

% Extract numbers from the output
numbers = [];
for i = 1:numel(output{1})
    line = output{1}{i};
    numbers_in_line = regexp(line, '[-+]?\d{1,4}(\.\d{4})?(e[-+]?\d{2})?', 'match');
    if ~isempty(numbers_in_line) 
        if abs(str2double(numbers_in_line))<10
            numbers = [numbers; str2double(numbers_in_line)];
        end
    end
end

if P.fdivfilt==1
    min(numbers(2:end))
    mean(numbers(2:end))
    max(numbers(2:end))
end
if P.fdivfilt==2
    min(numbers(6:5:end))
    mean(numbers(6:5:end))
    max(numbers(6:5:end))
end

if P.plotouts==1
    SMBplots(N,outtitle1,P.Glacier);
end
    
%% export geotifs
if P.exports==1
    writegeotiffs(N,THX,P)
end

disp('finished')
   
