% tendancy_dh.m: script to derive mean and filtered elevation change from a set of DEMs.
% Authors: Auguste Basset, Amaury Dehecq, Marin Kneib

cd C:\Users\kneibm\Documents\CAIRN\Accu_Argentiere\argentiere_pleiades_smb\code

%% Initialize
clear
close all

% path to useful codes
addpath('../UsefulCodes/topotoolbox')
addpath('../UsefulCodes/ImGRAFT')
addpath('../UsefulCodes/')

% folders
datafolder='../data/dh/DEMs/';
resultsfolder = '../output/dh_results/';

% list of DEMs
date = {'20210815','20200930','20200917','20200809','20190825',...
    '20180908','20180812','20170814','20160928','20160807', '20150830',...
    '20130920','20120819'};

% convert to datetime and calculate dt
date_t = datetime(date,"Format",'yyyyMMdd');

% load data (shapefiles)
Argentiere_shp = shaperead('../data/gis/outline_pleiades_2020-09-08.shp');
OffGl = shaperead('../data/gis/velocity_off_glacier.shp');

%% Load reference DEM
[refDEM,xdem,ydem,~] = geoimread([datafolder,'20160807_DEM_4m_shift_H-V_clip.tif']); % DEM with smallest extents
dgeot = geotiffinfo([datafolder,'20160807_DEM_4m_shift_H-V_clip.tif']);
key = dgeot.GeoTIFFTags.GeoKeyDirectoryTag;

[Xdem,Ydem] = meshgrid(xdem,ydem);

refDEM(refDEM<500) = NaN;

%% resample DEMs to same extents & store in 3D array

% Allocate
DEM3d = zeros(size(refDEM,1),size(refDEM,2),length(date));

% interpolate each DEM
for ii = 1:length(date)
    [DEMt,xt,yt,~] = geoimread([datafolder,date{ii},'_DEM_4m_shift_H-V_clip.tif']);
    DEMt(DEMt<500) = NaN;
    [Xt,Yt] = meshgrid(xt,yt);
    DEMt = interp2(Xt,Yt,double(DEMt),Xdem,Ydem,'bilinear');
    DEM3d(:,:,ii) = DEMt;
end


%% dh/dt for each pixel (linear regression)

% time
time = convertTo(date_t, 'posixtime'); % converts each element of time to the number of seconds that have elapsed since the epoch of January 1, 1970, 00:00:00 UTC

% initialize
tendancy = NaN(size(refDEM,1),size(refDEM,2));
r2 = NaN(size(refDEM,1),size(refDEM,2));

% Loop restricted on zone of interest
for ii = 1:size(refDEM,1)
    for jj = 1:size(refDEM,2)
        time_temp = time';
        pix_values = DEM3d(ii,jj,:);
        pix_values = pix_values(:);

        % remove Nans in time series
        time_temp(isnan(pix_values)) = [];
        pix_values(isnan(pix_values)) = [];

        % only do regression if less than 5 NaNs
        if length(pix_values) >= length(time)-5
            % only do regression if more than 5 years between first and
            % last data point
            if time_temp(1)-time_temp(end) >= 5*365*24*60*60
                % linear regression
                X = [ones(length(time_temp),1) time_temp];
                b = X\pix_values;
                yCalc = X*b;

                % calculate residuals and remove if too high
                condition = abs(yCalc - pix_values) > 3*std(yCalc - pix_values);
                time_temp(condition) = [];
                pix_values(condition) = [];

                % second linear regression
                X = [ones(length(time_temp),1) time_temp];
                b = X\pix_values;
                yCalc = X*b;

                % outputs
                tendancy(ii,jj) = b(2);
                r2(ii,jj) = 1 - sum((pix_values - yCalc).^2)/sum((pix_values - mean(pix_values)).^2);
            end 
        end
    end 
end

% convert m/s to m/yr
tendancy = tendancy*3600*24*365;
tendancy(tendancy<-10) = NaN;
tendancy(tendancy>5) = NaN;

%% mask

% rasterize mask
[Argentiere_r,~] = geotiffcrop_shp(Argentiere_shp,refDEM,dgeot);

% slope (to remove steep off-glacier)
avgDEM = nanmean(DEM3d,3); % average of all DEMs without accounting for NANs -> for slope calculation only 
DEM_gridobj = GRIDobj(xdem,ydem,avgDEM);
Slope = gradient8(DEM_gridobj,'deg');
slope_mask = Slope.Z;
slope_mask(slope_mask>40) = NaN;
slope_mask(~isnan(slope_mask)) = 1;
slope_mask(isnan(slope_mask)) = 0;
slope_mask = bwareaopen(slope_mask,10); % remove small patches

%% Plotting
cm = [1 0 0; 1 1 1; 0 0 1];                     % Basic Colormap
cmi = interp1([-10; 0; 10], cm, (-10:0.1:10));           % interpolated Colormap
cmi2 = interp1([-1; 0; 1], cm, (-1:0.1:1));

figure
imagesc(tendancy)
colormap(gca, cmi)
caxis([-10 10])

figure
imshow(r2)
colormap(gca, cmi2)
caxis([-10 10])

%% gap filling within glacier outlines
tendancy_gl = Argentiere_r.*tendancy; % for gap filling only consider dh on the glacier

tofill = isnan(tendancy_gl) & ~isnan(Argentiere_r); %masked results

tendancy2=nanmedfilt2(tendancy_gl,[3 3]);
tendancy_filling=griddata(Xdem(~isnan(tendancy_gl)),Ydem(~isnan(tendancy_gl)),tendancy2(~isnan(tendancy_gl)),Xdem(tofill),Ydem(tofill),'cubic');

tendancy_gl(tofill)=tendancy_filling;

% smooth data with 3x3 median filter
tendancy_gl = nanmedfilt2(tendancy_gl,[3 3]);

% fill holes
tendancy_filt = tendancy;
tendancy_filt(isnan(tendancy) & ~isnan(Argentiere_r)) = tendancy_gl(isnan(tendancy) & ~isnan(Argentiere_r));

%% uncertainties (mean, median and std of off-glacier terrain)
[~,OffGl_r] = geotiffcrop_shp(OffGl,tendancy,dgeot);
OffGl_r(Slope.Z>30) = NaN;
nanmean(OffGl_r(:))
nanstd(OffGl_r(:))

dh_mean = Argentiere_r.*tendancy_filt;
nanmean(dh_mean(:))
%% Output

geotiffwrite([resultsfolder,'tendancy.tif'], tendancy, dgeot.RefMatrix, 'GeoKeyDirectoryTag',key);
geotiffwrite([resultsfolder,'r2.tif'], r2, dgeot.RefMatrix, 'GeoKeyDirectoryTag',key);
geotiffwrite([resultsfolder,'tendancy_filt.tif'], tendancy_filt, dgeot.RefMatrix, 'GeoKeyDirectoryTag',key);


%% extrapolated DEM
refDEM2 = DEM3d(:,:,4); % reference DEM is from 09/08/2020 (less holes)
medDEM = refDEM2-tendancy_filt*days(datetime(2020,08,09)-datetime(2017,02,15))/365;

% gap filling
tofill = isnan(medDEM); %masked results
medDEM_filling=griddata(Xdem(~isnan(medDEM)),Ydem(~isnan(medDEM)),medDEM(~isnan(medDEM)),Xdem(tofill),Ydem(tofill),'cubic');
medDEM(tofill)=medDEM_filling;

geotiffwrite([resultsfolder,'meanDEM-2017_02_15.tif'], medDEM, dgeot.RefMatrix, 'GeoKeyDirectoryTag',key);


