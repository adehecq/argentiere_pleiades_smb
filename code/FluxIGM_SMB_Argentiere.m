%FluxIGM_SMB_Argentiere.m - Master script to estimate surface mass balance
%distribution for a glacier from inputs of IGM flux inversion and thinning
%inspired from the scripts by Miles et al, 2021.
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

% dh folder
dh_pleiades = 1;
dh_hugonnet = 0;

% thx constrained by obs or not
tcons = 1;
tuncons = 0;

%% settings
P.plotouts=1; %output plots or not
P.exports=0; %save geotiffs or not
P.exports_indiv = 0; % Save every single SMB geotiff or not
P.DX = 20; %resolution to run calculations at 
P.segdist=300; %effective linear distance between flowbands
P.Vfilter=0; %swtich to smooth velocity data or not
P.dhfilter=0; %switch to peform 2x 3-sigma outlier removal (from overall dataset - only if erroneous pixels are common)
P.FDIVfilter=0; %switch to apply simply Gaussian filter to thickness data (useful for field thickness measurements)
P.umult=0; %switch for column-average velocity [0.8-1] is physical range. '0' estimates based on FDIV dist. '2' estimates for each pixel. 
P.Vreproj=0; %0 if velocities are provided oriented in the correct coordinate system, 1 if they are in the source data projection, [2 if they are true north], [3 to determine from slope]
P.fdivfilt=1; %use VanTricht gradient filters (2), just flux filter (1) or not at al (0)
P.uncertainty=0; %0 for 'simple run without uncertainty, otherwise N runs; note that you will then need to setup perturbations within the N structure
P.Valongslope=0; %1 for velocity along slope, 0 for 'normal' (i.e. as in Van Tricht, 2021; Miles 2021) 

%% inputs     

% dh
if dh_pleiades
    DH.path = fullfile(homedir,'..\output\dh_results\tendancy_filt.tif');
    DH.mult=1; %scale to convert input units to m/a
    dh_suff = '_dhPleiades';
end
if dh_hugonnet
    DH.path = fullfile(homedir,'..\data\dh\Hugonnet_2010-2020_utm32N.tif');
    DH.mult=1; %scale to convert input units to m/a
    dh_suff = '_dhHugonnet';
end

% fdiv
if tcons
    FDIV.path = 'C:\Users\kneibm\Documents\CAIRN\Accu_Argentiere\IGM\ThxInversionGJ\2024_03_03_MK\full_50m\divflux.tif';
    FDIV_suff = '_fluxIGM_tcons';
end
if tuncons
    FDIV.path = 'C:\Users\kneibm\Documents\CAIRN\Accu_Argentiere\IGM\ThxInversionGJ\2024_03_04_MK\full_50m_nothxobs\divflux.tif';
    FDIV_suff = '_fluxIGM_tuncons';
end

% DEM (always the same)
DEM.path = fullfile(homedir,'..\output\dh_results\meanDEM-2017_02_15.tif');



%% initialization
addpath(genpath([homedir,'MilesSMB']))
addpath(genpath([homedir,'MilesSMB\cbrewer']))

if P.uncertainty > 0
    str = string(P.uncertainty);
    UNC_suff = ['_uncertainty_N=',str{:}];
else
    UNC_suff = '';
end

date = string(datetime("today","Format","uuuu-MM-dd"));
outdir = [resultsfolder,date{:},dh_suff,FDIV_suff,UNC_suff];
mkdir(outdir)

outtitle1 = [num2str(P.DX) 'mgrid'];
    
%% Load all input data, resample to common grid, etc
% load FDIV, determine LL bounding box
FDIV.data=geotiffread(FDIV.path); %read data
FDIV.info=geotiffinfo(FDIV.path); %read metadata
FDIV.R = FDIV.info.RefMatrix; %georeferencing matrix
[FDIV.xm,FDIV.ym] = pixcenters(FDIV.R,[FDIV.info.Height,FDIV.info.Width]);
[FDIV.xmG,FDIV.ymG] = meshgrid(FDIV.xm,FDIV.ym);
[FDIV.LatG,FDIV.LonG] = projinv(FDIV.info,FDIV.xmG(:),FDIV.ymG(:));

MASK0 = FDIV.data~=0;
        
%project glacier's bounding box into geographic coordinates (key to subset from larger datasets)
BBoxUTM = FDIV.info.BoundingBox;
BBoxLL = [min(FDIV.LonG(MASK0)) min(FDIV.LatG(MASK0));max(FDIV.LonG(MASK0)) max(FDIV.LatG(MASK0))];
    
% load, subset DEM
DEM.info=geotiffinfo(DEM.path); %read metadata
DEM.R = DEM.info.RefMatrix;
[DEM.PixelRegion,DEM.LatG,DEM.LonG,DEM.xmG,DEM.ymG] = subset_geo(DEM.info,BBoxLL); %subset DEM domain
DEM.data=imread(DEM.path,'PixelRegion',DEM.PixelRegion); %read subsetted DEM

% load, subset dH
DH.info=geotiffinfo(DH.path);
DH.R = DH.info.RefMatrix;
[DH.PixelRegion,DH.LatG,DH.LonG,DH.xmG,DH.ymG] = subset_geo(DH.info,BBoxLL);
DH.data=imread(DH.path,'PixelRegion',DH.PixelRegion).*DH.mult;

% clip to glacier outlines
outlines = shaperead('C:\Users\kneibm\Documents\CAIRN\Accu_Argentiere\argentiere_pleiades_smb\data\gis\outline_pleiades_2020-09-08.shp');
dgeot = geotiffinfo(FDIV.path);
[~,FDIV.data] = geotiffcrop_shp(outlines,FDIV.data,dgeot);

% compensate lost mass when cropping
FDIV.data = FDIV.data-nansum(FDIV.data(:))/sum(~isnan(FDIV.data(:)));

%% REPROJECT AND RESAMPLE ALL INPUTS (to FDIV coordinate system), derive MASK, etc
buffdist = 100; %expands the domain for subsetting around the glacier

%new coordinates based on the FDIV data, but at DX interval
DX = P.DX;
N.x3 = DX.*[(floor((FDIV.xm(1)-buffdist)/DX)):((ceil(FDIV.xm(end)+buffdist)/DX))];
N.y3 = DX.*[(ceil((FDIV.ym(1)+buffdist)/DX)):-1:(floor((FDIV.ym(end)-buffdist)/DX))];
[N.x3g,N.y3g] = meshgrid(N.x3,N.y3);

N.Rout = [0,-DX;DX,0;N.x3(1)-DX,N.y3(1)+DX]; %updated georeferencing matrix based on new domain
N.DX=DX;

N.FDIV = interp2(FDIV.xmG,FDIV.ymG,double(FDIV.data),N.x3g,N.y3g,'linear');
FDIV.data(isnan(FDIV.data))=0;%Nans result in contraction of the glacier domain after interpolation
N.MASK=N.FDIV>-100;

%resample DEM
[DEM.xN,DEM.yN] = projfwd(FDIV.info,DEM.LatG,DEM.LonG); %compute projected coordinates
N.DEM = griddata(DEM.xN(:),DEM.yN(:),double(DEM.data(:)),N.x3g(:),N.y3g(:),'natural');
N.DEM = reshape(N.DEM,size(N.x3g));

%resample dH
[DH.xN,DH.yN] = projfwd(FDIV.info,DH.LatG,DH.LonG); %compute projected coordinates
N.DH = griddata(DH.xN(:),DH.yN(:),double(DH.data(:)),N.x3g(:),N.y3g(:),'natural');
N.DH = reshape(N.DH,size(N.x3g));

N.FDIV(N.FDIV<-200) = NaN; % remove results from interpolation on the edges
N.fdivfilt=P.fdivfilt;
N.uncertainty=P.uncertainty;
N.Valongslope = P.Valongslope;
N.exports_indiv = P.exports_indiv;
N.Glacier = P.Glacier;
N.info.GeoTIFFTags.GeoKeyDirectoryTag = FDIV.info.GeoTIFFTags.GeoKeyDirectoryTag;


%% fill gaps in reprojected inputs if needed

%gap-fill dH based on elevation
N.DH0=N.DH; %unfilled values
N.DH(N.MASK)= index_nanfill(N.DH(N.MASK),N.DEM(N.MASK));
N.DH((N.MASK==0))=0;

%% zone segmentation

N.zones = uint16(segment_Gmask_slope2(N.DEM,N.MASK,P.DX,P.segdist)); %segments glacier mask and DEM into approximately uniformly-spaced elevation bands

%% Uncertainties on velocities and dh
N.sigdH = 0.07;
N.sigdensity_mixedzone = 0.15;
N.density_mixedzone = 0.75;
N.umult = 1;

%% SMB calculations
cd(outdir)

if P.uncertainty==0
    % density corrections
    N.Qdensity = 0.9; %900 kg m3 everywhere; 0.9 is actually the specific gravity
    N.Hdensity = 0.9.*N.MASK; %initial value everywhere of 900kg m3; 0.9 is actually the specific gravity
        
    %grid implementation of density correction
    ind1=(N.FDIV>0); 
    ind2=(N.DH>0);
    ind3=(abs(N.DH)>abs(N.FDIV));
    ind4 = N.DEM>nanmedian(N.DEM(N.MASK));
        
    N.Hdensity(ind1&~ind2)=0.9; %thinning and emergence = melt
    N.Hdensity(~ind1&ind2)=0.6; %thickening and submergence = acc
    N.Hdensity(ind1&ind2&ind3)=0.6; %emergence and thickening, more thickening - acc
    N.Hdensity(ind1&ind2&~ind3)=N.density_mixedzone; %emergence and thickening, more thickening - mixed
    N.Hdensity(~ind1&~ind2&ind3)=0.9; %submergence and thinning, more thinning - melt
    N.Hdensity(~ind1&~ind2&~ind3)=N.density_mixedzone; %submergence and thinning, less thinning - mixed
    N.Hdensity((ind4==0)&N.MASK)=0.9; % below median elevation
       
    % SMB
    N.SMB = N.Hdensity.*N.DH+N.Qdensity.*N.FDIV; %continuity equation. note that 'density' terms are actually specific gravity
    N.SMBz2= zonal_aggregate(N.zones,N.SMB); %aggregates values in the zone - simple mean
    
    %mask before plotting
    N.DH((N.MASK==0))=NaN;
    N.FDIV((N.MASK==0))=NaN;
    N.SMB((N.MASK==0))=NaN;
else
    %  nr random or uniform draws for random (sig density & dH) uncertainty
    nr = N.uncertainty;
    rns = [rand(nr,2) randn(nr,4)];
    prns=[cdf('Uniform',rns(:,1:2),0,1) cdf('Normal',rns(:,3:6),1,1)];
        
    % initialize vars
    tFDIV=zeros([size(N.DH),nr]);
    tSMB=zeros([size(N.DH),nr]);
    tHdens=zeros([size(N.DH),nr]);
    
    dx = mode(diff(N.x3));
    dy = mode(diff(N.y3));
    
    % loop through runs
    tic
    for irun=1:nr
        C=N;

        %	determine MC adjustments
        C.dHadd=icdf('normal',(prns(irun,5)),0,N.sigdH); %extra to determine random adjustments
        if N.sigdensity_mixedzone ~= 0 
            C.densadd=icdf('uniform',(prns(irun,2)),-N.sigdensity_mixedzone,N.sigdensity_mixedzone);
        else
            C.densadd=0;
        end

        if ~isnan(C.dHadd)
            C.DH=C.DH+C.dHadd;
        end
        if ~isnan(C.densadd)
            C.density_mixedzone=C.density_mixedzone+C.densadd;
        end

        % density corrections
        C.Qdensity = 0.9; %900 kg m3 everywhere; 0.9 is actually the specific gravity
        C.Hdensity = 0.9.*C.MASK; %initial value everywhere of 900kg m3; 0.9 is actually the specific gravity
            
        %grid implementation of density correction
        ind1=(C.FDIV>0); 
        ind2=(C.DH>0);
        ind3=(abs(C.DH)>abs(C.FDIV));
        ind4 = C.DEM>nanmedian(C.DEM(C.MASK));
            
        C.Hdensity(ind1&~ind2)=0.9; %thinning and emergence = melt
        C.Hdensity(~ind1&ind2)=0.6; %thickening and submergence = acc
        C.Hdensity(ind1&ind2&ind3)=0.6; %emergence and thickening, more thickening - acc
        C.Hdensity(ind1&ind2&~ind3)=C.density_mixedzone; %emergence and thickening, more thickening - mixed
        C.Hdensity(~ind1&~ind2&ind3)=0.9; %submergence and thinning, more thinning - melt
        C.Hdensity(~ind1&~ind2&~ind3)=C.density_mixedzone; %submergence and thinning, less thinning - mixed
        C.Hdensity((ind4==0)&C.MASK)=0.9; % below median elevation
           
        % SMB
        C.SMB = C.Hdensity.*C.DH+C.Qdensity.*C.FDIV; %continuity equation. note that 'density' terms are actually specific gravity
        C.SMBz2= zonal_aggregate(C.zones,C.SMB); %aggregates values in the zone - simple mean
        
        %mask before plotting
        C.DH((C.MASK==0))=NaN;
        C.FDIV((C.MASK==0))=NaN;
        C.SMB((C.MASK==0))=NaN;
        
        % index outputs into stack
        tFDIV(:,:,irun)=C.FDIV;
        tSMB(:,:,irun)=C.SMB;
        tHdensity(:,:,irun)=C.Hdensity;
        
        
        if C.exports_indiv == 1
            geotiffwrite([C.Glacier '_SMB_rho=',num2str(C.density_mixedzone),'dhadd=',num2str(C.dHadd),'.tif'],...
                C.SMB,C.Rout,'GeoKeyDirectoryTag',C.info.GeoTIFFTags.GeoKeyDirectoryTag,'TiffTags' , struct ( 'Compression' , Tiff.Compression.Deflate))
        end

        clear C
    end
    toc
    
    % postprocess MC stack into N
    N.SMB=nanmean(tSMB,3);
    N.FDIV=nanmean(tFDIV,3);
    N.Hdensity=nanmean(tHdensity,3);
    
    N.SMBu=nanstd(tSMB,[],3);
    N.FDIVu=nanstd(tFDIV,[],3);
    N.Hdensityu=nanstd(tHdensity,[],3);
    
    N.Qdensity = 0.9; %900 kg m3 everywhere; 0.9 is actually the specific gravity

    N.z2fdiv = zonal_aggregate(N.zones,N.FDIV); % aggregates values in the zone - simple mean excluding NaNs; same result as perimeter integration
    N.z2DH = zonal_aggregate(N.zones,N.DH); % aggregates values in the zone - simple mean
    N.SMBz2= zonal_aggregate(N.zones,N.SMB); %aggregates values in the zone - simple mean
    
    N.z2fdiv_unc = zonal_aggregate_v2(N.zones,N.FDIVu,'rssn'); % aggregates values in the zone - simple mean excluding NaNs; same result as perimeter integration
end
    
%% export geotifs
if FDIV.info.GeoTIFFTags.GeoKeyDirectoryTag.GTRasterTypeGeoKey==1 %Matlab geotiffexport off by 1 pixel
    N.Rout(3,1)=N.Rout(3,1)-N.DX;
    N.Rout(3,2)=N.Rout(3,2)+N.DX;
end
Glacier=P.Glacier;

% FDIV
geotiffwrite([Glacier '_FDIV.tif'],N.FDIV,N.Rout,'GeoKeyDirectoryTag',FDIV.info.GeoTIFFTags.GeoKeyDirectoryTag,'TiffTags' , struct ( 'Compression' , Tiff.Compression.Deflate))
geotiffwrite([Glacier '_FDIVu.tif'],N.FDIVu,N.Rout,'GeoKeyDirectoryTag',FDIV.info.GeoTIFFTags.GeoKeyDirectoryTag,'TiffTags' , struct ( 'Compression' , Tiff.Compression.Deflate))

% SMB
geotiffwrite([Glacier '_SMB.tif'],N.SMB,N.Rout,'GeoKeyDirectoryTag',FDIV.info.GeoTIFFTags.GeoKeyDirectoryTag,'TiffTags' , struct ( 'Compression' , Tiff.Compression.Deflate))
geotiffwrite([Glacier '_SMBu.tif'],N.SMBu,N.Rout,'GeoKeyDirectoryTag',FDIV.info.GeoTIFFTags.GeoKeyDirectoryTag,'TiffTags' , struct ( 'Compression' , Tiff.Compression.Deflate))

   

disp('finished')
   
