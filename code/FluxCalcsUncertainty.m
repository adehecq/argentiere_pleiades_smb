function N=FluxCalcsUncertainty(N) 
%FluxCalcsUncertainty.m - Function to estimate surface mass balance
% distribution for a glacier from inputs of ice thickness, thinning, and
% velocity, based on the continuity equation (see, e.g. Bisset et al, 2020: https://doi.org/10.3390/rs12101563)
% In addition to inputs for FluxCalcsSimple, requires N to contain variables of the following as single value or grid in variable units:
%       sigUmult - uncertainty of column-averaged velocity factor
%       sigU - uncertainty of u component of velocity (m/a)
%       sigV - uncertainty of v component of velocity (m/a)
%       sigTHX - uncertainty of ice thickness in terms of systematic bias (percent)
%       sigTHX2 - uncertainty of ice thickness in terms of random error (m)
%       sigDH - uncertainty of thinning (m/a)
% The code then runs through nr randomly-chosen parameter sets.
%
% Other m-files required: C2xyz.m, index_nanfill.m, remove_small_zones.m,
% segment_Gmask_slope2.m, subset_geo.m, through_fluxes.m, zonal_aggregate.m
%
% Author: Evan Miles
% Work address: Swiss Federal Research Institute WSL
% Email: evan.miles@wsl.ch
% Sep 2019; Last revision: 16-June-2020
% Adapted by Marin Kneib for SMB inversion on Argentière Glacier

%% baseline calculations for unperturbed case
    N.Smean=N.umult.*N.S;
    N.Umean=N.umult.*N.U;
    N.Vmean=N.umult.*N.V;
    
    %pixel-based flux magnitude
    N.FLUX = N.Smean.*N.THX(:,:,1);N.FLUX(N.FLUX<=0)=NaN;
    
%%  nr random or uniform draws for systematic (THXscale,umult,Vscale) and random (10m thx, sig m/a vel, sig density & dH) uncertainty
    nr = N.uncertainty;
    rns = [rand(nr,2) rand(nr,4)];
    prns=[cdf('Uniform',rns(:,1:6),0,1)];
        
% initialize vars
    tFDIV=zeros([size(N.DH),nr]);
    tSMB=zeros([size(N.DH),nr]);
    tHdens=zeros([size(N.DH),nr]);
    
    dx = mode(diff(N.x3));
    dy = mode(diff(N.y3));
    
%     loop through runs
    tic
    for irun=1:nr
        C=N;

        %	determine MC adjustments
        if N.sigUmult ~= 0
            C.umult=icdf('uniform',prns(irun,1),N.umult-N.sigUmult,N.umult+N.sigUmult);
        else
            C.umult = N.umult;
        end
        if N.sigU ~= 0
            C.Uadd=icdf('normal',(prns(irun,3)),0,N.sigU); %extra to determine random adjustments
            C.Vadd=icdf('normal',(prns(irun,4)),0,N.sigV); %extra to determine random adjustments
        else
            C.Uadd=0;
            C.Vadd=0;
        end
        if N.sigdH ~= 0
            C.dHadd=icdf('normal',(prns(irun,5)),0,N.sigdH); %extra to determine random adjustments
        else
            C.dHadd=0;
        end
        if N.sigdensity_mixedzone ~= 0
            C.densadd=icdf('uniform',(prns(irun,2)),-N.sigdensity_mixedzone,N.sigdensity_mixedzone);
        else
            C.densadd=0;
        end

        %   determine current inputs to flux calcs
        C.U=(C.U+C.Uadd);
        C.V=(C.V+C.Vadd);

        if size(C.THX,3)>0
            thx_idx = int32(1 + (size(C.THX,3)-1) .* rand(1));
        else
            thx_idx = 1;
        end
        C.THX=C.THX(:,:,thx_idx); % draw index in [1 100] range
        C.THX(C.THX<0) = 5;
        
        C.DH=C.DH+C.dHadd;
        C.density_mixedzone=C.density_mixedzone+C.densadd;
        
        % calculate fluxes and SMB
        C=FluxCalcsSimple(C); 
        
        % index outputs into stack
        tFDIV(:,:,irun)=C.FDIV;
        tSMB(:,:,irun)=C.SMB;
        tHdensity(:,:,irun)=C.Hdensity;
        
        
        if N.exports_indiv == 1
            geotiffwrite([C.Glacier '_SMB_rho=',num2str(C.density_mixedzone),...
                '_gamma=',num2str(C.umult),'_thx=',num2str(thx_idx),...
                '_uadd=',num2str(C.Uadd),'_vadd=',num2str(C.Vadd),'dhadd=',num2str(C.dHadd),'.tif'],...
                C.SMB,C.Rout,'GeoKeyDirectoryTag',C.info.GeoTIFFTags.GeoKeyDirectoryTag,'TiffTags' , struct ( 'Compression' , Tiff.Compression.Deflate))
        end

        clear C
    end
    toc
    %% postprocess MC stack into N
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
