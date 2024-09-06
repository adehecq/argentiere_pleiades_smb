function N=FluxCalcsSimple(N)
%FluxCalcsSimple.m - Function to estimate surface mass balance
% distribution for a glacier from inputs of ice thickness, thinning, and
% velocity, based on the continuity equation (see, e.g. Bisset et al, 2020: https://doi.org/10.3390/rs12101563)
%
%
% Other m-files required: C2xyz.m, index_nanfill.m, remove_small_zones.m,
% segment_Gmask_slope2.m, subset_geo.m, through_fluxes.m, zonal_aggregate.m
%
% Author: Evan Miles
% Work address: Swiss Federal Research Institute WSL
% Email: evan.miles@wsl.ch
% Sep 2019; Last revision: 16-June-2020
% Adapted by Marin Kneib for SMB inversion on Argentière Glacier
    
%% FLUX AND EMERGENCE CALCS 
    N.Smean=N.umult.*N.S;
    N.Umean=N.umult.*N.U;
    N.Vmean=N.umult.*N.V;

    %% calculate DEM gradient
    dx = mode(diff(N.x3));
    dy = mode(diff(N.y3));

    if N.Valongslope==1
        [N.gradX,N.gradY] = gradient(N.DEM,abs(dx));
        alphax = atan(N.gradX);
        alphay = atan(N.gradY);
        N.Umean = N.Umean./cos(alphax);
        N.Vmean = N.Vmean./cos(alphay);
        N.Umean = nanmedfilt2(N.Umean,3);
        N.Vmean = nanmedfilt2(N.Vmean,3);
        N.Smean = sqrt(N.Vmean.^2+N.Umean.^2);
    end

    
    %% calculate flux divergence per pixel
        
    %initialize
    N.FDIV = zeros(size(N.THX));
    N.FDIVx = zeros(size(N.THX));
    N.FDIVy = zeros(size(N.THX));

    %pixel-based flux magnitude
    N.FLUX = N.Smean.*N.THX;
    N.FLUX(N.FLUX<=0)=NaN;
        
    N.Umean(isnan(N.Umean))=0;
    N.Vmean(isnan(N.Vmean))=0;
    N.THX(isnan(N.THX))=0;
    
    %VanTricht2021 formulation of flux divergence - numerically equivalent to above
    [N.dUdx,~]=gradient(N.Umean,dx);%dx to normalize to pixels
    [~,N.dVdy]=gradient(N.Vmean,dy);N.dVdy=N.dVdy; 
    [N.dHdx,N.dHdy]=gradient(N.THX,dx,dy);

    N.dUdx((N.MASK)==0)=NaN;
    N.dVdy((N.MASK)==0)=NaN;
    N.dHdx((N.MASK)==0)=NaN;
    N.dHdy((N.MASK)==0)=NaN;
    
    N.FDIV = (N.Umean.*N.dHdx+N.Vmean.*N.dHdy+N.THX.*N.dUdx+N.THX.*N.dVdy); %VanTrich eq 5 
    if N.fdivfilt==2 %use spatially-filtered gradients
        N.FDIV0=N.FDIV;
        N.dUdx0=N.dUdx;N.dVdy0=N.dVdy;
        N.dUdx = flowfilter_varsig(N,N.dUdx);
        N.dVdy = flowfilter_varsig(N,N.dVdy);
        N.dHdx0=N.dHdx;N.dHdy0=N.dHdy;
        N.dHdx = flowfilter_varsig(N,N.dHdx);
        N.dHdy = flowfilter_varsig(N,N.dHdy);
        N.FDIV = (N.Umean.*N.dHdx+N.Vmean.*N.dHdy+N.THX.*N.dUdx+N.THX.*N.dVdy); %VanTrich eq 5 
    end
        
    %trim to mask
    N.FDIV(N.MASK==0)=NaN;
    
    
    %% filter FDIV 
    if N.fdivfilt>=1
        % filter with varying search distance, respecting NaNs
        N.FDIV0=N.FDIV;
        [N.FDIV,N.delFDIV] = flowfilter_varsig(N,N.FDIV0);
    end
    
    %% aggregate variables over zones
    N.z2fdiv = zonal_aggregate(N.zones,N.FDIV); % aggregates values in the zone - simple mean excluding NaNs; same result as perimeter integration
    N.z2DH = zonal_aggregate(N.zones,N.DH); % aggregates values in the zone - simple mean

    %% density corrections
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
   
     %% SMB
    
    N.SMB = N.Hdensity.*N.DH+N.Qdensity.*N.FDIV; %continuity equation. note that 'density' terms are actually specific gravity
    N.SMBz2= zonal_aggregate(N.zones,N.SMB); %aggregates values in the zone - simple mean

    %mask before plotting
    N.DH((N.MASK==0))=NaN;
    N.U((N.MASK==0))=NaN;
    N.V((N.MASK==0))=NaN;
    N.THX((N.MASK==0))=NaN;
    N.FDIV((N.MASK==0))=NaN;
    N.FDIVx((N.MASK==0))=NaN;
    N.FDIVy((N.MASK==0))=NaN;
    N.SMB((N.MASK==0))=NaN;
    N.zDH((N.MASK==0))=NaN;
    N.SMBz2((N.MASK==0))=NaN;
    N.z2DH((N.MASK==0))=NaN;
    N.z2fdiv((N.MASK==0))=NaN;