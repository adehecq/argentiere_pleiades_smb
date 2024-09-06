function [shp_raster,crop_raster] = geotiffcrop_shp(shapefile,raster,dgeot)
    %% This function 'crops' a geotiff to the extents of a shapefile - pixels outside shapefile are NaNs
    
    % INPUTS
    %   shapefile: shp file opened with matlab shaperead function
    %   raster: raster file used as template to rasterize shapefile (opened
    %   with geotiffread or ImGRAFT function geoimread for ex)
    %   dgeot: raster geo info file (opened with geotiffinfo)
    % OUTPUTS
    %   shp_raster: rasterized shapefile
    %   crop_raster: cropped raster (NaNs outside the shapefile bounds)
    
    % Author: Marin Kneib
    % Work address: Swiss Federal Research Institute WSL
    % Email: marin.kneib@wsl.ch
    
    shp_raster = zeros(size(raster,1),size(raster,2));
    for j = 1:length(shapefile)
        % Find NaN values
        nan_idx = find(isnan(shapefile(j).X));
        % main polygon defined by all x,y until first NaN
        X = shapefile(j).X(1:nan_idx(1,1)-1);
        Y = shapefile(j).Y(1:nan_idx(1,1)-1);
        [row,col] = map2pix(dgeot.RefMatrix,X,Y);
        row=row(~isnan(row));
        col=col(~isnan(col));
        shp_raster = shp_raster+poly2mask(col,row,size(raster,1),size(raster,2));
        % Add holes defined by coordinates after first NaN
        if length(nan_idx)>1
            for i = 2:length(nan_idx)
                X = shapefile(j).X(nan_idx(1,i-1)+1:nan_idx(1,i)-1);
                Y = shapefile(j).Y(nan_idx(1,i-1)+1:nan_idx(1,i)-1);
                [row,col] = map2pix(dgeot.RefMatrix,X,Y);
                row=row(~isnan(row));
                col=col(~isnan(col));
                shp_raster = shp_raster-poly2mask(col,row,size(raster,1),size(raster,2));
            end
        end
    end
    shp_raster(shp_raster==0)=NaN;
    shp_raster(~isnan(shp_raster))=1;
    crop_raster = double(shp_raster).*double(raster);
end