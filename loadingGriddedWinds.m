function [wind_full] = loadingGriddedWinds(yrs, product)

%%% Function loads in the full gridded wind dataset

% Input: yrs, an array of years, e.g. 2003:2024

% Output: wind_full, a structure with fields for u, v, wind_curl, lat, lon, time, etc. 

%%%%%%%%%%%%%%%%%%
%%% ERA5 Winds %%%
%%%%%%%%%%%%%%%%%%

if strcmp(product, 'ERA5')

%%% Formatting u-wind data
for ii = 1:length(yrs)
    nc_name = 'C:\Users\jak279\OneDrive - Yale University\Research\Data\2D Data Products\ERA5 Winds\era5_u_' + string(yrs(ii)) + '.nc';
    data_all{ii} = ncread(nc_name, 'u10'); 
    time_all{ii} = datetime(ncread(nc_name, 'valid_time') + datenum(yrs(ii), 1, 1), 'ConvertFrom', 'datenum', 'TimeZone', 'UTC'); 
end
wind_full.u_speed = data_all{1};
wind_full.datetime = time_all{1};
if length(data_all) > 1
    for i = 2:length(data_all)
        wind_full.u_speed = cat(3, wind_full.u_speed, data_all{i});
        wind_full.datetime = cat(1, wind_full.datetime, time_all{i});
    end
end
lati = ncread(nc_name, 'latitude');
loni = ncread(nc_name, 'longitude'); 
ind = loni < 0;
loni(ind) = loni(ind) + 360;
[wind_full.lat, wind_full.lon] = meshgrid(lati,loni); % transforms vectors lati and loni into arrays 

%%% Formatting v-wind data
for ii = 1:length(yrs)
    nc_name = 'C:\Users\jak279\OneDrive - Yale University\Research\Data\2D Data Products\ERA5 Winds\era5_v_' + string(yrs(ii)) + '.nc';
    data_all{ii} = ncread(nc_name, 'v10'); 
end
wind_full.v_speed = data_all{1};
if length(data_all) > 1
    for i = 2:length(data_all)
        wind_full.v_speed = cat(3, wind_full.v_speed, data_all{i});
    end
end

%%% Computing curl
rho0 = 1.25;
cd = 0.00125;
del = 1;
taux = rho0 .* cd .* wind_full.u_speed.^2;
tauy = rho0 .* cd .* wind_full.v_speed.^2;
wind_full.u_stress = rho0 .* cd .* wind_full.u_speed .* abs(wind_full.u_speed);
wind_full.v_stress = rho0 .* cd .* wind_full.v_speed .* abs(wind_full.v_speed);
curl_arr = NaN(size(wind_full.u_speed));
for i =  2:length(loni)-1
    for j = 2:length(lati)-1
        dtauy = tauy(i+del,j,:) - tauy(i-del,j,:);
        x_m = gc_dist(wind_full.lat(i+del, j), wind_full.lon(i+del,j), wind_full.lat(i-del,j), wind_full.lon(i-del,j)) * 1000;
        dtauy_dx = dtauy ./ x_m;

        dtaux = taux(i,j+del,:) - taux(i,j-del,:);
        y_m = gc_dist(wind_full.lat(i,j+del), wind_full.lon(i,j+del), wind_full.lat(i,j-del), wind_full.lon(i,j-del)) * 1000;
        dtaux_dy = dtaux ./ y_m;

        dtauy_dx = reshape(dtauy_dx,[],1);
        dtaux_dy = reshape(dtaux_dy,[],1);

        curl_arr(i,j,:) = dtauy_dx - dtaux_dy;
    end
end
wind_curl = -curl_arr;
wind_full.curl = wind_curl;

%%%%%%%%%%%%%%%%%%
%%% NCEP Winds %%%
%%%%%%%%%%%%%%%%%%

elseif strcmp(product, 'NCEP')

    %%% Formatting u-wind data
for ii = 1:length(yrs)
    nc_name = 'C:\Users\jak279\OneDrive - Yale University\Research\Data\2D Data Products\NCEP Wind\uwnd.sig995.' + string(yrs(ii)) + '.nc';
    data_all{ii} = ncread(nc_name, 'uwnd'); 
    time_all{ii} = hours(ncread(nc_name, 'time')) + datetime(1800, 1, 1, 0, 0, 0, 'TimeZone', 'UTC'); 
end
wind_full.u_speed = data_all{1};
wind_full.datetime = time_all{1};
if length(data_all) > 1
    for i = 2:length(data_all)
        wind_full.u_speed = cat(3, wind_full.u_speed, data_all{i});
        wind_full.datetime = cat(1, wind_full.datetime, time_all{i});
    end
end
lati = double(ncread(nc_name, 'lat'));
loni = double(ncread(nc_name, 'lon')); 
ind = loni < 0;
loni(ind) = loni(ind) + 360;
[wind_full.lat, wind_full.lon] = meshgrid(lati,loni); % transforms vectors lati and loni into arrays 

%%% Formatting v-wind data
for ii = 1:length(yrs)
    nc_name = 'C:\Users\jak279\OneDrive - Yale University\Research\Data\2D Data Products\NCEP Wind\vwnd.sig995.' + string(yrs(ii)) + '.nc';
    data_all{ii} = ncread(nc_name, 'vwnd'); 
    time_all{ii} = hours(ncread(nc_name, 'time')) + datetime(1800, 1, 1, 0, 0, 0, 'TimeZone', 'UTC'); 
end
wind_full.v_speed = data_all{1};
if length(data_all) > 1
    for i = 2:length(data_all)
        wind_full.v_speed = cat(3, wind_full.v_speed, data_all{i});
    end
end

%%% Computing curl
rho0 = 1.25;
cd = 0.00125;
del = 1;
taux = rho0 .* cd .* wind_full.u_speed.^2;
tauy = rho0 .* cd .* wind_full.v_speed.^2;
wind_full.u_stress = rho0 .* cd .* wind_full.u_speed .* abs(wind_full.u_speed);
wind_full.v_stress = rho0 .* cd .* wind_full.v_speed .* abs(wind_full.v_speed);
curl_arr = NaN(size(wind_full.u_speed));
for i =  2:length(loni)-1
    for j = 2:length(lati)-1
        dtauy = tauy(i+del,j,:) - tauy(i-del,j,:);
        x_m = gc_dist(wind_full.lat(i+del, j), wind_full.lon(i+del,j), wind_full.lat(i-del,j), wind_full.lon(i-del,j)) * 1000;
        dtauy_dx = dtauy ./ x_m;

        dtaux = taux(i,j+del,:) - taux(i,j-del,:);
        y_m = gc_dist(wind_full.lat(i,j+del), wind_full.lon(i,j+del), wind_full.lat(i,j-del), wind_full.lon(i,j-del)) * 1000;
        dtaux_dy = dtaux ./ y_m;

        dtauy_dx = reshape(dtauy_dx,[],1);
        dtaux_dy = reshape(dtaux_dy,[],1);

        curl_arr(i,j,:) = dtauy_dx - dtaux_dy;
    end
end
wind_curl = -curl_arr;
wind_full.curl = wind_curl;
end

end